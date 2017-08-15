#!/usr/bin/perl -w
# Script to extract counts from VCF file
#  and SNP positions on that chromosome
#
#

use File::Basename;
use Scalar::Util qw(looks_like_number);
use strict;

my $USAGE = "[VCF file][Order][Working Directory]\n";

# Read file provided from the commandline
unless( $ARGV[0] && $ARGV[1] && $ARGV[2])
{
  print $USAGE;
  exit;
}

# Read Commandline
my $vcf_file  = $ARGV[0];
my $order     = $ARGV[1];
my $work_dir  = $ARGV[2]; 

$vcf_file = fileparse($vcf_file);

# $_ is the default variable name
# $ keeps track of the current line number

# Set the character which will be used to indicate the end of a line
local $/ = "\n";

# Create Output file
my $output_file = "snvcounts_outputfile";
my $basepath = $work_dir . "/" . $output_file;
my $vcfpath  = $work_dir . "/" . $vcf_file;
my $medianpath = $work_dir . "/median_outputfile";

# Open filehandle for write access
open(my $fh, '>>', $basepath);

# Open filehandle for read access
open my $filehandle, '<', $vcfpath;

# Analysis will be set on chromsomes 1-22, X, and Y

my %valid_chrom;

for ( my $i = 1; $i <= 22; $i++)
{
  $valid_chrom{"chr$i"} = 1;
}

$valid_chrom{'chrX'} = 1;
$valid_chrom{'chrY'} = 1;


my @median_values = ();

my $header = "NA";
my @header_array = ();
my $header_length = "NA";
my $snpcount = 0;

# Loop through each line:
while (<$filehandle>)
{
  # The text of the line includeing the linebreak, is in the variable $_
  # Strip the linebreak character at the end
  chomp;

  # operate on the line
  if( substr($_, 0, 1) eq "#")
  {
    if ( substr($_, 1, 5) eq "CHROM")
    {
      if( $order eq 'TN' || $order eq 'NT')
      {
        print $fh "Chr\tPos\tTumorMutant\tTumorTotal\tNormalMutant\tNormalTotal\n";
      }
      else
      {
        print "unknown order: $order\n";
        exit;
      }
    }
  }
  else
  {
    my @line_array = split("\t", $_);
    my $chromosome = StandardName($line_array[0]);
    my $pos    = $line_array[1];
    my $id     = $line_array[2];
    my $ref    = $line_array[3];
    my $alt    = $line_array[4];
    my $format = $line_array[8];
    my $data1  = $line_array[9];
    my $data2 = 0;
    if ( scalar(@line_array) > 10)
    {
      $data2 = $line_array[10];
    }
    else
    {
      print "Do not have Paired Tumor/Normal Data\n";
      exit;
    }
   
    # This flag is used to ensure SNP ( not indel at location) and 
    # both Tumor and Normal data are present
    my $dataflag = 1;

    # Verify tumor and Normal data are both present
    if ( $data1 eq "." || $data2 eq ".")
    {
      $dataflag = 0;
    }

    # Verify that we have a SNP location
    if ( length($ref) > 1 || length($alt) > 1)
    {
      $dataflag = 0;
    }

    if ( $dataflag)
    {
      if (exists $valid_chrom{$chromosome})
      {
        my $fieldflag = 1;
        my $AD_exists = 0;
        my $RO_exists = 0;
        my $AO_exists = 0;
        
        # identify which tags are present in the format field
        my @field_array = split(":", $format);
        foreach my $element( @field_array)
        {
          if ( $element eq 'AD')
          {
            $AD_exists = 1;
          }
          if ( $element eq 'RO')
          {
            $RO_exists = 1;
          }
          if ( $element eq 'AO')
          {
            $AO_exists = 1;
          }
        }
        # identify values for present tags
        my @value1_array = split(":", $data1); 
        my %value1_hash;
        my $size = scalar(@field_array);
        for ( my $i = 0; $i < $size; $i++)
        {
          $value1_hash{$field_array[$i]} = $value1_array[$i];
        }
          
        my @value2_array = split(":", $data2);
        my %value2_hash;
        for ( my $i = 0; $i < $size; $i++)
        {
          $value2_hash{$field_array[$i]} = $value2_array[$i];
        }
        
        my @counts1 = ParseCounts(\%value1_hash, $AD_exists, $RO_exists, $AO_exists);
        my @counts2 = ParseCounts(\%value2_hash, $AD_exists, $RO_exists, $AO_exists);

        if ( $counts1[0] == -1 || $counts1[1] == -1)
        {
          $fieldflag = 0;
        } 
        if ( $counts2[0] == -1 || $counts2[1] == -1)
        {
          $fieldflag = 0;
        }
        if ($fieldflag)
        {
          $snpcount++;
          
          # The first case
          my $ref1_count = $counts1[0];
          my $mut1_count = $counts1[1];
          my $tot1_count = $ref1_count + $mut1_count;

          # The second case
          my $ref2_count = $counts2[0];
          my $mut2_count = $counts2[1];
          my $tot2_count = $ref2_count + $mut2_count;
          if( $order eq 'TN')
          {
            print $fh "$chromosome\t$pos\t$mut1_count\t$tot1_count\t$mut2_count\t$tot2_count\n";
            push ( @median_values, $tot2_count);
          }
          elsif( $order eq 'NT')
          {
            print $fh "$chromosome\t$pos\t$mut2_count\t$tot2_count\t$mut1_count\t$tot1_count\n";
            push ( @median_values, $tot1_count);
          }
        }
      }
    }
  }
  
  # error checking
  if (m/^ERROR/)
  {
    warn "Error on line $. - skipping rest of file";
    last;
  }
}

close $filehandle;
close $fh;

if (scalar(@median_values) > 0)
{
  my $median = median(@median_values);

  open(my $fh2, '>', $medianpath);
  print $fh2 "$median \n";
  close $fh2;
}

print "$snpcount\n";

exit;

sub StandardName
{
  my $chr = uc $_[0];
  my $string_length = length($chr);
  if ( $string_length < 3)
  {
    $chr = "chr" . $chr;
  }
  elsif ( $string_length == 4)
  {
    my $id = substr($chr, 3, 1);
    $chr = "chr" . $id;
  }
  elsif ( $string_length == 5)
  {
    my $id = substr($chr, 3, 2);
    $chr = "chr" . $id;
  }
  return $chr;
}

sub ParseCounts
{
  my ($hash_label_ref, $AD_flag, $RO_flag, $AO_flag) = @_;

  my $ref_count = -1;
  my $mut_count = -1;
  my @result = ();

  my $AD_incomplete = 0;
  my $RO_AO_incomplete = 0;

  if ($AD_flag)
  {
    my $ad_value = $hash_label_ref->{ 'AD'};
    if ( IsComma( $ad_value))
    {
      my @v_array = split(",", $ad_value);
      $mut_count = $v_array[1];
      if( looks_like_number($v_array[0]))
      {
        $ref_count = $v_array[0];
      }
      else
      {
        $AD_incomplete = 1;
      }
      if ( looks_like_number($v_array[1]))
      {
        $mut_count = $v_array[1];
      }
      else 
      {
        $AD_incomplete = 1;
      }
    }
    else
    {
      $AD_incomplete = 1;
    }
  }

  if ($AD_incomplete)
  {
    $ref_count = -1;
    $mut_count = -1;
  }
  
  if ( ($AD_incomplete && $RO_flag && $AO_flag) || (!$AD_flag && $RO_flag && $AO_flag))
  {
    if ( looks_like_number($hash_label_ref->{ 'RO'}))
    {
      $ref_count = $hash_label_ref->{ 'RO'};
    }
    else
    {
      $RO_AO_incomplete = 1;
    }
    if ( looks_like_number($hash_label_ref->{ 'AO'}))
    {
      $mut_count = $hash_label_ref->{ 'AO'};
    }
    else
    {
      $RO_AO_incomplete = 1;
    }
  }
  if ( $RO_AO_incomplete)
  {
    $ref_count = -1;
    $mut_count = -1;
  }
  
  push (@result, $ref_count);
  push (@result, $mut_count);
  return @result;
}

sub IsComma
{
  my $string = $_[0];
  my $length = length($string);
  for ( my $i = 0; $i < $length; $i++)
  {
    if (substr($string, $i, 1) eq ",")
    {
      return 1;
    }
  }
  return 0;
}

sub median
{
  my @vals = sort {$a <=> $b} @_;
  my $len = @vals;
  if($len%2) # then odd
  {
    return $vals[int($len/2)];
  }
  else
  {
    return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
  }
}
