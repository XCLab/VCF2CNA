#!/usr/bin/perl

# Script to determine filetype from input data

use strict;

my $USAGE = "[INPUT file]\n";

# Ensure file passed through commandline 
unless( $ARGV[0])
{
  print $USAGE;
  exit;
}

my $inputfile = $ARGV[0];
my $header = ReadHeader($inputfile);
my $type = ClassifyData($header);

print "$type\n";

exit;

sub ClassifyData
{
  my $header = $_[0];
  my $type = "unknown";
  
  my @line_array = split("\t", $header);
  my $line_elements = scalar(@line_array);
  if( substr($line_array[0], 0, 16) eq '##fileformat=VCF' && $line_elements == 1)
  {
    $type = "VCF";
    return $type;
  }

  my $high20_chrCol  = -1;
  my $high20_posCol  = -1;
  my $high20_typeCol = -1;
  my $high20_refCol  = -1;
  my $high20_altCol  = -1;
  my $high20_refCountCol = -1;
  my $high20_altCountCol = -1;

  my $maf_chrCol = -1;
  my $maf_posCol = -1;
  my $maf_typeCol = -1;
  my $maf_refCol = -1;
  my $maf_altCol = -1;
  my $maf_refCountCol = -1;
  my $maf_altCountCol = -1;

  for ( my $i = 0; $i < $line_elements; $i++)
  {
    if ($line_array[$i] eq 'Chr')
    {
      $high20_chrCol = 1;
    }
    elsif( $line_array[$i] eq 'Chromosome')
    {
      $maf_chrCol = 1;
    }
    elsif( $line_array[$i] eq 'Pos')
    {
      $high20_posCol = 1;
    }
    elsif( $line_array[$i] eq 'Start_Position' or $line_array[$i] eq 'Start_position')
    {
      $maf_posCol = 1;
    }
    elsif( $line_array[$i] eq 'Type')
    {
      $high20_typeCol = 1;
    }
    elsif( $line_array[$i] eq 'Variant_Type' or $line_array[$i] eq 'VariantType')
    {
      $maf_typeCol = 1;
    }
    elsif( $line_array[$i] eq 'Chr_Allele')
    {
      $high20_refCol = 1;
    }
    elsif( $line_array[$i] eq 'Tumor_ReadCount_Alt')
    {
      $maf_refCol = 1;
    }
    elsif( $line_array[$i] eq 'Alternative_Allele')
    {
      $high20_altCol = 1;
    }
    elsif( $line_array[$i] eq 'Tumor_ReadCount_Total')
    {
      $maf_altCol = 1;
    }
    elsif( $line_array[$i] eq 'reference_normal_count')
    {
      $high20_refCountCol= 1;
    }
    elsif( $line_array[$i] eq 'Normal_ReadCount_Alt')
    {
      $maf_refCountCol= 1;
    }
    elsif( $line_array[$i] eq 'alternative_normal_count')
    {
      $high20_altCountCol= 1;
    }
    elsif( $line_array[$i] eq 'Normal_ReadCount_Total')
    {
      $maf_altCountCol= 1;
    }
  }

  my $mafsum = $maf_chrCol + $maf_posCol + $maf_typeCol + $maf_refCol + $maf_altCol + $maf_refCountCol + $maf_altCountCol;
  my $h20sum = $high20_chrCol + $high20_posCol + $high20_typeCol + $high20_refCol + $high20_altCol + $high20_refCountCol + $high20_altCountCol;

  if ( $mafsum == 7)
  {
    $type = "MAF";
    return $type;
  }

  if ( $h20sum == 7)
  {
    $type = "HIGH20";
    return $type;
  }
  
  if ( $line_array[0] eq 'Chr'  && $line_array[1] eq 'Pos'  && $line_array[2] eq 'MinD' &&
       $line_array[3] eq 'TinD' && $line_array[4] eq 'MinN' && $line_array[5] eq 'TinN')
  {
    $type = "SJ_MAF";
    return $type;   
  }

  return $type;
}

sub ReadHeader
{
  my $file = $_[0];
 
  open( GET_FILE_DATA, $file);
  my $line = <GET_FILE_DATA>;
  chomp $line;
  close(GET_FILE_DATA);
  
  return $line;
}
