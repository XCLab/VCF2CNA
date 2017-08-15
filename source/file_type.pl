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
  if ( $line_array[0] eq 'TARGET_CASE_ID' && $line_elements == 40)
  {
    $type = "MAF";
    return $type;
  }
  if ( $line_array[0] eq 'NormalSample' && $line_elements == 45)
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
