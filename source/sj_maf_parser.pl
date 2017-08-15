#!/usr/bin/perl -w
# Script to compute median of NormalTotal values in snvcounts_outputfile

use File::Basename;
use Scalar::Util qw(looks_like_number);
use strict;

my $USAGE = "[snvcounts_outputfile]\n";

# Read file provided from the commandline
unless( $ARGV[0])
{
  print $USAGE;
  exit;
}

# Read Commandline
my $snv_counts  = $ARGV[0];

open( GET_FILE_DATA, $snv_counts);
my $header = <GET_FILE_DATA>;

my @median_values = ();

while( my $line = <GET_FILE_DATA>)
{
  chomp $line;
  my @line_array = split("\t", $line);
  my $chr = $line_array[0];
  my $pos = $line_array[1];
  my $tumor_mutant = $line_array[2];
  my $tumor_total  = $line_array[3];
  my $normal_mutant= $line_array[4];
  my $normal_total = $line_array[5];
  if (is_valid($chr))
  {
    push(@median_values, $normal_total);
  }
}
close( GET_FILE_DATA);

if ( scalar(@median_values) > 0)
{
  my $median = median(@median_values);
  print "$median\n";
}

exit;

sub median
{
  my @vals = sort {$a <=> $b} @_;
  my $len = @vals;
  if( $len%2) # then odd
  {
    return $vals[int($len/2)];
  }
  else
  {
    return ( $vals[int($len/2)-1] + $vals[int($len/2)] )/2;
  }
}

sub is_valid
{
  my $chr = $_[0];
  $chr =~ s/chr//g;
  if( $chr eq "X")
  {
    $chr = 23;
  }
  elsif( $chr eq "Y")
  {
    $chr = 24;
  }
  my $is_valid = 0;
  for ( my $i = 1; $i <= 24; $i++)
  {
    if ( $chr eq $i)
    {
      $is_valid = 1;
      return $is_valid;
    }
  }
  return $is_valid;
}
