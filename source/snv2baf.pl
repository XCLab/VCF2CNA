#!/usr/bin/perl

# This program will obtain parse the high20 file to extract BAF data
# This is a refactor of the java program

use strict;
use warnings;

my $USAGE = "[SNVCOUNTS_OUTPUT FILE]\n";

unless ( $ARGV[0])
{
  print $USAGE;
  exit;
}

my $snvcounts_file = $ARGV[0];

my $delta = 0.05;

print "Chromosome\tPosition\tBAF_G\tBAF_D\n";
open( GET_FILE_DATA, $snvcounts_file);
my $header = <GET_FILE_DATA>;

while( my $entry = <GET_FILE_DATA>)
{
  chomp $entry;
  my @Token = split("\t", $entry);
  my $chr              = $Token[0];
  my $pos              = $Token[1];
  my $tumor_mutant     = $Token[2];
  my $tumor_total      = $Token[3];
  my $normal_mutant    = $Token[4];
  my $normal_total     = $Token[5];

  my ($chrom, $is_valid) = chr_convert($chr);
  unless( $is_valid)
  {
    next;
  } 
  my $germ_a_total = $normal_total - $normal_mutant;
  my $germ_b_total = $normal_mutant;
  my $tumor_a_total = $tumor_total - $tumor_mutant;
  my $tumor_b_total = $tumor_mutant;
  
  my $ratioG = 0.0;
  my $ratioD = 0.0;
  
  if ( $germ_b_total > 0)
  {
    $ratioG = $germ_b_total / ( $germ_a_total + $germ_b_total);
  }
  
  if ( $tumor_b_total > 0)
  {
    $ratioD = $tumor_b_total / ( $tumor_a_total + $tumor_b_total);
  }

  if ( $germ_a_total < 10)
  {
    next;
  }

  if ( $ratioG >= ( 0.5 - $delta) && $ratioG <= ( 0.5 + $delta))
  {
    if ( $tumor_a_total + $tumor_b_total < 10)
    {
      next;
    }
    print "$chrom\t$pos\t$ratioG\t$ratioD\n"; 
  }
}
close( GET_FILE_DATA);

exit;

sub chr_convert
{
  my $chrom = $_[0];
  $chrom =~ s/chr//g;
  if ( $chrom eq "X")
  {
    $chrom = 23;
  }
  elsif( $chrom eq "Y")
  {
    $chrom = 24;
  }

  my $is_valid = 0;
  for ( my $i = 1; $i <=24; $i++)
  {
    if ( $chrom eq $i)
    {
      $is_valid = 1;
      return( $chrom, $is_valid);
    }
  }
  return ($chrom, $is_valid);
}

sub read_file
{
  my ( $file, $header_flag) = @_;
  my @read_file = ();
  my $header;
  open( GET_FILE_DATA, $file);
  if( $header_flag)
  {
    $header = <GET_FILE_DATA>;
  }
  while( my $line = <GET_FILE_DATA>)
  {
    chomp $line;
    push( @read_file, $line);
  }
  close( GET_FILE_DATA);
  return @read_file;
}
