#!/bin/perl
use strict;
use warnings;

my %spectra;
while(<STDIN>){
  chomp;
  if ($_ =~ m/^#/){next;}
  my @a = split(/\t/,$_);
  my $ref = $a[2];
  my $alt = $a[3];
  my @context = split(//,$a[4]);
  my $other_species = $a[5];
  my $ref_count = 0;
  my $alt_count = 0;
  my $called_count = 0;
  my @b = split(//,$other_species);
  my $max_called = $#b + 1;
  foreach my $other (@b){
    if ($other eq $ref){
      $ref_count++;
    }elsif ($other eq $alt){
      $alt_count++;
    }
    if ($other eq "-"){next;}
    if ($other eq "N"){next;}
    $called_count++;
  }
  my $ancestral;
  my $derived;
  my $percent_called = $called_count/$max_called;
  if ($percent_called < 0.2){next;}
  if (($alt_count > 0) and ($ref_count == 0)){
    $ancestral = $alt;
    $derived = $ref;
  }elsif (($alt_count == 0) and ($ref_count > 0)){
    $ancestral = $ref;
    $derived = $alt;
  }else{
    $ancestral = "?";
    $derived = "?";
    next;
  }
  $spectra{"$context[0]$ancestral$context[2]:$context[0]$derived$context[2]"}++;
}
print "site_type\tn";
foreach my $x (sort keys %spectra){
  print "\n$x\t$spectra{$x}";
}
