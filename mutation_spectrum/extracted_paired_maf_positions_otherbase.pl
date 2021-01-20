#!/bin/perl
use strict;
use warnings;

#pull out matching sites from two sample MAF. The first genome is aleutianus. In this case, we're using an aleutianus reference and extracting the base call from the other species. 
my $all_het_sites = $ARGV[0];
my $chosen_chr = $ARGV[1];
open(HET, "gunzip -c $all_het_sites |");

my %sites;
while(<HET>){
  chomp;
  if($_ =~ m/^$/){next;}
  my @a = split(/\t/,$_);
  my $chr = $a[0];
  my $pos = $a[1];
  if ($chr ne $chosen_chr){next;}
  $sites{$chr}{$pos}++;
}
close HET;
my %f;
$f{"A"} = "T";
$f{"T"} = "A";
$f{"C"} = "G";
$f{"G"} = "C";
$f{"-"} = "-";
$f{"N"} = "N";


my %stored_base;
my %stored_contig;
my %stored_position;
my %stored_direction;
my $sample_tracker = 1;

while(<STDIN>){
  chomp;
  if (($_ =~ m/^$/) or ($_ =~ m/^#/)){
    next;
  }
  if ($_ =~ m/^a/){
    foreach my $i (sort  {$a <=> $b} keys %stored_base){
      #Only look at second row for match
        if ($sites{$stored_contig{$i}{1}}{$stored_position{$i}{1}}){
          print "\n$stored_contig{$i}{1}\t$stored_position{$i}{1}\t";
          print "$stored_base{$i}{2}";
        }
    }
    #Reset stuff
    $sample_tracker = 1;
    undef(%stored_position);
    undef(%stored_base);
    undef(%stored_contig);
    undef(%stored_direction);
    next;
  }
  if ($_ =~ m/^s/){
    my @a = split(' ',$_);
    my $contig = $a[1];
    my $start = $a[2];
    my $length = $a[3];
    my $direction = $a[4];
    my $end = $a[5];
    my @seq = split(//,$a[6]);
    $stored_direction{$sample_tracker} = $direction;
    if ($direction eq "-"){    
      my $newstart = $end - $start;
      $start = $newstart;
    }
    my $current_pos = $start;
    foreach my $i (0..$#seq){
      if ($seq[$i] ne "-"){
        $current_pos++;
        $seq[$i] = uc($seq[$i]);
      }
      if ($direction eq "-"){
        $seq[$i] = $f{$seq[$i]};
      }
      $stored_position{$i}{$sample_tracker} = $current_pos;
      $stored_base{$i}{$sample_tracker} = $seq[$i];
      $stored_contig{$i}{$sample_tracker} = $contig;
    }
    $sample_tracker++;
  }
  
}


