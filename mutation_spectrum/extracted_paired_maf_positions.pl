#!/bin/perl
use strict;
use warnings;

#pull out matching sites from two sample MAF. The first genome is aleutianus, which we're trying to find the corresponding site in, and the second is the ragtag genome.
my $het_sites = $ARGV[0];
my $species = $ARGV[1]; #Added to start of contig name because ragtag scaffolds are all called the same.

open(HET, "cat $het_sites |");


print STDERR "Reading target sites..\n"; 
my %sites;
while(<HET>){
  chomp;
  my @a = split(/\t/,$_);
  my $chr = "$species.$a[2]";
  my $pos = $a[3];
  $sites{$chr}{$pos}++;
}
close HET;
print STDERR "Done reading target sites..\n";

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

my $progress_counter=1;
while(<STDIN>){
  chomp;
  if (($_ =~ m/^$/) or ($_ =~ m/^#/)){
    next;
  }
  if ($_ =~ m/^a/){
    foreach my $i (sort  {$a <=> $b} keys %stored_base){
      #Only look at second row for match
        if ($sites{$stored_contig{$i}{2}}{$stored_position{$i}{2}}){
          if (($stored_base{$i}{1} eq "-") or ($stored_base{$i}{1} eq "N")){next;}
          if (($stored_base{$i}{2} eq "-") or ($stored_base{$i}{2} eq "N")){next;}
          print "\n$stored_contig{$i}{1}\t$stored_position{$i}{1}";
          print "\t$stored_base{$i}{1}\t$stored_direction{1}\t";
          print "$stored_contig{$i}{2}\t$stored_position{$i}{2}\t";
          print "$stored_base{$i}{2}\t$stored_direction{2}";
        }
    }
    #Reset stuff
    $sample_tracker = 1;
    undef(%stored_position);
    undef(%stored_base);
    undef(%stored_contig);
    undef(%stored_direction);
    $progress_counter++;
#    if ($progress_counter % 100 == 0){ print STDERR "Progress $progress_counter\n";}
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
print STDERR "Done outputting sites"
