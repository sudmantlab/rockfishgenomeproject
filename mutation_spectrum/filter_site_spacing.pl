#!/bin/perl

use strict;
use warnings;

#This filters a het sites file to remove sites within 1 of another het.
my %line_storage;
my %pos_storage;
my $counter = 0;
while(<STDIN>){
  chomp;
  $line_storage{$counter}= $_;
  my @a = split(/\t/,$_);
  my $pos = $a[3];
  $pos_storage{$counter}= $pos;
  $counter++;
}
my $print_counter = 0;
my $skipped_sites;
foreach my $i (0..($counter-1)){
  #Check neighbours
  if ($pos_storage{$i - 1}){
    if (abs($pos_storage{$i - 1} - $pos_storage{$i}) == 1){
      $skipped_sites++;
      next;
    }
  }if ($pos_storage{$i + 1}){
    if (abs($pos_storage{$i + 1} - $pos_storage{$i}) == 1){
      $skipped_sites++;
      next;
    }
  }
  if($print_counter == 0){
    print "$line_storage{$i}";
    $print_counter++;
  }else{
    print "\n$line_storage{$i}";
  }
  
}
print STDERR "Skipped $skipped_sites sites due to spacing\n";
