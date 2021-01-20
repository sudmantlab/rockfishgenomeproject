#!/bin/perl

use strict;
use warnings;


my $min_exclude_length = $ARGV[0]; #Minimum length of excluded regions.

my $min_mappability = 1; #Minimum mappability to not put in excluded regions.

my $current_chr = "NA";
my $previous_start;
my $last_end;
while(<STDIN>){
  chomp;
  my @a = split(/\t/,$_);
  my $chr = $a[0];
  my $start = $a[1];
  my $end = $a[2];
  my $mappability = $a[3];
  if ($chr ne $current_chr){
    #If we've progressed into another chromosome.
    if ($last_end){
      my $length = $last_end - $previous_start;
      if ($length >= $min_exclude_length){
        print "$current_chr\t$previous_start\t$last_end\n";
      }
    }
    $current_chr = $chr;
    undef($previous_start);
    undef($last_end);
  }
  if ($mappability >= $min_mappability){
    #its mappable, so print out other windows or 
    if ($last_end){
      my $length = $last_end - $previous_start;
      if ($length >= $min_exclude_length){
        print "$current_chr\t$previous_start\t$last_end\n";
      }
    }
    undef($previous_start);
    undef($last_end);
  }else{
    #If its not mappable, record the start and end;
    unless(defined($previous_start)){
      $previous_start = $start;
    }
    $last_end =$end;  
  }
  
}

#Print out last window;
my $length = $last_end - $previous_start;
if ($length >= $min_exclude_length){
  print "$current_chr\t$previous_start\t$last_end";

}
