#!/bin/perl

use strict;
use warnings;


#This creates a bed file to extract the bases surrounding for the het sites
my $prefix = "Sebastes_aleutianus";
my $ref_seq = "/global/scratch/gregoryowens/sebastes/snp_calling/FALCON/Sebastes_aleutianus/referencegenome.fasta";

open my $bed, '>',  "$prefix.bed";
print STDERR "Creating BED file of sites ...\n";
my $counter = 0;
my %lines;
while(<STDIN>){
  chomp;
  if ($_ =~ m/^$/){next;}
  my @a = split(/\t/,$_);
  my $chr = $a[0];
  my $pos = $a[1];
  my $start = $pos-1;  #+1 for bed
  my $end = $pos; #+2 for bed
  print $bed "$chr\t$start\t$end\n";
  $lines{$counter} = $_;
  $counter++;

}
close $bed;

if ($counter == 0) {
  die("Invalid input file given.");
}
unless ( -e "$prefix.fasta"){

  print STDERR "Using bedtools to create fasta of surrounding sequence..\n";
  system("bedtools getfasta -fi $ref_seq -bed $prefix.bed -fo $prefix.fasta");
}

print STDERR "bedtools getfasta -fi $ref_seq -bed $prefix.bed -fo $prefix.fasta";
open FASTA, "$prefix.fasta";
my %seqs;
my $second_counter = 0;
while(<FASTA>){
  chomp;
  if ($second_counter % 2 == 0){
    $second_counter++;
    next;
  }else{
    my $n = ($second_counter-1) / 2;
    $seqs{$n} = $_;
    $second_counter++;
  }
}
foreach my $i (0..($counter - 1)){
  if($i == 0){
    print "$lines{$i}\t$seqs{$i}";

  }else{
    print "\n$lines{$i}\t$seqs{$i}";
  }
}
