#!/bin/perl
use strict;
use warnings;


my $sample = $ARGV[1];
my $species = $ARGV[0];

my $het_file = "/global/scratch2/rohitkolora/Rockfish/Data/Greg/snp_calling/aleutianus_scaffolded/$species/Illumina/map.$species.$sample.hetsites.filtered.txt";
my $site_match_file = "/global/scratch2/rohitkolora/Rockfish/Data/Greg/snp_calling/all_aleutianus_hetsites/$species/$species.$sample.hetsites.matched.txt";
my $species_list = "/global/scratch2/rohitkolora/Rockfish/Data/Greg/snp_calling/specieslist.txt";
#Load in a list of species
my @species;
open SPECIES, $species_list;
while(<SPECIES>){
  chomp;
  push(@species,$_);
}
close SPECIES;
print "# $sample $species\n";
print "#  Specieslist: ";
foreach my $chosen_species (@species){
  if ($species ne $chosen_species){
    print "$chosen_species ";
  }
}

my %context;
my %ref;
my %alt;
#Load in the alleles and surrounding bases
open HET, $het_file;
my %het_sites;
while(<HET>){
  chomp;
  if ($_ =~ m/^$/){next;}
  my @a = split(/\t/,$_);
  my $chr = "$species.$a[2]";
  my $pos = $a[3];
  my $ref = $a[4];
  my $alt = $a[5];
  my $context = uc($a[6]);
  $ref{$chr}{$pos} = $ref;
  $alt{$chr}{$pos} = $alt;
  $context{$chr}{$pos} = $context;
}
close HET;
print STDERR "Loaded all heterozygous sites\n";
#Find the aleutianus reference position
my %chr_match;
my %pos_match;
open REF, $site_match_file;
while(<REF>){
  chomp;
  if ($_ =~ m/^$/){next;}
  my @a = split(/\t/,$_);
  my $al_chr = $a[0];
  my $al_pos = $a[1];
  my $chr = $a[4];
  my $pos = $a[5];
  $chr_match{$al_chr}{$al_pos} = $chr;
  $pos_match{$al_chr}{$al_pos} = $pos;
}
print STDERR "Loaded all matching sites\n";
my %other;
foreach my $chosen_species (@species){
  if ($species eq $chosen_species){next;}
  my $other_species_file = "/global/scratch2/rohitkolora/Rockfish/Data/Greg/snp_calling/all_aleutianus_hetsites/$chosen_species/$chosen_species.fullset.txt";
  open OTHER_SPECIES, $other_species_file;
  my %loaded_bases;
  while(<OTHER_SPECIES>){
    chomp;
    if ($_ =~ m/^$/){next;}
    my @a = split(/\t/,$_);
    my $chr = $a[0];
    my $pos = $a[1];
    my $base = $a[2];
    if($chr_match{$chr}{$pos}){
      my $home_chr = $chr_match{$chr}{$pos};
      my $home_pos = $pos_match{$chr}{$pos};
      if($loaded_bases{$home_chr}{$home_pos}){
        $loaded_bases{$home_chr}{$home_pos} = "N";
      }else{
        $loaded_bases{$home_chr}{$home_pos} = $base;
      }
    }
  }
  close OTHER_SPECIES;
  foreach my $chr (sort keys %ref){
    foreach my $pos (sort keys %{$ref{$chr}}){
      if ($loaded_bases{$chr}{$pos}){
        $other{$chr}{$pos} .= $loaded_bases{$chr}{$pos};
      }else{
        $other{$chr}{$pos} .= "-";
      }
      
    }
  }
  print STDERR "Loaded all matching bases for $chosen_species\n";
}
foreach my $chr (sort keys %ref){
  foreach my $pos (sort {$a <=> $b} keys %{$ref{$chr}}){
    print "\n$chr\t$pos\t$ref{$chr}{$pos}\t";
    print "$alt{$chr}{$pos}\t$context{$chr}{$pos}\t";
    print "$other{$chr}{$pos}";
  }
}



