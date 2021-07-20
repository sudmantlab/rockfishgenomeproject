#!/usr/bin/env perl

my $num1 = $ARGV[0];
my $group = $ARGV[1];

# Go for a chi-square test and do multiple testing for the Branch model results
use Statistics::Distributions qw(chisqrprob);
use Statistics::Multtest qw(BY qvalue);
use Statistics::R;
 
#my $R = Statistics::R->new();


my $p_value = Statistics::Distributions::chisqrprob(1, $num1);
#my $q_value = qvalue($p_value);

print $group,"\t",$p_value,"\n";

# print $p_value,"\t",$q_value,"\n";


