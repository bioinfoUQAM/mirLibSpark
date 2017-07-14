#!/bin/perl
#Auteur : M.A. Remita
#Date : 11/04/2017

# The program run miRcheck
# counting position by 0 based for both start and stop

# pre-miRNA from miRBase:
# GUGGAGCUCCUAUCAUUCCAAUGAAGGGUCUACCGGAAGGGUUUGUGCAGCUGCUCGUUCAUGGUUCCCACUAUCCUAUCUCCAUAGAAAACGAGGAGAGAGGCCUGUGGUUUGCAUGACCGAGGAGCCGCUUCGAUCCCUCGCUGACCGCUGUUUGGAUUGAAGGGAGCUCUGCAU
# (((((((((((.(((.(((((...((((((...((..((((.(((.((.(((.((((.(((((....((((...(((.(((((...........))))).)))...))))....))))).)))).))).))..))).))))))..)))).))..))))).))).)))))))))))..
# start: 153	stop: 173
# line = perl eval_mircheck.pl "(((((((((((.(((.(((((...((((((...((..((((.(((.((.(((.((((.(((((....((((...(((.(((((...........))))).)))...))))....))))).)))).))).))..))).))))))..)))).))..))))).))).))))))))))).." 153 173 def
# output of mirCheck: 
# 3prime	1	173


use strict;
use warnings;
use lib '../lib/';
#use lib '/home/shared_dir/lib/';
use miRcheck;

my $usage = 'perl eval_mircheck.pl "((((((.((((((....).))))).)))))).........." 46 64 def\n';

die $usage."\n" if scalar @ARGV != 4;

my $fold= $ARGV[0];
my $miR_start= $ARGV[1];
my $miR_stop=$ARGV[2];
my $param=$ARGV[3];

my ($fback, $fback_start, $fback_stop);

($fback, $fback_start, $fback_stop) = miR_check($fold, $miR_start, $miR_stop, $param);

print "$fback\t$fback_start\t$fback_stop\n";


__END__
#$hairpin = "AAACUAAUGGACCACCUGUGUUUAGUAGUUUCAAUCAACUU";
#$fold = "((((((.((((((....).))))).))))))..........";
#$mirna = "UGUGUUUAGUAGUUUC";
#$miR_start = index ($hairpin,$mirna);				#46
#$miR_stop = $miR_start + length ($mirna) -1;			#64

