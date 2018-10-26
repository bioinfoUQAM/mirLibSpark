#!/usr/bin/perl

use warnings;
use strict;

use Getopt::Std;


my $usage =
"$0 blastparsefile (options)

With no options, this script simply parses a blastparsefile and outputs each line.
The options allow output only of lines with certain characteristics.
Matchfiles (which can be given as additional arguments) should contain
only ids separated by single newlines.

-a      Only output id of the queries fulfilling the criteria (non-redundant)
-b      Only output id of the subjects fulfilling the criteria (non-redundant)
-c real Only considers lines that have equal or better identity than designated by the real  
-d real Only considers lines that have equal or better bitscore than designated by the real
-e      Only considers lines where the alignment has full coverage and 100% identity
-f file Only considers lines where the query is contained in the matchfile
-g file Only considers lines where the query is not contained in the matchfile
-h file Only considers lines where the subject is contained in the matchfile
-i file Only considers lines where the subject is not contained in the matchfile
";

#imput
my $blastparsefile=shift or die $usage;

#options
my %options=();
getopts("abc:d:ef:g:h:i:",\%options);

#cut-offs
my $min_id=0;
my $min_bitscore=0;
my $min_query_cov=0;

#modify cut-offs according to options
if($options{c}){$min_id=$options{c};}
if($options{d}){$min_bitscore=$options{d};}
if($options{e}){$min_query_cov=1;$min_id=1;}

#global variables
my %queryhash;
my %subjecthash;
my %linehash;

#parse matchfiles
if($options{f}){%linehash=linestohash($options{f})};
if($options{g}){%linehash=linestohash($options{g})};
if($options{h}){%linehash=linestohash($options{h})};
if($options{i}){%linehash=linestohash($options{i})};


#parse blast_parse file
open(FILENAME, $blastparsefile) or die "Could not open $blastparsefile.\n";

while (my $line = <FILENAME>){
    if($line=~m/^(\S+)\s+(\S+)\s+(\d+)\.+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)\.+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/){
	
	my $query=$1;
	my $query_lng=$2;
	my $query_beg=$3;
	my $query_end=$4;
	my $subject=$5;
	my $subject_lng=$6;
	my $subject_beg=$7;
	my $subject_end=$8;
	my $evalue=$9;
	my $id=$10;
	my $bitscore=$11;
	my $strand=$12;

	#test format
	if($options{c} and !$id=~/^(\d|\.)+$/){print STDERR "problem with id percentage format\n"; exit;}
	if($options{d} and !$bitscore=~/^(\d|\.)+$/){print STDERR "problem with bitscore format\n"; exit;}

	#find query coverage
	$query_lng=~s/,//;
	my $query_cov=($query_end-$query_beg+1)/$query_lng;
	if($options{e} and $id eq 'n/a'){$id=1;}
	
	#test to see if relevant id is contained in matchfile 
	my $check_patterns="1";
	if($options{f} or $options{g}){$check_patterns=check_patterns($query);}
	if($options{h} or $options{i}){$check_patterns=check_patterns($subject);}

	#test to see if all requirements are fulfilled to print the line
	unless(($options{c} and $id<$min_id) or ($options{d} and $bitscore<$min_bitscore) or ($options{e} and ($query_cov<$min_query_cov or $id<$min_id)) or !$check_patterns){
	    unless($options{a} or $options{b}){
		print $line;
	    }else{
		$queryhash{$query}=1;
		$subjecthash{$subject}=1;
	    }
	}
    }
}

#if options -a or -b are used, print the relevant ids
my @queries=sort keys %queryhash;
my @subjects=sort keys %subjecthash;

if($options{a}){
    foreach my $query(@queries){
	print "$query\n";
    }
}

if($options{b}){
    foreach my $subject(@subjects){
	print "$subject\n";
    }
}



sub check_patterns{
    
    #see if the relevant id is contained in the matchfile
  
    my $id=shift;
        
    my $match=0;
    
    if($linehash{$id}){
	$match=1;
    }
    
    if(($match and ($options{f} or $options{h})) or (!$match and ($options{g} or $options{i}))){
	return 1;
    }
    
    return 0;
}




sub linestohash{

    #read matchfile into hash
    my($filename)=@_;
    
    my %linehash;
    
    open(FILENAME, $filename) or die "Could not open $filename.\n";
    
    while (my $line=<FILENAME>){
	chomp $line;
	$linehash{$line}=1;
    }
    
    return %linehash;
}
