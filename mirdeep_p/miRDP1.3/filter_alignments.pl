#!/usr/bin/perl

use warnings;
use strict;

use Getopt::Std;

my $usage =
"$0 blastparsefile

Parses a blastparsefile. As default, it prints out the lines where the query aligns perfectly
to the subject. With the -a option, it allows non-aligning nucleotides in the 3'end of the
query. With the -b option, it prints out the reads that align with the criteria given (in fasta
format). The reads are truncated with the number of nts needed to give a perfect alignment.

-a integer     Number of nucls in the 3'ends that are allowed mismatches
-b file_fasta  Print out the reads in fasta format
-c integer     Trim any sequences that are mapping more times than designated
";

#input
my $file_blast_parsed=shift or die $usage;

#options
my %options=();
getopts("a:b:c:",\%options);

#global variables
my %hash_fasta;
my %hash_align;
my %hash_align_best;


#parse blast_parsed input file
parse_file_blast_parsed($file_blast_parsed);

#print out fasta
if($options{b}){
    parse_file_fasta(\$options{b},\%hash_fasta);
    output_fasta();
}else{
#print out blast_parsed
    output_blast_parsed();
}



sub parse_file_blast_parsed{
    
    my $file_blast_parsed=shift;
    
    open (FILE_BLAST_PARSED, "<$file_blast_parsed") or die "can not open $file_blast_parsed\n";
    while (my $line=<FILE_BLAST_PARSED>){
	
	if($line=~/^(\S+)\s+(\S+)\s+(\d+)\.+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)\.+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(.+)$/){
	    
	    my $query=$1;
	    my $query_lng=$2;
	    my $query_beg=$3;
	    my $query_end=$4;
	    my $subject=$5;
	    my $subject_lng=$6;
	    my $subject_beg=$7;
	    my $subject_end=$8;
	    my $e_value=$9;
	    my $id=$10;
	    my $bitscore=$11;
	    my $other=$12;

	    #number of query nucleotides not covered in the alignment in the 5' and 3' end
	    my $trunc_beg=$query_beg-1;
	    my $trunc_end=$query_lng-$query_end;	

	    #number of uncovered nucleotides allowed in the 3' end
	    my $trunc_end_max=0;
	    if($options{a}){$trunc_end_max=$options{a};}

	    #test identity percentage
	    test_id($id);
	    
	    #full identity, all 5' nucleotides covered and limited number of 3' nucleotides uncovered 
	    if($id==1 and $trunc_beg==0 and $trunc_end<=$trunc_end_max){
	
		#read line into hash
		$hash_align{$query}{$trunc_end}.=$line;

		#register the least number of uncovered 3' nucleotides for the given query
		unless(defined($hash_align_best{$query}) and $hash_align_best{$query}<=$trunc_end){
		    $hash_align_best{$query}=$trunc_end;
		}
	    }
	}
    }
}


sub output_blast_parsed{

    #prints out the blast_parsed lines that pass the criteria

    my @queries=sort keys %hash_align_best;

    foreach my $query(@queries){
	if($options{c}){unless(test_query($query)){next;}}
	my $align_best=$hash_align_best{$query};
	my $lines_align_best=$hash_align{$query}{$align_best};
	print $lines_align_best;
    }
}



sub output_fasta{
  
    #prints out the queries that pass the criteria, in fasta format
    #At the end of each identifier, _t\d+ is added, with the integer
    #indicating how many nucleotides was omitted
  
    my @queries=sort keys %hash_align_best;
    
    foreach my $query(@queries){
	if($options{c}){unless(test_query($query)){next;}}
	my $align_best=$hash_align_best{$query};
	my $seq=$hash_fasta{$query};
	my $lng=length $seq;
	my $seq_align_best=substr($seq,0,$lng-$align_best);
	print ">$query";
	if($align_best){print "_t$align_best";}
	print "\n$seq_align_best\n";
    }
}


sub test_query{

    #makes sure that the query do not map too many times

    my $query=shift;

    my $align_best=$hash_align_best{$query};
    my $lines_align_best=$hash_align{$query}{$align_best};
    my $alignments=0;
    while($lines_align_best=~/\n/g){
	$alignments++;
    }
    if($alignments>$options{c}){
	return 0;
    }

    return 1;
}


sub parse_file_fasta{

    #read fasta file into a hash with the identifier as key and sequence as value

    my ($file,$hash) = @_;
    my ($id, $desc, $sequence) = ();

    open (FASTA, "<$$file") or die "can not open $$file\n";
    while (<FASTA>)
    {
        chomp;
        if (/^>(\S+)(.*)/)
	{
	    $id       = $1;
	    $desc     = $2;
	    $sequence = "";
	    while (<FASTA>){
                chomp;
                if (/^>(\S+)(.*)/){
		    $$hash{$id} = $sequence;
		    $id         = $1;
		    $desc       = $2;
		    $sequence   = "";
		    next;
                }
		$sequence .= $_;
            }
        }
	
    }
    $$hash{$id} = $sequence;
    close FASTA;
}



sub test_id{
  
    #test if the id percentage is in the correct format
  
    my $id=shift;

    unless($id=~/^(\d|\.)+$/){die "problem with id percentage format\n";}

    return;
}
