#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Std;


my $usage =
"$0 file_1 file_2

Prints out the ids of the annotations in file_1 that overlap with
annotations in file_2. Disregards strands.

-a first file is in gtf rather than blastparsed format
-b second file is in gtf rather than blastparsed format
";

#input
my $file_1=shift or die $usage;
my $file_2=shift or die $usage;

#options
my %options=();
getopts("ab",\%options);

#global variables
my %hash_1;
my %hash_2;
my %hash_local;
my %hash_overlap;

my $count;

#parse file 1, either in gtf or blast_parsed format
if($options{a}){
    parse_file_gtf($file_1);
}else{
    parse_file_blast_parsed($file_1);
}

#parse file 2, either in gtf or blast_parsed format
if($options{b}){
    parse_file_gtf($file_2);
}else{
    parse_file_blast_parsed($file_2);
}

find_overlap();

#print ids of overlapping annotations
print_keys(\%hash_overlap);


exit;





sub parse_file_blast_parsed{

    my($file_blast_parsed)=@_;

    open(FILENAME, $file_blast_parsed) or die "Could not open file $file_blast_parsed";
    
    while (my $line = <FILENAME>){
	if($line=~m/^(\S+)\s+(\S+)\s+(\d+)\.+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)\.+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(.+)$/){

	    my $query=$1;
	    my $query_lng=$2;
	    my $query_beg=$3;
	    my $query_end=$4;
	    my $subject=$5;
	    my $subject_lng=$6;
	    my $subject_beg=$7;
	    my $subject_end=$8;
	    my $e_value=$9;
	    my $pid=$10;
	    my $bitscore=$11;
	    my $other=$12;
	    
	    #reverse positions so that nucleotide 1 on the minus strand and nucleotide 2 on the plus
	    #strand base pairs. This makes it much easier to compare annotations on opposite strands.
	    my $strand=find_strand($other);
	    if($strand eq "-"){
		($subject_beg,$subject_end)=reverse_positions($subject_lng,$subject_beg,$subject_end);	
	    }

	    #read query, subject, begin, end into hash
	    if($file_blast_parsed eq $file_1){
		$hash_1{$subject}{$query}{$subject_beg}=$subject_end;
	    }else{
		$hash_2{$subject}{$query}{$subject_beg}=$subject_end;
	    }
	    $count++;
	    print STDERR "$file_blast_parsed $count\r";
	}
    }
    $count=0;
    print STDERR "\n\r";
    close FILENAME;
}





sub parse_file_gtf{

    my($file_gtf)=@_;

    open(FILE_GTF, $file_gtf) or die "Could not open $file_gtf.\n";
   
    while (my $line = <FILE_GTF>){

	if($line=~/^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+).+(\S+)/){
	    my $subject=$1;
	    my $source=$2;
	    my $type=$3;
	    my $subject_beg=$4;
	    my $subject_end=$5;
	    my $score=$6;
	    my $strand=$7;
	    my $query=$8;

	    #read query, subject, begin, end into hash
	    if($file_gtf eq $file_1){
		$hash_1{$subject}{$query}{$subject_beg}=$subject_end;
	    }else{
		$hash_2{$subject}{$query}{$subject_beg}=$subject_end;
	    }
	    $count++;
	    print STDERR "$file_gtf $count\r";

	}
    }
    $count=0;
    print STDERR "\n\r";
    close FILE_GTF;
}



sub find_overlap{

    #each subject (database sequence) is treated separately

    my @subjects=sort keys %hash_2;
    foreach my $subject(@subjects){

	#if no queries from file 1 map to the subject, no overlap is possible
	if(not defined $hash_1{$subject}){next;}

	#fill nucleotide hash
	fill_hash_local($subject);

	#resolve overlap
	find_overlap_local($subject);

	#initialize nucleotide hash
	%hash_local=();
    }
}






sub fill_hash_local{


    #fill all nucleotide positions that annotations in file 2 map to into a
    #hash with the position as key

    my ($subject)=@_;

    my @queries=sort keys %{$hash_2{$subject}};
    foreach my $query(@queries){
	
	my @subject_begs=sort keys %{$hash_2{$subject}{$query}};
	foreach my $subject_beg(@subject_begs){

	    print STDERR "Filling hash on $subject $count\r";
	    $count++;

	    my $subject_end=$hash_2{$subject}{$query}{$subject_beg};
	    for(my $pos=$subject_beg;$pos<=$subject_end;$pos++){
		$hash_local{$pos}=1;
	    }
	}
    }
    $count=0;
    print STDERR "\n\r";
}







sub find_overlap_local{

    my ($subject)=@_;

    my @queries=sort keys %{$hash_1{$subject}};
    foreach my $query(@queries){

	#if the query from file 1 is already overlapping with annotation
	#from file 2, then disregard
	if($hash_overlap{$query}){next;}
	
	#for each nucleotide position that the query map to, investigate
	#if annotation from file 2 map to it. If so, put into %hash_overlap
	#and exit the for loop.
	my @subject_begs=sort keys %{$hash_1{$subject}{$query}};
        OUTER: foreach my $subject_beg(@subject_begs){
	    
	    print STDERR "Finding overlap on $subject $count\r";
	    $count++;

	    my $subject_end=$hash_1{$subject}{$query}{$subject_beg};
	    INNER: for(my $pos=$subject_beg;$pos<=$subject_end;$pos++){
		if($hash_local{$pos}){
		    $hash_overlap{$query}=1;
		    last OUTER;
		}
	    }
	}
    }
    $count=0;
    print STDERR "\n\r";
}





sub print_keys{

    #print hash

    my $hash=shift;

    my @keys=sort keys %$hash;
    foreach my $key(@keys){
	print "$key\n";
    }
}




sub find_strand{

    #A subroutine to find the strand, parsing different blast formats

    my($other)=@_;

    my $strand="+";

    if($other=~/-/){
	$strand="-";
    }

    if($other=~/minus/i){
	$strand="-";
    }

    return($strand);
}


sub reverse_positions{

    #A subroutine to find positions relative to the minus strand

    my($length,$begin,$end)=@_;

    my $new_end=$length-$begin+1;
    my $new_beg=$length-$end+1;

    return($new_beg,$new_end);
}



sub rev{

    #reverses the order of nucleotides in a sequence

    my($sequence)=@_;

    my $rev=reverse $sequence;   

    return $rev;
}

sub com{

    #the complementary of a sequence

    my($sequence)=@_;

    $sequence=~tr/acgtuACGTU/TGCAATGCAA/;   
 
    return $sequence;
}

sub revcom{

    #reverse complement

    my($sequence)=@_;

    my $revcom=rev(com($sequence));

    return $revcom;
}


sub max2 {
    
    #max of two numbers
    
    my($a, $b) = @_;
    return ($a>$b ? $a : $b);
}

sub min2  {
    
    #min of two numbers
    
    my($a, $b) = @_;
    return ($a<$b ? $a : $b);
}
