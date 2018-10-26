#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Std;


my $usage =
"$0 file_fasta file_blast_parsed precursor_length

This script excised potential microRNA precursor sequences
from a genome using the positions of aligned reads as guidelines.
The fasta file given as input should be the genome in question,
and the file in blastparsed format should contain the alignments.
the precursor_length is the length of excised precursors. 
";


my $file_fasta=shift or die $usage;
my $file_blast_parsed=shift or die $usage;
my $precursor_length=shift or die $usage;

my %hash_fasta;
my %hash_align;

parse_file_fasta(\$file_fasta,\%hash_fasta);
parse_file_blast_parsed($file_blast_parsed);
excise();

exit;




sub parse_file_blast_parsed{

    my($file)=@_;

    open(FILENAME, $file) or die "Could not open file $file";
    
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
	    
	    my $strand=find_strand($other);
	    
	    unless($hash_fasta{$subject}){
		print STDERR "$subject not present in the genome fasta file\n";
		next;
	    }

	    insertfeature($query,$query_beg,$query_end,$query_lng,$subject,$subject_beg,$subject_end,$subject_lng,$strand);
	}
    }
    close FILENAME;
}




sub insertfeature{

    my($query,$query_beg,$query_end,$query_lng,$subject,$subject_beg,$subject_end,$subject_lng,$strand)=@_;
 
    my @prev_begs=sort {$a<=>$b} keys %{$hash_align{$subject}{$strand}};

    #testing for each previous query if the alignment is close enough to the present to fuse the two
    foreach my $prev_beg(@prev_begs){

	my $distance=$subject_beg-$prev_beg;

	if($distance>1000){next;}
	if($distance<-1000){last;}
	    
	my @prev_queries=sort keys %{$hash_align{$subject}{$strand}{$prev_beg}};
	    
	foreach my $prev_query(@prev_queries){

	    my $prev_end=$hash_align{$subject}{$strand}{$prev_beg}{$prev_query}{"subject_end"};
	    my $flank_beg=max2(1,$subject_beg-30);
	    my $flank_end=min2($subject_lng,$subject_end+30);
		
	    if(overlapping($flank_beg,$flank_end,$prev_beg,$prev_end)){
		my($new_beg,$new_end)=find_overlap($subject_beg,$subject_end,$prev_beg,$prev_end);
		$subject_beg=$new_beg;
		$subject_end=$new_end;
		    
		delete_feature($subject,$strand,$prev_beg,$prev_query);
	    }
	}
    }	

    #read new alignment into hash structure
    $hash_align{$subject}{$strand}{$subject_beg}{$query}{"subject_end"}=$subject_end;
    $hash_align{$subject}{$strand}{$subject_beg}{$query}{"subject_lng"}=$subject_lng;
    $hash_align{$subject}{$strand}{$subject_beg}{$query}{"query_beg"}=$query_beg;
    $hash_align{$subject}{$strand}{$subject_beg}{$query}{"query_end"}=$query_end;
    $hash_align{$subject}{$strand}{$subject_beg}{$query}{"query_lng"}=$query_lng;
}




sub parse_file_fasta{
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
		    $$hash{$id}  = $sequence;
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



sub excise{

    my $short_length = $precursor_length - 23;
    my @subjects=sort keys %hash_align;
    foreach my $subject(@subjects){

	my $count=0;
		
	my @strands=sort keys %{$hash_align{$subject}};
	foreach my $strand(@strands){

	    #use the genome contig or the reverse complement?
	    my $seq=();
	    if($strand eq "+"){
		$seq=$hash_fasta{$subject};
	    }else{
		$seq=revcom($hash_fasta{$subject});
	    }
	    my $seq_lng=length($seq);
    
	    my @subject_begs=sort {$a<=>$b} keys %{$hash_align{$subject}{$strand}};
	    foreach my $subject_beg(@subject_begs){
		
		my @queries=sort keys %{$hash_align{$subject}{$strand}{$subject_beg}};
		foreach my $query(@queries){
		    
		    my $subject_end=$hash_align{$subject}{$strand}{$subject_beg}{$query}{"subject_end"};

		    #longer alignments should be excised as a single potential precursor, shorter as two
		    if(($subject_end-$subject_beg)>30){
			excise_position(\$subject,\$seq,\$seq_lng,\$strand,\($subject_beg-22),\($subject_end+22),\$count);
			$count++;
		    }else{	
			excise_position(\$subject,\$seq,\$seq_lng,\$strand,\($subject_beg-22),\($subject_beg+$short_length),\$count);
			$count++;
			excise_position(\$subject,\$seq,\$seq_lng,\$strand,\($subject_end-$short_length),\($subject_end+22),\$count);
			$count++;
		    }
		}
	    }
	}
    }
    return;
}




sub excise_position{

    my($subject,$seq,$seq_lng,$strand,$excise_beg_old,$excise_end_old,$count)=@_;

    my $length_account=$precursor_length + 30;
    my $excise_beg=max2(1,$$excise_beg_old);
    my $excise_end=min2($$seq_lng,$$excise_end_old);

    my $excise_lng=$excise_end-$excise_beg+1;
    if($length_account<$excise_lng){return;}
    my $seq_sub=substr($$seq,$excise_beg-1,$excise_lng);

    print ">$$subject\_$$count strand:$$strand excise_beg:$excise_beg excise_end:$excise_end\n$seq_sub\n";
    return;
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

    my($sequence)=@_;

    my $rev=reverse $sequence;   

    return $rev;
}

sub com{

    my($sequence)=@_;

    $sequence=~tr/acgtuACGTU/TGCAATGCAA/;   
 
    return $sequence;
}

sub revcom{

    my($sequence)=@_;

    my $revcom=rev(com($sequence));

    return $revcom;
}




sub delete_feature{

    #Deletes a feature in the hash structure, removing stubs afterwards.
    my($subject,$strand,$subject_beg,$query)=@_;

    test_feature($subject,$strand,$subject_beg,$query);
 
    my $begins=scalar(keys %{$hash_align{$subject}{$strand}{$subject_beg}});

    #If this feature is the last of the query, remove the stub.
    if($begins==1){
	delete($hash_align{$subject}{$strand}{$subject_beg});
    }else{
	delete($hash_align{$subject}{$strand}{$subject_beg}{$query});
    }
}


sub test_feature{

    my($subject,$strand,$subject_beg,$query)=@_;

    unless($hash_align{$subject}{$strand}{$subject_beg}{$query}){
	die "The feature does not exist\n";
    }
}


sub find_overlap{

    #What is the overlap between the two stretches?
    my($begin1,$end1,$begin2,$end2)=@_;

    test_begin_end($begin1,$end1,$begin2,$end2);

    my $begin=min2($begin1,$begin2);
    my $end=max2($end1,$end2);

    return($begin,$end);
}



sub overlapping{

    #Do the two stretches overlap or is one contained in the other?
    my($begin1,$end1,$begin2,$end2)=@_;

    test_begin_end($begin1,$end1,$begin2,$end2);
    
    if(($begin1<=$begin2 and $begin2<=$end1+1) or ($begin1<=$end2+1 and $end2<=$end1) or contained($begin1,$end1,$begin2,$end2) or contained($begin2,$end2,$begin1,$end1)){
	return 1;
    }else{
	return 0;
    }
}



sub contained{

    #Is the stretch defined by the first positions contained in the stretch defined by the second?
    my($begin1,$end1,$begin2,$end2)=@_;

    test_begin_end($begin1,$end1,$begin2,$end2);

    if($begin2<=$begin1 and $end1<=$end2){
	return 1;
    }else{
	return 0;
    }
}


sub test_begin_end{

    #Are the beginposition numerically smaller than the endposition for each pair?
    my($begin1,$end1,$begin2,$end2)=@_;

    unless($begin1<$end1 and $begin2<$end2){
	die "Begin positions must be numerically smaller than endpositions\n";
    }
}


sub max2 {
        my($a, $b) = @_;
        return ($a>$b ? $a : $b);
}

sub min2  {
        my($a, $b) = @_;
        return ($a<$b ? $a : $b);
}
