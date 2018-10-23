#!/usr/bin/perl -w

use strict;
use Getopt::Std;


######################################### USAGE ################################

my $usage =
"$0 annotated_miR precursor retrieved_miR

This script is used to compare annotated miRs and retrieved miRs achieved by miRDeep.
The file, annotated_miR, consists of annotated miRs info, including miR position at genome, strand direction, etc.
Its format likes the file ath_miR_info.
The file, retrieved_miR, is the result of retrieved or predicted miRs achieved by miRDeep.
The file, precursors, is the fasta file achieved by excise_candidate.pl, one script in miRDeep package.
 
";


####################################### INPUT FILES ############################


my $annotated_miR=shift or die $usage;
my $precursor=shift or die $usage;
my $retrieved_miR=shift or die $usage;

##################################### GLOBAL VARIBALES #########################


my %MIR_ID=();
my %MIR_chr=();
my %MIR_B=();
my %MIR_E=();
my %miR_chr=();
my %miR_dir=();
my %miR_B=();
my %miR_E=();


my %miRNA_strand=();
my %precursor_bg=();
my %precursor_end=();
my %precursor_chr=();


my %predicted_miRNA_chr=();
my %predicted_miRNA_bg=();
my %predicted_miRNA_end=();
my %predicted_star_bg=();
my %predicted_star_end=();



######################################### MAIN ################################# 

#extract known miRNA information
extract_miRNA_info($annotated_miR);


#parse excised precursor information
parse_precursor($precursor);


#extract predicted miRNA information
parse_miRDeep_result($retrieved_miR);


#print $predicted_miRNA_chr{"osa-MIR1317_2"}."\t".$predicted_miRNA_bg{"osa-MIR1317_2"}."\t".$predicted_miRNA_end{"osa-MIR1317_2"}."\t".$predicted_star_bg{"osa-MIR1317_2"}."\t".$predicted_star_end{"osa-MIR1317_2"}."\n";

#compare known and predicted miRNA
compare_position();


#
exit;


############################################ SUBROUTINES #######################


sub extract_miRNA_info{
	#obtain miRNA precursor and mature miRNA position at chromosome
	my $file=shift @_;
	open(FILE,"$file")||die"can not open the miRNA_info file:$!\n";
	

	while(<FILE>){
		chomp;
		my @tmp=split("\t",$_);
		
		my $precursor_length = $tmp[6] - $tmp[5] + 1;
		
		$MIR_chr{$tmp[1]} = $tmp[3];		
		$MIR_B{$tmp[1]} = $tmp[5];
		$MIR_E{$tmp[1]} = $tmp[6];
		$MIR_ID{$tmp[0]} = $tmp[1];
		$miR_chr{$tmp[0]} = $tmp[3];
		$miR_dir{$tmp[0]} = $tmp[4];
		
		
		if($tmp[4] eq "+"){
			$miR_B{$tmp[0]} = $tmp[5] + $tmp[7];
			$miR_E{$tmp[0]} = $tmp[5] + $tmp[8];
			
			
			}
		if($tmp[4] eq "-"){
			$miR_B{$tmp[0]} = $tmp[6] - $tmp[8];
			$miR_E{$tmp[0]} = $tmp[6] - $tmp[7];
			
			}

		}
		
		close FILE;
}



sub parse_precursor{
	#extract the precursor position information in the miRNA extension sequences
	my $file = shift @_;
	open(FILE,"$file")||die"can not open the precursor file:$!\n";


	while(<FILE>){
		chomp;
		if($_=~/^>/){
			my @tmp=split(" ",$_);
			my $tmp_precursor_id;
			
			if($tmp[0]=~/>(\S+)_[0-9]+/){
				$tmp_precursor_id=$1;
				
					my $tmp_miRNA_id = substr($tmp[0],1,);
										
					my @tmp1=split(":",$tmp[1]);
					my @tmp2=split(":",$tmp[2]);
					my @tmp3=split(":",$tmp[3]);
					
					$miRNA_strand{$tmp_miRNA_id}=$tmp1[1];
					$precursor_chr{$tmp_miRNA_id}=$MIR_chr{$tmp_precursor_id};
					
					if($tmp1[1] eq "+"){
						$precursor_bg{$tmp_miRNA_id} = $MIR_B{$tmp_precursor_id} - 500 + $tmp2[1];
						$precursor_end{$tmp_miRNA_id} = $MIR_B{$tmp_precursor_id} - 500 + $tmp3[1];
						
						}
					if($tmp1[1] eq "-"){
						$precursor_bg{$tmp_miRNA_id} = $MIR_E{$tmp_precursor_id} + 500 - $tmp3[1];
						$precursor_end{$tmp_miRNA_id} = $MIR_E{$tmp_precursor_id} + 500 - $tmp2[1];
						
						}

				}
		}

	}
	close FILE;
	
}


sub parse_miRDeep_result{
	#extract mature miRNA position predicted by miRDeep
	my $file=shift @_;
	

	my $predicted_miRNA_id=();
	my $predicted_miRNA_strand=();
	my $predicted_miRNA_bg=();
	my $predicted_miRNA_end=();
	my $predicted_star_bg=();
	my $predicted_star_end=();

	
	$/="\n\n\n";
	
	open(FILE,"$file")||die"can not open the miRDeep result file:$!\n";
	while(<FILE>){
		chomp;
		
		if($_=~/mature_beg\s+([0-9]+)/){
			$predicted_miRNA_bg=$1;
			}
		if($_=~/mature_end\s+([0-9]+)/){
			$predicted_miRNA_end=$1;
			}
		if($_=~/pri_id\s+(\S+)\n/){
			$predicted_miRNA_id=$1;
			}
		if($_=~/mature_strand\s+(\S)\n/){
			$predicted_miRNA_strand=$1;			
			}
		if($_=~/star_beg\s+([0-9]+)/){
			$predicted_star_bg=$1;
			}
		if($_=~/star_end\s+([0-9]+)/){
			$predicted_star_end=$1;
			}
		
		$predicted_miRNA_chr{$predicted_miRNA_id} = $precursor_chr{$predicted_miRNA_id};
		
		if($miRNA_strand{$predicted_miRNA_id} eq "+"){
			$predicted_miRNA_bg{$predicted_miRNA_id} = $precursor_bg{$predicted_miRNA_id} + $predicted_miRNA_bg;
			$predicted_miRNA_end{$predicted_miRNA_id} = $precursor_bg{$predicted_miRNA_id} + $predicted_miRNA_end;
			$predicted_star_bg{$predicted_miRNA_id} = $precursor_bg{$predicted_miRNA_id} + $predicted_star_bg;
			$predicted_star_end{$predicted_miRNA_id} = $precursor_bg{$predicted_miRNA_id} + $predicted_star_end;
			}
		if($miRNA_strand{$predicted_miRNA_id} eq "-"){
			$predicted_miRNA_bg{$predicted_miRNA_id} = $precursor_end{$predicted_miRNA_id} - $predicted_miRNA_end;
			$predicted_miRNA_end{$predicted_miRNA_id} = $precursor_end{$predicted_miRNA_id} - $predicted_miRNA_bg;
			$predicted_star_bg{$predicted_miRNA_id} = $precursor_end{$predicted_miRNA_id} - $predicted_star_end;
			$predicted_star_end{$predicted_miRNA_id} = $precursor_end{$predicted_miRNA_id} - $predicted_star_bg;
					
			}
		
		
		}
	close FILE;
	
}



sub compare_position{
	my $key=();
	my $id=();
	foreach $key (sort(keys %predicted_miRNA_bg)){
		foreach $id (sort(keys %miR_chr)){
			my $position_shift1 = $predicted_miRNA_end{$key} - $miR_B{$id};
			my $position_shift2 = $miR_E{$id} - $predicted_miRNA_bg{$key};
			my $position_star_shift1 = $predicted_star_end{$key} - $miR_B{$id};
			my $position_star_shift2 = $miR_E{$id} - $predicted_star_bg{$key};
			if($miR_chr{$id} eq $predicted_miRNA_chr{$key}){
				if(($position_shift1 >= 15)and ($position_shift2 >= 15)){
					print $key."\t".$id."\t".$predicted_miRNA_bg{$key}."\t".$predicted_miRNA_end{$key}."\t".$miR_B{$id}."\t".$miR_E{$id}."\t".$miRNA_strand{$key}."\t".$MIR_ID{$id}."\t"."normal"."\n";
							
				}
				if(($position_star_shift1 >= 15) and ($position_star_shift2 >= 15)){
					print $key."\t".$id."\t".$predicted_star_bg{$key}."\t".$predicted_star_end{$key}."\t".$miR_B{$id}."\t".$miR_E{$id}."\t".$miRNA_strand{$key}."\t".$MIR_ID{$id}."\t"."*"."\n";
				
				}
						
			}
			}
				
		}

}


