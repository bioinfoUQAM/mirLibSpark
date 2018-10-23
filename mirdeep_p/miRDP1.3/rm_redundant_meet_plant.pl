#!/usr/bin/perl -w

use strict;
use Getopt::Std;


######################################### USAGE ################################

my $usage =
"$0 chr_length precursor retrieved_miR nr_prediction plant_criteria

This script is used to remove the identical MIR retrieved by miRDeep.
The file, chr_length, is the info on the length of chromosomes.
The file, precursors, is the fasta file achieved by excise_candidate.pl, one script in miRDeep package.
The file, retrieved_miR, is the result of retrieved or predicted miRs achieved by miRDeep.
The file, nr_prediction, is the output name of prediction items which do not include redundant ones
The file, plant_criteria, is the output name of prediction items which meet the criteria of plant miRNA 
";


####################################### INPUT FILES ############################


my $chr_length=shift or die $usage;
my $precursor=shift or die $usage;
my $retrieved_miR=shift or die $usage;
my $nr_prediction=shift or die $usage;
my $meet_criteria=shift or die $usage;


##################################### GLOBAL VARIBALES #########################


my %chr_length=();


my %miRNA_strand=();
my %precursor_bg=();
my %precursor_end=();
my %precursor_chr=();


my %miRNA_m_id=();
my %miRNA_bg=();
my %miRNA_end=();
my %miRNA_seq=();
my %miRNA_struct=();
my %miRNA_arm=();

my %star_bg=();
my %star_end=();
my %star_seq=();
my %star_struct=();

my %pre_id=();
my %pre_seq=();
my %pre_bg=();
my %pre_end=();

my %item=();

my @tmp_XX=();
######################################### MAIN ################################# 

#extract known miRNA information
parse_chr($chr_length);

#print $chr_length{chr01}."\n";

#parse excised precursor information
parse_precursor($precursor);

#print $miRNA_strand{"chr05_112931"}."\t".$precursor_bg{"chr05_112931"}."\t".$precursor_end{"chr05_112931"}."\t".$precursor_chr{"chr05_112931"}."\n";

#extract predicted miRNA information
parse_miRDeep_result($retrieved_miR);

#print $miRNA_chr{"chr03_70799"}."\t".$miRNA_bg{"chr03_70799"}."\t".$miRNA_end{"chr03_70799"}."\t".$star_bg{"chr03_70799"}."\t".$star_end{"chr03_70799"}."\n";
#print $miRNA_chr{"osa-MIR1317_2"}."\t".$miRNA_bg{"osa-MIR1317_2"}."\t".$miRNA_end{"osa-MIR1317_2"}."\t".$star_bg{"osa-MIR1317_2"}."\t".$star_end{"osa-MIR1317_2"}."\n";
#print %pre_seq;
#compare known and predicted miRNA
remove_identical();

#print items which do not include the redundant ones
print_out($nr_prediction);

#filter predicted items which not meet the plant miRNA/miRNA* duplex 
meet_plant_criteria();

#print items which not only are non-redundant and meet the plant criteria
print_out($meet_criteria);

#
exit;


############################################ SUBROUTINES #######################

sub parse_chr{
	my $file=shift @_;
	open(FILE,"$file")||die"can not open the chr_length file:$!\n";
	
	while(<FILE>){
		chomp;
		my @tmp=split("\t",$_);
		$chr_length{$tmp[0]}=$tmp[1];	
		
		}
		
	close FILE;
	
	}

#

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
					$precursor_chr{$tmp_miRNA_id}=$tmp_precursor_id;
					
					if($tmp1[1] eq "+"){
						$precursor_bg{$tmp_miRNA_id} =  $tmp2[1];
						$precursor_end{$tmp_miRNA_id} =  $tmp3[1];
						
						}
					if($tmp1[1] eq "-"){
						$precursor_bg{$tmp_miRNA_id} = $chr_length{$tmp_precursor_id}  - $tmp3[1];
						$precursor_end{$tmp_miRNA_id} = $chr_length{$tmp_precursor_id}  - $tmp2[1];
						
						}

				}
		}

	}
	close FILE;
	
}




sub parse_miRDeep_result{
	#extract mature miRNA position predicted by miRDeep
	my $file=shift @_;
	
	my $id=();
	my $item=();
	
	
	my $miRNA_m_id=();
	my $miRNA_seq=();
	my $miRNA_bg=();
	my $miRNA_end=();
	my $miRNA_struct=();
	my $miRNA_arm=();
	
	my $star_bg=();
	my $star_end=();
	my $star_seq=();
	my $star_struct=();
	
	my $pre_id=();
	my $pre_bg=();
	my $pre_end=();
	my $pre_seq=();
	
	
	$/="\n\n\n";
	
	open(FILE,"$file")||die"can not open the miRDeep result file:$!\n";
	while(<FILE>){
		chomp;
		$id++;
	
		if($_=~/mature_arm\s+(\S+)\n/){
			$miRNA_arm=$1;
			}
		if($_=~/mature_beg\s+(\S+)\n/){
			$miRNA_bg=$1;
			}
		if($_=~/mature_end\s+(\S+)/){
			$miRNA_end=$1;
			}
		if($_=~/mature_query\s+(\S+)\n/){
			$miRNA_m_id=$1;
			}
		if($_=~/mature_seq\s+(\S+)\n/){
			$miRNA_seq=$1;
			}
		if($_=~/mature_struct\s+(\S+)\n/){
			$miRNA_struct=$1;
			}
		if($_=~/pre_seq\s+(\S+)\n/){
			$pre_seq=$1;			
			}
		if($_=~/pri_id\s+(\S+)\n/){
			$pre_id=$1;
			}
		if($_=~/star_beg\s+(\S+)/){
			$star_bg=$1;
			}
		if($_=~/star_end\s+(\S+)/){
			$star_end=$1;
			}
		if($_=~/star_seq\s+(\S+)\n/){
			$star_seq=$1;
			}
		if($_=~/star_struct\s+(\S+)\n/){
			$star_struct=$1;
			}
		
		$miRNA_m_id{$id}=$miRNA_m_id;
		$miRNA_bg{$id}=$miRNA_bg;
		$miRNA_end{$id}=$miRNA_end;
		$miRNA_arm{$id}=$miRNA_arm;
		$miRNA_seq{$id}=$miRNA_seq;
		$miRNA_struct{$id}=$miRNA_struct;
		
		$star_bg{$id}=$star_bg;
		$star_end{$id}=$star_end;
		$star_seq{$id}=$star_seq;
		$star_struct{$id}=$star_struct;
		
		$pre_seq{$id}=$pre_seq;
		$pre_id{$id}=$pre_id;
				
		
		my $item=$_;
		$item{$id} = $item;
		
	
		if($miRNA_strand{$pre_id} eq "+"){
			$miRNA_bg{$id} = $precursor_bg{$pre_id} + $miRNA_bg;
			$miRNA_end{$id} = $precursor_bg{$pre_id} + $miRNA_end;
			$star_bg{$id} = $precursor_bg{$pre_id} + $star_bg;
			$star_end{$id} = $precursor_bg{$pre_id} + $star_end;
			
			if($miRNA_arm eq "first"){
				$pre_bg=$miRNA_bg{$id};
				$pre_end=$star_end{$id};
				
				}
			if($miRNA_arm eq "second"){
				$pre_bg=$star_bg{$id};
				$pre_end=$miRNA_end{$id};
				
				}
			
			}
		if($miRNA_strand{$pre_id} eq "-"){
			$miRNA_bg{$id} = $precursor_end{$pre_id} - $miRNA_end;
			$miRNA_end{$id} = $precursor_end{$pre_id} - $miRNA_bg;
			$star_bg{$id} = $precursor_end{$pre_id} - $star_end;
			$star_end{$id} = $precursor_end{$pre_id} - $star_bg;
					
			if($miRNA_arm eq "first"){
				$pre_bg=$star_bg{$id};
				$pre_end=$miRNA_end{$id};
								
				}
			if($miRNA_arm eq "second"){
				$pre_bg=$miRNA_bg{$id};
				$pre_end=$star_end{$id};
				
				}
								
			}
		
		$pre_bg{$id}=$pre_bg;
		$pre_end{$id}=$pre_end;
		
		}
	close FILE;
	
}

#


sub remove_identical{
	
	my @MIR=sort(keys %pre_seq);
		
	foreach my $key1 (@MIR){

			foreach my $key2 (@MIR){
				if(($key1 != $key2) and ($precursor_chr{$pre_id{$key1}} eq $precursor_chr{$pre_id{$key2}}) and ($miRNA_bg{$key1} eq $miRNA_bg{$key2})){
					delete $pre_seq{$key2};
					}
				}
			@MIR=sort(keys %pre_seq);
		}

}



sub meet_plant_criteria{
	
	my @MIR1=sort(keys %pre_seq);

	foreach my $key (@MIR1){
		my $mature_stem=substr($miRNA_struct{$key},0,-2);
		my $star_stem=substr($star_struct{$key},0,-2);
		my @mature_mismatch=$mature_stem =~/\./g;
		my @star_mismatch=$star_stem =~/\./g;

		if((@mature_mismatch >4) or (@star_mismatch >4) or (abs(@mature_mismatch - @star_mismatch)>2)){
			
			delete $pre_seq{$key};
			}
					
		}
		
	}
	
	


sub print_out{
	my $file=shift @_;
  open(OUT,">$file")||die"can not open the output file:$!\n";
	my @MIR2=sort(keys %pre_seq);
		
	foreach my $key (@MIR2){
		print OUT $precursor_chr{$pre_id{$key}}."\t".$miRNA_strand{$pre_id{$key}}."\t".$miRNA_m_id{$key}."\t".$pre_id{$key}."\t".$miRNA_bg{$key}.'..'.$miRNA_end{$key}."\t".$pre_bg{$key}.'..'.$pre_end{$key}."\t".$miRNA_seq{$key}."\t".$pre_seq{$key}."\n";
		}
	close OUT;
	}
	
	
	