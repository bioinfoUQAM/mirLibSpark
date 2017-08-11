#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Std;


my $usage =
"$0 file_fasta file_precursors_gff

This script extended annotated microRNA precursor sequences, 
increasing 500 bp at both ends from genomic or BAC sequences.
The fasta file given as input should be the genome sequences 
in fasta format, and the file in precursors gff form should 
be General Feature Format which can be download from miRBase,
like ath.gff. 
";


my $file_fasta=shift or die $usage;
my $file_precursors_gff=shift or die $usage;

my %hash_fasta;
my %hash_chr;
my %hash_strand;
my %hash_start;
my %hash_end;

parse_file_fasta(\$file_fasta,\%hash_fasta);
parse_file_precursors_gff($file_precursors_gff);
excise();

exit;


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



sub parse_file_precursors_gff{

    my($file)=@_;

    open(FILENAME, $file) or die "Could not open the gff file $file";
    
    while (my $line = <FILENAME>){
          
          #the file format like: chr1	.	miRNA	28500	28706	.	+	.	ACC="MI0005394"; ID="ath-MIR838";
			
			if($line=~m/^(\S+)\s+\.\s+\S+\s+(\d+)\s+(\d+)\s+\.\s+(\S+)\s+\.\s+\S+\s+ID\=\"(\S+)\"\;/){
			
	    my $miR_id=$5;
	    	    
	    $hash_chr{$miR_id}=$1;
	    $hash_start{$miR_id}=$2;
	    $hash_end{$miR_id}=$3;
	    $hash_strand{$miR_id}=$4;
	    
	    }
   }
    close FILENAME;
}



sub excise{
			
		foreach my $miR_id (sort keys %hash_chr){
			
			my $chr_length=length $hash_fasta{$hash_chr{$miR_id}};

			$hash_start{$miR_id} = $hash_start{$miR_id} - 500;
			my $start = $hash_start{$miR_id} - 1;
			if($hash_start{$miR_id} <= 0){				
				$hash_start{$miR_id} = 0;
				$start = 0;				
				}
				
			$hash_end{$miR_id} = $hash_end{$miR_id} + 500;
			if($hash_end{$miR_id} >= $chr_length){
				$hash_end{$miR_id}=$chr_length;
				}
			
			my $extended_length=$hash_end{$miR_id} - $hash_start{$miR_id};
			
			my $seq=substr($hash_fasta{$hash_chr{$miR_id}},$start,$extended_length);
			
			if($hash_strand{$miR_id} eq '-'){
				$seq = RC_seq($seq);
				}
			
			print ">".$miR_id."\n";
			for ( my $pos = 0 ; $pos < length($seq) ; $pos += 60 ) {
        print substr($seq, $pos, 60),"\n";
    	}

			}
	
	}



sub RC_seq{
		
		my($seq)=@_;
		$seq=~tr/ATCGNatcgn/TAGCNtagcn/;
		$seq = reverse $seq;
		
		return $seq;
	
	}