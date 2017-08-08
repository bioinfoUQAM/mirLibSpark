#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Std;



################################# MIRDEEP #################################################

################################## USAGE ##################################################


my $usage=
"$0 file_signature file_structure

This is the core algorithm of miRDeep. It takes as input a file in blastparsed format with
information on the positions of reads aligned to potential precursor sequences (signature).
It also takes as input an RNAfold output file, giving information on the sequence, structure
and mimimum free energy of the potential precursor sequences.

Extra arguments can be given. -s specifies a fastafile containing the known mature miRNA
sequences that should be considered for conservation purposes. -t prints out the potential
precursor sequences that do _not_ exceed the cut-off (default prints out the sequences that
exceeds the cut-off). -u gives limited output, that is only the ids of the potential precursors
that exceed the cut-off. -v varies the cut-off. -x is a sensitive option for Sanger sequences
obtained through conventional cloning. -z consider the number of base pairings in the lower
stems (this option is not well tested).

-h print this usage
-s fasta file with known miRNAs
-t print filtered
-u limited output (only ids)
-v cut-off (default 1)
-x sensitive option for Sanger sequences
-y use Randfold
-z consider Drosha processing
";





############################################################################################

################################### INPUT ##################################################


#signature file in blast_parsed format
my $file_blast_parsed=shift or die $usage;

#structure file outputted from RNAfold
my $file_struct=shift or die $usage;

#options
my %options=();
getopts("hs:tuv:xyz",\%options);






#############################################################################################

############################# GLOBAL VARIABLES ##############################################


#parameters
my $nucleus_lng=7;

my $score_star=3.9;
my $score_star_not=-1.3;
my $score_nucleus=3;
my $score_nucleus_not=-0.6;
my $score_randfold=1.6;
my $score_randfold_not=-2.2;
my $score_intercept=0.3;
my @scores_stem=(-3.1,-2.3,-2.2,-1.6,-1.5,0.1,0.6,0.8,0.9,0.9,0);
my $score_min=1;
if($options{v}){$score_min=$options{v};}
if($options{x}){$score_min=-5;}

my $e=2.718281828;

#hashes
my %hash_desc;
my %hash_seq;
my %hash_struct;
my %hash_mfe;
my %hash_nuclei;
my %hash_mirs;
my %hash_query;
my %hash_comp;
my %hash_bp;

#other variables
my $subject_old;
my $message_filter;
my $message_score;
my $lines;
my $out_of_bound;




##############################################################################################

################################  MAIN  ###################################################### 


#print help if that option is used
if($options{h}){die $usage;}

#parse structure file outputted from RNAfold
parse_file_struct($file_struct);

#if conservation is scored, the fasta file of known miRNA sequences is parsed
if($options{s}){create_hash_nuclei($options{s})};

#parse signature file in blast_parsed format and resolve each potential precursor
parse_file_blast_parsed($file_blast_parsed);

exit;




##############################################################################################

############################## SUBROUTINES ###################################################



sub parse_file_blast_parsed{

#    read through the signature blastparsed file, fills up a hash with information on queries
#    (deep sequences) mapping to the current subject (potential precursor) and resolve each
#    potential precursor in turn
 
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
	    my $pid=$10;
	    my $bitscore=$11;
	    my $other=$12;
	    
	    #if the new line concerns a new subject (potential precursor) then the old subject must be resolved
	    if($subject_old and $subject_old ne $subject){
		resolve_potential_precursor();
	    }

	    #resolve the strand
	    my $strand=find_strand($other);

	    #resolve the number of reads that the deep sequence represents
	    my $freq=find_freq($query);

	    #read information of the query (deep sequence) into hash
	    $hash_query{$query}{"subject_beg"}=$subject_beg;
	    $hash_query{$query}{"subject_end"}=$subject_end;
	    $hash_query{$query}{"strand"}=$strand;
	    $hash_query{$query}{"freq"}=$freq;

	    #save the signature information
	    $lines.=$line;

	    $subject_old=$subject;
	}
    }
    resolve_potential_precursor();
}

sub resolve_potential_precursor{
    
#    dissects the potential precursor in parts by filling hashes, and tests if it passes the
#    initial filter and the scoring filter

#    binary variable whether the potential precursor is still viable
    my $ret=1;

    fill_structure();
    
    fill_pri();

    fill_mature();
   
    fill_star();

    fill_loop();

    fill_lower_flanks();

#    do_test_assemble();

#    this is the actual classification
    unless(pass_filtering_initial() and pass_threshold_score()){$ret=0;}

    print_results($ret);
    
    reset_variables();
    
    return;
    
}



sub print_results{

    my $ret=shift;

#    print out if the precursor is accepted and accepted precursors should be printed out
#    or if the potential precursor is discarded and discarded potential precursors should
#    be printed out
	
    if((!$options{t} and $ret) or ($options{t} and !$ret)){
	#full output
	unless($options{u}){
	    if($message_filter){print $message_filter;}
	    if($message_score){print $message_score;}
	    print_hash_comp();
	    print $lines."\n\n";
	    return;
	}
	#limited output (only ids)
	my $id=$hash_comp{"pri_id"};
	print "$id\n";
    }    
}







sub pass_threshold_score{

#    this is the scoring

    #minimum free energy of the potential precursor
    my $score_mfe=score_mfe($hash_comp{"pri_mfe"});

    #count of reads that map in accordance with Dicer processing
    my $score_freq=score_freq($hash_comp{"freq"});

    #basic score
    my $score=$score_mfe+$score_freq;

    #scoring of conserved nucleus/seed (optional)
    if($options{s}){

	#if the nucleus is conserved
	if(test_nucleus_conservation()){

	    #nucleus from position 2-8
	    my $nucleus=substr($hash_comp{"mature_seq"},1,$nucleus_lng);

	    #resolve DNA/RNA ambiguities
	    $nucleus=~tr/[T]/[U]/;

	    #print score contribution
	    score_s("score_nucleus\t$score_nucleus");

	    #print the ids of known miRNAs with same nucleus
	    score_s("$hash_mirs{$nucleus}");

	    #add to score
	    $score+=$score_nucleus;

	#if the nucleus is not conserved
	}else{
	    #print (negative) score contribution
	    score_s("score_nucleus\t$score_nucleus_not");

	    #add (negative) score contribution
	    $score+=$score_nucleus_not;
	}
    }
   
    #if the majority of potential star reads fall as expected from Dicer processing
    if($hash_comp{"star_read"}){
	score_s("score_star\t$score_star");
	$score+=$score_star;
    }else{
	score_s("score_star\t$score_star_not");
	$score+=$score_star_not;
    }

    #score lower stems for potential for Drosha recognition (highly optional)
    if($options{z}){
	my $stem_bp=$hash_comp{"stem_bp"};
	my $score_stem=$scores_stem[$stem_bp];
	$score+=$score_stem;
	score_s("score_stem\t$score_stem");
    }


    $score+=$score_intercept;

    #score for randfold (optional)
    if($options{y}){

#	only calculate randfold value if it can make the difference between the potential precursor
#	being accepted or discarded
	if($score+$score_randfold>=$score_min and $score+$score_randfold_not<=$score_min){

	    #randfold value<0.05
	    if(test_randfold()){$score+=$score_randfold;score_s("score_randfold\t$score_randfold");}

	    #randfold value>0.05
	    else{$score+=$score_randfold_not;score_s("score_randfold\t$score_randfold_not");}
	}
    }

    #round off values to one decimal
    my $round_mfe=round($score_mfe*10)/10;
    my $round_freq=round($score_freq*10)/10;
    my $round=round($score*10)/10;

    #print scores
    score_s("score_mfe\t$round_mfe\nscore_freq\t$round_freq\nscore\t$round");
 
    #return 1 if the potential precursor is accepted, return 0 if discarded
    unless($score>=$score_min){return 0;}
    return 1;
}

sub test_randfold{

    #print sequence to temporary file, test randfold value, return 1 or 0

    print_file("pri_seq.fa",">pri_seq\n".$hash_comp{"pri_seq"});

    my $p_value=`randfold -s pri_seq.fa 999 | cut -f 3`;

    system "rm pri_seq.fa";

    if($p_value<=0.05){return 1;}

    return 0;
}


sub print_file{

    #print string to file

    my($file,$string)=@_;

    open(FILE, ">$file");
    print FILE "$string";
    close FILE;
}


sub test_nucleus_conservation{

    #test if nucleus is identical to nucleus from known miRNA, return 1 or 0

    my $nucleus=substr($hash_comp{"mature_seq"},1,$nucleus_lng);
    $nucleus=~tr/[T]/[U]/;
    if($hash_nuclei{$nucleus}){return 1;}
    
    return 0;
}



sub pass_filtering_initial{

    #test if the structure forms a plausible hairpin
    unless(pass_filtering_structure()){filter_p("structure problem"); return 0;}

    #test if >90% of reads map to the hairpin in consistence with Dicer processing
    unless(pass_filtering_signature()){filter_p("signature problem");return 0;}

    return 1;

}


sub pass_filtering_signature{

    #number of reads that map in consistence with Dicer processing
    my $consistent=0;

    #number of reads that map inconsistent with Dicer processing
    my $inconsistent=0;
   
#   number of potential star reads map in good consistence with Drosha/Dicer processing
#   (3' overhangs relative to mature product)
    my $star_perfect=0;

#   number of potential star reads that do not map in good consistence with 3' overhang
    my $star_fuzzy=0;
 

    #sort queries (deep sequences) by their position on the hairpin
    my @queries=sort {$hash_query{$a}{"subject_beg"} <=> $hash_query{$b}{"subject_beg"}} keys %hash_query;

    foreach my $query(@queries){

	#number of reads that the deep sequence represents
	unless(defined($hash_query{$query}{"freq"})){next;}
	my $query_freq=$hash_query{$query}{"freq"};

	#test which Dicer product (if any) the deep sequence corresponds to
	my $product=test_query($query);

	#if the deep sequence corresponds to a Dicer product, add to the 'consistent' variable
	if($product){$consistent+=$query_freq;}

	#if the deep sequence do not correspond to a Dicer product, add to the 'inconsistent' variable
	else{$inconsistent+=$query_freq;}

	#test a potential star sequence has good 3' overhang
	if($product eq "star"){
	    if(test_star($query)){$star_perfect+=$query_freq;}
	    else{$star_fuzzy+=$query_freq;}
	}
    }

#   if the majority of potential star sequences map in good accordance with 3' overhang
#    score for the presence of star evidence
    if($star_perfect>$star_fuzzy){$hash_comp{"star_read"}=1;}

    #total number of reads mapping to the hairpin
    my $freq=$consistent+$inconsistent;
    $hash_comp{"freq"}=$freq;
    unless($freq>0){filter_s("read frequency too low"); return 0;}

    #unless >90% of the reads map in consistence with Dicer processing, the hairpin is discarded
    my $inconsistent_fraction=$inconsistent/($inconsistent+$consistent);
    unless($inconsistent_fraction<=0.1){filter_p("inconsistent\t$inconsistent\nconsistent\t$consistent"); return 0;}

    #the hairpin is retained
    return 1;
}

sub test_star{

    #test if a deep sequence maps in good consistence with 3' overhang

    my $query=shift;

    #5' begin and 3' end positions
    my $beg=$hash_query{$query}{"subject_beg"};
    my $end=$hash_query{$query}{"subject_end"};

    #the difference between observed and expected begin positions must be 0 or 1
    my $offset=$beg-$hash_comp{"star_beg"};
    if($offset==0 or $offset==1 or $offset==-1){return 1;}

    return 0;
}



sub test_query{

    #test if deep sequence maps in consistence with Dicer processing
    
    my $query=shift;

    #begin, end, strand and read count
    my $beg=$hash_query{$query}{"subject_beg"};
    my $end=$hash_query{$query}{"subject_end"};
    my $strand=$hash_query{$query}{"strand"};
    my $freq=$hash_query{$query}{"freq"};

    #should not be on the minus strand (although this has in fact anecdotally been observed for known miRNAs)
    if($strand eq '-'){return 0;}

    #the deep sequence is allowed to stretch 2 nt beyond the expected 5' end
    my $fuzz_beg=2;
    #the deep sequence is allowed to stretch 5 nt beyond the expected 3' end
    my $fuzz_end=5;

    #if in accordance with Dicer processing, return the type of Dicer product 
    if(contained($beg,$end,$hash_comp{"mature_beg"}-$fuzz_beg,$hash_comp{"mature_end"}+$fuzz_end)){return "mature";}
    if(contained($beg,$end,$hash_comp{"star_beg"}-$fuzz_beg,$hash_comp{"star_end"}+$fuzz_end)){return "star";}
    if(contained($beg,$end,$hash_comp{"loop_beg"}-$fuzz_beg,$hash_comp{"loop_end"}+$fuzz_end)){return "loop";}
  
    #if not in accordance, return 0
    return 0;
}


sub pass_filtering_structure{

    #The potential precursor must form a hairpin with miRNA precursor-like characteristics

    #return value
    my $ret=1;

    #potential mature, star, loop and lower flank parts must be identifiable
    unless(test_components()){return 0;}

    #no bifurcations
    unless(no_bifurcations_precursor()){$ret=0;}

    #minimum 14 base pairings in duplex
    unless(bp_duplex()>=14){$ret=0; filter_s("too few pairings in duplex");}

    #not more than 6 nt difference between mature and star length
    unless(-6<diff_lng() and diff_lng()<6){$ret=0; filter_s("too big difference between mature and star length") }

    return $ret;
}






sub test_components{

    #tests whether potential mature, star, loop and lower flank parts are identifiable

    unless($hash_comp{"mature_struct"}){
	filter_s("no mature");
	return 0;
    }

    unless($hash_comp{"star_struct"}){
	filter_s("no star");
	return 0;
    }

    unless($hash_comp{"loop_struct"}){
	filter_s("no loop");
   	return 0;
    }

    unless($hash_comp{"flank_first_struct"}){
	filter_s("no flanks");
   	return 0;
    }
     
    unless($hash_comp{"flank_second_struct"}){
	filter_s("no flanks");
    	return 0;
    }
    return 1;
}





sub no_bifurcations_precursor{

    #tests whether there are bifurcations in the hairpin

    #assembles the potential precursor sequence and structure from the expected Dicer products
    #this is the expected biological precursor, in contrast with 'pri_seq' that includes
    #some genomic flanks on both sides

    my $pre_struct;
    my $pre_seq;
    if($hash_comp{"mature_arm"} eq "first"){
	$pre_struct.=$hash_comp{"mature_struct"}.$hash_comp{"loop_struct"}.$hash_comp{"star_struct"};
	$pre_seq.=$hash_comp{"mature_seq"}.$hash_comp{"loop_seq"}.$hash_comp{"star_seq"};
    }else{
	$pre_struct.=$hash_comp{"star_struct"}.$hash_comp{"loop_struct"}.$hash_comp{"mature_struct"};
	$pre_seq.=$hash_comp{"star_seq"}.$hash_comp{"loop_seq"}.$hash_comp{"mature_seq"};
    }

    #read into hash
    $hash_comp{"pre_struct"}=$pre_struct;
    $hash_comp{"pre_seq"}=$pre_seq;

    #simple pattern matching checks for bifurcations
    unless($pre_struct=~/^((\.|\()+..(\.|\))+)$/){
	filter_s("bifurcation in precursor");
	return 0;
    }

    return 1;
}

sub bp_precursor{
 
    #total number of bps in the precursor
 
    my $pre_struct=$hash_comp{"pre_struct"};

    #simple pattern matching
    my $pre_bps=0;
    while($pre_struct=~/\(/g){
	$pre_bps++;
    }
    return $pre_bps;
}


sub bp_duplex{

    #total number of bps in the duplex

    my $duplex_bps=0;
    my $mature_struct=$hash_comp{"mature_struct"};

    #simple pattern matching
    while($mature_struct=~/(\(|\))/g){
	$duplex_bps++;
    }
    return $duplex_bps;
}

sub diff_lng{

    #find difference between mature and star lengths

    my $mature_lng=length $hash_comp{"mature_struct"};
    my $star_lng=length $hash_comp{"star_struct"};
    my $diff_lng=$mature_lng-$star_lng;
    return $diff_lng;
}



sub do_test_assemble{

#    not currently used, tests if the 'pri_struct' as assembled from the parts (Dicer products, lower flanks)
#    is identical to 'pri_struct' before disassembly into parts

    my $assemble_struct;

    if($hash_comp{"flank_first_struct"} and $hash_comp{"mature_struct"} and $hash_comp{"loop_struct"} and $hash_comp{"star_struct"} and $hash_comp{"flank_second_struct"}){
	if($hash_comp{"mature_arm"} eq "first"){
	    $assemble_struct.=$hash_comp{"flank_first_struct"}.$hash_comp{"mature_struct"}.$hash_comp{"loop_struct"}.$hash_comp{"star_struct"}.$hash_comp{"flank_second_struct"};
	}else{
	    $assemble_struct.=$hash_comp{"flank_first_struct"}.$hash_comp{"star_struct"}.$hash_comp{"loop_struct"}.$hash_comp{"mature_struct"}.$hash_comp{"flank_second_struct"};
	}
	unless($assemble_struct eq $hash_comp{"pri_struct"}){
	    $hash_comp{"test_assemble"}=$assemble_struct;
	    print_hash_comp();
	}
    }
    return;
 }



sub fill_structure{

    #reads the dot bracket structure into the 'bp' hash where each key and value are basepaired

    my $struct=$hash_struct{$subject_old};
    my $lng=length $struct;

    #local stack for keeping track of basepairings
    my @bps;

    for(my $pos=1;$pos<=$lng;$pos++){
	my $struct_pos=excise_struct($struct,$pos,$pos,"+");

	if($struct_pos eq "("){
	    push(@bps,$pos);
	}

	if($struct_pos eq ")"){
	    my $pos_prev=pop(@bps);
	    $hash_bp{$pos_prev}=$pos;
	    $hash_bp{$pos}=$pos_prev;
	}
    }
    return;
}



sub fill_star{

    #fills specifics on the expected star strand into 'comp' hash ('component' hash)
    
    #if the mature sequence is not plausible, don't look for the star arm
    my $mature_arm=$hash_comp{"mature_arm"};
    unless($mature_arm){$hash_comp{"star_arm"}=0; return;}
 
    #if the star sequence is not plausible, don't fill into the hash
    my($star_beg,$star_end)=find_star();
    my $star_arm=arm_star($star_beg,$star_end);
    unless($star_arm){return;}

    #excise expected star sequence and structure
    my $star_seq=excise_seq($hash_comp{"pri_seq"},$star_beg,$star_end,"+");
    my $star_struct=excise_seq($hash_comp{"pri_struct"},$star_beg,$star_end,"+");

    #fill into hash
    $hash_comp{"star_beg"}=$star_beg;
    $hash_comp{"star_end"}=$star_end;
    $hash_comp{"star_seq"}=$star_seq;
    $hash_comp{"star_struct"}=$star_struct;
    $hash_comp{"star_arm"}=$star_arm;

    return;
}


sub find_star{

    #uses the 'bp' hash to find the expected star begin and end positions from the mature positions

    #the -2 is for the overhang
    my $mature_beg=$hash_comp{"mature_beg"};
    my $mature_end=$hash_comp{"mature_end"}-2;
    my $mature_lng=$mature_end-$mature_beg+1;

    #in some cases, the last nucleotide of the mature sequence does not form a base pair,
    #and therefore does not basepair with the first nucleotide of the star sequence.
    #In this case, the algorithm searches for the last nucleotide of the mature sequence
    #to form a base pair. The offset is the number of nucleotides searched through.
    my $offset_star_beg=0;
    my $offset_beg=0;

    #the offset should not be longer than the length of the mature sequence, then it
    #means that the mature sequence does not form any base pairs
    while(!$offset_star_beg and $offset_beg<$mature_lng){
	if($hash_bp{$mature_end-$offset_beg}){
	    $offset_star_beg=$hash_bp{$mature_end-$offset_beg};
	}else{
	    $offset_beg++;
	}
    }
    #when defining the beginning of the star sequence, compensate for the offset
    my $star_beg=$offset_star_beg-$offset_beg;

    #same as above
    my $offset_star_end=0;
    my $offset_end=0;
    while(!$offset_star_end and $offset_end<$mature_lng){
	if($hash_bp{$mature_beg+$offset_end}){
	    $offset_star_end=$hash_bp{$mature_beg+$offset_end};
	}else{
	    $offset_end++;
	}
    }
    #the +2 is for the overhang
    my $star_end=$offset_star_end+$offset_end+2;

    return($star_beg,$star_end);
}


sub fill_pri{

    #fills basic specifics on the precursor into the 'comp' hash
    
    my $seq=$hash_seq{$subject_old};
    my $struct=$hash_struct{$subject_old};
    my $mfe=$hash_mfe{$subject_old};
    my $length=length $seq;
    
    $hash_comp{"pri_id"}=$subject_old;
    $hash_comp{"pri_seq"}=$seq;
    $hash_comp{"pri_struct"}=$struct;
    $hash_comp{"pri_mfe"}=$mfe;
    $hash_comp{"pri_beg"}=1;
    $hash_comp{"pri_end"}=$length;
    
    return;
}


sub fill_mature{

    #fills specifics on the mature sequence into the 'comp' hash

    my $mature_query=find_mature_query();
    my($mature_beg,$mature_end)=find_positions_query($mature_query);
    my $mature_strand=find_strand_query($mature_query);
    my $mature_seq=excise_seq($hash_comp{"pri_seq"},$mature_beg,$mature_end,$mature_strand);
    my $mature_struct=excise_struct($hash_comp{"pri_struct"},$mature_beg,$mature_end,$mature_strand);
    my $mature_arm=arm_mature($mature_beg,$mature_end,$mature_strand);

    $hash_comp{"mature_query"}=$mature_query;
    $hash_comp{"mature_beg"}=$mature_beg;
    $hash_comp{"mature_end"}=$mature_end;
    $hash_comp{"mature_strand"}=$mature_strand;
    $hash_comp{"mature_struct"}=$mature_struct;
    $hash_comp{"mature_seq"}=$mature_seq;
    $hash_comp{"mature_arm"}=$mature_arm;

    return;
}



sub fill_loop{

    #fills specifics on the loop sequence into the 'comp' hash

    #unless both mature and star sequences are plausible, do not look for the loop
    unless($hash_comp{"mature_arm"} and $hash_comp{"star_arm"}){return;}

    my $loop_beg;
    my $loop_end;

    #defining the begin and end positions of the loop from the mature and star positions
    #excision depends on whether the mature or star sequence is 5' of the loop ('first')
    if($hash_comp{"mature_arm"} eq "first"){
	$loop_beg=$hash_comp{"mature_end"}+1;
    }else{
	$loop_end=$hash_comp{"mature_beg"}-1;
    }
    
    if($hash_comp{"star_arm"} eq "first"){
	$loop_beg=$hash_comp{"star_end"}+1;
    }else{
	$loop_end=$hash_comp{"star_beg"}-1;
    }

    #unless the positions are plausible, do not fill into hash
    unless(test_loop($loop_beg,$loop_end)){return;}

    my $loop_seq=excise_seq($hash_comp{"pri_seq"},$loop_beg,$loop_end,"+");
    my $loop_struct=excise_struct($hash_comp{"pri_struct"},$loop_beg,$loop_end,"+");

    $hash_comp{"loop_beg"}=$loop_beg;
    $hash_comp{"loop_end"}=$loop_end;
    $hash_comp{"loop_seq"}=$loop_seq;
    $hash_comp{"loop_struct"}=$loop_struct;

    return;
}


sub fill_lower_flanks{

    #fills specifics on the lower flanks and unpaired strands into the 'comp' hash

    #unless both mature and star sequences are plausible, do not look for the flanks
    unless($hash_comp{"mature_arm"} and $hash_comp{"star_arm"}){return;}

    my $flank_first_end;
    my $flank_second_beg;

    #defining the begin and end positions of the flanks from the mature and star positions
    #excision depends on whether the mature or star sequence is 5' in the potenitial precursor ('first')
    if($hash_comp{"mature_arm"} eq "first"){
	$flank_first_end=$hash_comp{"mature_beg"}-1;
    }else{
	$flank_second_beg=$hash_comp{"mature_end"}+1;
    }
    
    if($hash_comp{"star_arm"} eq "first"){
	$flank_first_end=$hash_comp{"star_beg"}-1;
    }else{
	$flank_second_beg=$hash_comp{"star_end"}+1;
    }   

    #unless the positions are plausible, do not fill into hash
    unless(test_flanks($flank_first_end,$flank_second_beg)){return;}

    $hash_comp{"flank_first_end"}=$flank_first_end;
    $hash_comp{"flank_second_beg"}=$flank_second_beg;
    $hash_comp{"flank_first_seq"}=excise_seq($hash_comp{"pri_seq"},$hash_comp{"pri_beg"},$hash_comp{"flank_first_end"},"+");
    $hash_comp{"flank_second_seq"}=excise_seq($hash_comp{"pri_seq"},$hash_comp{"flank_second_beg"},$hash_comp{"pri_end"},"+");
    $hash_comp{"flank_first_struct"}=excise_struct($hash_comp{"pri_struct"},$hash_comp{"pri_beg"},$hash_comp{"flank_first_end"},"+");
    $hash_comp{"flank_second_struct"}=excise_struct($hash_comp{"pri_struct"},$hash_comp{"flank_second_beg"},$hash_comp{"pri_end"},"+");

    if($options{z}){
	fill_stems_drosha();
    }

    return;
}


sub fill_stems_drosha{

    #scores the number of base pairings formed by the first ten nt of the lower stems
    #in general, the more stems, the higher the score contribution
    #warning: this options has not been thoroughly tested

    my $flank_first_struct=$hash_comp{"flank_first_struct"};
    my $flank_second_struct=$hash_comp{"flank_second_struct"};
    
    my $stem_first=substr($flank_first_struct,-10);
    my $stem_second=substr($flank_second_struct,0,10);
    
    my $stem_bp_first=0;
    my $stem_bp_second=0;

    #find base pairings by simple pattern matching
    while($stem_first=~/\(/g){
	$stem_bp_first++;
    }
    
    while($stem_second=~/\)/g){
	$stem_bp_second++;
    }
    
    my $stem_bp=min2($stem_bp_first,$stem_bp_second);
    
    $hash_comp{"stem_first"}=$stem_first;
    $hash_comp{"stem_second"}=$stem_second;
    $hash_comp{"stem_bp_first"}=$stem_bp_first;
    $hash_comp{"stem_bp_second"}=$stem_bp_second;
    $hash_comp{"stem_bp"}=$stem_bp;
    
    return;
}




sub arm_mature{
 
    #tests whether the mature sequence is in the 5' ('first') or 3' ('second') arm of the potential precursor

    my ($beg,$end,$strand)=@_;
 
    #mature and star sequences should alway be on plus strand
    if($strand eq "-"){return 0;}

    #there should be no bifurcations and minimum one base pairing
    my $struct=excise_seq($hash_comp{"pri_struct"},$beg,$end,$strand);
    if(defined($struct) and $struct=~/^(\(|\.)+$/ and $struct=~/\(/){
	return "first";
    }elsif(defined($struct) and $struct=~/^(\)|\.)+$/ and $struct=~/\)/){
	return "second";
    }
    return 0;
}


sub arm_star{

    #tests whether the star sequence is in the 5' ('first') or 3' ('second') arm of the potential precursor

    my ($beg,$end)=@_;

    #unless the begin and end positions are plausible, test negative
    unless($beg>0 and $beg<=$hash_comp{"pri_end"} and $end>0 and $end<=$hash_comp{"pri_end"} and $beg<=$end){return 0;}

    #no overlap between the mature and the star sequence
    if($hash_comp{"mature_arm"} eq "first"){
	($hash_comp{"mature_end"}<$beg) or return 0;
    }elsif($hash_comp{"mature_arm"} eq "second"){
	($end<$hash_comp{"mature_beg"}) or return 0;
    }

    #there should be no bifurcations and minimum one base pairing
    my $struct=excise_seq($hash_comp{"pri_struct"},$beg,$end,"+");
    if($struct=~/^(\(|\.)+$/ and $struct=~/\(/){
	return "first";
    }elsif($struct=~/^(\)|\.)+$/ and $struct=~/\)/){
	return "second";
    }
    return 0;
}


sub test_loop{

    #tests the loop positions

    my ($beg,$end)=@_;

    #unless the begin and end positions are plausible, test negative
    unless($beg>0 and $beg<=$hash_comp{"pri_end"} and $end>0 and $end<=$hash_comp{"pri_end"} and $beg<=$end){return 0;}

    return 1;
}


sub test_flanks{

    #tests the positions of the lower flanks

    my ($beg,$end)=@_;

    #unless the begin and end positions are plausible, test negative
    unless($beg>0 and $beg<=$hash_comp{"pri_end"} and $end>0 and $end<=$hash_comp{"pri_end"} and $beg<=$end){return 0;}

    return 1;
}


sub comp{

    #subroutine to retrive from the 'comp' hash

    my $type=shift;
    my $component=$hash_comp{$type};
    return $component;
}


sub find_strand_query{

    #subroutine to find the strand for a given query

    my $query=shift;
    my $strand=$hash_query{$query}{"strand"};
    return $strand;
}


sub find_positions_query{

    #subroutine to find the begin and end positions for a given query

    my $query=shift;
    my $beg=$hash_query{$query}{"subject_beg"};
    my $end=$hash_query{$query}{"subject_end"};
    return ($beg,$end);
}



sub find_mature_query{

    #finds the query with the highest frequency of reads and returns it
    #is used to determine the positions of the potential mature sequence

    my @queries=sort {$hash_query{$b}{"freq"} <=> $hash_query{$a}{"freq"}} keys %hash_query;
    my $mature_query=$queries[0];
    return $mature_query;
}




sub reset_variables{

    #resets the hashes for the next potential precursor

    %hash_query=();
    %hash_comp=();
    %hash_bp=();

    $message_filter=();
    $message_score=();
    $lines=();

    return;
}



sub excise_seq{

    #excise sub sequence from the potential precursor

    my($seq,$beg,$end,$strand)=@_;

    #begin can be equal to end if only one nucleotide is excised
    unless($beg<=$end){print STDERR "begin can not be smaller than end for $subject_old\n";exit;}

    #rarely, permuted combinations of signature and structure cause out of bound excision errors.
    #this happens once appr. every two thousand combinations
    unless($beg<=length($seq)){$out_of_bound++;return 0;}

    #if on the minus strand, the reverse complement should be excised
    if($strand eq "-"){$seq=revcom($seq);}

    #the blast parsed format is 1-indexed, substr is 0-indexed
    my $sub_seq=substr($seq,$beg-1,$end-$beg+1);

    return $sub_seq;

}

sub excise_struct{

    #excise sub structure

    my($struct,$beg,$end,$strand)=@_;
    my $lng=length $struct;

    #begin can be equal to end if only one nucleotide is excised
    unless($beg<=$end){print STDERR "begin can not be smaller than end for $subject_old\n";exit;}

    #rarely, permuted combinations of signature and structure cause out of bound excision errors.
    #this happens once appr. every two thousand combinations
    unless($beg<=length($struct)){return 0;}

    #if excising relative to minus strand, positions are reversed
    if($strand eq "-"){($beg,$end)=rev_pos($beg,$end,$lng);}

    #the blast parsed format is 1-indexed, substr is 0-indexed
    my $sub_struct=substr($struct,$beg-1,$end-$beg+1);
 
    return $sub_struct;
}


sub create_hash_nuclei{
 
    #parses a fasta file with sequences of known miRNAs considered for conservation purposes
    #reads the nuclei into a hash

    my ($file) = @_;
    my ($id, $desc, $sequence, $nucleus) = ();

    open (FASTA, "<$file") or die "can not open $file\n";
    while (<FASTA>)
    {
        chomp;
        if (/^>(\S+)(.*)/)
	{
	    $id       = $1;
	    $desc     = $2;
	    $sequence = "";
	    $nucleus  = "";
	    while (<FASTA>){
                chomp;
                if (/^>(\S+)(.*)/){
		    $nucleus                = substr($sequence,1,$nucleus_lng);
		    $nucleus                =~ tr/[T]/[U]/;
		    $hash_mirs{$nucleus}   .="$id\t";
		    $hash_nuclei{$nucleus} += 1;

		    $id               = $1;
		    $desc             = $2;
		    $sequence         = "";
		    $nucleus          = "";
		    next;
                }
		$sequence .= $_;
            }
        }
    }
    $nucleus                = substr($sequence,1,$nucleus_lng);
    $nucleus                =~ tr/[T]/[U]/;
    $hash_mirs{$nucleus}   .="$id\t";
    $hash_nuclei{$nucleus} += 1;
    close FASTA;
}
    

sub parse_file_struct{
 
    #parses the output from RNAfoldand reads it into hashes

    my($file) = @_;
    my($id,$desc,$seq,$struct,$mfe) = ();

    open (FILE_STRUCT, "<$file") or die "can not open $file\n";
    while (<FILE_STRUCT>)
    {
        chomp;
        if (/^>(\S+)\s*(.*)/)
	{
	    $id          = $1;
	    $desc        = $2;
	    $seq         = "";
	    $struct      = "";
	    $mfe         = "";
	    while (<FILE_STRUCT>){
                chomp;
                if (/^>(\S+)\s*(.*)/){
		    $hash_desc{$id}   = $desc;
		    $hash_seq{$id}    = $seq;
		    $hash_struct{$id} = $struct;
		    $hash_mfe{$id}    = $mfe;

		    $id          = $1;
		    $desc        = $2;
		    $seq         = "";
		    $struct      = "";
		    $mfe         = "";

		    next;
                }
		if(/^\w/){
		    tr/uU/tT/;
		    $seq .= $_;
		}if(/((\.|\(|\))+)/){
		    $struct .=$1;
		}
		if(/\((\s*-\d+\.\d+)\)/){
		    $mfe = $1;
		}
	    
	    }
        }
    }

    $hash_desc{$id}        = $desc;
    $hash_seq{$id}         = $seq;
    $hash_struct{$id}      = $struct;
    $hash_mfe{$id}         = $mfe;

    close FILE_STRUCT;
    return;
}


sub score_s{

    #this score message is appended to the end of the string of score messages outputted for the potential precursor

    my $message=shift;
    $message_score.=$message."\n";;
    return;
}



sub score_p{

   #this score message is appended to the beginning of the string of score messages outputted for the potential precursor

    my $message=shift;
    $message_score=$message."\n".$message_score;
    return;
}



sub filter_s{

    #this filtering message is appended to the end of the string of filtering messages outputted for the potential precursor

    my $message=shift;
    $message_filter.=$message."\n";
    return;
}


sub filter_p{

    #this filtering message is appended to the beginning of the string of filtering messages outputted for the potential precursor

    my $message=shift;
    if(defined $message_filter){$message_filter=$message."\n".$message_filter;}
    else{$message_filter=$message."\n";}
    return;
}

    
sub find_freq{

    #finds the frequency of a given read query from its id.

    my($query)=@_;

    if($query=~/x(\d+)/){
	my $freq=$1;
	return $freq;
    }else{
	print STDERR "Problem with read format\n";
	return 0;
    }
}


sub print_hash_comp{

    #prints the 'comp' hash

    my @keys=sort keys %hash_comp;
    foreach my $key(@keys){
	my $value=$hash_comp{$key};
	print "$key  \t$value\n";
    }
}



sub print_hash_bp{

    #prints the 'bp' hash

    my @keys=sort {$a<=>$b} keys %hash_bp;
    foreach my $key(@keys){
	my $value=$hash_bp{$key};
	print "$key\t$value\n";
    }
    print "\n";
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


sub contained{

    #Is the stretch defined by the first positions contained in the stretch defined by the second?

    my($beg1,$end1,$beg2,$end2)=@_;

    testbeginend($beg1,$end1,$beg2,$end2);

    if($beg2<=$beg1 and $end1<=$end2){
	return 1;
    }else{
	return 0;
    }
}


sub testbeginend{

    #Are the beginposition numerically smaller than the endposition for each pair?

    my($begin1,$end1,$begin2,$end2)=@_;

    unless($begin1<=$end1 and $begin2<=$end2){
	print STDERR "beg can not be larger than end for $subject_old\n";
	exit;
    }
}


sub rev_pos{

#   The blast_parsed format always uses positions that are relative to the 5' of the given strand
#   This means that for a sequence of length n, the first nucleotide on the minus strand base pairs with
#   the n't nucleotide on the plus strand

#   This subroutine reverses the begin and end positions of positions of the minus strand so that they
#   are relative to the 5' end of the plus strand	
   
    my($beg,$end,$lng)=@_;
    
    my $new_end=$lng-$beg+1;
    my $new_beg=$lng-$end+1;
    
    return($new_beg,$new_end);
}

sub round {

    #rounds to nearest integer
   
    my($number) = shift;
    return int($number + .5);
    
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



sub score_freq{

#   scores the count of reads that map to the potential precursor
#   Assumes geometric distribution as described in methods section of manuscript

    my $freq=shift;

    #parameters of known precursors and background hairpins
    my $parameter_test=0.999;
    my $parameter_control=0.6;

    #log_odds calculated directly to avoid underflow
    my $intercept=log((1-$parameter_test)/(1-$parameter_control));
    my $slope=log($parameter_test/$parameter_control);
    my $log_odds=$slope*$freq+$intercept;

    #if no strong evidence for 3' overhangs, limit the score contribution to 0
    unless($options{x} or $hash_comp{"star_read"}){$log_odds=min2($log_odds,0);}
   
    return $log_odds;
}



sub score_mfe{

#   scores the minimum free energy in kCal/mol of the potential precursor
#   Assumes Gumbel distribution as described in methods section of manuscript

    my $mfe=shift;

    #numerical value, minimum 1
    my $mfe_adj=max2(1,-$mfe);

    #parameters of known precursors and background hairpins, scale and location
    my $prob_test=prob_gumbel_discretized($mfe_adj,5.5,32);
    my $prob_background=prob_gumbel_discretized($mfe_adj,4.8,23);
    my $odds=();
    my $log_odds=();
    if (($prob_background == 0)or($prob_test == 0)){
    	$log_odds=5.4;
    	}
    if($prob_background != 0){
    	$odds=$prob_test/$prob_background;
    	$log_odds=log($odds);
    	} 

    return $log_odds;
}



sub prob_gumbel_discretized{

#   discretized Gumbel distribution, probabilities within windows of 1 kCal/mol
#   uses the subroutine that calculates the cdf to find the probabilities

    my ($var,$scale,$location)=@_;

    my $bound_lower=$var-0.5;
    my $bound_upper=$var+0.5;

    my $cdf_lower=cdf_gumbel($bound_lower,$scale,$location);
    my $cdf_upper=cdf_gumbel($bound_upper,$scale,$location);

    my $prob=$cdf_upper-$cdf_lower;

    return $prob;
}


sub cdf_gumbel{

#   calculates the cumulative distribution function of the Gumbel distribution

    my ($var,$scale,$location)=@_;

    my $cdf=$e**(-($e**(-($var-$location)/$scale)));

    return $cdf;
}

