package miRcheck;
use Exporter();
@ISA = qw(Exporter);
@EXPORT = qw(&miR_check &RNAfold &get_subseq &hash_FASTA &max &min);


sub miR_check{
    #checks potential miRNAs against a set of requirements
    #returns orientation of miRNA for those that pass all requirments
    #returns first requirment failed if not
    my ($fold,$beg,$end,$param) = @_;
    my (@pairs) =  map_pairs($fold);
    my (@mir,@long_mir,@mir_star,@long_mir_star,$nonpair,$orientation,$x,$a,$b,$mir_check) = ();
    my($mir_star,$star_start,$star_stop) = ();
    my($fback_start,$fback_stop) = (0,0);

    if ($param eq "def"){
        # #parameters for miRNA checking (default)
        $MAX_UNPAIR = 6;          # max # unpaired bases in putative mir
        $MAX_STAR_UNPAIR = 6;     # max # unpaired bases in putative mir*
        $MAX_SIZEDIFFERENCE = 3;  # max size difference between mir and mir*
        $MAX_MIR_GAP = 3;         # longest acceptable run of unpaired bases in mir
        $MAX_STAR_GAP = 3;        # longest acceptable run of unpaired bases in mir*
        $MIN_FBACK_SIZE = 54;     # shortest acceptable length of folback including mir and mir*
        $MAX_MIR_AS_BULGE = 2;    # maximum total # assymetrically unpaired bases in mir
        $MIN_UNPAIR = 1;          # minimum number of unpair bases in acceptable extended mirs/mir*s
        $BP_EXTENSION = 2;        # number of nt to extend mirs and mir*s
        $PRINT = 0;               # print diagnostic output to screen
    }elsif($param eq "mey"){
        # #parameters for miRNA checking (selon Jones-Rhoades et david bartel 2004, Meyers et al 2008) 
        $MAX_UNPAIR = 4;          # max # unpaired bases in putative mir									#6 par defaut
        $MAX_STAR_UNPAIR = 5;     # max # unpaired bases in putative mir*							#6 par defaut
        $MAX_SIZEDIFFERENCE = 3;  # max size difference between mir and mir*
        $MAX_MIR_GAP = 2;         # longest acceptable run of unpaired bases in mir						#3 par defaut
        $MAX_STAR_GAP = 3;        # longest acceptable run of unpaired bases in mir*
        $MIN_FBACK_SIZE = 60;     # shortest acceptable length of folback including mir and mir*	#54 par defaut
        $MAX_MIR_AS_BULGE = 2;    # maximum total # assymetrically unpaired bases in mir			#2 par defaut
        $MIN_UNPAIR = 1;          # minimum number of unpair bases in acceptable extended mirs/mir*s
        $BP_EXTENSION = 2;        # number of nt to extend mirs and mir*s
        $PRINT = 0;               # print diagnostic output to screen    
    }elsif($param eq "m18"){
        # #parameters for miRNA checking (Meyers et al 2018) 
        $MAX_UNPAIR = 4;          # max # unpaired bases in putative mir									#6 par defaut
        $MAX_STAR_UNPAIR = 5;     # max # unpaired bases in putative mir*							#6 par defaut
        $MAX_SIZEDIFFERENCE = 3;  # max size difference between mir and mir*
        $MAX_MIR_GAP = 2;         # longest acceptable run of unpaired bases in mir						#3 par defaut
        $MAX_STAR_GAP = 3;        # longest acceptable run of unpaired bases in mir*
        $MIN_FBACK_SIZE = 60;     # shortest acceptable length of folback including mir and mir*	#54 par defaut
        $MAX_MIR_AS_BULGE = 2;    # maximum total # assymetrically unpaired bases in mir			#2 par defaut
        $MIN_UNPAIR = 1;          # minimum number of unpair bases in acceptable extended mirs/mir*s
        $BP_EXTENSION = 2;        # number of nt to extend mirs and mir*s
        $PRINT = 0;               # print diagnostic output to screen    
    }else{
        die "Error!! you have to specify the parameters def or mey or m18\n";
    }

    while (@_)
    {
	$thisarg = shift @_;
	if ($thisarg eq "-unpair" ) {$MAX_UNPAIR=shift @_;} 
	if ($thisarg eq "-star_unpair" ) {$MAX_STAR_UNPAIR=shift @_;}
	if ($thisarg eq "-size_diff" ) {$MAX_SIZEDIFFERENCE=shift @_;}
	if ($thisarg eq "-mir_bulge") {$MAX_MIR_GAP=shift @_;}
	if ($thisarg eq "-star_bulge") {$MAX_STAR_GAP=shift @_;}
	if ($thisarg eq "-fback_min") {$MIN_FBACK_SIZE=shift @_;}
	if ($thisarg eq "-ass") {$MAX_MIR_AS_BULGE=shift @_;}
	if ($thisarg eq "-min_unpair") {$MIN_UNPAIR=shift @_;}
	if ($thisarg eq "-bp_ext") {$BP_EXTENSION=shift @_;}
	if ($thisarg eq "-print") {$PRINT =1;}
    }

    @mir = @pairs[$beg..$end];

    #check extent, pattern, and direction of pairing in miRNA
    $orientation = check_mir_pairing(\@mir,$beg);
    if (not($orientation =~ /prime/)){return($orientation,0,0);}

    #define miRNAstar, here meaning the sequence paired to miRNA in foldback
    ($mir_star,$star_start,$star_stop) = define_mir_star(\@mir,\@pairs,$orientation);

    #check size of foldback
    if ($orientation eq  "5prime"){
        $fbacksize = abs($beg - $star_start);
	$fback_start = $beg;
	$fback_stop = $star_start;
    }
    else {
	$fbacksize = abs($end-$star_stop);
	$fback_start = $star_stop;
        $fback_stop = $end;
    }
    if ($fbacksize < $MIN_FBACK_SIZE){return("FBACK_SHORT",0,0);}
   
    #check for assymmetric non-pairing in miRNA
    $mir_assymetry = count_assymetry(@mir);    
    if ($mir_assymetry > $MAX_MIR_AS_BULGE){return("MIR_ASSYMETRY",0,0);}

    #make extended mir and mir*s.  Check that one pair meets pairing requirements
    #the point of this is to ensure that the pairing extend in at least
    #one direction beyond the miRNA itself.
    $orientation = '';
    foreach $add_on (0..$BP_EXTENSION){
	if ($PRINT){ print "add on: $add_on\n";}
	if ($beg - $add_on < 0){next;} #check that miRNA is not too near begining
	if ($end+$BP_EXTENSION-$add_on > $#pairs){next;} #check that miRNA is not too near end
	@long_mir = @pairs[$beg - $add_on..$end+$BP_EXTENSION-$add_on];
	if ($PRINT){print "@long_mir\n";}
	$mir_check =  check_mir_pairing(\@long_mir,$beg - $add_on);
	if (not($mir_check =~ /prime/)){next;}
	if ($PRINT){print "mir pass\n";}

	#define long miR*
	($long_mir_star,$a,$b) = define_mir_star(\@long_mir,\@pairs,$mir_check);
	$star_check =  check_star_pairing($long_mir_star,\@long_mir);
	if ($PRINT){print "@$long_mir_star\n";}
	if (not($star_check =~ /good/)){next;}
	if ($PRINT){print "star pass\n";}
	#check that mir pairing is not too perfect
        if (not ((min_imperfection(@long_mir) >= $MIN_UNPAIR) or (min_imperfection(@$long_mir_star) >= $MIN_UNPAIR)) ){next;}
	if ($PRINT){print "imperfection pass\n";}
	$orientation = $mir_check; # this "long miR" and "long miR star" pair passed pairing tests
    }
    if (not($orientation)){return("FAILED_MIR_EXTENSION",0,0);}

    #at least one long mir passed all tests, return orientation
    return($orientation,$fback_start,$fback_stop);
}

sub min_imperfection{
    #returns # of unpaired bases
    my (@a) = @_;
    my ($a,$count) = ();
    $count = 0;
    foreach $a (@a){if ($a eq '-'){++$count;}}
    return($count);
}

sub count_assymetry{
    #returns the amount assymmetry in the mir/mir* pairing
    #assymetry = sum of (#unpaired mir - #unpaired star) 
    #over all non-pairing portions
    #note that this only discriminates against extra bases on one side of hairpin
    my (@mir) = @_;
    my ($ass,$temp,$in_bulge,$base,$last) = ();
    $ass =0;
    $in_bulge = 0;
    $last = 0;
    foreach $base (@mir){
	if (($in_bulge) and (not($base eq "-"))){#paired base next to unpaired base
	    if (not($last)){$in_bulge = 0;} #skip if first paired base
	    else{ 
		$diff = abs($base-$last); #find difference in partners of this base and last paired base
		$temp = $in_bulge - ($diff -1); #compare #unpaired bases to offset in pairs
		if ($temp > 0) {$ass += $temp;}
		$in_bulge = 0;
	    }
	}
	if ($base eq "-"){$in_bulge += 1;} # keep track of consecutive unpaired bases
	else {$last = $base;} #remember partner of last paired base
    }	    
    return($ass);
}

sub check_mir_pairing{
    my ($mir,$pos) = @_;
    my (@mir) = @{$mir};
    my ($nonpair,$orientation) = ();

    #check amount of pairing in miR
    $nonpair = 0;
    foreach $b (@mir){
        if ($b eq '-'){++$nonpair;}
    }
    if ($PRINT){ print "$nonpair of $MAX_UNPAIR nonpairs\n";}
    if ($nonpair > $MAX_UNPAIR){return("NONPAIRING");}

    #check direction of pairing in miR
    $orientation = '';
    my ($left_pairing,$right_pairing) = ();
    foreach $b (0 .. $#mir){
        if ($mir[$b] eq '-'){next;}
        if ($mir[$b] > $pos + $b) {++$left_pairing;}
        if ($mir[$b] < $pos + $b) {++$right_pairing;}
    }
    if (($left_pairing) and ($right_pairing)){return("MIXED_PAIRING");}
    elsif ($right_pairing){$orientation = "3prime";}
    elsif ($left_pairing){$orientation = "5prime";}

    #check longest unpaired gap
    if (biggest_gap(@mir) > $MAX_MIR_GAP){return("MIR_BULGE");}
    return($orientation);
}

sub check_star_pairing{
    my ($star,$mir) = @_;
    my (@star) = @{$star};
    my (@mir) = @{$mir};
    my ($size_difference,$nonpair, $b) = ();
    $size_difference = @star - @mir;
    if ($size_difference >= $MAX_SIZEDIFFERENCE){return("STAR_TOO_LONG");}
    #check pairing in miRNA*
    $nonpair = 0;
    foreach $b (@star){
        if ($b eq '-'){++$nonpair;}
    }
    if ($nonpair > $MAX_STAR_UNPAIR){return("STAR_NONPAIRING");}
    #check consecutive unpaired in mir*
    if (biggest_gap(@star) > $MAX_STAR_GAP){return("STAR_BULGE");}
    return("good");
}

sub define_mir_star{
    #defines the stretch of sequence pairing to the miRNA, counting equal length
    #unpaired sequence at the ends
    my ($mir,$pairs,$orientation) = @_;
    my (@mir) = @{$mir};
    my (@pairs) = @{$pairs};
    my ($star_start,$star_stop,$star_length,$x,@miRNA_star) = ();
    $x = 0;
    until ($star_start){
        if (not($mir[$x] eq "-")){
            if ($orientation eq  "5prime"){
                $star_start = $mir[$x]+$x;
            }
            else {$star_start = $mir[$x]+$x;
	      }
        }
        ++$x;
    }
    $x = $#mir;
    until ($star_stop){
        if (not($mir[$x] eq '-')){
            if ($orientation eq  "5prime"){
                $star_stop = $mir[$x] - ($#mir-$x);
            }
            else {$star_stop = $mir[$x] -  ($#mir-$x);
	      }
        }
        --$x;
    }
    $star_stop = max(0,$star_stop);
    $star_start = min($star_start,$#pairs);
    @miRNA_star = @pairs[$star_stop..$star_start];
    return(\@miRNA_star,$star_start,$star_stop);
}

sub biggest_gap{
    #returns largest run of unpaired bases
    my (@pairs) = @_;
    my ($p,$tmp,$high);
    $tmp = 0;
    $high = 0;
    foreach $p (@pairs){
	if ($p eq '-'){
	    ++$tmp;
	    if ($tmp > $high){$high = $tmp;}
	}
	else {$tmp =0;}
    }
    return($high);
}

sub RNAfold{
    #subroutine for handling RNAfold from with a perl script.
    #requires RNAfold to be installed
    #returns stucture in parentheses format and predicted energy of structure
    my ($name,$seq,$ID) = @_;
    if (not($ID)){$ID = '1';}
    my ($fold,$energy,@lines) = ();
    open(TMP,">RNASEQ.tmp.$ID");
    print TMP ">$name\n$seq\n";
    close(TMP);
    system("touch RNAFOLD.tmp.$ID; rm RNAFOLD.tmp.$ID");
    system("/home2/mjonesrh/exe/linuxexe/RNAfold_original < RNASEQ.tmp.$ID >RNAFOLD.tmp.$ID");
    open(TMP,"RNAFOLD.tmp.$ID");
    @lines = <TMP>;
    ($fold,$energy) = split(" ",$lines[2]);
    $energy =~ s/\(//g;  $energy =~ s/\)//g;
    return($fold,$energy);
}    

sub get_subseq{
    my($rseq,$x,$y) = @_;
    my($subseq) = ();
    if ($y >= $x){
        $beg = max(1,$x);
        $end = min($y,length($rseq));
        $subseq = substr($rseq,$beg-1,$end-$beg+1);
    }
    else {
        $beg = min(length($rseq),$x);
        $end = max($y,1);

        $subseq = substr($rseq,$end-1,$beg-$end+1);
        $subseq = revcomp($subseq);
    }
    return($subseq);
}

sub map_pairs{
    #input: RNAfold structure
    #output: array containing coordinate of BP, dash if unpaired
    my ($structure) = @_;
    my (@s,@left,@pairs,$x,$pairing) = ();
    @s = split("",$structure);
    foreach $x (0..$#s ){
        $pairs[$x] = '-';
    }
    foreach $x (0..$#s ){
        if ($s[$x] eq "\("){push(@left,$x);}
        if ($s[$x] eq "\)"){
            $pairing = pop(@left);
            $pairs[$x] = $pairing;
            $pairs[$pairing] = $x;
        }
    }
    return(@pairs);
}

sub min {  
    #returns the minimum value of an array
    my (@a) = @_;
    my ($min) = ($a[0]);
    foreach $a (@a){
        if ($a < $min){$min = $a;}
    }
    return($min);
}

sub max {
    #returns the maximum value of an array
    my (@a) = @_;
    my ($max) = ($a[0]);
    foreach $a (@a){
        if ($a > $max){$max = $a;}
    }
    return($max);
}

sub hash_FASTA{
    my($file) = @_;
    my(%hash) = ();
    open(HASHFILE,$file);
    while (<HASHFILE>) {
        if (/^\s+$/){next;}
        elsif (/>(\S+)/){$name = $1;}
        elsif (/(\S+)/){
            if (exists($hash{$name})) { $hash{$name} .= $1;}
            else {  $hash{$name} = $1;}
#           print "$name\n$hash{$name}\n";
        }
    }
    return(\%hash);
}

sub revcomp{
    my ($seq) = @_;
    $seq = reverse($seq);
    $seq =~ tr/Uu/Tt/;
    $seq =~ tr/ACGTacgt/TGCAtgca/;
    return($seq);
}
