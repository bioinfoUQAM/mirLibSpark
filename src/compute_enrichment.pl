#!/usr/bin/perl
#Auteur : Amine M. Remita
#Date : 12/07/2018

# Le programme crée un heatmap d'enrechissement des entites (GO, KEGG, etc.)


# Chao-Jung Wu, 
# update: 2018-09-30
# mask lines 88, 89 so that it does not plot heatmaps. 
# And lease from the the dependency of R and RColorBrewer, gplots. 
# I will have to find a way to plot the heatmaps, such as python seaborne.

use v5.10.1;
use strict;
use warnings;
use File::Path qw(make_path);
use File::Basename;

# get rid of "Smartmatch is experimental at.." warning
no if $] >= 5.017011, warnings => 'experimental::smartmatch';

my $USAGE = "
\n*******************************************************************************
perl compute_enrichment.pl [inputfile] [output_path] [1|0]

inputfile   : TSV file containing the background and the conditions files and their names
output_path : Directory of the output files
[1|0]       : 1 to compute the upper pvalues, 0 to compute lower pvalues
*******************************************************************************\n\n";

my $cmd_R = "Rscript --no-save";
#
#######################################################################
# Dimensions des figures
my $_width = "3000";  #3600 #pathway
my $_height = "3900"; #1600 #library
# Margins
my $_margin_row = "100";#"40" #pathway name margin
my $_margin_col = "40";#30 #library name margin
# Axis label font size
my $cexRow = "6"; #3.5 #library name
my $cexCol = "10"; #4   #pathway name
my $notecex = "6";		# cell note font size
# Legend
my $legend_size = "2.5"; #2.5
# Key
my $keysize = "0.5";
# Palettes
#https://moderndata.plot.ly/create-colorful-graphs-in-r-with-rcolorbrewer-and-plotly/
my $palette = "PuRd";   # Reds, Oranges, Purples, Blues, BuGn, PuRd
#pvalue threshold for heatmap colors
my $pv_threshold = "0.05";
#######################################################################


###########################
# Arguments
my $inputfile   = $ARGV[0];
my $outputdir   = $ARGV[1];
my $pvalue_type = $ARGV[2];          # 1 : upper pvalue, 0 : lower pvalue
die "$USAGE" if scalar @ARGV != 3;
###########################

# Add a trailing slash to the output path
$outputdir .= "/" if $outputdir !~ /.+\/$/;
# Verify if the path exists or not
unless (-e $outputdir){
  make_path($outputdir, {chmod => 0755}) or die "Erreur, impossible de creer $outputdir \n$USAGE";
}

# 
unless ($pvalue_type ~~ [0, 1]){die "non valid parameters for upper/lower pvalue : 1 or 0 !!\n$USAGE";}
my $pval_str = $pvalue_type == 1 ? "upper":"lower";

# Get the background and conditions' files
my $GOSFILE = getBackgroundFileName($inputfile);           # background
my ($tfiles, $tCases) = getConditionFileNames($inputfile);

my %goDesc = getAnnotationDescription($GOSFILE);

my $DISTFLE = $outputdir."enrichment_dist";
my $PVALFLE = $outputdir."enrichment_pval_".$pval_str;

my %unirefs = &getUnirefList(\@$tfiles,\@$tCases);

###############################################################################
# Compute frequencies and pvalues
&computeEnrichment (\$GOSFILE,\%unirefs, $tCases,\$DISTFLE);  # $tCases is a already a reference
&computePvalues ($GOSFILE, $DISTFLE, $PVALFLE, $pvalue_type);
###############################################################################

# Plot the heatmaps ###########################################################
#&createPvalueHeatmap ($DISTFLE,$PVALFLE);
#&createPvalueHeatmapTrunc ($DISTFLE,$PVALFLE);

## SVG ########################################################################
# &createPvalueHeatmapSVG ($DISTFLE,$PVALFLE);
# &createPvalueHeatmapTruncSVG ($DISTFLE,$PVALFLE);
###############################################################################

print STDOUT "Fin normal du programme";

#Fonction :
##################################
#
# PREMIERE PARTIE DU PROGRAMME
###############################################################################
###############################################################################
#
sub getAnnotationDescription{
  my %hash = ();

  open my $IN,$_[0] or die $!;
  while (my $line = <$IN>){
    chomp $line;
    my @splt = split /\t/, $line;
    $hash{$splt[1]} = $splt[2];
  }
  close $IN;
  return %hash;
}
#
sub getUnirefList{
  my @caseFiles = @{$_[0]};
  my @cases = @{$_[1]};
  my %unis = ();
  
  for (my $i=0;$i<=($#cases);$i++){
    @{$unis{$cases[$i]}} = ();
    open my $inf , $caseFiles[$i] or die $caseFiles[$i];
    while (my $line = <$inf>){
      chomp $line;
      my $elem = join "\t", (split /\t/, $line)[0..1];
      push @{$unis{$cases[$i]}} , $elem;
    }
    close $inf;
  }
  return %unis;
}
#
sub computeEnrichment {                   # cette fonction est nommée &getEnrichment dans les scripts originaux
  my ($ref1,$ref2, $ref3, $ref4) = @_;
  my $bgFile = $$ref1;
  my %UnirefList = %$ref2;                # Lib -> @unirefs
  my @tabCases = @$ref3;                  # Liste of sorted conditions
  my $outF = $$ref4.".csv";
  
  my %goSlimDescp = ();                   # GO ---> Description
  my %goFreqClus = ();	                  # key ---> GOs  ---> frequency

  my ($elem, $goSlim, $goslim_desc,$uniref);
  print "Calcul de l'enrechissement \n";
  
  open IN, $bgFile or die $!;

  # ATGTTTTCGTCTAGCAATATA	GO:0048544	recognition of pollen
  while (my $line = <IN>){
    chomp $line;
    my @splt = split /\t/, $line;
    $elem = $splt[0] ."\t". $splt[1];
    $goSlim = $splt[1];
    $goslim_desc = $splt[1];

    unless (defined $goSlimDescp{$goSlim}){     # je recupère la description du goSlim
      $goSlimDescp{$goSlim} = $goslim_desc;
    }
    
    # cuillette des données
    foreach my $case (keys %UnirefList){
      if ($elem ~~ @{$UnirefList{$case}}){      # si l'uniref est ciblé par le groupe de mirnas exprimé dans cette condition
        ${$goFreqClus{$case}}{$goSlim} ++;
      }
    }
  }
  close IN;
  
  my @slimGroups = sort {$goSlimDescp{$a} cmp $goSlimDescp{$b}} keys %goSlimDescp;
  
  # impression des données
  print "Impression des donnees\n";
  open (DIST, ">".$outF) or die $!;
  
  # impression de l'entête
  print DIST "Clusters;";
  foreach my $slim (@slimGroups){
    print DIST $goSlimDescp{$slim};
    if ($slim ne $slimGroups[$#slimGroups]){print DIST ";";}
  }
  print DIST "\n";
  #impression de la matrice
  foreach my $case (@tabCases){
    print DIST $case.";"; 
    foreach my $slim (@slimGroups){
      if (defined ${$goFreqClus{$case}}{$slim}){
        print DIST ${$goFreqClus{$case}}{$slim};
      }else{
        print DIST "0";
      }
      if ($slim ne $slimGroups[$#slimGroups]){print DIST ";";}
    }
    print DIST "\n";
  }
  # print DIST "\n";
  close DIST;
}
# DEUXIEME PARTIE DU PROGRAMME
###############################################################################
###############################################################################
sub computePvalues {
  my $bgFile = $_[0];
  my $distF = $_[1].".csv";
  my $pvalF = $_[2].".csv";
  my $pval_type = $_[3];
  my %pvalues = ();
  
  print "Calcul des pvalues \n";
  
  # ************** One background for all conditions  *************************************************
  # ***************************************************************************************************
  my $N = &getPopSize ($bgFile);											# Population size of C : 180659, Population size of F : 201651, Population size of P : 370828
  # # # print $N."\n\n";
  my %bipartitions = &getBipartSize ($bgFile);							# bipartition  ---> size
  my $nbGroup = scalar keys %bipartitions;
  # foreach my $bipar (keys %bipartitions) {print $bipar."\t".$bipartitions{$bipar}."\n"};
  # ***************************************************************************************************
  
  # #Calcul de la frequences d'enrechissement pour chaque cluster pour chaque bpg, à partir du fichier $distF
  my %freqClusters = &getFreqClusters ($distF);							# cluster ---> @ freq
  
  # foreach my $clust (sort keys %freqClusters){
    # print $clust."\t".$freqClusters{$clust}[$nbGroup]."\n";
  # }
  
  # recuperer les etiquette dans le meme ordre:
  my @caseClusters = &getCaseClsuters ($distF);
  
  # Calcul de la pvalue
  # hypergeom( 300, 700, 100, 40 );
  
  foreach my $clust (@caseClusters){									#** il faut respecter l'ordre **#
    for (my $i = 0; $i <= ($nbGroup - 1); $i ++){
      # $m est le nombre de boules tirées désirées
      # $n est le nombre totale de boules tirées
      # $M est le nombre de boules désirés dans l'espace
      # $N est le nombre totale des boules dans l'espace
      my $m = ${$freqClusters{$clust}}[$i];
      my $n = ${$freqClusters{$clust}}[$nbGroup];
      my $M = $bipartitions{(${$freqClusters{"Clusters"}}[$i])};			# (${$freqClusters{"Clusters"}}[$i]) ==> le nom de la bipartition
      # print $M."\n";
      # my $pvalue =  hypergeom( 300, 700, 100, 40 );
      # my $pvalue =  hypergeom( $M, ($N-$M), $n, $m );
      my $pvalue;
      if($pval_type == 1){
        $pvalue =  &pvalue_upper( $M, ($N-$M), $n, $m );
      }else{
        $pvalue =  &pvalue_lower( $M, ($N-$M), $n, $m );
      }
      ${$pvalues{$clust}}[$i] = &roundFunction ($pvalue);
      # ${$pvalues{$clust}}[$i] = $pvalue;
    }
  }
  # # impression
  open (OUT, ">".$pvalF) or die $!;
  print OUT "Clusters;";

  for (my $j = 0 ; $j<= ($nbGroup - 1);$j++){
    unless (defined $goDesc{${$freqClusters{"Clusters"}}[$j]}){
      print (${$freqClusters{"Clusters"}}[$j]);
      print "\n";
    }
    
    print OUT $goDesc{${$freqClusters{"Clusters"}}[$j]};
    if ($j != ($nbGroup - 1)){print OUT ";"}
  }
  print OUT "\n";
  
  foreach my $clust (@caseClusters){
    print OUT $clust.";";
    # foreach my $pval (@{$pvalues{$clust}}){print OUT $pval; if ($pval ne ${$pvalues{$clust}}[$#{$pvalues{$clust}}]){print OUT ","}}
    for (my $j = 0 ; $j<= $#{$pvalues{$clust}};$j++){
      print OUT ${$pvalues{$clust}}[$j]; 
      if ($j != $#{$pvalues{$clust}}){print OUT ";"}
    }
    print OUT "\n";
  }
  
  close OUT;
}
#
sub getPopSize {
  my $size = 0;
  
  open (IN,$_[0]) or die $!;
  while (my $line = <IN>){
    $size ++;
  }
  close IN;
  print "Population size : $size\n";
  return $size;
}
# # C	GO:0000228	nuclear chromosome	O15392	GO:0000228	nuclear chromosome
# A0A075B6K6	GO:0002250	P
sub getBipartSize {
  my %hash = ();
  my $goSlim;
  
  open (IN,$_[0]) or die $!;
  while (my $line = <IN>){
    $goSlim = (split /\t/ ,$line)[1];
    $hash{$goSlim}++;
  }
  close IN;
  return %hash;
}
#
sub getCaseClsuters {
  my @tab = ();
  open (IN, $_[0]) or die $!;
  my $line = <IN>;
  while ($line = <IN>){
    chomp $line;
    my @splt = split (/;/, $line);
    push (@tab,$splt[0]);
  }
  close IN;
  return @tab;
}
#
sub getFreqClusters {
  my %hash = ();
  open (IN, $_[0]) or die $!;
  
  while (my $line = <IN>){
    chomp $line;
    my @splt = split (/;/, $line);
    #Clusters;reproduction;developmental maturation;DNA metabolic process;mitotic nuclear division;biosynthetic process;autophagy;neurological system process;protein folding;protein targeting;homeostatic process;cell differentiation;signal transduction;cellular amino acid metabolic process;generation of precursor metabolites and energy;cellular component assembly;macromolecular complex assembly;mitochondrion organization;cell adhesion;cellular nitrogen compound metabolic process;nucleobase-containing compound catabolic process;extracellular matrix organization;embryo development;cytoskeleton organization;translation;anatomical structure formation involved in morphogenesis;biological_process;protein complex assembly;lipid metabolic process;chromosome segregation;response to stress;small molecule metabolic process;vesicle-mediated transport;symbiosis, encompassing mutualism through parasitism;immune system process;cell proliferation;cell motility;plasma membrane organization;carbohydrate metabolic process;transport;mRNA processing;cell morphogenesis;transmembrane transport;chromosome organization;cell junction organization;cell-cell signaling;cell death;membrane organization;aging;cellular protein modification process;nucleocytoplasmic transport;cell division;locomotion;protein maturation;circulatory system process;anatomical structure development;growth;catabolic process;cell cycle
    #Conserved_gene;2;0;2;2;48;2;46;1;1;2;27;78;2;0;3;1;2;0;52;8;0;4;4;1;3;334;3;26;1;21;18;8;1;24;21;2;1;0;34;0;2;8;4;0;2;20;3;0;16;2;2;1;2;0;21;2;12;7
    
    my $total = 0;
    for (my $i=1;$i<=$#splt;$i++){
      push (@{$hash{$splt[0]}}, $splt[$i]);							# ajouter la fréquence pour chaque slim
      if ($splt[0] ne "Clusters"){$total += $splt[$i];}				# eliminer la premiere ligne dans le calcul (legende)
    }
    push (@{$hash{$splt[0]}}, $total);									# ajouter la frequence = tout le cluster
  }
  close IN;
  return %hash;
}
# http://hosho.ees.hokudai.ac.jp/~kubo/Rdoc/library/gplots/html/heatmap.2.html
sub createPvalueHeatmap {
  print "in createPvalueHeatmap\n";
  my $DF = $_[0].".csv";
  my $PF = $_[1].".csv";
  my $img = $_[1].".png";
  # creation du fichier script R temporaire
  my $scrFile = $outputdir."temp.r";
  
  my %bipartitions = &getBipartSize ($GOSFILE);								# bipartition  ---> size
  my $nbGroup = scalar keys %bipartitions;
  my $nbPlusOne = $nbGroup + 1;
  
  open (RS, ">".$scrFile) or die $!;
  print RS "library (\"RColorBrewer\")\n",
      "library (\"gplots\")\n",
      
      "fr <- read.csv(\"$DF\",sep=\";\")\n",	 # ,header = TRUE
      "pv <- read.csv(\"$PF\",sep=\";\")\n",
      
      "row.names(fr) <- fr\$Clusters\n",
      "row.names(pv) <- pv\$Clusters\n",
      "fr <- fr[,2:$nbPlusOne]\n",
      "pv <- pv[,2:$nbPlusOne]\n",
      
      "fr_matrix <- data.matrix(fr)\n",
      "pv_matrix <- data.matrix(pv)\n",
      "png(\"$img\", $_width, $_height)\n",
      
      # "myCol <- rev (c(colorRampPalette(brewer.pal(5,\"$palette\"))(4)))\n",
      "myCol <- rev (c(colorRampPalette(brewer.pal(3,\"$palette\"))(2)))\n",
      
			# "myBreaks <- c(0, 0.00001, 0.001, 0.005 , 1)\n",
      #"myBreaks <- c(0, 0.00005, 0.0005, 0.005, 0.05, 1)\n",
      "myBreaks <- c(0, $pv_threshold, 1)\n",

      
      "fr_heatmap <- heatmap.2(pv_matrix, Rowv=NA, Colv=NA, dendrogram=\'none\', col = myCol, scale=\"none\", cexRow=$cexRow, cexCol= $cexCol, margins=c($_margin_row,$_margin_col), key=TRUE ,trace=\"none\", symkey=FALSE, denscol=\"red\", srtCol=45, density.info=\"density\", keysize=$keysize ,rowsep=1:dim(fr_matrix)[1], colsep=1:$nbPlusOne,  sepcolor=\"black\", sepwidth=c(0.008,0.008), cellnote = fr_matrix , notecex=$notecex, notecol=\"black\", breaks = myBreaks)\n", 
  
      # "legend(.0001,0.94, fill = myCol, legend = c(\"0 - 10e-5\", \"10e-5 - 10e-3\", \"10e-3 - .005\", \".005 - 1\"))\n",
      # "legend(.0001,0.94, fill = myCol, legend = c(\"0 - 5e-5\", \"5e-5 - 5e-4\", \"5e-4 - 5e-3\", \"5e-3 - .05\", \".05 - 1\"))\n",
      "legend(.0001,0.80, fill = myCol, legend = c(\"0 - $pv_threshold\", \"$pv_threshold - 1\"), cex=$legend_size)\n",
      
      "";
      
  close RS;
  # lancer le script R
  my $cmd = $cmd_R." ".$scrFile;
  system ($cmd);
  
  # unlink $scrFile;
}
#
sub createPvalueHeatmapTrunc {
  my $DF = $_[0].".csv";
  my $PF = $_[1].".csv";
  my $img = $_[1]."_trunc.png";
  
  # creation du fichier script R temporaire
  my $scrFile = $outputdir."temp.r";
  my %bipartitions = &getBipartSize ($GOSFILE);								# bipartition  ---> size
  my $nbGroup = scalar keys %bipartitions;
  my $nbPlusOne = $nbGroup + 1;
  
  open (RS, ">".$scrFile) or die $!;
  print RS "library (\"RColorBrewer\")\n",
      "library (\"gplots\")\n",
      
      "infFunc <- function (x){\n",
# 			"  return (x>0.005)\n",
      "  return (x>0.05)\n",
      "}\n",
      
      "fr <- read.csv(\"$DF\",sep=\";\")\n",	 # ,header = TRUE
      "pv <- read.csv(\"$PF\",sep=\";\")\n",
      
      "row.names(fr) <- fr\$Clusters\n",
      "row.names(pv) <- pv\$Clusters\n",
      "fr <- fr[,2:$nbPlusOne]\n",
      "pv <- pv[,2:$nbPlusOne]\n",
      
      "fr_matrix <- data.matrix(fr)\n",
      "pv_matrix <- data.matrix(pv)\n",
      
      "pv1 <- pv_matrix[,!apply(infFunc(pv_matrix), 2, all)]\n",
      "fr1 <- fr_matrix[,!apply(infFunc(pv_matrix), 2, all)]\n",
      
      "pv2 <- pv1[!apply(infFunc(pv_matrix), 1, all),]\n",
      "fr2 <- fr1[!apply(infFunc(pv_matrix), 1, all),]\n",
      
      "png(\"$img\", $_width, $_height)\n",
      
      # "myCol <- rev (c(colorRampPalette(brewer.pal(5,\"$palette\"))(4)))\n",
      "myCol <- rev (c(colorRampPalette(brewer.pal(3,\"$palette\"))(2)))\n",
      
			# "myBreaks <- c(0, 0.00001, 0.001, 0.005 , 1)\n",
      # "myBreaks <- c(0, 0.00005, 0.0005, 0.005, 0.05, 1)\n",
      "myBreaks <- c(0, $pv_threshold, 1)\n",
      
      "fr_heatmap <- heatmap.2(pv2, Rowv=NA, Colv=NA, dendrogram=\'none\', col = myCol, scale=\"none\", cexRow=$cexRow, cexCol= $cexCol, margins=c($_margin_row,$_margin_col), key=TRUE ,trace=\"none\", srtCol=90, symkey=FALSE, denscol=\"red\", density.info=\"density\", keysize=$keysize ,rowsep=1:dim(fr2)[1], colsep=1:$nbPlusOne,  sepcolor=\"black\", sepwidth=c(0.008,0.008), cellnote = fr2 , notecex=$notecex, notecol=\"black\",breaks = myBreaks)\n",
      
      # "legend(.0001,0.80, fill = myCol, legend = c(\"0 - 10e-5\", \"10e-5 - 10e-3\", \"10e-3 - .005\", \".005 - 1\"))\n",
      # "legend(.0001,0.80, fill = myCol, legend = c(\"0 - 5e-5\", \"5e-5 - 5e-4\", \"5e-4 - 5e-3\", \"5e-3 - .05\", \".05 - 1\"))\n",
      "legend(.0001, 0.80, fill = myCol, legend = c(\"0 - $pv_threshold\", \"$pv_threshold - 1\"), cex=$legend_size)\n",
      
      "";
      
  close RS;
  # lancer le script R
  my $cmd = $cmd_R." ".$scrFile;
  system ($cmd);
  
  # unlink $scrFile;
}
#
sub createPvalueHeatmapSVG {
  my $DF = $_[0].".csv";
  my $PF = $_[1].".csv";
  my $img = $_[1].".svg";
  # creation du fichier script R temporaire
  my $scrFile = $outputdir."temp.r";
  
  my %bipartitions = &getBipartSize ($GOSFILE);								# bipartition  ---> size
  my $nbGroup = scalar keys %bipartitions;
  my $nbPlusOne = $nbGroup + 1;
  
  open (RS, ">".$scrFile) or die $!;
  print RS "library (\"RColorBrewer\")\n",
      "library (\"gplots\")\n",
      
      "fr <- read.csv(\"$DF\",sep=\";\")\n",	 # ,header = TRUE
      "pv <- read.csv(\"$PF\",sep=\";\")\n",
      
      "row.names(fr) <- fr\$Clusters\n",
      "row.names(pv) <- pv\$Clusters\n",
      "fr <- fr[,2:$nbPlusOne]\n",
      "pv <- pv[,2:$nbPlusOne]\n",
      
      "fr_matrix <- data.matrix(fr)\n",
      "pv_matrix <- data.matrix(pv)\n",
      
      # "png(\"$img\", $_width, $_height)\n",
      "svg(\"$img\", 16.5,12)\n",
      
      "myCol <- rev (c(colorRampPalette(brewer.pal(5,\"$palette\"))(4)))\n",
      "myBreaks <- c(0, 0.00001, 0.001, 0.005 , 1)\n",
      
      "fr_heatmap <- heatmap.2(pv_matrix, Rowv=NA, Colv=NA, dendrogram=\'none\', col = myCol, scale=\"none\", cexRow=1.2, cexCol= 1.2, margins=c($_margin_row,$_margin_col), key=TRUE ,trace=\"none\", symkey=FALSE, denscol=\"red\", density.info=\"density\", keysize=0.55 ,rowsep=1:dim(fr_matrix)[1], colsep=1:$nbPlusOne,  sepcolor=\"black\", sepwidth=c(0.008,0.008), cellnote = fr_matrix , notecex=$notecex, notecol=\"black\",breaks = myBreaks)\n",
      "legend(.000001,0.94, fill = myCol, legend = c(\"0 - 10e-5\", \"10e-5 - 10e-3\", \"10e-3 - .005\", \".005 - 1\"), cex = 0.65)\n",
      
      "";
      
  close RS;
  # lancer le script R
  my $cmd = $cmd_R." ".$scrFile;
  system ($cmd);
  
  # unlink $scrFile;
}
#
sub createPvalueHeatmapTruncSVG {
  my $DF = $_[0].".csv";
  my $PF = $_[1].".csv";
  my $img = $_[1]."_trunc.svg";
  # creation du fichier script R temporaire
  my $scrFile = $outputdir."temp.r";
  
  my %bipartitions = &getBipartSize ($GOSFILE);								# bipartition  ---> size
  my $nbGroup = scalar keys %bipartitions;
  my $nbPlusOne = $nbGroup + 1;
  
  open (RS, ">".$scrFile) or die $!;
  print RS "library (\"RColorBrewer\")\n",
      "library (\"gplots\")\n",
      
      "infFunc <- function (x){\n",
      "  return (x>0.005)\n",
      "}\n",
      
      "fr <- read.csv(\"$DF\",sep=\";\")\n",	 # ,header = TRUE
      "pv <- read.csv(\"$PF\",sep=\";\")\n",
      
      "row.names(fr) <- fr\$Clusters\n",
      "row.names(pv) <- pv\$Clusters\n",
      "fr <- fr[,2:$nbPlusOne]\n",
      "pv <- pv[,2:$nbPlusOne]\n",
      
      "fr_matrix <- data.matrix(fr)\n",
      "pv_matrix <- data.matrix(pv)\n",
      
      "pv1 <- pv_matrix[,!apply(infFunc(pv_matrix), 2, all)]\n",
      "fr1 <- fr_matrix[,!apply(infFunc(pv_matrix), 2, all)]\n",
      
      "pv2 <- pv1[!apply(infFunc(pv_matrix), 1, all),]\n",
      "fr2 <- fr1[!apply(infFunc(pv_matrix), 1, all),]\n",
      
      "svg(\"$img\", 16.5,12)\n",					# avec les deg mirnas
      
      "myCol <- rev (c(colorRampPalette(brewer.pal(5,\"$palette\"))(4)))\n",#   \"orange\" 
      "myBreaks <- c(0, 0.00001, 0.001, 0.005 , 1)\n",
      
      "fr_heatmap <- heatmap.2(pv2, Rowv=NA, Colv=NA, dendrogram=\'none\', col = myCol, scale=\"none\", cexRow=1.2, cexCol= 1.2, margins=c($_margin_row,$_margin_col), key=TRUE ,trace=\"none\", symkey=FALSE, denscol=\"red\", density.info=\"density\", keysize=0.55 ,rowsep=1:dim(fr2)[1], colsep=1:$nbPlusOne,  sepcolor=\"black\", sepwidth=c(0.008,0.008), cellnote = fr2 , notecex=$notecex, notecol=\"black\",breaks = myBreaks)\n",
      "legend(.000001,0.94, fill = myCol, legend = c(\"0 - 10e-5\", \"10e-5 - 10e-3\", \"10e-3 - .005\", \".005 - 1\"), cex = 0.65)\n",
      
      "";
      
  close RS;
  # lancer le script R
  my $cmd = $cmd_R." ".$scrFile;
  system ($cmd);
  
  # unlink $scrFile;
}
# print hypergeom( 300, 700, 100, 40 ); 
sub hypergeom {
  # There are m "bad" and n "good" balls in an urn.
  # Pick N of them. The probability of i or more successful selections:
  # (m!n!N!(m+n-N)!)/(i!(n-i)!(m+i-N)!(N-i)!(m+n)!)
  my ($n, $m, $N, $i) = @_;

  my $loghyp1 = logfact($m) +logfact($n)+logfact($N)+logfact($m+$n-$N);
  my $loghyp2 = logfact($i)+logfact($n-$i)+logfact($m+$i-$N)+logfact($N-$i)+logfact($m+$n);
  return exp($loghyp1 - $loghyp2);
}
#
sub logfact {
  my $x = shift;
  my $ser = (   1.000000000190015
        + 76.18009172947146   / ($x + 2)
        - 86.50532032941677   / ($x + 3)
        + 24.01409824083091   / ($x + 4)
        - 1.231739572450155   / ($x + 5)
        + 0.12086509738661e-2 / ($x + 6)
        - 0.5395239384953e-5  / ($x + 7) );
  my $tmp = $x + 6.5;
  ($x + 1.5) * log($tmp) - $tmp + log(2.5066282746310005 * $ser / ($x+1));
}
sub roundFunction {
  my $rounded = sprintf("%12e", $_[0]);
  # if ($rounded == "100.00"){$rounded = "100";}
  return $rounded;
}
# my $pvalue =  &pvalue_upper( $M, ($N-$M), $m, $n );
sub pvalue_upper {
  my ($n, $m, $N, $i) = @_;
  my $pval;
  
  for (my $j = $i; $j <= (&min($n,$N)); $j++){
    $pval += hypergeom ($n, $m, $N, $j);
  }
  
  return $pval;
}
#
sub pvalue_lower {
  my ($n, $m, $N, $i) = @_;
  my $pval=0;
  
  for (my $j = 1; $j <= $i; $j++){                 # test for under-representation
    $pval += hypergeom ($n, $m, $N, $j);
  }
  
  return $pval;
}
#
sub min {
  my ($a, $b) = @_;
  
  if ($a < $b){
    return $a;
  }else{
    return $b;
  }
}
#
sub getBackgroundFileName {
  my $infile = $_[0];
  my $name;
  my $dirname = dirname($infile) . "/";
  
  open my $IN, $infile or die $!;
  
  while (my $line = <$IN>){
    chomp $line;
    my @splt = split /\t/ ,$line;
    if ($splt[1]=~ /^background$/i){
      $name = $dirname.$splt[0];
      # print $name . "\n";
    }
  }
  close $IN;
  
  die "The program cannot find the name of the background file in $infile \n" if not defined $name;
  
  return $name;
}
#
sub getConditionFileNames {
  my $infile = $_[0];
  my @files = ();
  my @cases = ();
  my $dirname = dirname($infile) . "/";
  
  open my $IN, $infile or die $!;
  
  while (my $line = <$IN>){
    chomp $line;
    my @splt = split /\t/ ,$line;
    if ($splt[1] !~ /^background$/i){
      push @files, $dirname.$splt[0];
      # print $dirname.$splt[0] . "\n";
      push @cases, $splt[1];
    }
  }
  close $IN;
  
  return (\@files, \@cases);
}
__END__
