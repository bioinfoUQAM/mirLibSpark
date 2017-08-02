'''
$ bowtie-build -f a_thaliana_t10.fa bowtie-index/TAIR10_genome
$ bowtie -a -v 0 bowtie-index/TAIR10_genome -f high_conf_mature_ath_uniq_collapsed.fa >high_conf_mature_ath_uniq_collapsed.aln
$ perl convert_bowtie_to_blast.pl high_conf_mature_ath_uniq_collapsed.aln high_conf_mature_ath_uniq_collapsed.fa a_thaliana_t10.fa >high_conf_mature_ath_uniq_collapsed.bst
$ perl filter_alignments.pl high_conf_mature_ath_uniq_collapsed.bst -c 15 >high_conf_mature_ath_uniq_collapsed_filter15.bst
$ cp ../New_M/ncRNA_CDS.gff .
$ perl overlap.pl high_conf_mature_ath_uniq_collapsed_filter15.bst ncRNA_CDS.gff -b >id_overlap_ncRNA_CDS
	stdout: high_conf_mature_ath_uniq_collapsed_filter15.bst 154
			ncRNA_CDS.gff 181359
$ perl alignedselected.pl high_conf_mature_ath_uniq_collapsed_filter15.bst -g id_overlap_ncRNA_CDS >high_conf_mature_ath_uniq_collapsed_filter15_ncRNA_CDS.bst
$ perl filter_alignments.pl high_conf_mature_ath_uniq_collapsed_filter15_ncRNA_CDS.bst -b high_conf_mature_ath_uniq_collapsed.fa > high_conf_mature_ath_uniq_collapsed_filtered.fa
$ perl excise_candidate.pl a_thaliana_t10.fa high_conf_mature_ath_uniq_collapsed_filter15_ncRNA_CDS.bst 250 >high_conf_mature_ath_uniq_collapsed_precursors.fa
$ cat high_conf_mature_ath_uniq_collapsed_precursors.fa | RNAfold --noPS > high_conf_mature_ath_uniq_collapsed_structures
$ bowtie-build -f high_conf_mature_ath_uniq_collapsed_precursors.fa bowtie-index/high_conf_mature_ath_uniq_collapsed_precursors
$ bowtie -a -v 0 bowtie-index/high_conf_mature_ath_uniq_collapsed_precursors -f high_conf_mature_ath_uniq_collapsed_filtered.fa > high_conf_mature_ath_uniq_collapsed_precursors.aln
$ perl convert_bowtie_to_blast.pl high_conf_mature_ath_uniq_collapsed_precursors.aln high_conf_mature_ath_uniq_collapsed_filtered.fa high_conf_mature_ath_uniq_collapsed_precursors.fa >high_conf_mature_ath_uniq_collapsed_precursors.bst
$ sort +3 -25 high_conf_mature_ath_uniq_collapsed_precursors.bst >high_conf_mature_ath_uniq_collapsed_signatures
$ perl miRDP.pl high_conf_mature_ath_uniq_collapsed_signatures high_conf_mature_ath_uniq_collapsed_structures > high_conf_mature_ath_uniq_collapsed_predictions
$ cp ../New_M/chromosome_length .
$ perl rm_redundant_meet_plant.pl chromosome_length high_conf_mature_ath_uniq_collapsed_precursors.fa high_conf_mature_ath_uniq_collapsed_predictions high_conf_mature_ath_uniq_collapsed_nr_prediction high_conf_mature_ath_uniq_collapsed_filter_P_prediction
'''
import os

def init_tair10 ():
  #= [~/workspace/miRNA_predictor/mirDeep-P/JulieTest]
  os.system('cp ../miRDP1.3/* .')
  os.system('bowtie-build -f ../a_thaliana_t10.fa bowtie-index/TAIR10_genome')
  os.system('cp ../New_M/ncRNA_CDS.gff .')
  os.system('cp ../New_M/chromosome_length .')

def run_mirdp_new (infile):
  inBase = infile[:-3]
  os.system('bowtie -a -v 0 bowtie-index/TAIR10_genome -f ' + infile + '>' + inBase + '.aln')
  os.system('perl convert_bowtie_to_blast.pl ' + inBase + '.aln ' + inBase + '.fa ../a_thaliana_t10.fa >' + inBase + '.bst')
  os.system('perl filter_alignments.pl ' + inBase + '.bst -c 15 >' + inBase + '_filter15.bst')
  os.system('perl overlap.pl ' + inBase + '_filter15.bst ncRNA_CDS.gff -b >' + inBase + '_id_overlap_ncRNA_CDS')
  os.system('perl alignedselected.pl ' + inBase + '_filter15.bst -g ' + inBase + '_id_overlap_ncRNA_CDS >' + inBase + '_filter15_ncRNA_CDS.bst')
  os.system('perl filter_alignments.pl ' + inBase + '_filter15_ncRNA_CDS.bst -b ' + infile + ' > ' + inBase + '_filtered.fa')
  os.system('perl excise_candidate.pl ../a_thaliana_t10.fa ' + inBase + '_filter15_ncRNA_CDS.bst 250 >' + inBase + '_precursors.fa')
  os.system('cat ' + inBase + '_precursors.fa | RNAfold --noPS > ' + inBase + '_structures')
  os.system('bowtie-build -f ' + inBase + '_precursors.fa bowtie-index/' + inBase + '_precursors')
  os.system('bowtie -a -v 0 bowtie-index/' + inBase + '_precursors -f ' + inBase + '_filtered.fa > ' + inBase + '_precursors.aln')
  os.system('perl convert_bowtie_to_blast.pl ' + inBase + '_precursors.aln ' + inBase + '_filtered.fa ' + inBase + '_precursors.fa >' + inBase + '_precursors.bst')
  os.system('sort +3 -25 ' + inBase + '_precursors.bst >' + inBase + '_signatures')
  os.system('perl miRDP.pl ' + inBase + '_signatures ' + inBase + '_structures > ' + inBase + '_predictions')
  os.system('perl rm_redundant_meet_plant.pl chromosome_length ' + inBase + '_precursors.fa ' + inBase + '_predictions ' + inBase + '_nr_prediction ' + inBase + '_filter_P_prediction')

  
infile = 'high_conf_mature_ath_uniq_collapsed.fa'
run_mirdp_known (infile)