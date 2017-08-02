'''
cp /home/cloudera/workspace/mirLibHadoop/test/Run_mirdeep_P.py .

'''
import os

def init_bowtie ():
  #= [~/workspace/miRNA_predictor/mirDeep-P/JulieTest]
  os.system('cp ../miRDP1.3/* .')
  os.system('cp ../a_thaliana_t10.fa .')
  os.system('mkdir bowtie-index')
  #os.system('bowtie-build -f a_thaliana_t10.fa bowtie-index/TAIR10_genome >/dev/null')
  os.system('cp ../bowtie-index/* bowtie-index/')
  os.system('cp ../New_M/ncRNA_CDS.gff .')
  os.system('cp ../New_M/chromosome_length .')
  #
  os.system('cp ../Known_M/ath.gff .')
  os.system('cp ../Known_M/annotated_miRNA_extended.fa .')
  #os.system('bowtie-build -f annotated_miRNA_extended.fa bowtie-index/annotated_miRNA_extended >/dev/null')
  

def run_mirdp_new (infile):
  inBase = infile[:-3]
  os.system('bowtie -a -v 0 bowtie-index/TAIR10_genome -f ' + infile + '>' + inBase + '.aln')
  os.system('perl convert_bowtie_to_blast.pl ' + inBase + '.aln ' + inBase + '.fa a_thaliana_t10.fa >' + inBase + '.bst')
  os.system('perl filter_alignments.pl ' + inBase + '.bst -c 15 >' + inBase + '_filter15.bst')
  os.system('perl overlap.pl ' + inBase + '_filter15.bst ncRNA_CDS.gff -b >' + inBase + '_id_overlap_ncRNA_CDS')
  os.system('perl alignedselected.pl ' + inBase + '_filter15.bst -g ' + inBase + '_id_overlap_ncRNA_CDS >' + inBase + '_filter15_ncRNA_CDS.bst')
  os.system('perl filter_alignments.pl ' + inBase + '_filter15_ncRNA_CDS.bst -b ' + infile + ' > ' + inBase + '_filtered.fa')
  os.system('perl excise_candidate.pl a_thaliana_t10.fa ' + inBase + '_filter15_ncRNA_CDS.bst 250 >' + inBase + '_precursors.fa')
  os.system('cat ' + inBase + '_precursors.fa | RNAfold --noPS > ' + inBase + '_structures')
  os.system('bowtie-build -f ' + inBase + '_precursors.fa bowtie-index/' + inBase + '_precursors >/dev/null')
  os.system('bowtie -a -v 0 bowtie-index/' + inBase + '_precursors -f ' + inBase + '_filtered.fa > ' + inBase + '_precursors.aln')
  os.system('perl convert_bowtie_to_blast.pl ' + inBase + '_precursors.aln ' + inBase + '_filtered.fa ' + inBase + '_precursors.fa >' + inBase + '_precursors.bst')
  os.system('sort +3 -25 ' + inBase + '_precursors.bst >' + inBase + '_signatures')
  os.system('perl miRDP.pl ' + inBase + '_signatures ' + inBase + '_structures > ' + inBase + '_predictions')
  #os.system('perl rm_redundant_meet_plant.pl chromosome_length ' + inBase + '_precursors.fa ' + inBase + '_predictions ' + inBase + '_nr_prediction ' + inBase + '_filter_P_prediction')

def run_mirdp_known (infile):
  inBase = infile[:-3]
  os.system('bowtie -a -v 0 bowtie-index/annotated_miRNA_extended -f high_conf_mature_ath_uniq_collapsed.fa >high_conf_mature_ath_uniq_collapsed_extended.aln')
  os.system('perl convert_bowtie_to_blast.pl ' + inBase + '_extended.aln ' + inBase + '.fa annotated_miRNA_extended.fa  > ' + inBase + '_extended.bst')
  os.system('perl excise_candidate.pl annotated_miRNA_extended.fa ' + inBase + '_extended.bst 250 >precursors_250.fa')
  os.system('bowtie-build -f precursors_250.fa bowtie-index/precursors_250 >/dev/null')
  os.system('cat precursors_250.fa|RNAfold --noPS >precursors_250_structure')
  os.system('bowtie -a -v 0 bowtie-index/precursors_250 -f ' + inBase + '.fa >' + inBase + '_250.aln')
  os.system('perl convert_bowtie_to_blast.pl ' + inBase + '_250.aln ' + inBase + '.fa precursors_250.fa >' + inBase + '_250.bst')
  os.system('sort +3 -25 ' + inBase + '_250.bst >' + inBase + '_250_signature')
  os.system('perl miRDP.pl ' + inBase + '_250_signature precursors_250_structure >' + inBase + '_250_prediction')





os.system('cp ../high_conf_mature_ath_uniq_collapsed.fa .')
infile = 'high_conf_mature_ath_uniq_collapsed.fa'
init_bowtie ()
run_mirdp_new (infile) #= takes 15 secs, validates 75 unique mirna
run_mirdp_known (infile)

print('mode of predicting new: ')
os.system('cut -f1 high_conf_mature_ath_uniq_collapsed_predictions | grep \'seq_\' | sort | uniq | wc -l')
print('mode of predicting known: ')
os.system('cut -f1 high_conf_mature_ath_uniq_collapsed_250_prediction | grep \'seq_\' | sort | uniq | wc -l')
