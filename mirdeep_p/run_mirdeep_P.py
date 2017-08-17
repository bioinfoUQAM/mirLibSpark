'''
time python run_mirdeep_P.py
'''
import os

options = ['tair9', 'tair10']


def init_mirdeep_p (op):
  os.system('cp miRDP1.3/* .')
  #os.system('mkdir bowtie-index')
  if op == options[0]:
    genome = 'TAIR9_genome.fa'
    rep = 'tair9/'
    f_annotated = 'annotated_miRNA_extended.fa'
    #os.system('bowtie-build -f ' + genome + ' bowtie-index/' + genome[:-3] + ' >/dev/null') 
    os.system('cp ' + rep + 'ath.gff .')#= miRBase v15
    os.system('cp ' + rep + 'annotated_miRNA_extended.fa .')
    os.system('cp genome/' + genome + ' .')
    os.system('cp ' + rep + 'ncRNA_CDS.gff .')#
    os.system('cp ' + rep + 'chromosome_length .')#
    #os.system('perl fetch_extended_precursors.pl ' + genome + ' ath.gff >annotated_miRNA_extended.fa')
    #os.system('bowtie-build -f annotated_miRNA_extended.fa bowtie-index/annotated_miRNA_extended.fa >/dev/null')
  elif op == options[1]:
    genome = 'a_thaliana_t10.fa'
    rep = 'tair10/'
    f_annotated = 'annotated_miRNA_v21_extended.fa'
    #os.system('bowtie-build -f ' + genome + ' bowtie-index/' + genome[:-3] + ' >/dev/null') 
    os.system('cp ' + rep + 'ath.gff3.edited .')#= miRBase v21
    #os.system('cp ' + rep + 'annotated_miRNA_v21_extended.fa .')
    os.system('cp genome/' + genome + ' .')
    os.system('cp ' + rep + 'ncRNA_CDS.gff .')#
    os.system('cp ' + rep + 'chromosome_length .')#
    os.system('perl fetch_extended_precursors.pl ' + genome + ' ath.gff3.edited >annotated_miRNA_v21_extended.fa')
    os.system('bowtie-build -f annotated_miRNA_v21_extended.fa bowtie-index/annotated_miRNA_v21_extended.fa >/dev/null')
  return genome, f_annotated

def run_mirdp_new (infile):
  os.system('bowtie -a -v 0 bowtie-index/' + genome[:-3] + ' -f ' + infile + '>indata.aln')
  os.system('perl convert_bowtie_to_blast.pl indata.aln ' + infile + ' ' + genome + ' >indata.bst')
  os.system('perl filter_alignments.pl indata.bst -c 15 >indata_filter15.bst')
  os.system('perl overlap.pl indata_filter15.bst ncRNA_CDS.gff -b >indata_id_overlap_ncRNA_CDS')
  os.system('perl alignedselected.pl indata_filter15.bst -g indata_id_overlap_ncRNA_CDS >indata_filter15_ncRNA_CDS.bst')
  os.system('perl filter_alignments.pl indata_filter15_ncRNA_CDS.bst -b ' + infile + ' > indata_filtered.fa')
  os.system('perl excise_candidate.pl ' + genome + ' indata_filter15_ncRNA_CDS.bst 250 >indata_precursors.fa')
  os.system('cat indata_precursors.fa | RNAfold --noPS > indata_structures')
  os.system('bowtie-build -f indata_precursors.fa bowtie-index/indata_precursors >/dev/null')
  os.system('bowtie -a -v 0 bowtie-index/indata_precursors -f indata_filtered.fa > indata_precursors.aln')
  os.system('perl convert_bowtie_to_blast.pl indata_precursors.aln indata_filtered.fa indata_precursors.fa >indata_precursors.bst')
  os.system('sort +3 -25 indata_precursors.bst >indata_signatures')
  print('begin miRDP')
  os.system('perl miRDP.pl indata_signatures indata_structures > indata_predictions')
  #os.system('perl rm_redundant_meet_plant.pl chromosome_length indata_precursors.fa indata_predictions indata_nr_prediction indata_filter_P_prediction')

def run_mirdp_known (infile, f_annotated):
  os.system('bowtie -a -v 0 bowtie-index/' + f_annotated + ' -f ' + infile + ' >indata.aln')
  os.system('perl convert_bowtie_to_blast.pl indata.aln ' + infile + ' ' + f_annotated + ' > indata_extended.bst')
  os.system('perl excise_candidate.pl ' + f_annotated + ' indata_extended.bst 250 >precursors_250.fa')
  os.system('bowtie-build -f precursors_250.fa bowtie-index/precursors_250 >/dev/null')
  os.system('cat precursors_250.fa|RNAfold --noPS >precursors_250_structure')
  os.system('bowtie -a -v 0 bowtie-index/precursors_250 -f ' + infile + ' >indata_250.aln')
  os.system('perl convert_bowtie_to_blast.pl indata_250.aln ' + infile + ' precursors_250.fa >indata_250.bst')
  os.system('sort +3 -25 indata_250.bst >indata_250_signature')
  os.system('perl miRDP.pl indata_250_signature precursors_250_structure >indata_250_prediction')

infile = 'high_conf_mature_ath_uniq_collapsed.fa'
#infile = '100_collapsed.fa'
os.system('cp input_storage/' + infile + ' .')
genome, f_annotated = init_mirdeep_p ('tair9')
#run_mirdp_new (infile) #= takes 15 secs, validates 75 unique mirna
run_mirdp_known (infile, f_annotated) #= also takes 15 secs, validates 75 unique mirna

#print('mode of predicting new: ')
#os.system('cut -f1 indata_predictions | grep \'seq_\' | sort | uniq | wc -l')
print('mode of predicting known: ')
os.system('cut -f1 indata_250_prediction | grep \'seq_\' | sort | uniq | wc -l')

