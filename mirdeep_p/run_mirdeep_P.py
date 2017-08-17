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
    os.system('perl fetch_extended_precursors.pl ' + genome + ' ath.gff >annotated_miRNA_extended.fa')
    os.system('bowtie-build -f annotated_miRNA_extended.fa bowtie-index/annotated_miRNA_extended.fa >/dev/null')
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
  inBase = infile[:-3]
  os.system('bowtie -a -v 0 bowtie-index/' + genome[:-3] + ' -f ' + infile + '>' + inBase + '.aln')
  os.system('perl convert_bowtie_to_blast.pl ' + inBase + '.aln ' + inBase + '.fa ' + genome + ' >' + inBase + '.bst')
  os.system('perl filter_alignments.pl ' + inBase + '.bst -c 15 >' + inBase + '_filter15.bst')
  os.system('perl overlap.pl ' + inBase + '_filter15.bst ncRNA_CDS.gff -b >' + inBase + '_id_overlap_ncRNA_CDS')
  os.system('perl alignedselected.pl ' + inBase + '_filter15.bst -g ' + inBase + '_id_overlap_ncRNA_CDS >' + inBase + '_filter15_ncRNA_CDS.bst')
  os.system('perl filter_alignments.pl ' + inBase + '_filter15_ncRNA_CDS.bst -b ' + infile + ' > ' + inBase + '_filtered.fa')
  os.system('perl excise_candidate.pl ' + genome + ' ' + inBase + '_filter15_ncRNA_CDS.bst 250 >' + inBase + '_precursors.fa')
  os.system('cat ' + inBase + '_precursors.fa | RNAfold --noPS > ' + inBase + '_structures')
  os.system('bowtie-build -f ' + inBase + '_precursors.fa bowtie-index/' + inBase + '_precursors >/dev/null')
  os.system('bowtie -a -v 0 bowtie-index/' + inBase + '_precursors -f ' + inBase + '_filtered.fa > ' + inBase + '_precursors.aln')
  os.system('perl convert_bowtie_to_blast.pl ' + inBase + '_precursors.aln ' + inBase + '_filtered.fa ' + inBase + '_precursors.fa >' + inBase + '_precursors.bst')
  os.system('sort +3 -25 ' + inBase + '_precursors.bst >' + inBase + '_signatures')
  print('begin miRDP')
  os.system('perl miRDP.pl ' + inBase + '_signatures ' + inBase + '_structures > ' + inBase + '_predictions')
  #os.system('perl rm_redundant_meet_plant.pl chromosome_length ' + inBase + '_precursors.fa ' + inBase + '_predictions ' + inBase + '_nr_prediction ' + inBase + '_filter_P_prediction')

def run_mirdp_known (infile, f_annotated):
  inBase = infile[:-3]
  os.system('bowtie -a -v 0 bowtie-index/' + f_annotated + ' -f ' + infile + ' >' + inBase + '.aln')
  os.system('perl convert_bowtie_to_blast.pl ' + inBase + '.aln ' + inBase + '.fa ' + f_annotated + ' > ' + inBase + '_extended.bst')
  os.system('perl excise_candidate.pl ' + f_annotated + inBase + '_extended.bst 250 >precursors_250.fa')
  os.system('bowtie-build -f precursors_250.fa bowtie-index/precursors_250 >/dev/null')
  os.system('cat precursors_250.fa|RNAfold --noPS >precursors_250_structure')
  os.system('bowtie -a -v 0 bowtie-index/precursors_250 -f ' + inBase + '.fa >' + inBase + '_250.aln')
  os.system('perl convert_bowtie_to_blast.pl ' + inBase + '_250.aln ' + inBase + '.fa precursors_250.fa >' + inBase + '_250.bst')
  os.system('sort +3 -25 ' + inBase + '_250.bst >' + inBase + '_250_signature')
  os.system('perl miRDP.pl ' + inBase + '_250_signature precursors_250_structure >' + inBase + '_250_prediction')


infile = 'high_conf_mature_ath_uniq_collapsed.fa'
#infile = '100_collapsed.fa'
os.system('cp input_storage/' + infile + ' .')
genome, f_annotated = init_mirdeep_p ('tair10')
run_mirdp_new (infile) #= takes 15 secs, validates 75 unique mirna
#run_mirdp_known (infile, f_annotated) #= also takes 15 secs, validates 75 unique mirna

inBase = infile[:-3]
print('mode of predicting new: ')
os.system('cut -f1 ' + inBase + '_predictions | grep \'seq_\' | sort | uniq | wc -l')
#print('mode of predicting known: ')
#os.system('cut -f1 ' + inBase + '_250_prediction | grep \'seq_\' | sort | uniq | wc -l')

