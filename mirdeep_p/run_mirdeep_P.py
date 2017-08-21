'''
time python run_mirdeep_P.py

author: Chao-Jung Wu
date: 2017-08-
'''
import os
import time
from os import listdir

options = ['tair9', 'tair10']

def convert_raw_to_fasta500 (infile, rep_input):
  inBasename = os.path.splitext(infile)[0]
  outfile = '../tmp/' + inBasename + '.fa'
  fh_out = open (outfile, 'w')
  nb = 0
  with open (rep_input + infile, 'r') as fh:
    for i in fh:
      seq = i.split('\t')[0]
      freq = i.rstrip('\n').split('\t')[1]
      header = '>seq_' + str(nb) + '_x' + freq
      print >> fh_out, header, '\n', seq
      nb += 1
  fh_out.close()
  return outfile

def init_mirdeep_p (op):
  os.system('cp miRDP1.3/* .')
  os.system('mkdir bowtie-index 2>/dev/null')
  if op == options[0]:
    genome = 'TAIR9_genome.fa'
    rep = 'tair9/'
    f_annotated = 'annotated_miRNA_extended.fa'
    #os.system('bowtie-build -f ' + genome + ' bowtie-index/' + genome[:-3] + ' >/dev/null') 
    os.system('cp ' + rep + 'ath.gff .')#= miRBase v15
    #os.system('cp ' + rep + 'annotated_miRNA_extended.fa .')
    os.system('cp genome/' + genome + ' .')
    os.system('cp ' + rep + 'ncRNA_CDS.gff .')#
    os.system('cp ' + rep + 'chromosome_length .')#
    os.system('perl fetch_extended_precursors.pl ' + genome + ' ath.gff3.edited >' + f_annotated)
    os.system('bowtie-build -f annotated_miRNA_extended.fa bowtie-index/annotated_miRNA_extended.fa >/dev/null')
  elif op == options[1]:
    genome = 'a_thaliana_t10.fa'
    rep = 'tair10/'
    f_annotated = 'annotated_miRNA_v21_extended.fa'
    #os.system('bowtie-build -f ' + genome + ' bowtie-index/' + genome[:-3] + ' >/dev/null') 
    os.system('cp ../dbs/bowtie_index/*.ebwt bowtie-index/')
    os.system('cp ' + rep + 'ath.gff3.edited .')#= miRBase v21
    os.system('cp ' + rep + 'annotated_miRNA_v21_extended.fa .')
    #os.system('cp genome/' + genome + ' .')
    os.system('cat ../dbs/ATH/Genome/TAIR10_*.fas > ' + genome)#########
    os.system('cp ' + rep + 'ncRNA_CDS.gff .')#
    os.system('cp ' + rep + 'chromosome_length .')#
    #os.system('perl fetch_extended_precursors.pl ' + genome + ' ath.gff3.edited >' + f_annotated)
    os.system('bowtie-build -f annotated_miRNA_v21_extended.fa bowtie-index/annotated_miRNA_v21_extended.fa >/dev/null')
  return genome, f_annotated

def run_mirdp_new (infile):
  os.system('bowtie -a -v 0 bowtie-index/' + genome[:-3] + ' -f ' + infile + '>indata.aln 2>/dev/null')
  os.system('perl convert_bowtie_to_blast.pl indata.aln ' + infile + ' ' + genome + ' >indata.bst 2>/dev/null')
  os.system('perl filter_alignments.pl indata.bst -c 15 >indata_filter15.bst 2>/dev/null')
  os.system('perl overlap.pl indata_filter15.bst ncRNA_CDS.gff -b >indata_id_overlap_ncRNA_CDS 2>/dev/null')
  os.system('perl alignedselected.pl indata_filter15.bst -g indata_id_overlap_ncRNA_CDS >indata_filter15_ncRNA_CDS.bst 2>/dev/null')
  os.system('perl filter_alignments.pl indata_filter15_ncRNA_CDS.bst -b ' + infile + ' > indata_filtered.fa 2>/dev/null')
  os.system('perl excise_candidate.pl ' + genome + ' indata_filter15_ncRNA_CDS.bst 250 >indata_precursors.fa 2>/dev/null')
  os.system('cat indata_precursors.fa | RNAfold --noPS > indata_structures')
  os.system('bowtie-build -f indata_precursors.fa bowtie-index/indata_precursors >/dev/null 2>/dev/null')
  os.system('bowtie -a -v 0 bowtie-index/indata_precursors -f indata_filtered.fa > indata_precursors.aln 2>/dev/null')
  os.system('perl convert_bowtie_to_blast.pl indata_precursors.aln indata_filtered.fa indata_precursors.fa >indata_precursors.bst 2>/dev/null')
  os.system('sort +3 -25 indata_precursors.bst >indata_signatures')
  os.system('perl miRDP.pl indata_signatures indata_structures > result_new_predictions')
  #os.system('perl rm_redundant_meet_plant.pl chromosome_length indata_precursors.fa indata_predictions indata_nr_prediction indata_filter_P_prediction')

def run_mirdp_known (infile, f_annotated):
  os.system('bowtie -a -v 0 bowtie-index/' + f_annotated + ' -f ' + infile + ' >indata.aln 2>/dev/null')
  os.system('perl convert_bowtie_to_blast.pl indata.aln ' + infile + ' ' + f_annotated + ' > indata_extended.bst 2>/dev/null')
  os.system('perl excise_candidate.pl ' + f_annotated + ' indata_extended.bst 250 >precursors_250.fa 2>/dev/null')
  os.system('bowtie-build -f precursors_250.fa bowtie-index/precursors_250 >/dev/null 2>/dev/null')
  os.system('cat precursors_250.fa|RNAfold --noPS >precursors_250_structure')
  os.system('bowtie -a -v 0 bowtie-index/precursors_250 -f ' + infile + ' >indata_250.aln 2>/dev/null')
  os.system('perl convert_bowtie_to_blast.pl indata_250.aln ' + infile + ' precursors_250.fa >indata_250.bst 2>/dev/null')
  os.system('sort +3 -25 indata_250.bst >indata_250_signature')
  os.system('perl miRDP.pl indata_250_signature precursors_250_structure > result_known_250_predictions')


def process_one_file (infile):
  start = time.time()

  run_mirdp_new (infile)
  run_mirdp_known (infile, f_annotated)

  end = time.time()

  duration = end - start
  duration = format(duration, '.0f')
  print 'miRDP duration: ', infile, duration, 'sec'

  #print('mode of predicting new: ')
  #os.system('cut -f1 result_new_predictions | grep \'seq_\' | sort | uniq | wc -l')
  #print('mode of predicting known: ')
  #os.system('cut -f1 result_known_250_predictions | grep \'seq_\' | sort | uniq | wc -l')
  #print('combinded prediction: ')
  os.system('cat result_* > result_combinded.txt')
  os.system('cut -f1 result_combinded.txt | grep \'seq_\' | sort | uniq | wc -l')
  os.system('rm -f indata* *.gff* *.pl precursors_250_structure *.fa chromosome_length result_combinded.txt bowtie-index/*precursors*')

def run_miRDP ():
  #rep_input = '/home/cloudera/Desktop/mirLibHadoop/input/'
  rep_input = '../input/'
  #rep_input = '/home/cjwu/demo/srna2/'
  infiles = [f for f in listdir(rep_input) if os.path.isfile(os.path.join(rep_input, f))]
  totaltime = 0
  for infile in infiles:
    print 'start processing', infile
    infile = convert_raw_to_fasta500 (infile, rep_input)
    start = time.time()

    run_mirdp_new (infile)
    run_mirdp_known (infile, f_annotated)

    end = time.time()

    duration = end - start
    totaltime += duration
    duration = format(duration, '.0f')
    print 'miRDP duration: ', infile, duration, 'sec'

    os.system('cat result_* > result_combinded.txt')
    os.system('cut -f1 result_combinded.txt | grep \'seq_\' | sort | uniq | wc -l')
    os.system('rm -f indata* *.gff* *.pl precursors_250_structure *.fa chromosome_length result_combinded.txt bowtie-index/*precursors*')
  
  totaltime = format(totaltime, '.0f')
  print 'miRDP duration: ', totaltime, 'sec'

genome, f_annotated = init_mirdeep_p ('tair10')
run_miRDP ()


