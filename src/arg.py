'''
program: arg.py
author: Chao-Jung Wu
date: 2018-10-30
version: 0.00.01
'''
from __future__ import print_function
import argparse
import os
import sys


def find_project_path ():
  import os
  cwd = os.getcwd()
  project_path = cwd.split('/mirLibSpark')[0] + '/mirLibSpark'
  return project_path

def getOpt (parser): 
    project_path = find_project_path ()
    #
    parser.add_argument('--dummy', action='store_true')
    parser.add_argument('--message')
    parser.add_argument('--project_path', default = project_path)
    parser.add_argument('--input_path')
    parser.add_argument('--output_path')
    parser.add_argument('--species', default = 'ath',\
                         choices=['ath', 'wheat', 'corn', 'rice', 'potato', 'brome', 'wheatD'], \
                         help='ath: 	Arabidopsis_thaliana.TAIR10;\
                               wheat: 	Triticum_aestivum.IWGSC;\
                               corn: 	Zea_mays.AGPv4;\
                               rice: 	Oryza_sativa.IRGSP-1.0;\
                               potato: Solanum_tuberosum.SolTub_3.0;\
                               brome: 	Brachypodium_distachyon.Brachypodium_distachyon_v3.0;\
                               wheatD: Aegilops_tauschii.ASM34733v1.\
                               Please use provided script to construct the dbs folder for selected species from ensembl-release40.')
    parser.add_argument('--input_type', default='w', choices=['raw', 'w', 'reads','r', 'fasta', 'a', 'fastq', 'q'])
    parser.add_argument('--adapter', default='no', help='example = TGGAATTCTCGGGTGCCAAGGAACTC')
    parser.add_argument('--bowtie_index_suffix')
    parser.add_argument('--genome_path')
    parser.add_argument('--b_index_path')
    parser.add_argument('--known_non_file', help='Only ath is provided.')
    parser.add_argument('--chromosomes', choices=['All', 'split'], help='Only ath has options. Other species have only one default choice. Genomes larger than 1G are splitted by chromosomes, otherwise all choromsomes are in one file.')
    parser.add_argument('--target_file')
    parser.add_argument('--perform_differnatial_analysis', action='store_true')
    parser.add_argument('--diffguide_file')
    parser.add_argument('--perform_KEGGpathways_enrichment_analysis', action='store_true')
    parser.add_argument('--gene_vs_pathway_file')
    parser.add_argument('--pathway_description_file')
    parser.add_argument('--mirdup_model', default='thaliana.model',\
                         choices=['Viridiplantae.model', 'thaliana.model'])
    parser.add_argument('--limit_s_freq', default='10')
    parser.add_argument('--limit_m_freq', default='100')
    parser.add_argument('--limit_len', default='18')
    parser.add_argument('--limit_nbLoc', default='15')
    parser.add_argument('--pri_l_flank', default='500')
    parser.add_argument('--pri_r_flank', default='200')
    parser.add_argument('--pre_flank', default='10')
    parser.add_argument('--temperature', default='25', help='Celsius')
    parser.add_argument('--mirdup_limit', default='0.98')
    parser.add_argument('--mcheck_param', default='def')
    parser.add_argument('--Max_Score_cutoff', default='170')
    parser.add_argument('--Max_Energy_cutoff', default='-15')
    parser.add_argument('--Gap_Penalty', default='-15')
    parser.add_argument('--nbTargets', default='100')
    parser.add_argument('--sc_partition', default='64')
    parser.add_argument('--sc_mstrmemory', default='20g')
    parser.add_argument('--sc_execmemory', default='20g')
    parser.add_argument('--sc_master', default='local[*]')
    #
    args = parser.parse_args()
    #
    args.project_path = args.project_path.rstrip('/')
    if args.input_path == None: args.input_path = args.project_path + '/input/'
    if args.output_path == None: args.output_path = args.project_path + '/output/'
    #
    if args.species == 'ath': 
      bowtie_index_suffix = 'ATH_TAIR10'
      filename1 = 'Arabidopsis_thaliana.TAIR10.cdna.all.fa'
      filename2 = 'ath_gene_vs_pathway.txt'
      filename3 = 'ath_pathway_description.txt'
      chromosomes = '1,2,3,4,5,Mt,Pt'
      if args.chromosomes == None: args.chromosomes = 'All'
      elif args.chromosomes == 'split': args.chromosomes = chromosomes
    elif args.species == 'wheat': 
      bowtie_index_suffix = 'WHEAT_IWGSC'
      filename1 = 'mRNA_contigs_smaller.fasta'
      filename2 = 'ko_gene_vs_pathway.txt'
      filename3 = 'ko_pathway_description.txt'
      chromosome = '1A,1B,1D,2A,2B,2D,3A,3B,3D,4A,4B,4D,5A,5B,5D,6A,6B,6D,7A,7B,7D'
      args.chromosomes = chromosomes
    elif args.species == 'corn': 
      bowtie_index_suffix = 'CORN_AGPv4'
      filename1 = 'Zea_mays.AGPv4.cdna.all.fa.gz'
      filename2 = 'zma_gene_vs_pathway.txt'
      filename3 = 'zma_pathway_description.txt'
      chromosome = '1,2,3,4,5,6,7,8,9,10,Mt,Pt'
      args.chromosomes = chromosomes
    elif args.species == 'rice': 
      bowtie_index_suffix = 'RICE_IRGSP_1'
      filename1 = 'Oryza_sativa.IRGSP-1.0.cdna.all.fa'
      filename2 = 'osa_gene_vs_pathway.txt'
      filename3 = 'osa_pathway_description.txt'
      chromosome = '1,2,3,4,5,6,7,8,9,10,11,12,Mt,Pt'
      args.chromosomes = 'All'
    elif args.species == 'potato': 
      bowtie_index_suffix = 'POTATO_SolTub_3'
      filename1 = 'Solanum_tuberosum.SolTub_3.0.cdna.all.fa'
      filename2 = 'sot_gene_vs_pathway.txt'
      filename3 = 'sot_pathway_description.txt'
      chromosome = '1,2,3,4,5,6,7,8,9,10,11,12'
      args.chromosomes = 'All'
    elif args.species == 'brome': 
      bowtie_index_suffix = 'STIFFBROME'
      filename1 = 'Brachypodium_distachyon.Brachypodium_distachyon_v3.0.cdna.all.fa'
      filename2 = 'bdi_gene_vs_pathway.txt'
      filename3 = 'bdi_pathway_description.txt'
      chromosome = '1,2,3,4,5'
      args.chromosomes = 'All'
    elif args.species == 'wheatD': 
      bowtie_index_suffix = 'WHEAT_D'
      filename1 = 'Aegilops_tauschii.ASM34733v1.cdna.all.fa'
      filename2 = 'ats_gene_vs_pathway.txt'
      filename3 = 'ats_pathway_description.txt'
      chromosome = 'toplevel'
      args.chromosomes = 'All'
    #
    if args.bowtie_index_suffix == None: args.bowtie_index_suffix = bowtie_index_suffix
    #
    key = '/dbs/' + args.bowtie_index_suffix + '/Genome/'
    if args.genome_path == None: args.genome_path = args.project_path + key   
    # 
    key = '/dbs/' + args.bowtie_index_suffix + '/bowtie_index/'
    if args.b_index_path == None: args.b_index_path = args.project_path + key   
    #= only ath is provided for know_non miRNA list
    key = '/dbs/' + 'ATH_TAIR10' + '/TAIR10_ncRNA_CDS.gff'
    if args.known_non_file == None: args.known_non_file = args.project_path + key   
    #
    key = '/dbs/' + args.bowtie_index_suffix + '/' + filename1
    if args.target_file == None: args.target_file = args.project_path + key
    #
    if args.perform_differnatial_analysis == True:
      args.perform_differnatial_analysis = 'yes' 
      if args.diffguide_file == None: 
        sys.stderr.write('diffguide_file is required for perform_differnatial_analysis.\n\
                          Exit the program.')
        sys.exit()
    #
    if args.perform_KEGGpathways_enrichment_analysis == True:
      args.perform_KEGGpathways_enrichment_analysis = 'yes'
      if args.perform_differnatial_analysis == 'no': 
        sys.stderr.write('perform_differnatial_analysis is required for perform_KEGGpathways_enrichment_analysis.\n\
                          Exit the program.')
        sys.exit()
    #
    key = '/dbs/' + args.bowtie_index_suffix + '/' + filename2
    if args.gene_vs_pathway_file == None: args.gene_vs_pathway_file = args.project_path + key
    #
    key = '/dbs/' + args.bowtie_index_suffix + '/' + filename3
    if args.pathway_description_file == None: args.pathway_description_file = args.project_path + key
    #
    if args.input_type == 'raw': args.input_type = 'w'
    elif args.input_type == 'reads': args.input_type = 'r'
    elif args.input_type == 'fasta': args.input_type = 'a'
    elif args.input_type == 'fastq': args.input_type = 'q'
    #
    #= store args in dict
    paramDict = vars(args)
    #= add additional parameters in dict
    paramDict['sc_appname'] = 'mirLibSpark'
    if args.dummy == True:
      for k, v in sorted(paramDict.items()): print(k, ': ', v)
      print('============================================================\n')
      sys.stderr.write('Display registered parameters.\n\
                        Exit the program.')
      sys.exit()
    return paramDict


def init_params ():
  
  parser = argparse.ArgumentParser()
  paramDict = getOpt(parser)
  return paramDict





if __name__ == '__main__' :
  paramDict = init_params ()
  for k in sorted(paramDict): print(k, ':', paramDict[k])




