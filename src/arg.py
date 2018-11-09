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
  #project_path = cwd.split('/mirLibSpark')[0] + '/mirLibSpark'
  project_path = cwd.split('/src')[0].split('/workdir')[0]
  return project_path

def getOpt (parser): 
    project_path = find_project_path ()
    #
    parser.add_argument('--dummy', action='store_true')
    parser.add_argument('--reporting', default='0', action='store_true', help='report the number of instances passing each of the parameters. Activating this option will result in significant prolonged execution time')
    parser.add_argument('--message', default = 'None')
    parser.add_argument('--project_path', default = project_path)
    parser.add_argument('--input_path')
    parser.add_argument('--output_path')
    parser.add_argument('--species', default = 'ath',\
                         choices=['ath', 'wheat', 'corn', 'rice', 'potato', 'brome', 'wheatD'],
                         help='ath: 	Arabidopsis_thaliana.TAIR10;\
                               wheat: 	Triticum_aestivum.IWGSC;\
                               corn: 	Zea_mays.AGPv4;\
                               rice: 	Oryza_sativa.IRGSP-1.0;\
                               potato: Solanum_tuberosum.SolTub_3.0;\
                               brome: 	Brachypodium_distachyon.Brachypodium_distachyon_v3.0;\
                               wheatD: Aegilops_tauschii.ASM34733v1.\
                               Please use provided script to construct the dbs folder for selected species from ensembl-release40.')
    parser.add_argument('--input_type', default='w', choices=['raw', 'w', 'reads','r', 'fasta', 'a', 'fastq', 'q'])
    parser.add_argument('--adapter', default='none', help='example = TGGAATTCTCGGGTGCCAAGGAACTC')
    parser.add_argument('--bowtie_index_prefix')
    parser.add_argument('--genome_path')
    parser.add_argument('--b_index_path')
    parser.add_argument('--known_non_file', help='Only ath is provided.')
    parser.add_argument('--chromosomes', choices=['All', 'split'], help='Only ath has options. Other species have only one default choice. Genomes larger than 1G are splitted by chromosomes, otherwise all choromsomes are in one file.')
    parser.add_argument('--target_file')
    parser.add_argument('--perform_differnatial_analysis', default='False', action='store_true')
    parser.add_argument('--diffguide_file', default='none')
    parser.add_argument('--perform_KEGGpathways_enrichment_analysis', action='store_true')
    parser.add_argument('--gene_vs_pathway_file')
    parser.add_argument('--pathway_description_file')
    parser.add_argument('--mirdup_model', default='thaliana.model',\
                         choices=['Viridiplantae.model', 'thaliana.model'])
    parser.add_argument('--limit_s_freq', default='10')
    parser.add_argument('--limit_m_freq', default='100')
    parser.add_argument('--limit_len', default='19')
    parser.add_argument('--miRNA_len_upperlimit', default='24')
    parser.add_argument('--miRNA_len_lowerlimit', default='21')
    parser.add_argument('--premirna_max_len', default='300')
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
    if args.reporting == True: args.reporting = '1'
    #
    args.project_path = args.project_path.rstrip('/')
    if args.input_path == None: args.input_path = args.project_path + '/input/'
    else: args.input_path = args.input_path.rstrip('/') + '/'
    if args.output_path == None: args.output_path = args.project_path + '/output/'
    else: args.output_path = args.output_path.rstrip('/') + '/'
    #
    if args.species == 'ath': 
      bowtie_index_prefix = 'ATH_TAIR10'
      filename1 = 'Arabidopsis_thaliana.TAIR10.cdna.all.fa'
      filename2 = 'ath_gene_vs_pathway.txt'
      filename3 = 'ath_pathway_description.txt'
      chromosomes = '1,2,3,4,5,Mt,Pt'
      if args.chromosomes == None: args.chromosomes = 'All'
      elif args.chromosomes == 'split': args.chromosomes = chromosomes
    elif args.species == 'wheat': 
      bowtie_index_prefix = 'WHEAT_IWGSC'
      filename1 = 'mRNA_contigs_smaller.fasta'
      filename2 = 'ko_gene_vs_pathway.txt'
      filename3 = 'ko_pathway_description.txt'
      chromosomes = '1A,1B,1D,2A,2B,2D,3A,3B,3D,4A,4B,4D,5A,5B,5D,6A,6B,6D,7A,7B,7D'
      args.chromosomes = chromosomes
    elif args.species == 'corn': 
      bowtie_index_prefix = 'CORN_AGPv4'
      filename1 = 'Zea_mays.AGPv4.cdna.all.fa.gz'
      filename2 = 'zma_gene_vs_pathway.txt'
      filename3 = 'zma_pathway_description.txt'
      chromosomes = '1,2,3,4,5,6,7,8,9,10,Mt,Pt'
      args.chromosomes = chromosomes
    elif args.species == 'rice': 
      bowtie_index_prefix = 'RICE_IRGSP_1'
      filename1 = 'Oryza_sativa.IRGSP-1.0.cdna.all.fa'
      filename2 = 'osa_gene_vs_pathway.txt'
      filename3 = 'osa_pathway_description.txt'
      chromosome = '1,2,3,4,5,6,7,8,9,10,11,12,Mt,Pt'
      args.chromosomes = 'All'
    elif args.species == 'potato': 
      bowtie_index_prefix = 'POTATO_SolTub_3'
      filename1 = 'Solanum_tuberosum.SolTub_3.0.cdna.all.fa'
      filename2 = 'sot_gene_vs_pathway.txt'
      filename3 = 'sot_pathway_description.txt'
      chromosomes = '1,2,3,4,5,6,7,8,9,10,11,12'
      args.chromosomes = 'All'
    elif args.species == 'brome': 
      bowtie_index_prefix = 'STIFFBROME'
      filename1 = 'Brachypodium_distachyon.Brachypodium_distachyon_v3.0.cdna.all.fa'
      filename2 = 'bdi_gene_vs_pathway.txt'
      filename3 = 'bdi_pathway_description.txt'
      chromosomes = '1,2,3,4,5'
      args.chromosomes = 'All'
    elif args.species == 'wheatD': 
      bowtie_index_suffix = 'WHEAT_D'
      bowtie_index_prefix = 'WHEAT_D'
      filename1 = 'Aegilops_tauschii.ASM34733v1.cdna.all.fa'
      filename2 = 'ats_gene_vs_pathway.txt'
      filename3 = 'ats_pathway_description.txt'
      chromosomes = 'toplevel'
      args.chromosomes = 'All'
    #
    if args.bowtie_index_prefix == None: args.bowtie_index_prefix = bowtie_index_prefix
    #
    key = '/dbs/' + args.bowtie_index_prefix + '/Genome/'
    if args.genome_path == None: args.genome_path = args.project_path + key   
    # 
    key = '/dbs/' + args.bowtie_index_prefix + '/bowtie_index/'
    if args.b_index_path == None: args.b_index_path = args.project_path + key   
    #= only ath is provided for know_non miRNA list
    key = '/dbs/' + 'ATH_TAIR10' + '/TAIR10_ncRNA_CDS.gff'
    if args.known_non_file == None: args.known_non_file = args.project_path + key   
    #
    key = '/dbs/' + args.bowtie_index_prefix + '/' + filename1
    if args.target_file == None: args.target_file = args.project_path + key
    #
    if args.perform_differnatial_analysis == False: args.perform_differnatial_analysis = 'False'
    if args.perform_differnatial_analysis == True:
      args.perform_differnatial_analysis = 'yes' 
      if args.diffguide_file == 'none': 
        sys.stderr.write('diffguide_file is required for perform_differnatial_analysis.\n\
                          Exit the program.')
        sys.exit()
    #
    if args.perform_KEGGpathways_enrichment_analysis == False: args.perform_KEGGpathways_enrichment_analysis = 'False'
    if args.perform_KEGGpathways_enrichment_analysis == True:
      args.perform_KEGGpathways_enrichment_analysis = 'yes'
      if args.perform_differnatial_analysis == 'False': 
        sys.stderr.write('perform_differnatial_analysis is required for perform_KEGGpathways_enrichment_analysis.\n\
                          Exit the program.')
        sys.exit()
    #
    key = '/dbs/' + args.bowtie_index_prefix + '/' + filename2
    if args.gene_vs_pathway_file == None: args.gene_vs_pathway_file = args.project_path + key
    #
    key = '/dbs/' + args.bowtie_index_prefix + '/' + filename3
    if args.pathway_description_file == None: args.pathway_description_file = args.project_path + key
    #
    if args.input_type == 'w': args.input_type = 'raw'
    elif args.input_type == 'r': args.input_type = 'reads'
    elif args.input_type == 'a': args.input_type = 'fasta'
    elif args.input_type == 'q': args.input_type = 'fastq'
    #
    #= store args in dict
    paramDict = vars(args)
    #= add additional parameters in dict
    paramDict['sc_appname'] = 'mirLibSpark'
    #
    if args.dummy == False: paramDict['dummy'] = 'False'
    elif args.dummy == True:
      paramDict['dummy'] = 'True'
      for k, v in sorted(paramDict.items()): print(k + ': ' +v)
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



