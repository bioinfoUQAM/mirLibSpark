import argparse
import utils as ut

def getOpt (paramDict):
    parser = argparse.ArgumentParser()
    parser.add_argument('--message')
    parser.add_argument('--project_path')
    parser.add_argument('--input_path')
    parser.add_argument('--output_path')



    parser.add_argument('--species', choices=['ath', 'wheat', 'corn', 'rice', 'potato', 'brome', 'wheatD'])
    parser.add_argument('--bowtie_index_suffix')



    parser.add_argument('--genome_path')
    parser.add_argument('--b_index_path')
    parser.add_argument('--known_non_file')
    parser.add_argument('--target_file')
    parser.add_argument('--perform_differnatial_analysis', choices=['yes', 'no'])
    parser.add_argument('--diffguide_file')
    parser.add_argument('--perform_KEGGpathways_enrichment_analysis', choices=['yes', 'no'])
    parser.add_argument('--gene_vs_pathway_file')
    parser.add_argument('--pathway_description_file')
    parser.add_argument('--input_type', choices=['raw', 'w', 'reads','r', 'fasta', 'a', 'fastq', 'q'])
    parser.add_argument('--adapter')
    parser.add_argument('--mirdup_model', choices=['Viridiplantae.model', 'thaliana.model'])
    parser.add_argument('--limit_s_freq')
    parser.add_argument('--limit_m_freq')
    parser.add_argument('--limit_len')
    parser.add_argument('--limit_nbLoc')
    parser.add_argument('--pri_l_flank')
    parser.add_argument('--pri_r_flank')
    parser.add_argument('--pre_flank')
    parser.add_argument('--temperature')
    parser.add_argument('--mirdup_limit')
    parser.add_argument('--mcheck_param')
    parser.add_argument('--Max_Score_cutoff')
    parser.add_argument('--Max_Energy_cutoff')
    parser.add_argument('--Gap_Penalty')
    parser.add_argument('--nbTargets')
    parser.add_argument('--sc_partition')
    parser.add_argument('--sc_mstrmemory')
    parser.add_argument('--sc_execmemory')
    parser.add_argument('--sc_master')
    #parser.add_argument('--sc_execcores')
    

 
    
    args = parser.parse_args()

    if args.message == None: args.species = paramDict['message']
    if args.project_path == None: args.project_path = paramDict['project_path'] 
    if args.project_path == None: args.project_path = paramDict['project_path'] 
    if args.project_path == None: args.project_path = paramDict['project_path'] 
    if args.project_path == None: args.project_path = paramDict['project_path'] 
    #
    if args.species == 'ath': 
      bowtie_index_suffix = 'ATH_TAIR10'
      filename1 = 'Arabidopsis_thaliana.TAIR10.cdna.all.fa'
      filename2 = 'ath_gene_vs_pathway.txt'
      filename3 = 'ath_pathway_description.txt'
    elif args.species == 'wheat': 
      bowtie_index_suffix = 'WHEAT_IWGSC'
      filename1 = 'mRNA_contigs_smaller.fasta'
      filename2 = 'ko_gene_vs_pathway.txt'
      filename3 = 'ko_pathway_description.txt'
    elif args.species == 'corn': 
      bowtie_index_suffix = 'CORN_AGPv4'
      filename1 = 'Zea_mays.AGPv4.cdna.all.fa.gz'
      filename2 = 'zma_gene_vs_pathway.txt'
      filename3 = 'zma_pathway_description.txt'
    elif args.species == 'rice': 
      bowtie_index_suffix = 'RICE_IRGSP_1'
      filename1 = 'Oryza_sativa.IRGSP-1.0.cdna.all.fa'
      filename2 = 'osa_gene_vs_pathway.txt'
      filename3 = 'osa_pathway_description.txt'
    elif args.species == 'potato': 
      bowtie_index_suffix = 'POTATO_SolTub_3'
      filename1 = 'Solanum_tuberosum.SolTub_3.0.cdna.all.fa'
      filename2 = 'sot_gene_vs_pathway.txt'
      filename3 = 'sot_pathway_description.txt'
    elif args.species == 'brome': 
      bowtie_index_suffix = 'STIFFBROME'
      filename1 = 'Brachypodium_distachyon.Brachypodium_distachyon_v3.0.cdna.all.fa'
      filename2 = 'bdi_gene_vs_pathway.txt'
      filename3 = 'bdi_pathway_description.txt'
    elif args.species == 'wheatD': 
      bowtie_index_suffix = 'WHEAT_D'
      filename1 = 'Aegilops_tauschii.ASM34733v1.cdna.all.fa'
      filename2 = 'ats_gene_vs_pathway.txt'
      filename3 = 'ats_pathway_description.txt'


    if args.bowtie_index_suffix == None: args.bowtie_index_suffix = bowtie_index_suffix
    #
    key = 'dbs/' + args.bowtie_index_suffix + '/Genome/'
    if args.genome_path == None: args.genome_path = args.project_path + key   
    # 
    key = 'dbs/' + args.bowtie_index_suffix + '/bowtie_index/'
    if args.b_index_path == None: args.b_index_path = args.project_path + key   
    #
    key = 'dbs/' + 'ATH_TAIR10' + '/TAIR10_ncRNA_CDS.gff'
    if args.known_non_file == None: args.known_non_file = args.project_path + key   
    #
    key = 'dbs/' + args.bowtie_index_suffix + '/' + filename1
    if args.target_file == None: args.target_file = args.project_path + key
    #
    if args.perform_differnatial_analysis == None: args.perform_differnatial_analysis = 'no'
    if args.perform_differnatial_analysis == 'yes': 
      if args.diffguide_file == None: print('ERROR, to do sys.exit()')
    #
    if args.perform_KEGGpathways_enrichment_analysis == None: args.perform_differnatial_analysis = 'no'
    if args.perform_KEGGpathways_enrichment_analysis == 'yes': 
      if args.perform_differnatial_analysis == 'no': print('ERROR, to do sys.exit()')
    #
    key = 'dbs/' + args.bowtie_index_suffix + '/' + filename2
    if args.gene_vs_pathway_file == None: args.gene_vs_pathway_file = args.project_path + key
    #
    key = 'dbs/' + args.bowtie_index_suffix + '/' + filename3
    if args.pathway_description_file == None: args.pathway_description_file = args.project_path + key
    #
    if args.input_type == 'raw': args.input_type = 'w'
    elif args.input_type == 'reads': args.input_type = 'r'
    elif args.input_type == 'fasta': args.input_type = 'a'
    elif args.input_type == 'fastq': args.input_type = 'q'
    #
    if args.mirdup_model == None: args.mirdup_model = 'Viridiplantae.model'
    if args.limit_s_freq == None: args.limit_s_freq = '10'
    if args.limit_s_freq == None: args.limit_s_freq = '10'
    



    





    a = args.bowtie_index_suffix

    print(a)




paramfile = '../getoptParame.txt'
paramDict = ut.readParam (paramfile)

#for k in sorted(paramDict): print(k, '=', paramDict[k])


getOpt (paramDict)




'''
ath: 	Arabidopsis_thaliana.TAIR10
wheat: 	Triticum_aestivum.IWGSC	
corn: 	Zea_mays.AGPv4	
rice: 	Oryza_sativa.IRGSP-1.0
potato: Solanum_tuberosum.SolTub_3.0
brome: 	Brachypodium_distachyon.Brachypodium_distachyon_v3.0
wheatD: Aegilops_tauschii.ASM34733v1



ATH_TAIR10  POTATO_SolTub_3  WHEAT_D
CORN_AGPv4  RICE_IRGSP_1     WHEAT_inhouse
miRBase_21  STIFFBROME	     WHEAT_IWGSC
'''

