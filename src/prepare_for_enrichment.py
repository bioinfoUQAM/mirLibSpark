'''
This is the main to create the background files for enrichment analysis, so it is a preperation step.

Chao-Jung Wu
2018-07-13

update: 2018-09-25


This is a refactoring version 2018-09-25
WIP


How Res1 was created?

'''

import os

## global variables
#
#cleanup = 'TBD 0 No_GO No_GOBP No_GOMF No_GOCC No_Pathway'.split()
#
diffKey = 'UP DOWN'.split()
rep_output = '../output/'
if not os.path.exists(rep_output): os.makedirs(rep_output)####180925


#= ref files
Res2 = '../dbs/TAE/GOKEGG/Res2_miRNA_GOKEGG_functions.txt'
ref3 = '../dbs/TAE/GOKEGG/ref3_GOcode_description.txt'
ref2 = '../dbs/TAE/GOKEGG/ref2_Pathwaycode_description.txt'
#= diff files
rep_diff_y2010 = '../enrichment_training/dif_files/y2010/'
#rep_diff_y2010 = '../dif_files/y2010/'
#rep_diff_y2013 = '../dif_files/y2013/'
infiles  = [rep_diff_y2010 + f for f in os.listdir(rep_diff_y2010) if os.path.isfile(os.path.join(rep_diff_y2010, f))]
#infiles += [rep_diff_y2013 + f for f in os.listdir(rep_diff_y2013) if os.path.isfile(os.path.join(rep_diff_y2013, f))]
#
#= DATA
infile = Res2
with open (infile, 'r') as fh: DATA = [x.rstrip('\n').split('\t') for x in fh.readlines()]
ids = DATA[0][1:]
field_to_look = [ x for x in range (1, len(ids)+ 1) ]
#


#= store GO and Pathway description in the same dictionary
#= {ko04978: Mineral absorption	Organismal Systems}, {GO:0006351: transcription, DNA-templated}
d_GOorPathway_desc = {}
infile = ref3 #=GO
with open (infile, 'r') as fh: data = [x.rstrip('\n').split('\t') for x in fh.readlines()]
for i in data: d_GOorPathway_desc[i[0]] = i[1]
infile = ref2 #= KEGG pathways
with open (infile, 'r') as fh: data = [x.rstrip('\n').split('\t') for x in fh.readlines()]
for i in data: d_GOorPathway_desc[i[0]] = i[1]


#########################################################################
def __write_namecodefile (folder, ID):
  outfile = folder + '/namecode.txt'
  with open (outfile, 'w') as fh:
    print('background_' + ID + '\tbackground' , file = fh)
    print('Norstar_dif_2_1_1013\tNorstar_210' , file = fh)
    print('Norstar_dif_3_1_1014\tNorstar_310' , file = fh)
    print('Norstar_dif_4_1_1015\tNorstar_410' , file = fh)
    print('Norstar_dif_5_1_1016\tNorstar_510' , file = fh)
    #print('Norstar_dif_2_1_1313\tNorstar_213' , file = fh)
    #print('Norstar_dif_3_1_1314\tNorstar_313' , file = fh)
    #print('Norstar_dif_4_1_1315\tNorstar_413' , file = fh)
    #print('Norstar_dif_5_1_1316\tNorstar_513' , file = fh)
    print('Iso8S_dif_2_1_1001\tIso8S_210' , file = fh)
    print('Iso8S_dif_3_1_1002\tIso8S_310' , file = fh)
    print('Iso8S_dif_4_1_1003\tIso8S_410' , file = fh)
    print('Iso8S_dif_5_1_1004\tIso8S_510' , file = fh)
    #print('Iso8S_dif_2_1_1301\tIso8S_213' , file = fh)
    #print('Iso8S_dif_3_1_1302\tIso8S_313' , file = fh)
    #print('Iso8S_dif_4_1_1303\tIso8S_413' , file = fh)
    #print('Iso8S_dif_5_1_1304\tIso8S_513' , file = fh)
    print('Manitou_dif_2_1_1009\tManitou_210' , file = fh)
    print('Manitou_dif_3_1_1010\tManitou_310' , file = fh)
    print('Manitou_dif_4_1_1011\tManitou_410' , file = fh)
    print('Manitou_dif_5_1_1012\tManitou_510' , file = fh)
    #print('Manitou_dif_2_1_1309\tManitou_213' , file = fh)
    #print('Manitou_dif_3_1_1310\tManitou_313' , file = fh)
    #print('Manitou_dif_4_1_1311\tManitou_413' , file = fh)
    #print('Manitou_dif_5_1_1312\tManitou_513' , file = fh)
    print('Iso11W_dif_2_1_1005\tIso11W_210' , file = fh)
    print('Iso11W_dif_3_1_1006\tIso11W_310' , file = fh)
    print('Iso11W_dif_4_1_1007\tIso11W_410' , file = fh)
    print('Iso11W_dif_5_1_1008\tIso11W_510' , file = fh)
    #print('Iso11W_dif_2_1_1305\tIso11W_213' , file = fh)
    #print('Iso11W_dif_3_1_1306\tIso11W_313' , file = fh)
    #print('Iso11W_dif_4_1_1307\tIso11W_413' , file = fh)
    #print('Iso11W_dif_5_1_1308\tIso11W_513' , file = fh)

def __create_background_of_this_field (outfile, field, DATA):
  fh_out = open (outfile, 'w')
  for i in DATA[1:]:
    seq = i[0]
    lookat = i[int(field)]
    if 'No_' in lookat or 'TBD' in lookat or len(lookat) < 2: continue
    lookat = lookat.split(',')
    for j in lookat: 
      #if 'A) tail shortening' in j : j = 'GO:0060213'
      if j == '': desc = '0'
      else:
        desc = d_GOorPathway_desc[j]
      items = '\t'.join([seq, j, desc])
      print (items, file=fh_out, flush=True)
  fh_out.close()

#= create background files in each feature folder
count = 0
d_field_background = {}
for field in field_to_look:
  ID = '_'.join(ids[count].split(' '))
  folder = rep_output + ID
  if not os.path.exists(folder): os.makedirs(folder)
  __write_namecodefile (folder, ID)
  outfile = folder + '/background_' + ID 
  __create_background_of_this_field (outfile, field, DATA)
  count += 1
#############################################################################  
 


 
############################################################################## 
d_seq_FuncAnnot = {}
for i in DATA[1:]:
  seq = i[0]
  FuncAnnot = i[1:16]
  d_seq_FuncAnnot[seq] = '#'.join(FuncAnnot)

def __one_seq_annotation (outfile, lookat):
  if 'No_' in lookat or 'TBD' in lookat or len(lookat) < 2: return
  fh_out = open (outfile, 'a')
  lookat = lookat.split(',')
  for j in lookat: 
    #if 'A) tail shortening' in j : j = 'GO:0060213'
    print (seq + '\t' + j, file=fh_out, flush=True)
  fh_out.close()  




def __purge (outfile):
  ''' 
  if the first time open this file, empty the content of this file, because the next opening is by 'a' ajouter
  '''
  fh_out = open (outfile, 'w')
  fh_out.close()
  
purge = []
for infile in infiles:
  namecode = infile.split('/')[-1][:-4]
  if infile.endswith('.swp') or infile.endswith('~'): continue
  with open (infile, 'r') as fh: data = [x.rstrip('\n').split('\t') for x in fh.readlines()]
  for i in data:
    seq = i[0]
    Diff_exp = i[7]
    if Diff_exp in diffKey: 
      infos = d_seq_FuncAnnot[seq]
      for index in range(len(ids)):
        ID = '_'.join(ids[index].split(' '))
        rep = rep_output + ID + '/'
        enrich_rep = rep_output + ID + '/output_comput_enrich/'
        if not os.path.exists(enrich_rep): os.makedirs(enrich_rep)
        outfile = rep + namecode
        if outfile not in purge:
          purge.append(outfile)
          __purge (outfile)
        seq_info = infos.split('#')[index]  
        __one_seq_annotation (outfile, seq_info)  
##########################################################################################  
