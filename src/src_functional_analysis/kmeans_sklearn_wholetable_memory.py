'''
author: Dylan (initial script)
date: 2018-01-12

the approach to see the whole table of occurence is not discriminant
'''
import csv
import numpy as np
from sklearn.cluster import KMeans


infile = 'csv127x859.csv' # to be grouped
#infile = 'occ_of_each_miRNA.txt'
outfile = infile + 'cluster_result.txt'
cluster_file = infile + 'cluster.txt'
confidence_file = 'confidence.txt'
occSum_file = 'occ_of_each_miRNA.txt'
k=30 #k >=2

max_occ = 120



def kmeans_cluster (k):
  list_cluster = []
  X = []

  with open(infile) as f:
    reader = list(csv.reader(f))

  for r in reader: X.append(r)

  kmeans = KMeans(n_clusters=k, random_state=0).fit(X)

  for i, label in enumerate(kmeans.labels_):
    list_cluster.append(label)
  return list_cluster


def confidence (infile):
  list_confidence = []
  with open (infile, 'r') as fh:
    fh.readline()
    for i in fh:
      confidence = i[0] 
      list_confidence.append(confidence)
  return list_confidence

def parse_cluster_file (k, list_cluster, list_confidence, position_max, outfile):
  Continue = 'True'
  fh_out = open (outfile, 'a')

  dict_master = {}
  count = 0
  max_nbN = 0
  for i in list_cluster:
      cluster = i
      catagory = list_confidence[count]
      
      if position_max == count: Constitutive = cluster
      count += 1
      if cluster not in dict_master.keys():
        dict_master[cluster] = [0, 0, 0] # [H, L, N]
      if catagory == 'H': dict_master[cluster][0] += 1
      elif catagory == 'L': dict_master[cluster][1] += 1
      elif catagory == 'N': dict_master[cluster][2] += 1
      else: print(catagory)
      N = dict_master[cluster][2]
      if N > max_nbN: max_nbN = N
  print('\nnbCluster = ' + str(k), file=fh_out, flush=True)
  for k, v in sorted(dict_master.items()):
    if v[2] == max_nbN:
      Basal = k
      print(str(k) + str(v) + '* basal group', file=fh_out, flush=True)
      continue
    print(str(k)+ str(v), file=fh_out, flush=True)
    if v[0] == 0 and v[1] == 0: 
      Continue = 'False'
  fh_out.close()
  return Continue, Basal, Constitutive

def max_occ_pos (infile):
  occ = []
  with open (infile, 'r') as fh:
    for i in fh:
      i = i.rstrip('\n')
      occ.append(int(i))
  max_occ = max(occ)
  
  position_max = occ.index(max_occ)
  return position_max 

##
## RUN
##

def run():

  position_max =  max_occ_pos (occSum_file)

  fh_out = open (outfile, 'w')
  print('Cluster_id [nbH, nbL, nbN]', file=fh_out, flush=True)
  fh_out.close()
  for i in range(2, k+1):
    list_cluster = kmeans_cluster (i)
    list_confidence = confidence (confidence_file)
    Continue, Basal, Constitutive = parse_cluster_file (i, list_cluster, list_confidence, position_max, outfile)
    if Continue == 'False': break
  
  nbCluster = i-1
  list_cluster = kmeans_cluster (nbCluster)
  Continue, Basal, Constitutive = parse_cluster_file (i, list_cluster, list_confidence, position_max, outfile)

  print('We divided the data into ' + str(nbCluster) + ' clusters. nbCluster = ' + str(nbCluster))
  print('See file: ' + cluster_file )
  print('The basal expression cluster is cluster #' + str(Basal))
  print('The constitutive expression cluster is cluster #' + str(Constitutive))
  print('The cluster assignment of each prediction element is in: ' + outfile )

  return nbCluster, int(Basal), int(Constitutive)




print(run())
#kmeans_cluster (10)


