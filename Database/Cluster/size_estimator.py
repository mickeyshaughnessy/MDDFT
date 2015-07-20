#!/usr/bin/python

"""
This script estimates the required database size and cluster size required to achieve
a given level of accuracy from a DFT model. How dense must the points in a database be to give reliable interpolation (from a given DFT system)?

input: a file called DB_VASP.dat
output: a text file output_SizeEstimator
usage: ./Size_estimator.py 

psuedocode:
	read in whole database
	for each database size
		select N search clusters
		for each search cluster
			find the closest match(es)
			estimate the force and compute the Delta_F
		print database size report - average error, average RMSD distance
		

"""
import fileinput
import sys
import numpy as np
import math
import array
import random
from operator import itemgetter

##############  global data ########################################
N_clusters = 1000  #overall size of database
N_pairs = N_clusters*(N_clusters-1)/2
N_neighbs = 6 #size of clusters
N_points = 20 
#DBsizes = [10,100,1000,10000,100000,1000000]
DBsizes = list(xrange(20))
del DBsizes[0]
Query_cluster = []
Clusters = [] # list, id --> cluster
Pair_distances = [] # list, pair --> distance
max_d = math.sqrt(3.*N_neighbs)
big = max_d*10

###################################################################
def computeRMSD(a, b):
  a = np.array(a)
  b = np.array(b)
  n_vec = np.shape(a)[0]
  correlation_matrix = np.dot(np.transpose(a), b)
  v,s,w = np.linalg.svd(correlation_matrix)
  is_reflection = (np.linalg.det(v) * np.linalg.det(w)) < 0.0
  if is_reflection: s[-1] = - s[-1]
  E0 = sum(sum(a * a)) + sum(sum(b * b))
  rmsd_sq = (E0 - 2.0*sum(s)) / float(n_vec)
  rmsd_sq = max([rmsd_sq, 0.0])
  return np.sqrt(rmsd_sq)

###################################################################
def compute_distances(ref_ids,query_clust):
  d_min = big
  id_min = -1
  for ref in ref_ids:
    d = computeRMSD(Clusters[ref],query_clust)
    if (d < d_min):
      id_min = ref
      d_min = d
  return [id_min,d_min]

###################################################################
def get_new_refs(query_clust, d_min):
  r = []
  for i in range(len(Pair_distances)):
    p = pair(i)
    d = Pair_distances[i]
    if   ((p[1] == min_id) and d < 2*d_min):  r.append(p[0])
    elif ((p[0] == min_id) and d < 2*d_min):  r.append(p[1])
  unique_refs = uniquify(r)
  return unique_refs

###################################################################
def initial_refs(type="random"):
  refs0 = []
  if (type == "random"):
    for j in range(0,N_refs):
      refs0.append(random.randint(0,N_clusters-1))
  elif (type == "most_distant"):
    refs0 = most_distant(N_refs)
  elif (type == "uniformly_distant"):
    refs0 = winnow(N_refs)
  return refs0


#####################################################################
## MAIN
#####################################################################
print "> reading in database"
DB = open('DB_VASP.dat', 'r')
FullClusters=[]
Clusters=[]
rawdat =[]
DBlinecount = 0
for line in DB.readlines():
	rawdat.append(line)
np.random.shuffle(rawdat)
print "> full database shuffled"
print "> linecount:"
for line in rawdat:
	DBlinecount+=1
	if (DBlinecount%100 == 1):
		print DBlinecount	
	words = line.split()
	cluster = map(float, words[5:])
	xyz = []
	xs = []
	ys = []
	zs = []
	count = 0
	for coor in cluster:
		count+=1
		if (count%3 == 1): #x
			xs.append(coor)
		elif (count%3 == 2): #x
			ys.append(coor)
		elif (count%3 == 0): #x
			zs.append(coor)
	xyz.append(xs)
	xyz.append(ys)
	xyz.append(zs)
	Clusters.append(xyz)
	fullcluster = map(float, words)
	FullClusters.append(fullcluster)
#### Database read in ######
print "> read database"

	
for size in DBsizes:   #Loop over fractions of the total DFT dataset to scan
	if (size < DBlinecount):
		cur_Clusters = Clusters[0:size]	
		avgDF=0
		Search_clusters = map(int,size*np.random.rand(int(size*0.1))) #vector of random integers between 1 and size, of length 10% of size
		print "Search targets: %s" % (Search_clusters)
		for search in Search_clusters:	
			count = 0
			dmin = 100000
			for clust in cur_Clusters:				
				count+=1		
				d = computeRMSD(cur_Clusters[search],clust)
				if (d < dmin and count != search and d != 0):
					min_id = count
					dmin = d
					print "> dmin= %s, min_id= %s, search_id= %s " % (dmin, min_id, search) 
			DeltaFx = FullClusters[search][2]-FullClusters[min_id][2] 
			DeltaFy = FullClusters[search][3]-FullClusters[min_id][3] 
			DeltaFz = FullClusters[search][4]-FullClusters[min_id][4]
			magDF = math.sqrt(DeltaFx**2+DeltaFy**2+DeltaFz**2)
			print "> dmin= %s, min_id= %s, search_id= %s DeltaF = %s" % (dmin, min_id, search, magDF) 
			avgDF+=magDF
		avgDF = avgDF/(size*0.1)
		print "> DB Size: %s, Avg DeltaF = %s" % (size, avgDF) 
