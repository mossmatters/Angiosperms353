#!/usr/bin/env python3

#This script implements the clustering techniques of my SequenceClusters jupyter notebook.

import sys

#These lines prevent the rocket ship from opening even when we're not doing anything graphical
#import matplotlib
#matplotlib.use("Agg")

from skbio import TabularMSA, DNA, DistanceMatrix
from skbio.sequence.distance import hamming
import pandas as pd
import numpy as np
#import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from scipy.spatial import distance
import random,os,itertools

def gap_dectector(sequence_column):
	'''Returns the number of gap characters in a column of a sequence matrix'''
	#In skbio, DNA sequences are stored as bytecode, (b'A') so need to convert back to strings
	try:
		return sequence_column.value_counts()[b"-"]
	except KeyError:
		return 0

def remove_gapped_columns(msa,site_threshold=0.95):
	msa_dict = msa.to_dict()
	msa_df = pd.DataFrame(msa_dict)
	gapped_columns = msa_df.apply(gap_dectector ,axis=1)
	nogaps_df = msa_df[gapped_columns < len(msa_df.columns) * site_threshold]
	nogap_seqs = [DNA(nogaps_df[i].str.decode("utf-8").str.cat(), metadata = {"id":i}) for i in nogaps_df]
	msa_nogap = TabularMSA(nogap_seqs)
	return msa_nogap
	
def get_reduced_alignment(fasta_filename,angio_1kp_ids,site_threshold=0.95,sample_threshold = 0.5,write_alignment=True):
	"""Given a 1KP alignment that contains all sequences, return one with only angiosperms 
	that has been trimmed to remove sites with fewer than X sequences, and to remove
	sequences that have less than Y percentage of the remaining alignment length."""
	msa = TabularMSA.read(fasta_filename, constructor=DNA)
	seqs_to_keep = []
	for seq in msa:
		if seq.metadata["id"] in angio_1kp_ids:
			seqs_to_keep.append(seq)
		
	angio_msa = TabularMSA(seqs_to_keep)        
	angio_msa.reassign_index(minter="id")
	sys.stderr.write("After removing non-angiosperms: {}\n".format(angio_msa.shape))

	angio_msa_nogap = remove_gapped_columns(angio_msa)

	sys.stderr.write("After removing gappy sites: {}\n".format(angio_msa_nogap.shape))

	seqs_to_keep = []
	for seq in angio_msa_nogap:
		num_gaps = len([x for x in seq.gaps() if x])
		if num_gaps < angio_msa_nogap.shape[1] * sample_threshold:
			seqs_to_keep.append(seq)
		
	angio_msa_nogap_noshort = TabularMSA(seqs_to_keep)
	angio_msa_nogap_noshort.reassign_index(minter="id")
	if write_alignment:
		angio_msa_nogap_noshort.write("onekp_only_angios_degapped/{}.onlyangios.noshort.fasta".format(gene))
	sys.stderr.write("After removing gappy sequences: {}\n".format(angio_msa_nogap_noshort.shape))
	return angio_msa_nogap

def p_distance(seq1,seq2):
	'''Modified hamming distance to include only non-gap sites'''
	from skbio.sequence import Sequence
	from numpy import isnan
	myseq1 = str(seq1)
	myseq2 = str(seq2)

	degapped1 = []
	degapped2 = []

	for i in range(len(myseq1)):
		if myseq1[i] != "-":
			if myseq2[i] != "-":
				degapped1.append(myseq1[i])
				degapped2.append(myseq2[i])
	degapped1 = "".join(degapped1)
	degapped2 = "".join(degapped2)

	#sys.stderr.write(degapped1)
	#sys.stderr.write(degapped2)

	hamming_dist = hamming(Sequence(degapped1),Sequence(degapped2))
	#sys.stderr.write(hamming_dist)
	if isnan(hamming_dist):
		#sys.stderr.write(seq1.metadata["id"], seq2.metadata["id"])
		return 0.0
	else:
		return hamming_dist

def kMedoids(D, k, tmax=100):
	'''Code from: https://github.com/letiantian/kmedoids'''
	# determine dimensions of distance matrix D
	m, n = D.shape

	# randomly initialize an array of k medoid indices
	M = np.sort(np.random.choice(n, k))

	# create a copy of the array of medoid indices
	Mnew = np.copy(M)

	# initialize a dictionary to represent clusters
	C = {}
	for t in range(tmax):
		# determine clusters, i. e. arrays of data indices
		J = np.argmin(D[:,M], axis=1)
		for kappa in range(k):
			C[kappa] = np.where(J==kappa)[0]
		# update cluster medoids
		for kappa in range(k):
			J = np.mean(D[np.ix_(C[kappa],C[kappa])],axis=1)
			j = np.argmin(J)
			Mnew[kappa] = C[kappa][j]
		np.sort(Mnew)
		# check for convergence
		if np.array_equal(M, Mnew):
			break
		M = np.copy(Mnew)
	else:
		# final update of cluster memberships
		J = np.argmin(D[:,M], axis=1)
		for kappa in range(k):
			C[kappa] = np.where(J==kappa)[0]

	# return results
	return M, C


def reduced_alignment(fasta_filename,medoids,gene):
	genomes = ["Ambtr_v1.0.27","Orysa_v7.0","Arath_TAIR10"]
	medoids = medoids + genomes
	msa = TabularMSA.read(fasta_filename,format = 'fasta',constructor=DNA)
	seqs_to_keep = []
	#print(medoids)
	for seq in msa:
		#print(seq.metadata["id"])
		if seq.metadata["id"] in medoids:
			seqs_to_keep.append(seq)
	medoid_msa = TabularMSA(seqs_to_keep)        
	medoid_msa.reassign_index(minter="id")
	medoid_msa_nogap = remove_gapped_columns(medoid_msa)
	medoid_msa_nogap.write("medoid_alignments/{}_medoids.fasta".format(gene))
	sys.stderr.write("Shape of degapped medoid alignment: {}\n".format(medoid_msa_nogap.shape))

gene = sys.argv[1]
sys.stderr.write(gene+'\n')
fasta_filename = "genes/{}/FNA2AA-upp-masked.fasta".format(gene)
angiosperm_id_fn = "1kp_angio_codes.txt"
angio_1kp_ids = set([x.rstrip() for x in open(angiosperm_id_fn)])
distance_matrix_fn = "onekp_only_angios_pdistance/{}_angio_P_dm.csv".format(gene)
degapped_alignment_fn = "onekp_only_angios_degapped/{}.onlyangios.noshort.fasta".format(gene)

#Read in sequence alignment and trim three ways: remove non angiosperms, remove gappy sites, remove gappy sequences


if os.path.isfile(degapped_alignment_fn):
	angio_msa_nogap_noshort = TabularMSA.read(degapped_alignment_fn,constructor=DNA)
	sys.stderr.write("Read in degapped alignment: {}\n".format(angio_msa_nogap_noshort.shape))
else: 
	angio_msa_nogap_noshort = get_reduced_alignment("genes/{}/FNA2AA-upp-masked.fasta".format(gene),angio_1kp_ids)



if os.path.isfile(distance_matrix_fn):
	p_dm = DistanceMatrix.read(distance_matrix_fn)    
	p_dm_df = p_dm.to_data_frame()
	sys.stderr.write("Read in pre-determined distance matrix!\n")     
else:
	p_dm = DistanceMatrix.from_iterable(angio_msa_nogap_noshort,metric=p_distance,key="id")
	p_dm_df = p_dm.to_data_frame()
	p_dm_df.to_csv("onekp_only_angios_pdistance/{}_angio_p_dm.csv".format(gene))

    

# Cluster sequences

divergent_seqs_medoids = []
runs = {}
best_run = len(p_dm_df)
best_run_idx = (6,0)
for k,i in itertools.product(range(6,16),range(100)):
	try:
		medoids,membership = kMedoids(p_dm,k)
		medoid_dist = p_dm_df[p_dm_df.ix[medoids].index].apply(min,1)
		num_over_25 = len(medoid_dist[medoid_dist > 0.25])
#			divergent_seqs_medoids.append((k,num_over_25))
		runs[(k,i)] = (medoids,membership,medoid_dist)
	except ValueError:
#			divergent_seqs_medoids.append((k,np.nan))
		num_over_25 = np.nan
		runs[(k,i)] = (0,0,0)
	if num_over_25 < best_run:
		best_run = num_over_25
		best_run_idx = (k,i)
		print(k,i,best_run)
	medoids,membership,medoid_dist = runs[best_run_idx]
# 	best_num_over_30 = len(medoid_dist[medoid_dist > 0.30])
# 	divergent_seqs_medoids.append(best_num_over_30)
	if num_over_25 == 0:
		break
	if k > 9 and num_over_25 < len(p_dm_df) * 0.03:
		break
	
sys.stderr.write("{}\t{}\t{}\n".format(gene,p_dm_df.shape[1],len(medoids)))        
         

#sys.stdout.write("\nNumber of Distances > 30%: {}\n".format())

best_medoids = p_dm_df.ix[medoids].index
with open("best_medoids.txt",'a') as outfile:
	outfile.write("{}\t{}\t{}\t{}\t{}\n".format(gene,p_dm_df.shape[1],len(medoids),len(medoid_dist[medoid_dist > 0.30]),"\t".join(best_medoids)))

reduced_alignment(fasta_filename,list(best_medoids),gene)

