from dataset import preprocessing
import collections
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats
import time
from sklearn.cluster import KMeans
import datetime

t1=time.perf_counter()

ids,seqs=preprocessing.sequences('spikeprot0904/spikeprot0904.fasta') #all spike proteins

proteins,vectors=preprocessing.prot_vec('protVec_100d_3grams.csv') #proteins and vectors
    
df3=pd.read_csv('metadata_processed.tsv',sep='\t') #subset of spike proteins

aligned_sequences="aligned_sequences.fasta"

ids2,seqs2=preprocessing.sequences(aligned_sequences)

sequence_vectors_dir="sequence_vectors/"

dates_all=[]
ids_all=[]
seqs_all=[]

cnt3=0
for i in range(len(ids)):
    try:
        dates_all.append(ids[i].split("|")[2]) #collecting all dates
        ids_all.append(ids[i].split("|")[3]) #collecting all IDs
        seqs_all.append(seqs[i]) #collecting all sequences
        cnt3=cnt3+1
    except:
        continue

cnt1=0

dates_subset=[]
ids_subset=[]
seqs_subset=[]


for i in range(len(ids_all)):
    if ids_all[i] in df3['gisaid_epi_isl'].values:
        ids_subset.append(ids_all[i])
        dates_subset.append(dates_all[i])
        seqs_subset.append(seqs_all[i])
        cnt1=cnt1+1

print(cnt1)

seqs_subset=[]
cnt2=0

for i in range(len(ids_subset)):
    if ids2[i]==ids_subset[i]:
        seqs_subset.append(seqs2[i])
        cnt2=cnt2+1

print(cnt2)

cnt=0

for i in range(len(ids_subset)):
    if ids_subset[i] in df3['gisaid_epi_isl'].values:
        cnt=cnt+1
        print(cnt)
        m=preprocessing.sequence_to_vec(seqs_subset[i],proteins,vectors) #converting full sequence to vector
        f=open(sequence_vectors_dir+str(dates_subset[i])+".txt",'a')
        f.write(str(seqs_subset[i])+"\t"+str(m)+"\n")
        f.close()


print(cnt)

t2=time.perf_counter()

total_time=t2-t1

print("Time : ",total_time)
