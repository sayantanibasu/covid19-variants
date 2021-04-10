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

ids,seqs=preprocessing.sequences('spikeprot0407/spikeprot0407.fasta') #all spike proteins

proteins,vectors=preprocessing.prot_vec('protVec_100d_3grams.csv') #proteins and vectors
    
df3=pd.read_csv('nextstrain_ncov_global_metadata.tsv',sep='\t') #subset of spike proteins

unaligned_sequences="unaligned_sequences.fasta"

#aligned_sequences="aligned_sequences.fasta"

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

f2=open(unaligned_sequences, "w")
for i in range(len(ids_all)):
    if ids_all[i] in df3['gisaid_epi_isl'].values:
        f2.write(">"+ids_all[i]+"\n"+seqs_all[i] +"\n")
        ids_subset.append(ids_all[i])
        dates_subset.append(dates_all[i])
        cnt1=cnt1+1
f2.close()

print(cnt1)

_,seqs=preprocessing.sequences(unaligned_sequences) #unaligned protein sequences

seqs_subset=[]

print("len sequences")
print(len(seqs))
print(len(dates_subset))

for i in range(len(ids_subset)):
        seqs_all.append(seqs[i]) #collecting all sequences

cnt=0

print(seqs_all)
date_sequence_vectors=[]
seq_sequence_vectors=[]
sequence_vectors=[]

for i in range(len(ids_subset)):
    if ids_subset[i] in df3['gisaid_epi_isl'].values:
        m=preprocessing.sequence_to_vec(seqs_all[i],proteins,vectors) #converting full sequence to vector
        sequence_vectors.append(m)
        date_sequence_vectors.append(dates_subset[i]) #collecting corresponding date
        seq_sequence_vectors.append(seqs_all[i]) #collecting corresponding seqs
        cnt=cnt+1
        print(cnt,dates_subset[i])

dates_unique=list(df3['Collection Data'].unique()) #collecting unique dates for clustering

cnt=0

for i in range(len(dates_unique)):
    print(dates_unique[i])
    for j in range(len(date_sequence_vectors)):
        if date_sequence_vectors[j]==datetime.datetime.strptime(dates_unique[i],'%m/%d/%y').strftime('%Y-%m-%d'):
            cnt=cnt+1
            f=open(sequence_vectors_dir+str(date_sequence_vectors[j])+".txt",'a')
            f.write(str(seq_sequence_vectors[j])+"\t"+str(sequence_vectors[j])+"\n")
            f.close()

t2=time.perf_counter()

total_time=t2-t1

print("Time : ",total_time)

print(cnt)
