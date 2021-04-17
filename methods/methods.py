#!/usr/bin/env python -W ignore::DeprecationWarning
import os
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import MinMaxScaler
from sklearn.covariance import EllipticEnvelope
from sklearn.ensemble import IsolationForest
from sklearn.neighbors import LocalOutlierFactor
from sklearn.svm import OneClassSVM
from ast import literal_eval
import warnings

def principal_components(all_vectors_file):
    warnings.filterwarnings("ignore")
    f=open(all_vectors_file,'r')
    lines=f.readlines()
    sequences=[]
    vectors=[]
    for line in lines: #get every sequence and vector
        sequence=line.split("\t")[0]
        vector=line.split("\t")[1]
        sequences.append(sequence)
        vector=literal_eval(vector)
        vectors.append(vector)
    f.close()
    scale=MinMaxScaler()
    vectors=scale.fit_transform(vectors)
    print(vectors)
    model=PCA(n_components=0.95,svd_solver='full') #PCA for 95% explained variance
    model.fit(vectors)
    pca_vectors=model.transform(vectors)
    return pca_vectors

def novelty(all_vectors_file,labels_file,variant_name):
    warnings.filterwarnings("ignore")
    f=open(all_vectors_file,'r')
    lines=f.readlines()
    sequences=[]
    vectors=[]
    for line in lines: #get every sequence and vector
        sequence=line.split("\t")[0]
        vector=line.split("\t")[1]
        sequences.append(sequence)
        vector=literal_eval(vector)
        vectors.append(vector)
    f.close()
    #print("Before PCA")
    #vectors=principal_components(all_vectors_file)
    #print("After PCA")
    f2=open(labels_file,'r')
    lines=f2.readlines()
    labels=[]
    for line in lines: #get every label
        label=line
        label=label.strip('\n')
        labels.append(label)
    f2.close()
    split=labels.index(variant_name)
    train_x=vectors[:split]
    train_y=labels[:split]
    test_x=vectors[split:]
    test_y=labels[split:]
    variant_test_pos=[]
    for i in range(len(test_y)):
        if test_y[i]==variant_name:
            variant_test_pos.append(i)
    model1=EllipticEnvelope()
    model1.fit(train_x)
    pred_y_model1=model1.predict(test_x)
    model2=IsolationForest()
    model2.fit(train_x)
    pred_y_model2=model2.predict(test_x)
    model3=LocalOutlierFactor(novelty=True)
    model3.fit(train_x)
    pred_y_model3=model3.predict(test_x)
    model4=OneClassSVM()
    model4.fit(train_x)
    pred_y_model4=model4.predict(test_x)
    cnt1=0
    cnt2=0
    cnt3=0
    cnt4=0
    total_variant_name=len(variant_test_pos)
    print(total_variant_name)
    print(len(vectors))
    print(len(labels))
    print(len(test_y))
    print(len(pred_y_model1))
    print(len(pred_y_model2))
    print(len(pred_y_model3))
    print(len(pred_y_model4))
    #print(variant_test_pos)
    for i in range(total_variant_name):
        if pred_y_model1[variant_test_pos[i]]==-1:    #detects novelty
            cnt1=cnt1+1
        if pred_y_model2[variant_test_pos[i]]==-1:    #detects novelty
            cnt2=cnt2+1
        if pred_y_model3[variant_test_pos[i]]==-1:    #detects novelty
            cnt3=cnt3+1
        if pred_y_model4[variant_test_pos[i]]==-1:    #detects novelty
            cnt4=cnt4+1
    print(cnt1/total_variant_name)
    print(cnt2/total_variant_name)
    print(cnt3/total_variant_name)
    print(cnt4/total_variant_name)
