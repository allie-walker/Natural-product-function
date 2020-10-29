# -*- coding: utf-8 -*-
"""
Created on Fri Aug 17 14:14:27 2018

@author: allis
"""

import numpy as np
from sklearn import linear_model
import random

def readFeatureMatrix(filename):
    in_file = open(filename,'r')
    feature_matrix = []
    for line in in_file:
        if "," not in line:
            continue
        entries = line.split(",")
        feature_matrix.append([])
        for i in range(0,len(entries)-1):
            feature_matrix[len(feature_matrix)-1].append(int(entries[i]))
    in_file.close()
    return np.array(feature_matrix)

def readFeatureMatrixFloat(filename):
    in_file = open(filename,'r')
    feature_matrix = []
    for line in in_file:
        if "," not in line:
            continue
        entries = line.split(",")
        feature_matrix.append([])
        for i in range(0,len(entries)-1):
            feature_matrix[len(feature_matrix)-1].append(float(entries[i]))
    in_file.close()
    return np.array(feature_matrix)
            
def readFeatureList(filename):
    in_file = open(filename,'r')
    feature_list = []
    for line in in_file:
        feature_list.append(line.replace("\n","").replace("\r",""))
    return feature_list

def readClassesMatrix(filename):
    class_matrix = []
    in_file = open(filename,'r')
    for line in in_file:
        if len(line)>0:
            class_matrix.append(int(line))
    in_file.close()
    return np.array(class_matrix)     
       
def getNotUnknownIndices(is_unknown):
    is_not_unknown_indices = []
    for i in range(0, len(is_unknown)):
        if not is_unknown[i]:
            is_not_unknown_indices.append(i)
    return is_not_unknown_indices
    
def writeROC(probs, validation_y):
    outfile_prob= open("classifier_results/logistic_antibiotic_prob_vs_tp3.txt",'w')
    outfile_roc= open("classifier_results/logistic_antibiotic_roc3.txt",'w')
    outfile_raw_count= open("classifier_results/logistic_antibiotic_raw_count3.txt",'w')
    is_class_prob = probs[:,1]
    for i in range(0, 21):
        cutoff_prob = i*.05
        classifications = (is_class_prob>cutoff_prob).astype(int)
        #print classifications
        true_positives = classifications*validation_y
        false_positives = (1-validation_y)*classifications
        tp_rate = (true_positives.sum()*1.0)/validation_y.sum()
        fp_rate = (false_positives.sum()*1.0)/(1-validation_y).sum()
        outfile_prob.write(str(1-i)+","+str(tp_rate)+"\n")
        outfile_roc.write(str(fp_rate)+","+str(tp_rate)+"\n")
        outfile_raw_count.write(str(false_positives.sum()) +"," + str(true_positives.sum()) + "\n")
    outfile_prob.close()
    outfile_roc.close()
    outfile_raw_count.close()
    
def readClusterList(filename):
    infile = open(filename,'r')
    cluster_list = []
    for line in infile:
        cluster_list.append(line.replace("\n",""))
    infile.close()
    return cluster_list

def getFeatureName(i, pfam_list, card_list, smCOG_list):
    merged_list = pfam_list+card_list+smCOG_list
    return merged_list[i]
    if i < len(pfam_list):
         #print(pfam_list[-3])
         return pfam_list[i]
    elif i < len(pfam_list)+len(card_list):
        return card_list[i-len(pfam_list)]
    else:
        return smCOG_list[i-len(pfam_list)-len(card_list)]
                              
#def normalizeFeatures(feature_matrix, method):
    #if method == "MIN
    #return