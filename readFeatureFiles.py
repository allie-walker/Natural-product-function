# -*- coding: utf-8 -*-
"""
Created on Fri Aug 17 14:14:27 2018

@author: allison
"""

import numpy as np
from sklearn import linear_model
import random
import os
import sys


def getFeatureFilesList(dir_name):
    """Get a list of all the feature types in the feature directory based on csv files
    """
    try:
        file_list = os.listdir(dir_name)
    except:
        sys.exit("Feature directory not found")
    
    #get all csv
    csv_files = [f for f in file_list if ".csv" in f]
    
    if len(csv_files) == 0:
        sys.exit("No csv files found in directory")
   
    feature_type_list = [f[0:f.find(".csv")] for f in csv_files]
    
    return feature_type_list

def readFeatures(dir_name, feature_type_list):
    """Reads csv files from directory in order of feature_type_list
    """
    i = 0
    for f in feature_type_list:
        in_file = f + ".csv"
        current_matrix = readFeatureMatrix(dir_name + in_file)
        if i == 0:
            features = current_matrix
            i +=1
        else:
            features = np.concatenate((features, current_matrix), axis=1)
    return features

def readFeatureNames(dir_name, feature_type_list):
    """Reads feature names from text files in directory in order of feature_type_list
    """
    i = 0
    for f in feature_type_list:
        in_file = f + "_list.txt"
        try:
            current_list = readFeatureList(dir_name + in_file)
        except:
            sys.exit("No matching label file found for " + f + ".csv, this label file should be called " + in_file)
        if i == 0:
            feature_list = current_list
        else:
            feature_list = feature_list + current_list
        i += 1
    return feature_list

def readFeaturesByType(dir_name, feature_type_list):
    """Reads features and returns a dictionary matching type to list of features
    """
    featureTypeDict = {}
    for f in feature_type_list:
        in_file = f + "_list.txt"
        try:
            current_list = readFeatureList(dir_name + in_file)
        except:
            sys.exit("No matching label file found for " + f + ".csv, this label file should be called " + in_file)
        featureTypeDict[f] = current_list
    return featureTypeDict

#TODO: write this method
def readClassifications(dir_name):
    return


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