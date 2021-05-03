# -*- coding: utf-8 -*-
"""
Created on Thu Mar 11 12:06:04 2021

@author: allis
"""

import numpy as np
import random
from sklearn.linear_model import SGDClassifier
from sklearn import preprocessing
import readFeatureFiles

#classification to rank features for
classification = "antibacterial"
#set random seed so results are consistent
random.seed(1)
#include features from SSN 
include_SSN = True

#parameters for classifiers
log_params = {}
log_params["antibacterial"] = {"l1_ratio":.05,"alpha":.01}
log_params["antigrampos"] = {"l1_ratio":.001,"alpha":.001}
log_params["antigramneg"] = {"l1_ratio":.05,"alpha":.001}
log_params["antieuk"] = {"l1_ratio":.001,"alpha":.001}
log_params["antifungal"] = {"l1_ratio":.0001,"alpha":.01}
log_params["cytotoxic_antitumor"] = {"l1_ratio":.001,"alpha":.001}


#read in features
pfam_features = readFeatureFiles.readFeatureMatrix("feature_matrices/PFAM5.csv")
card_features = readFeatureFiles.readFeatureMatrix("feature_matrices/CARD5_genes.csv")
smCOG_features = readFeatureFiles.readFeatureMatrix("feature_matrices/SMCOG5.csv")
SSN_features = readFeatureFiles.readFeatureMatrix("feature_matrices/SSN.csv")
CDS_features = readFeatureFiles.readFeatureMatrix("feature_matrices/CDS_motifs5.csv")


pk_consensus_features = readFeatureFiles.readFeatureMatrix("feature_matrices/pk_nrp_consensus5.csv")

#read in list of feature lables
pfam_list = readFeatureFiles.readFeatureList("feature_matrices/pfam_list5.txt")
card_list = readFeatureFiles.readFeatureList("feature_matrices/CARD5_gene_list.txt")
smCOG_list = readFeatureFiles.readFeatureList("feature_matrices/SMCOG_list5.txt")
SSN_list = readFeatureFiles.readFeatureList("feature_matrices/SSN_list.txt")
CDS_list = readFeatureFiles.readFeatureList("feature_matrices/CDS_motifs_list5.txt")

pk_consensus_list = readFeatureFiles.readFeatureList("feature_matrices/pk_nrp_consensus_list5.txt")

#concatenate features
features = np.concatenate((pfam_features, card_features), axis=1)
features = np.concatenate((features,  smCOG_features), axis=1)
features = np.concatenate((features,  CDS_features), axis=1)
if include_SSN:
    features = np.concatenate((features,  SSN_features), axis=1)
features = np.concatenate((features,  pk_consensus_features), axis=1)


feature_list = pfam_list +card_list+smCOG_list+CDS_list+SSN_list+pk_consensus_list


is_antibacterial = readFeatureFiles.readClassesMatrix("feature_matrices/is_antibacterial.csv")
is_antifungal = readFeatureFiles.readClassesMatrix("feature_matrices/is_antifungal.csv")
is_cytotoxic = readFeatureFiles.readClassesMatrix("feature_matrices/is_cytotoxic.csv")
is_unknown = readFeatureFiles.readClassesMatrix("feature_matrices/is_unknown.csv")
targets_gram_pos = readFeatureFiles.readClassesMatrix("feature_matrices/targets_gram_pos.csv")
targets_gram_neg = readFeatureFiles.readClassesMatrix("feature_matrices/targets_gram_neg.csv")
full_cluster_list = readFeatureFiles.readClusterList("feature_matrices/cluster_list_CARD.txt")
is_not_unknown_indices = readFeatureFiles.getNotUnknownIndices(is_unknown)
target_unannotated = is_antibacterial*((targets_gram_pos+targets_gram_neg)<1)
is_not_unknown_indices_gram =  readFeatureFiles.getNotUnknownIndices(is_unknown + target_unannotated)

is_antibacterial = (is_antibacterial >= 1).astype(int)
is_antieuk = ((is_antifungal + is_cytotoxic)>=1).astype(int)
is_gram_pos = (targets_gram_pos >= 1).astype(int)
is_gram_neg = (targets_gram_neg >= 1).astype(int)

y_vars = []
if classification == "antibacterial":
    y_vars = is_antibacterial
    y_vars = y_vars[is_not_unknown_indices]
    features = features[is_not_unknown_indices,:]
    
if classification == "antieuk":
    y_vars = is_antieuk
    y_vars = y_vars[is_not_unknown_indices]
    features = features[is_not_unknown_indices,:]

    
if classification == "antifungal":
    y_vars = (is_antifungal >= 1).astype(int)
    y_vars = y_vars[is_not_unknown_indices]
    features = features[is_not_unknown_indices,:]
    
if classification == "cytotoxic_antitumor":
    y_vars = (is_cytotoxic >= 1).astype(int)
    y_vars = y_vars[is_not_unknown_indices]
    features = features[is_not_unknown_indices,:]
    
if classification == "antigramneg":
    y_vars = is_gram_neg
    y_vars = y_vars[is_not_unknown_indices_gram]
    features = features[is_not_unknown_indices_gram,:]
    
if classification == "antigrampos":
    y_vars = is_gram_pos
    y_vars = y_vars[is_not_unknown_indices_gram]
    features = features[is_not_unknown_indices_gram,:]

#train logistic regression classifier    
min_max_scaler = preprocessing.MinMaxScaler()
scaled_training_x = min_max_scaler.fit_transform(features)
opt_log_params = log_params[classification]
log_classifier = SGDClassifier(loss='log',penalty='elasticnet',max_iter=100,alpha=opt_log_params["alpha"],l1_ratio=opt_log_params["l1_ratio"],tol=None)
log_classifier.fit(scaled_training_x, y_vars)

#get and sort coefficients and write results to file
log_coeff = log_classifier.coef_[0]
index_log_coeff = np.argsort(-1*log_coeff)
output_log = open("feature_scores/as5_"+classification+"_log_coeff.txt",'w')
for i in index_log_coeff:
    output_log.write(feature_list[i].replace(",","_")+",")
    output_log.write(str(log_coeff[i])+"\n")
output_log.close()

#print top and bottom 10 classifiers and coefficients
for i in range(0, 10):
    print(feature_list[index_log_coeff[i]] + ": " + str(log_coeff[index_log_coeff[i]]))
print("")                                                                
for i in range(0, 10):
    print(feature_list[index_log_coeff[len(index_log_coeff)-(i+1)]] + ": " + str(log_coeff[index_log_coeff[len(index_log_coeff)-(i+1)]]))                                                      