# -*- coding: utf-8 -*-
"""
Created on Thu Mar 11 12:06:04 2021

@author: Allison Walker
Rank features for logistic regression
"""

import numpy as np
import random
from sklearn.linear_model import SGDClassifier
from sklearn import preprocessing
import readFeatureFiles
import argparse
import os
import sys
import joblib

parser = argparse.ArgumentParser()
parser.add_argument('model_name',help='file containing the antismash results for the cluster in a genbank file')
args = parser.parse_args()

model_name = args.model_name
data_path = os.path.dirname(sys.argv[0]) + "/"


try:
    feature_dir = data_path + "feature_matrices/" + model_name  + "/features/"   
    feature_type_list = readFeatureFiles.getFeatureFilesList(feature_dir)
    feature_type_list_lower = [f.lower() for f in feature_type_list]
    feature_type_list_sort_ind = [i[0] for i in sorted(enumerate(feature_type_list_lower), key=lambda x:x[1])]
    feature_type_list_sorted = []
    for i in feature_type_list_sort_ind:
        feature_type_list_sorted.append(feature_type_list[i])
    feature_type_list = feature_type_list_sorted
    training_features = readFeatureFiles.readFeatures(feature_dir, feature_type_list)
    feature_list = readFeatureFiles.readFeatureNames(feature_dir, feature_type_list)
    feature_list_by_type = readFeatureFiles.readFeaturesByType(feature_dir, feature_type_list)
except:
    print("ERROR: did not find file containing training data, make sure your model and associated feature data exists and is in the correct location")
    
try:
    svm_bacterial, tree_bacterial, log_bacterial = joblib.load(data_path+"trained_models/" + model_name + "_antibacterial.sav")
    
    svm_antieuk, tree_antieuk, log_antieuk = joblib.load(data_path+"trained_models/" + model_name + "_antieuk.sav")
    
    svm_antifungal, tree_antifungal, log_antifungal = joblib.load(data_path+"trained_models/" + model_name + "_antifungal.sav")
    
    svm_antitumor, tree_antitumor, log_antitumor = joblib.load(data_path+"trained_models/" + model_name + "_cytotoxic_antitumor.sav")
    
    svm_antigramneg, tree_antigramneg, log_antigramneg = joblib.load(data_path+"trained_models/" + model_name + "_antigramneg.sav")

    svm_antigrampos, tree_antigrampos, log_antigrampos = joblib.load(data_path+"trained_models/" + model_name + "_antigrampos.sav")
except:
   print("ERROR: could not find pretrained model, make sure all data files are in correct location and your model exists")
   exit() 


#get and sort coefficients and write results to file
log_coeffs = {}
log_coeffs["antibacterial"] = log_bacterial.coef_[0]
log_coeffs["antieuk"] = log_antieuk.coef_[0]
log_coeffs["antifungal"] = log_antifungal.coef_[0]
log_coeffs["antitumor"] = log_antitumor.coef_[0]
log_coeffs["antigramneg"] = log_antigramneg.coef_[0]
log_coeffs["antigrampos"] = log_antigrampos.coef_[0]


for classifier in log_coeffs:
    log_coeff = log_coeffs[classifier]
    index_log_coeff = np.argsort(-1*log_coeff)
    output_log = open("feature_scores/" + model_name + "_"+classifier+"_log_coeff.txt",'w')
    for i in index_log_coeff:
        output_log.write(feature_list[i].replace(",","_")+",")
        output_log.write(str(log_coeff[i])+"\n")
    output_log.close()
