# -*- coding: utf-8 -*-
"""
Created on Tue Jun  3 16:10:12 2025

@author: Allison Walker
Ranks genes in a BGC by the total contribution of the features within those genes
"""
#TODO: add support for antiSMASH 6+

import argparse
from Bio import SeqIO
import readFeatureFiles
import joblib
import numpy as np
from Bio import SeqIO
import readInputFiles

parser = argparse.ArgumentParser()
parser.add_argument('antismash_file',help='path to genbank for antiMSASH output')
parser.add_argument('antismash_version',type=int,choices=[5,6,7,8], help='path to genbank for antiMSASH output')
parser.add_argument('model_name',help='Model name, there should be a correponding pretrained model in the trained_models directory and features in the feature_matrices directory')
parser.add_argument('-r','--rgi_file',help='path to txt for RGI output, needed if you are using a model with RGI features')
parser.add_argument('-rv','--rgi_version',type=int,choices=[None,5,6], help='version of rgi used to generate antismash input file, supported versions are 3 and 5, default 5') 
parser.add_argument('-o','--outdir',default= "./",help="output directory")
args = parser.parse_args()

antismash_file = args.antismash_file
antismash_version = args.antismash_version
model_name = args.model_name
rgi_file = args.rgi_file
rgi_version = args.rgi_version
outdir = args.outdir
#read feature list, models and coefficients
try:
    feature_dir = "feature_matrices/" + model_name  + "/features/"   
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
    exit()
    
try:
    svm_bacterial, tree_bacterial, log_bacterial = joblib.load("trained_models/" + model_name + "_antibacterial.sav")
    svm_antieuk, tree_antieuk, log_antieuk = joblib.load("trained_models/" + model_name + "_antieuk.sav")
    svm_antifungal, tree_antifungal, log_antifungal = joblib.load("trained_models/" + model_name + "_antifungal.sav")
    svm_antitumor, tree_antitumor, log_antitumor = joblib.load("trained_models/" + model_name + "_cytotoxic_antitumor.sav")
    svm_antigramneg, tree_antigramneg, log_antigramneg = joblib.load("trained_models/" + model_name + "_antigramneg.sav")
    svm_antigrampos, tree_antigrampos, log_antigrampos = joblib.load("trained_models/" + model_name + "_antigrampos.sav")
except:
   print("ERROR: could not find pretrained model, make sure all data files are in correct location and your model exists")
   exit() 
   
log_coeffs = {}
log_coeffs["antibacterial"] = log_bacterial["log"].coef_[0]
log_coeffs["antieuk"] = log_antieuk["log"].coef_[0]
log_coeffs["antifungal"] = log_antifungal["log"].coef_[0]
log_coeffs["antitumor"] = log_antitumor["log"].coef_[0]
log_coeffs["antigramneg"] = log_antigramneg["log"].coef_[0]
log_coeffs["antigrampos"] = log_antigrampos["log"].coef_[0]

feat_coeffs = {}
for classifier in log_coeffs:
    feat_coeffs[classifier] = {}
    log_coeff = log_coeffs[classifier]
    index_log_coeff = np.argsort(-1*log_coeff)
    for i in index_log_coeff:
        feat_coeffs[classifier][feature_list[i].replace(",","_")] = log_coeff[i]

#read in antiSMASH file
try:
    record = SeqIO.read(open(antismash_file, 'rU'),"genbank")
except:
    print("ERROR: error reading antismash output file")
    exit()
as_features = record.features

#check if RGI file is needed and read in 
rgi_infile = None
if rgi_file is not None and rgi_version != 0:
    try:
        rgi_infile = open(rgi_file, 'r')
    except:
        print("ERROR: error reading rgi output file")
        exit() 
if "CARD_genes" in feature_type_list and rgi_file == None:
    print("model uses CARD features, provide an RGI file")
    exit()
(gene_data, features) = readInputFiles.readInputFilesLocs(as_features, antismash_version, rgi_infile, rgi_version, training_features, feature_type_list, feature_list_by_type, "feature_matrices/"+model_name)

#rank genes
coef_totals = {}
for classfier in log_coeffs:
    coef_totals[classifier] = {}
for g in gene_data:
    gene_coef_totals = {}
    for classifier in log_coeffs:
        gene_coef_totals[classifier] = 0.0
    start = g["start"]
    end = g["end"]
    feature_types = []
    feature_names = []
    feature_coeffs = []
    for feat_type in features:
        if feat_type not in feature_type_list:
            print("ERROR: feature " + feat_type + " missing from feature list")
            exit()
        
        for x in features[feat_type]:
            x_start = x["start"]
            x_end = x["end"]
            if not (x_start >= start and x_start < end and x_end <= end):
                continue
            for classifier in log_coeffs:
                if x["description"] in feat_coeffs[classifier]:
                    gene_coef_totals[classifier] += feat_coeffs[classifier][x["description"]]
                    
    for classifier in log_coeffs:
        if classifier not in coef_totals:
            coef_totals[classifier] = {}
        coef_totals[classifier][g["locus_tag"][0] + "_" + g["product"][0]] = gene_coef_totals[classifier]
        
outfile = open(outdir + "/" + antismash_file[antismash_file.rfind("/"):antismash_file.rfind(".")] + "_ranked_genes.txt",'w')
for classifier in coef_totals:
    outfile.write(classifier + "\n")
    sorted_genes = sorted(coef_totals[classifier].items(), key=lambda item: item[1], reverse=True)
    for gene, score in sorted_genes:
        outfile.write(gene + "," + str(score) + "\n")
outfile.close()               
    