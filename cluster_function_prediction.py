#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 11:38:35 2019

@author: Allison Walker
"""
#TODO: add error handling for reading of files
#TODO: warning for not finding any features

import argparse
import cluster_function_prediction_tools as tools
import os, sys
from Bio import SeqIO
import SSN_tools
import readFeatureFiles
import numpy as np
import readInputFiles
import joblib


#read arguments given by user
parser = argparse.ArgumentParser()
parser.add_argument('antismash_results',help='file containing the antismash results for the cluster in a genbank file')
parser.add_argument('rgi_results',help='file containing the rgi results for the cluster')
parser.add_argument('--output', help='set directory to write predictions to, default write to current directory') 
parser.add_argument('--seed', help='random seed to use for training classifiers',type=int) 
parser.add_argument('--no_SSN', help="deprecated: use v1 for SSN features", nargs='?', default=True, const=True)
parser.add_argument('--blastp_path', help="path to blastp executable, only neeeded if using SSN, default is blastp")
parser.add_argument('--write_features', help='set directory to write features to, default do not write features') 
parser.add_argument('--antismash_version', help='version of antismash used to generate antismash input file, supported versions are 4 and 5, enter 0 for not using an rgi file, defualt 0') 
parser.add_argument('--rgi_version', help='version of rgi used to generate antismash input file, supported versions are 3 and 5, default 5') 
parser.add_argument('--webserver_output', help="output more basic machine readable format", nargs='?', default=False, const=True)
parser.add_argument('--database_version', help='version of database, default 1') 


#TODO: remove ssn and all related code

args = parser.parse_args()

data_path = os.path.dirname(sys.argv[0]) + "/"

if args.write_features == None:
    write_features = False
    feature_dir = ""
else:
    write_features = True
    write_feature_dir = args.write_features

if args.seed == None:
    seed = 0
else:
    seed = args.seed
    
if args.blastp_path == None:
    blastp_path = "blastp"
else:
    blastp_path = args.blastp_path

antismash_infilename = args.antismash_results
rgi_infilename = args.rgi_results
no_SSN = args.no_SSN

if args.output == None:
    out_directory = "./"
else:
    out_directory = args.output
    
if args.rgi_version == "5":
    rgi_version = 5
elif args.rgi_version == "3":
    rgi_version = 3
elif args.rgi_version == None:
    rgi_version = 0
else:
    print("please enter a valid rgi version, program currently accepts output from versions 3 and 5")
    exit()

antismash_version = 5    
if args.antismash_version == "5":
    antismash_version = 5
elif args.antismash_version == "4":
    antismash_version = 4
elif args.antismash_version == None:
    antismash_version = 5
else:
    print("please enter a valid antismash version, program currently accepts output from versions 4 and 5")
    exit()
    
database_version = 1
#TODO: add parsing for different database versions

    
#check validity of files and directories given by user
if not tools.checkIfFileExists(antismash_infilename, "antismash") or not tools.checkIfFileExists(rgi_infilename, "rgi"):
    exit()
if not os.path.isdir(out_directory):
    print("The given out directory does not exist, please enter a valid directory")
    exit()
if not os.access(out_directory, os.W_OK):
    print("You do not have permission to write to the given output directory, please use a different directory")
    exit()    

#figure out appropriate model name
if antismash_version == 4 and rgi_version == 3 and database_version == 1:
    model_name = "antismash4rgi3"
elif antismash_version == 4 and rgi_version == 5 and database_version == 1:
    model_name = "antismash4rgi5"
elif antismash_version == 4 and rgi_version == 0 and database_version == 1:
    model_name = "antismash4"
elif antismash_version == 5 and rgi_version == 3 and database_version == 1:
    model_name = "antismash5rgi3"
elif antismash_version == 5 and rgi_version == 5 and database_version == 1:
    model_name = "antismash5rgi5"
elif antismash_version == 5 and rgi_version == 0 and database_version == 1:
    model_name = "antismash5"
elif antismash_version == 6 and rgi_version == 3 and database_version == 1:
    model_name = "antismash6rgi3"
elif antismash_version == 6 and rgi_version == 5 and database_version == 1:
    model_name = "antismash6rgi5"
elif antismash_version == 6 and rgi_version == 0 and database_version == 1:
    model_name = "antismash6"
else:
    #TODO: throw error
    print("options not compatible with this version")
    exit()
try:    
    feature_dir = data_path + "feature_matrices/" + model_name  + "/features/"
    feature_type_list = readFeatureFiles.getFeatureFilesList(feature_dir)
    training_features = readFeatureFiles.readFeatures(feature_dir, feature_type_list)
    feature_list = readFeatureFiles.readFeatureNames(feature_dir, feature_type_list)
except:
    print("did not find file containing training data, please keep script located in directory downloaded from github")


#read the list of features
#TODO: automate this?
"""try:    
    training_SSN_features = readFeatureFiles.readFeatureMatrix(data_path+"feature_matrices/SSN.csv")
    if antismash_version == 4:  
        training_pfam_features = readFeatureFiles.readFeatureMatrix(data_path+"feature_matrices/PFAM.csv")
        training_smCOG_features = readFeatureFiles.readFeatureMatrix(data_path+"feature_matrices/SMCOG.csv")
        #SSN_calc_features = readFeatureFiles.readFeatureMatrixFloat("gene_feature_matrices/test_compounds_SSN.csv")
        training_CDS_features = readFeatureFiles.readFeatureMatrix(data_path+"feature_matrices/CDS_motifs.csv")
        
        training_pks_nrps_type_features = readFeatureFiles.readFeatureMatrix(data_path+"feature_matrices/pks_nrps_type.csv")
        training_pk_signature_features = readFeatureFiles.readFeatureMatrix(data_path+"feature_matrices/pk_signature.csv")
        training_pk_minowa_features = readFeatureFiles.readFeatureMatrix(data_path+"feature_matrices/pk_minowa.csv")
        training_pk_consensus_features = readFeatureFiles.readFeatureMatrix(data_path+"feature_matrices/pk_consensus.csv")
        
        training_nrp_stachelhaus_features = readFeatureFiles.readFeatureMatrix(data_path+"feature_matrices/nrp_stachelhaus.csv")
        training_nrp_nrpspredictor_features = readFeatureFiles.readFeatureMatrix(data_path+"feature_matrices/nrp_nrpspredictor.csv")
        training_nrp_pHMM_features = readFeatureFiles.readFeatureMatrix(data_path+"feature_matrices/nrp_pHMM.csv")
        training_nrp_predicat_features = readFeatureFiles.readFeatureMatrix(data_path+"feature_matrices/nrp_predicat.csv")
        training_nrp_sandpuma_features = readFeatureFiles.readFeatureMatrix(data_path+"feature_matrices/nrp_sandpuma.csv")
    elif antismash_version == 5:
        training_pfam_features = readFeatureFiles.readFeatureMatrix(data_path+"feature_matrices/PFAM5.csv")
        training_smCOG_features = readFeatureFiles.readFeatureMatrix(data_path+"feature_matrices/SMCOG5.csv")
        training_CDS_features = readFeatureFiles.readFeatureMatrix(data_path+"feature_matrices/CDS_motifs5.csv")        
        training_pk_consensus_features = readFeatureFiles.readFeatureMatrix(data_path+"feature_matrices/pk_nrp_consensus5.csv")

        
    if rgi_version == 3:
        training_card_features = readFeatureFiles.readFeatureMatrix(data_path+"feature_matrices/CARD_gene.csv")        
        used_resistance_genes_list = readFeatureFiles.readFeatureList(data_path+"feature_matrices/CARD_gene_list.txt")
    elif rgi_version == 5:
        training_card_features = readFeatureFiles.readFeatureMatrix(data_path+"feature_matrices/CARD5_genes.csv")
        used_resistance_genes_list = readFeatureFiles.readFeatureList(data_path+"feature_matrices/CARD5_gene_list.txt")
    
    is_antibacterial = readFeatureFiles.readClassesMatrix(data_path+"feature_matrices/is_antibacterial.csv")
    is_antifungal = readFeatureFiles.readClassesMatrix(data_path+"feature_matrices/is_antifungal.csv")
    is_cytotoxic = readFeatureFiles.readClassesMatrix(data_path+"feature_matrices/is_cytotoxic.csv")
    is_unknown = readFeatureFiles.readClassesMatrix(data_path+"feature_matrices/is_unknown.csv")
    targets_gram_pos = readFeatureFiles.readClassesMatrix(data_path+"feature_matrices/targets_gram_pos.csv")
    targets_gram_neg = readFeatureFiles.readClassesMatrix(data_path+"feature_matrices/targets_gram_neg.csv")
    full_cluster_list = readFeatureFiles.readClusterList(data_path+"feature_matrices/cluster_list_CARD.txt")        
except:
    print("did not find file containing training data, please keep script located in directory downloaded from github")
    exit()"""


#read the antismash input file

try:
    record = SeqIO.read(open(antismash_infilename, 'rU'),"genbank")
except:
    print("error reading antismash output file")
    exit()
as_features = record.features
try:
    rgi_infile = open(rgi_infilename, 'r')
except:
   print("error reading rgi output file")
   exit() 

"""
#make the feature matrices for the cluster
training_features = np.concatenate((training_pfam_features, training_card_features), axis=1)
training_features = np.concatenate((training_features,  training_smCOG_features), axis=1)
training_features = np.concatenate((training_features,  training_CDS_features), axis=1)
if not no_SSN:
    training_features = np.concatenate((training_features,  training_SSN_features), axis=1)
if antismash_version == 4:
    training_features = np.concatenate((training_features,  training_pks_nrps_type_features), axis=1)
    training_features = np.concatenate((training_features,  training_pk_signature_features), axis=1)
    training_features = np.concatenate((training_features,  training_pk_minowa_features), axis=1)
    training_features = np.concatenate((training_features,  training_pk_consensus_features), axis=1)
    training_features = np.concatenate((training_features,  training_nrp_stachelhaus_features), axis=1)
    training_features = np.concatenate((training_features,  training_nrp_nrpspredictor_features), axis=1)
    training_features = np.concatenate((training_features,  training_nrp_pHMM_features), axis=1)
    training_features = np.concatenate((training_features,  training_nrp_predicat_features), axis=1)
    training_features = np.concatenate((training_features,  training_nrp_sandpuma_features), axis=1)
else:
    training_features = np.concatenate((training_features,  training_pk_consensus_features), axis=1)

original_training_features = training_features"""


    



cluster_name = antismash_infilename

if "/" in cluster_name:
    cluster_name = cluster_name[cluster_name.rfind("/")+1:len(cluster_name)]
cluster_name = cluster_name[0:cluster_name.find(".gbk")] 
if not no_SSN:   
    print("SSNs no lonter supported in this version, please use version 1 for SSNs features")
    exit()
else:
    test_features = readInputFiles.readInputFiles(as_features, antismash_version, rgi_infile, rgi_version, training_features, feature_dir, [])


#TODO: fix this for webserver
if write_features:
    test_features_out =open(write_feature_dir + "/" + antismash_infilename[antismash_infilename.rfind("/"):antismash_infilename.rfind(".")]+".csv",'w')
    for f in test_features[0]:
        test_features_out.write(str(f)+",")
    test_features_out.close()


#Load appropriate pretrained model
try:
    svm_bacterial, tree_bacterial, log_bacterial = joblib.load(data_path+"trained_models/" + model_name + "_antibacterial.sav")
    svm_bacterial_prob = svm_bacterial.predict_proba(test_features)
    tree_bacterial_prob = tree_bacterial.predict_proba(test_features)
    log_bacterial_prob = log_bacterial.predict_proba(test_features)

    svm_antieuk, tree_antieuk, log_antieuk = joblib.load(data_path+"trained_models/" + model_name + "_antieuk.sav")
    svm_antieuk_prob = svm_antieuk.predict_proba(test_features)
    tree_antieuk_prob = tree_antieuk.predict_proba(test_features)
    log_antieuk_prob = log_antieuk.predict_proba(test_features)

    svm_antifungal, tree_antifungal, log_antifungal = joblib.load(data_path+"trained_models/" + model_name + "_antifungal.sav")
    svm_antifungal_prob = svm_antifungal.predict_proba(test_features)
    tree_antifungal_prob = tree_antifungal.predict_proba(test_features)
    log_antifungal_prob = log_antifungal.predict_proba(test_features)

    svm_antitumor, tree_antitumor, log_antitumor = joblib.load(data_path+"trained_models/" + model_name + "_cytotoxic_antitumor.sav")
    svm_antitumor_prob = svm_antitumor.predict_proba(test_features)
    tree_antitumor_prob = tree_antitumor.predict_proba(test_features)
    log_antitumor_prob = log_antitumor.predict_proba(test_features)
    
    svm_antigramneg, tree_antigramneg, log_antigramneg = joblib.load(data_path+"trained_models/" + model_name + "_antigramneg.sav")
    svm_antigramneg_prob = svm_antigramneg.predict_proba(test_features)
    tree_antigramneg_prob = tree_antigramneg.predict_proba(test_features)
    log_antigramneg_prob = log_antigramneg.predict_proba(test_features)

    svm_antigrampos, tree_antigrampos, log_antigrampos = joblib.load(data_path+"trained_models/" + model_name + "_antigrampos.sav")
    svm_antigrampos_prob = svm_antigrampos.predict_proba(test_features)
    tree_antigrampos_prob = tree_antigrampos.predict_proba(test_features)
    log_antigrampos_prob = log_antigrampos.predict_proba(test_features)
except:
   print("could not find pretrained model, make sure all data files are in correct location")
   exit() 

#TODO: make more machine readable for webserver
#print the results
print("probabilities of antibacterial activity:")
print("tree classifier: " + str(tree_bacterial_prob[0,1])  + " logistic regression classifier: " + str(log_bacterial_prob[0,1]) + " svm classifier: " + str(svm_bacterial_prob[0,1]))
print("probabilities of anti-gram positive activity:")
print("tree classifier: " + str(tree_antigrampos_prob[0,1])  + " logistic regression classifier: " + str(log_antigrampos_prob[0,1]) + " svm classifier: " + str(svm_antigrampos_prob[0,1]))
print("probabilities of anti-gram negative activity:")
print("tree classifier: " + str(tree_antigramneg_prob[0,1])  + " logistic regression classifier: " + str(log_antigramneg_prob[0,1]) + " svm classifier: " + str(svm_antigramneg_prob[0,1]))
print("probabilities of antifungal or antitumor or cytotoxic activity:")
print("tree classifier: " + str(tree_antieuk_prob[0,1])  + " logistic regression classifier: " + str(log_antieuk_prob[0,1]) + " svm classifier: " + str(svm_antieuk_prob[0,1]))
print("probabilities of antifungal activity:")
print("tree classifier: " + str(tree_antifungal_prob[0,1])  + " logistic regression classifier: " + str(log_antifungal_prob[0,1]) + " svm classifier: " + str(svm_antifungal_prob[0,1]))
print("probabilities of antitumor or cytotoxic activity:")
print("tree classifier: " + str(tree_antitumor_prob[0,1])  + " logistic regression classifier: " + str(log_antitumor_prob[0,1]) + " svm classifier: " + str(svm_antitumor_prob[0,1]))
                      
#write output
#TODO: remove for webserver?

try:
    if out_directory == "":
        outfile = open(cluster_name + ".txt",'w')
    else:
        outfile = open(out_directory + "/" +cluster_name + ".txt",'w')
except:
    print("couldn't open output file, please provide an output directory that can be written to")
    exit()

tools.writeProbabilitiesToFile(outfile, "antibacterial", tree_bacterial_prob, log_bacterial_prob, svm_bacterial_prob)
tools.writeProbabilitiesToFile(outfile, "anti-gram positive", tree_antigrampos_prob, log_antigrampos_prob, svm_antigrampos_prob)
tools.writeProbabilitiesToFile(outfile, "anti-gram negative", tree_antigramneg_prob, log_antigramneg_prob, svm_antigramneg_prob)
tools.writeProbabilitiesToFile(outfile, "antifugnal or antitumor or cytotoxic", tree_antieuk_prob, log_antieuk_prob, svm_antieuk_prob)
tools.writeProbabilitiesToFile(outfile, "antifungal", tree_antifungal_prob, log_antifungal_prob, svm_antifungal_prob)
tools.writeProbabilitiesToFile(outfile, "antitumor or cytotoxic", tree_antitumor_prob, log_antitumor_prob, svm_antitumor_prob)
outfile.close()