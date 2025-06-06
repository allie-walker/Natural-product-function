#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 11:38:35 2019

@author: Allison Walker
"""
#TODO: add error handling for reading of files
#TODO: warning for not finding any features
#TODO: handle not having rgi file

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
parser.add_argument('--rgi_results',help='file containing the rgi results for the cluster')
parser.add_argument('--output', help='set directory to write predictions to, default write to current directory',default="./") 
parser.add_argument('--seed', help='random seed to use for training classifiers',type=int, default=0) 
parser.add_argument('--no_SSN', help="deprecated: use v1 for SSN features", nargs='?', default=True, const=True)
parser.add_argument('--blastp_path', help="path to blastp executable, only neeeded if using SSN, default is blastp")
parser.add_argument('--write_features', help='set directory to write features to, default do not write features') 
parser.add_argument('--antismash_version', help='version of antismash used to generate antismash input file, supported versions are 4-8, defualt 8') 
parser.add_argument('--rgi_version', help='version of rgi used to generate antismash input file, supported versions are 3 and 5, default 5') 
parser.add_argument('--webserver_output', help="output more basic machine readable format", nargs='?', default=False, const=True)
parser.add_argument('--on_webserver', help="change data dir path for webserver", nargs='?', default=False, const=True)
parser.add_argument('--database_version', help='version of database, default 1') 


#TODO: remove ssn and all related code

args = parser.parse_args()
if not args.on_webserver:
    data_path = os.path.dirname(sys.argv[0]) + "/"
else:
    data_path = "./"

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
 
if  args.rgi_version == None and rgi_infilename == None:
    rgi_version = 0
elif args.rgi_version == "5":
    rgi_version = 5
elif args.rgi_version == "3":
    rgi_version = 3    
elif args.rgi_version == "6":
    rgi_version = 6
else:
    print("ERROR: please enter a valid rgi version, program currently accepts output from versions 3 and 5")
    exit()

   
if args.antismash_version == None:
    antismash_version = 8
else:
    try:
        antismash_version = int(args.antismash_version)
    except:
        print("ERROR: enter a number for the antiSMASH version")
if antismash_version < 4 or antismash_version > 8:
    print("ERROR: please enter a valid antismash version, program currently accepts output from versions 4-8")
    exit()
    
database_version = 1
#TODO: add parsing for different database versions

  
#check validity of files and directories given by user
if not tools.checkIfFileExists(antismash_infilename, "antismash"):
    exit()
if rgi_infilename is not None and rgi_version != 0 and not tools.checkIfFileExists(rgi_infilename, "rgi"):
    exit()
if not os.path.isdir(out_directory):
    print("ERROR: The given out directory does not exist, please enter a valid directory")
    exit()
if not os.access(out_directory, os.W_OK):
    print("ERROR: You do not have permission to write to the given output directory, please use a different directory")
    exit()    

#figure out appropriate model name
if rgi_version == 0:
    model_name = "antismash" + str(antismash_version)
else:
    model_name = "antismash" + str(antismash_version) + "rgi" + str(rgi_version)
if not os.path.exists( data_path + "trained_models/" + model_name + "_antibacterial.sav"):
    #TODO: throw error
    print("ERROR: options not compatible with this version")
    exit()
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
    print("ERROR: did not find file containing training data, please keep script located in directory downloaded from github")



#read the list of features
#TODO: automate this?

#read the antismash input file

try:
    record = SeqIO.read(open(antismash_infilename, 'rU'),"genbank")
except:
    print("ERROR: error reading antismash output file")
    exit()
as_features = record.features
rgi_infile = None
if rgi_infilename is not None and rgi_version != 0:
    try:
        rgi_infile = open(rgi_infilename, 'r')
    except:
        print("ERROR: error reading rgi output file")
        exit() 



cluster_name = antismash_infilename

if "/" in cluster_name:
    cluster_name = cluster_name[cluster_name.rfind("/")+1:len(cluster_name)]
cluster_name = cluster_name[0:cluster_name.find(".gbk")] 
if not no_SSN:   
    print("ERROR: SSNs no longer supported in this version, please use version 1 for SSNs features")
    exit()
else:
    #TODO: need to fix this for antismash6
    test_features = readInputFiles.readInputFiles(as_features, antismash_version, rgi_infile, rgi_version, training_features, feature_type_list, feature_list_by_type, data_path+ "feature_matrices/"+model_name)


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
   print("ERROR: could not find pretrained model, make sure all data files are in correct location")
   exit() 

#TODO: make more machine readable for webserver
#print the results
if not args.webserver_output:
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
    tools.writeProbabilitiesToFile(outfile, "antifungal or antitumor or cytotoxic", tree_antieuk_prob, log_antieuk_prob, svm_antieuk_prob)
    tools.writeProbabilitiesToFile(outfile, "antifungal", tree_antifungal_prob, log_antifungal_prob, svm_antifungal_prob)
    tools.writeProbabilitiesToFile(outfile, "antitumor or cytotoxic", tree_antitumor_prob, log_antitumor_prob, svm_antitumor_prob)
    outfile.close()
else:
    #will have format Success: \n
    print("Success:")
    #ab: tree, log, svm
    print(str(tree_bacterial_prob[0,1])  + "," + str(log_bacterial_prob[0,1]) + "," + str(svm_bacterial_prob[0,1]))
    #agrampos: tree, log, svm
    print(str(tree_antigrampos_prob[0,1])  + "," + str(log_antigrampos_prob[0,1]) + "," + str(svm_antigrampos_prob[0,1]))
    #antigramneg: tree, log, svm
    print(str(tree_antigramneg_prob[0,1])  + "," + str(log_antigramneg_prob[0,1]) + "," + str(svm_antigramneg_prob[0,1]))
    #antieuk: tree, log, svm
    print(str(tree_antieuk_prob[0,1])  + "," + str(log_antieuk_prob[0,1]) + "," + str(svm_antieuk_prob[0,1]))
    #antifungal: tree, log, svm
    print(str(tree_antifungal_prob[0,1])  + "," + str(log_antifungal_prob[0,1]) + "," + str(svm_antifungal_prob[0,1]))
    #antitumor: tree, log, svm
    print(str(tree_antitumor_prob[0,1])  + "," + str(log_antitumor_prob[0,1]) + " ," + str(svm_antitumor_prob[0,1]))
   