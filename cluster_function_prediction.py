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

SSN_pfam_names = ["Thiolase, N-terminal domain","ABC transporter","Acyl transferase domain","AAA domain",
                 "ABC-2 family transporter protein","Acyl-CoA dehydrogenase, C-terminal domain","Acyl-CoA dehydrogenase, N-terminal domain",
                 "Alcohol dehydrogenase GroES-like domain","Alpha/beta hydrolase family","Aminotransferase class I and II",
                 "Beta-ketoacyl synthase, C-terminal domain","Beta-ketoacyl synthase, N-terminal domain","Cytochrome P450","DegT/DnrJ/EryC1/StrS aminotransferase family",
                 "Enoyl-(Acyl carrier protein) reductase","Erythronolide synthase docking","FAD binding domain","Glycosyl transferase family 2",
                 "Glycosyltransferase family 28 N-terminal domain","Glycosyl transferases group 1","Glycosyltransferase like family 2","Glyoxalase/Bleomycin resistance protein/Dioxygenase superfamily",
                 "KR domain","Lanthionine synthetase C-like protein",
                 "Major Facilitator Superfamily","Methyltransferase small domain","Methyltransferase domain",
                 "NAD dependent epimerase/dehydratase family","NDP-hexose 2,3-dehydratase",
                 "O-methyltransferase","Oxidoreductase family, C-terminal alpha/beta domain","Oxidoreductase family, NAD-binding Rossmann fold",
                 "Phosphopantetheine attachment site","Polyketide cyclase / dehydrase and lipid transport","Polyketide synthase dehydratase",
                 "Protein of unknown function (DUF1205)",
                 "short chain dehydrogenase","SnoaL-like domain","SpaB C-terminal domain",
                 "Sugar (and other) transporter","transcriptional_regulatory_protein,_c_terminal_domains","Thioesterase superfamily","ubiE/COQ5 methyltransferase family","UDP-glucoronosyl and UDP-glucosyl transferase","YcaO-like family",
                 "Zinc-binding dehydrogenase","pyridine_nucleotide-disulphide_oxidoreductase"]

#read arguments given by user
parser = argparse.ArgumentParser()
parser.add_argument('antismash_results',help='file containing the antismash results for the cluster in a genbank file')
parser.add_argument('rgi_results',help='file containing the rgi results for the cluster')
parser.add_argument('--output', help='set directory to write predictions to, default write to current directory') 
parser.add_argument('--seed', help='random seed to use for training classifiers',type=int) 
parser.add_argument('--no_SSN', help="don't use pfam subfamilies in classification, program will run faster with only small impact on accuracy (default: use sub-PFAMs)", nargs='?', default=False, const=True)
parser.add_argument('--blastp_path', help="path to blastp executable, only neeeded if using SSN, default is blastp")
parser.add_argument('--write_features', help='set directory to write features to, default do not write features') 
parser.add_argument('--antismash_version', help='version of antismash used to generate antismash input file, supported versions are 4 and 5, defualt 5') 
parser.add_argument('--rgi_version', help='version of rgi used to generate antismash input file, supported versions are 3 and 5, default 5') 


args = parser.parse_args()

data_path = os.path.dirname(sys.argv[0]) + "/"

if args.write_features == None:
    write_features = False
    feature_dir = ""
else:
    write_features = True
    feature_dir = args.write_features

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
    rgi_version = 5
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

    
#check validity of files and directories given by user
if not tools.checkIfFileExists(antismash_infilename, "antismash") or not tools.checkIfFileExists(rgi_infilename, "rgi"):
    exit()
if not os.path.isdir(out_directory):
    print("The given out directory does not exist, please enter a valid directory")
    exit()
if not os.access(out_directory, os.W_OK):
    print("You do not have permission to write to the given output directory, please use a different directory")
    exit()    

#read the list of features
try:    
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
    exit()


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

original_training_features = training_features


    
#read SSN features
SSN_list = readFeatureFiles.readFeatureList(data_path+"feature_matrices/SSN_list.txt")
for i in range(0, len(SSN_list)):
    SSN_list[i] = SSN_list[i].replace("\r","")
    SSN_list[i] = SSN_list[i].replace("\n","")

included_SSN_clusters =  {}
if not no_SSN:
    for pfam_name in SSN_list:
        base = pfam_name[0:pfam_name.rfind("_")]
        if base not in included_SSN_clusters:
            included_SSN_clusters[base] = []
        numbering = pfam_name[pfam_name.rfind("_")+1:len(pfam_name)].replace("\r","")
        included_SSN_clusters[base].append(numbering.replace("\n",""))    

cluster_name = antismash_infilename

if "/" in cluster_name:
    cluster_name = cluster_name[cluster_name.rfind("/")+1:len(cluster_name)]
cluster_name = cluster_name[0:cluster_name.find(".gbk")] 
if not no_SSN:   
    test_SSN_feature_matrix= SSN_tools.generateSSNFeatureMatrix([antismash_infilename], SSN_pfam_names, SSN_list, included_SSN_clusters, blastp_path,cluster_name, data_path) 
    test_features = readInputFiles.readInputFiles(as_features, antismash_version, rgi_infile, rgi_version, training_features, data_path, test_SSN_feature_matrix)
else:
    test_features = readInputFiles.readInputFiles(as_features, antismash_version, rgi_infile, rgi_version, training_features, data_path, [])


if write_features:
    test_features_out =open(feature_dir + "/" + antismash_infilename[antismash_infilename.rfind("/"):antismash_infilename.rfind(".")]+".csv",'w')
    for f in test_features[0]:
        test_features_out.write(str(f)+",")
    test_features_out.close()

#do classifications
is_not_unknown_indices = readFeatureFiles.getNotUnknownIndices(is_unknown)
target_unannotated = is_antibacterial*((targets_gram_pos+targets_gram_neg)<1)
is_not_unknown_indices_gram =  readFeatureFiles.getNotUnknownIndices(is_unknown + target_unannotated)

is_antibacterial = (is_antibacterial >= 1).astype(int)
is_antieuk = ((is_antifungal + is_cytotoxic)>=1).astype(int)
is_gram_pos = (targets_gram_pos >= 1).astype(int)
is_gram_neg = (targets_gram_neg >= 1).astype(int)

training_features = training_features[is_not_unknown_indices,:]

y_vars = []
y_vars = is_antibacterial
y_vars = y_vars[is_not_unknown_indices]



#antibacterial predictions
svm_params = {"kernel":'rbf',"C":10,"gamma":0.1}
tree_params = {"depth":100,"n":50}
log_params = {"l1_ratio":.05,"alpha":.01}

tree_bacterial_prob = tools.treePrediction(training_features, y_vars, test_features, tree_params, seed)
log_bacterial_prob = tools.logPrediction(training_features, y_vars, test_features, log_params, seed)
svm_bacterial_prob = tools.svmPrediction(training_features, y_vars, test_features, svm_params, seed)

#antieuk predictions
y_vars = is_antieuk
y_vars = y_vars[is_not_unknown_indices]
svm_params = {"kernel":'linear',"C":.1}
tree_params = {"depth":None,"n":25}
log_params = {"l1_ratio":.001,"alpha":.001}

tree_antieuk_prob = tools.treePrediction(training_features, y_vars, test_features, tree_params, seed)
log_antieuk_prob = tools.logPrediction(training_features, y_vars, test_features, log_params, seed)
svm_antieuk_prob = tools.svmPrediction(training_features, y_vars, test_features, svm_params, seed)


#antifungal predictions
y_vars = (is_antifungal >= 1).astype(int)
y_vars = y_vars[is_not_unknown_indices]
svm_params = {"kernel":'rbf',"C":10,"gamma":0.1}
tree_params = {"depth":50,"n":50}
log_params = {"l1_ratio":.0001,"alpha":.01}

tree_antifungal_prob = tools.treePrediction(training_features, y_vars, test_features, tree_params, seed)
log_antifungal_prob = tools.logPrediction(training_features, y_vars, test_features, log_params, seed)
svm_antifungal_prob = tools.svmPrediction(training_features, y_vars, test_features, svm_params, seed)

#cytotox and antitumor predictons
y_vars = (is_cytotoxic >= 1).astype(int)
y_vars = y_vars[is_not_unknown_indices]
svm_params = {"kernel":'rbf',"C":100,"gamma":0.1}
tree_params = {"depth":50,"n":100}
log_params = {"l1_ratio":.001,"alpha":.001}

tree_antitumor_prob = tools.treePrediction(training_features, y_vars, test_features, tree_params, seed)
log_antitumor_prob = tools.logPrediction(training_features, y_vars, test_features, log_params, seed)
svm_antitumor_prob = tools.svmPrediction(training_features, y_vars, test_features, svm_params, seed)

#antigram negative predictions
y_vars = is_gram_neg
y_vars = y_vars[is_not_unknown_indices_gram]
training_features = original_training_features[is_not_unknown_indices_gram,:]
svm_params = {"kernel":'rbf',"C":10,"gamma":0.01}
tree_params = {"depth":100,"n":25}
log_params = {"l1_ratio":.05,"alpha":.001}

tree_antigramneg_prob = tools.treePrediction(training_features, y_vars, test_features, tree_params, seed)
log_antigramneg_prob = tools.logPrediction(training_features, y_vars, test_features, log_params, seed)
svm_antigramneg_prob = tools.svmPrediction(training_features, y_vars, test_features, svm_params, seed)

#antigram positive predictions
y_vars = is_gram_pos
y_vars = y_vars[is_not_unknown_indices_gram]
svm_params = {"kernel":'rbf',"C":10,"gamma":.01}
tree_params = {"depth":100,"n":50}
log_params = {"l1_ratio":.001,"alpha":.001}

tree_antigrampos_prob = tools.treePrediction(training_features, y_vars, test_features, tree_params, seed)
log_antigrampos_prob = tools.logPrediction(training_features, y_vars, test_features, log_params, seed)
svm_antigrampos_prob = tools.svmPrediction(training_features, y_vars, test_features, svm_params, seed)

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