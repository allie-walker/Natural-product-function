#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 11:38:35 2019

@author: allis
"""


import argparse
import cluster_function_prediction_tools as tools
import os, sys
from Bio import SeqIO
import SSN_tools
import readFeatureFiles
import numpy as np

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

#TODO: add arguments for writing features
parser.add_argument('antismash_results',help='file containing the antismash results for the cluster in a genbank file')
parser.add_argument('rgi_results',help='file containing the rgi results for the cluster')
parser.add_argument('--output', help='set directory to write predictions to, default write to current directory') 
parser.add_argument('--seed', help='random seed to use for training classifiers',type=int) 
parser.add_argument('--no_SNN', help="don't use pfam subfamilies in classification, program will run faster with only small impact on accuracy (default: use sub-PFAMs)", nargs='?', default=False, const=True)
args = parser.parse_args()

data_path = os.path.dirname(sys.argv[0]) + "/"
#TODO: FIX no_SSN option!!
#TODO: blastp path

if args.seed == None:
    seed = 0
else:
    seed = args.seed

antismash_infilename = args.antismash_results
rgi_infilename = args.rgi_results

if args.output == None:
    out_directory = ""
else:
    out_directory = args.output
    
    
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
    used_pfam_list = readFeatureFiles.readFeatureList(data_path+"gene_feature_matrices/pfam_list.txt")
    used_CDS_list = readFeatureFiles.readFeatureList(data_path+"gene_feature_matrices/CDS_motifs_list.txt")
    used_smCOG_list = readFeatureFiles.readFeatureList(data_path+"gene_feature_matrices/smCOG_list.txt")
    used_pks_nrps_type_list = readFeatureFiles.readFeatureList(data_path+"gene_feature_matrices/pks_nrps_type_list.txt")
    used_pk_signature_list = readFeatureFiles.readFeatureList(data_path+"gene_feature_matrices/pk_signature_list.txt")
    used_pk_minowa_list = readFeatureFiles.readFeatureList(data_path+"gene_feature_matrices/pk_minowa_list.txt")
    used_pk_consensus_list = readFeatureFiles.readFeatureList(data_path+"gene_feature_matrices/pk_consensus_list.txt")
    
    used_nrp_stachelhaus_list = readFeatureFiles.readFeatureList(data_path+"gene_feature_matrices/nrp_stachelhaus_list.txt")
    used_nrp_nrps_predictor_list = readFeatureFiles.readFeatureList(data_path+"gene_feature_matrices/nrp_nrpspredictor_list.txt")
    used_nrp_pHMM_list = readFeatureFiles.readFeatureList(data_path+"gene_feature_matrices/nrp_pHMM_list.txt")
    used_nrp_predicat_list = readFeatureFiles.readFeatureList(data_path+"gene_feature_matrices/nrp_predicat_list.txt")
    used_nrp_sandpuma_list = readFeatureFiles.readFeatureList(data_path+"gene_feature_matrices/nrp_sandpuma_list.txt")
    
    used_resistance_genes_list = readFeatureFiles.readFeatureList(data_path+"gene_feature_matrices/CARD_gene_list.txt")
    
    is_antibacterial = readFeatureFiles.readClassesMatrix(data_path+"gene_feature_matrices/is_antibacterial.csv")
    is_antifungal = readFeatureFiles.readClassesMatrix(data_path+"gene_feature_matrices/is_antifungal.csv")
    is_cytotoxic = readFeatureFiles.readClassesMatrix(data_path+"gene_feature_matrices/is_cytotoxic.csv")
    is_unknown = readFeatureFiles.readClassesMatrix(data_path+"gene_feature_matrices/is_unknown.csv")
    targets_gram_pos = readFeatureFiles.readClassesMatrix(data_path+"gene_feature_matrices/targets_gram_pos.csv")
    targets_gram_neg = readFeatureFiles.readClassesMatrix(data_path+"gene_feature_matrices/targets_gram_neg.csv")
    full_cluster_list = readFeatureFiles.readClusterList(data_path+"gene_feature_matrices/cluster_list.txt")
    
    training_pfam_features = readFeatureFiles.readFeatureMatrix(data_path+"gene_feature_matrices/PFAM.csv")
    training_card_features = readFeatureFiles.readFeatureMatrix(data_path+"gene_feature_matrices/CARD_genes.csv")
    training_smCOG_features = readFeatureFiles.readFeatureMatrix(data_path+"gene_feature_matrices/SMCOG.csv")
    training_SSN_features = readFeatureFiles.readFeatureMatrix(data_path+"gene_feature_matrices/SSN.csv")
    #SSN_calc_features = readFeatureFiles.readFeatureMatrixFloat("gene_feature_matrices/test_compounds_SSN.csv")
    training_CDS_features = readFeatureFiles.readFeatureMatrix(data_path+"gene_feature_matrices/CDS_motifs.csv")
    
    training_pks_nrps_type_features = readFeatureFiles.readFeatureMatrix(data_path+"gene_feature_matrices/pks_nrps_type.csv")
    training_pk_signature_features = readFeatureFiles.readFeatureMatrix(data_path+"gene_feature_matrices/pk_signature.csv")
    training_pk_minowa_features = readFeatureFiles.readFeatureMatrix(data_path+"gene_feature_matrices/pk_minowa.csv")
    training_pk_consensus_features = readFeatureFiles.readFeatureMatrix(data_path+"gene_feature_matrices/pk_consensus.csv")

    training_nrp_stachelhaus_features = readFeatureFiles.readFeatureMatrix(data_path+"gene_feature_matrices/nrp_stachelhaus.csv")
    training_nrp_nrpspredictor_features = readFeatureFiles.readFeatureMatrix(data_path+"gene_feature_matrices/nrp_nrpspredictor.csv")
    training_nrp_pHMM_features = readFeatureFiles.readFeatureMatrix(data_path+"gene_feature_matrices/nrp_pHMM.csv")
    training_nrp_predicat_features = readFeatureFiles.readFeatureMatrix(data_path+"gene_feature_matrices/nrp_predicat.csv")
    training_nrp_sandpuma_features = readFeatureFiles.readFeatureMatrix(data_path+"gene_feature_matrices/nrp_sandpuma.csv")
except:
    print("did not find file containing training data")
    exit()
#read the antismash input file
try:
    record = SeqIO.read(open(antismash_infilename, 'rU'),"genbank")
except:
    print("error reading antismash output file")
    exit()
features = record.features
score_cutoff = 20
CDS_motif_list = []
smCOG_list = []
pfam_list = []
resistance_genes_list = []

CDS_motifs = {}
smCOGs = {}
pfam_counts = {}
pks_nrp_subtypes = {}
pk_monomers_signature = {}
pk_monomers_minowa = {}
pk_monomers_consensus = {}
nrp_monomers_stachelhaus = {}
nrp_monomers_nrps_predictor= {}
nrp_monomers_pHMM = {}
nrp_monomers_predicat = {}
nrp_monomers_sandpuma= {}
subtype = ""
pks_signature = ""
minowa = ""
consensus = ""
stachelhaus = ""
nrpspredictor = ""
pHMM = ""
predicat = ""
sandpuma = ""
for feature in features:
    subtype = ""
    pks_signature = ""
    minowa = ""
    consensus = ""
    stachelhaus = ""
    nrpspredictor = ""
    pHMM = ""
    predicat = ""
    sandpuma = ""
    if feature.type == "CDS" and "sec_met" in feature.qualifiers:
        for f in feature.qualifiers["sec_met"]:
            if "NRPS/PKS subtype" in f:
                subtype = f.split(":")[1]
            if "Substrate specificity predictions" in f:
                if "PKS_AT" in f:
                    predictions = f.split(":")[4]
                    pks_signature = predictions.split(",")[0].split()[0]
                    minowa = predictions.split(",")[1].split()[0]
                    consensus = predictions.split(",")[2].split()[0]
                if "AMP-binding" in f:
                    predictions = f.split(":")[5]
                    stachelhaus = predictions.split(",")[0].split()[0]
                    nrpspredictor = predictions.split(",")[1].split()[0]
                    pHMM = predictions.split(",")[2].split()[0]
                    predicat = predictions.split(",")[3].split()[0]
                    sandpuma = predictions.split(",")[4].split()[0]
    if subtype != "":
        if subtype not in pks_nrp_subtypes:
            pks_nrp_subtypes[subtype] = 0
        pks_nrp_subtypes[subtype] += 1
    if pks_signature != "" and "no_call" not in pks_signature and "N/A" not in pks_signature and "n/a" not in pks_signature:
        if pks_signature not in pk_monomers_signature:
            pk_monomers_signature[pks_signature] = 0
        pk_monomers_signature[pks_signature] += 1
    if minowa != "" and "no_call" not in minowa and "N/A" not in minowa and "n/a" not in minowa:
        if minowa not in pk_monomers_minowa:
            pk_monomers_minowa[minowa] = 0
        pk_monomers_minowa[minowa] += 1
    if consensus != "" and "no_call" not in consensus and "N/A" not in consensus and "n/a" not in consensus:
        if consensus not in pk_monomers_consensus:
            pk_monomers_consensus[consensus] = 0
        pk_monomers_consensus[consensus] += 1
    if stachelhaus != "" and "no_call" not in stachelhaus and "N/A" not in stachelhaus and "n/a" not in stachelhaus:
        if stachelhaus not in nrp_monomers_stachelhaus:
            nrp_monomers_stachelhaus[stachelhaus] = 0
        nrp_monomers_stachelhaus[stachelhaus] += 1
    if nrpspredictor != "" and "no_call" not in nrpspredictor and "N/A" not in nrpspredictor and "n/a" not in nrpspredictor:
        if nrpspredictor not in nrp_monomers_nrps_predictor:
            nrp_monomers_nrps_predictor[nrpspredictor] = 0
        nrp_monomers_nrps_predictor[nrpspredictor] += 1
    if pHMM != "" and "no_call" not in pHMM and "N/A" not in pHMM and "n/a" not in pHMM:
        if pHMM not in nrp_monomers_pHMM:
            nrp_monomers_pHMM[pHMM] = 0
        nrp_monomers_pHMM[pHMM] += 1
    if predicat != "" and "no_call" not in predicat and "N/A" not in predicat and "n/a" not in predicat:
        if predicat not in nrp_monomers_predicat:
            nrp_monomers_predicat[predicat] = 0
        nrp_monomers_predicat[predicat] += 1
    if sandpuma != "" and "no_call" not in sandpuma and "N/A" not in sandpuma and "n/a" not in sandpuma:
        if sandpuma not in nrp_monomers_sandpuma:
            nrp_monomers_sandpuma[sandpuma] = 0
        nrp_monomers_sandpuma[sandpuma] += 1
        
    if feature.type == "CDS_motif":
        note_text = feature.qualifiers['note'][0]
        if "(" not in note_text:
            continue
        motif_name = note_text[0:note_text.index("(")-1]
        if motif_name not in CDS_motif_list:
            CDS_motif_list.append(motif_name)
        if motif_name not in CDS_motifs:
            CDS_motifs[motif_name] = 0
        CDS_motifs[motif_name] += 1
    elif feature.type == "CDS":
        if "note" in feature.qualifiers:
            for note in feature.qualifiers["note"]:
                #print note
                if "smCOG" in note:
                    #print note
                    if ":" not in note or "(" not in note:
                        continue
                    smCOG_type = note[note.index(":")+2:note.index("(")-1]
                    if smCOG_type not in smCOG_list:
                        smCOG_list.append(smCOG_type)
                    if smCOG_type not in smCOGs:
                        smCOGs[smCOG_type] = 0
                    smCOGs[smCOG_type] += 1
    elif feature.type == "PFAM_domain":
        score = float(feature.qualifiers["score"][0])
        if score <score_cutoff:
            continue
        domain_description = feature.qualifiers["description"][0]
        if domain_description not in pfam_list:
            pfam_list.append(domain_description)
        if domain_description not in pfam_counts:
            pfam_counts[domain_description] = 0
        pfam_counts[domain_description] += 1

#read resistance features
in_file = open(rgi_infilename, 'r')
e_value_threshold = 0.1
resistance_genes = {}
for line in in_file:
    if "ORF_ID" in line:
        continue
    entries = line.split("\t")
    e_value = float(entries[7])
    if e_value > e_value_threshold:
        continue
    best_hit = entries[8]
    hit_names = entries[11]
    if best_hit not in resistance_genes_list:
        resistance_genes_list.append(best_hit)
    if best_hit not in resistance_genes:
        resistance_genes[best_hit] = 0
    resistance_genes[best_hit] += 1
in_file.close()
    
#read SSN features
SSN_list = readFeatureFiles.readFeatureList(data_path+"gene_feature_matrices/SSN_list.txt")
for i in range(0, len(SSN_list)):
    SSN_list[i] = SSN_list[i].replace("\r","")
    SSN_list[i] = SSN_list[i].replace("\n","")

included_SSN_clusters =  {}
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
test_SSN_feature_matrix = SSN_tools.generateSSNFeatureMatrix([antismash_infilename], SSN_pfam_names, SSN_list, included_SSN_clusters, "blastp",cluster_name, data_path) 

#make the feature matrices for the cluster
training_features = np.concatenate((training_pfam_features, training_card_features), axis=1)
training_features = np.concatenate((training_features,  training_smCOG_features), axis=1)
training_features = np.concatenate((training_features,  training_CDS_features), axis=1)
training_features = np.concatenate((training_features,  training_SSN_features), axis=1)
training_features = np.concatenate((training_features,  training_pks_nrps_type_features), axis=1)
training_features = np.concatenate((training_features,  training_pk_signature_features), axis=1)
training_features = np.concatenate((training_features,  training_pk_minowa_features), axis=1)
training_features = np.concatenate((training_features,  training_pk_consensus_features), axis=1)
training_features = np.concatenate((training_features,  training_nrp_stachelhaus_features), axis=1)
training_features = np.concatenate((training_features,  training_nrp_nrpspredictor_features), axis=1)
training_features = np.concatenate((training_features,  training_nrp_pHMM_features), axis=1)
training_features = np.concatenate((training_features,  training_nrp_predicat_features), axis=1)
training_features = np.concatenate((training_features,  training_nrp_sandpuma_features), axis=1)
original_training_features = training_features


test_features = np.zeros((1, training_features.shape[1]))
i = 0
(test_features, i) = tools.addToFeatureMatrix(test_features, i, pfam_counts, used_pfam_list)
(test_features, i) = tools.addToFeatureMatrix(test_features, i, resistance_genes, used_resistance_genes_list)
(test_features, i) = tools.addToFeatureMatrix(test_features, i, smCOGs, used_smCOG_list)
(test_features, i) = tools.addToFeatureMatrix(test_features, i, CDS_motifs, used_CDS_list)


test_features[0,i:i+test_SSN_feature_matrix.shape[1]] = test_SSN_feature_matrix
i += test_SSN_feature_matrix.shape[1]



(test_features, i) = tools.addToFeatureMatrix(test_features, i, pks_nrp_subtypes, used_pks_nrps_type_list)
(test_features, i) = tools.addToFeatureMatrix(test_features, i, pk_monomers_signature, used_pk_signature_list)
(test_features, i) = tools.addToFeatureMatrix(test_features, i, pk_monomers_minowa, used_pk_minowa_list)
(test_features, i) = tools.addToFeatureMatrix(test_features, i, pk_monomers_consensus, used_pk_consensus_list)
(test_features, i) = tools.addToFeatureMatrix(test_features, i, nrp_monomers_stachelhaus, used_nrp_stachelhaus_list)
(test_features, i) = tools.addToFeatureMatrix(test_features, i, nrp_monomers_nrps_predictor, used_nrp_nrps_predictor_list)
(test_features, i) = tools.addToFeatureMatrix(test_features, i, nrp_monomers_pHMM, used_nrp_pHMM_list)
(test_features, i) = tools.addToFeatureMatrix(test_features, i, nrp_monomers_predicat, used_nrp_predicat_list)
(test_features, i) = tools.addToFeatureMatrix(test_features, i, nrp_monomers_sandpuma, used_nrp_sandpuma_list)

test_features_out =open("features/" + antismash_infilename[antismash_infilename.rfind("/"):antismash_infilename.rfind(".")]+".csv",'w')
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
