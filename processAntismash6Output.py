# -*- coding: utf-8 -*-
"""
Created on Thu Apr 24 12:27:14 2025

@author: Allison Walker
Scripts for getting feature matrices for antiSMASH6-8 output
"""

import os
from Bio import SeqIO
import readInputFiles
import argparse
import subprocess

#parameters
parser = argparse.ArgumentParser(description='Extracts features from antiSMASH outputs version 6-8 and RGI versions 3 and 5')
parser.add_argument('antismash_input_dir', type=str, default='directory with antiSMASH genbanks')
parser.add_argument('ouptut_dir', type=str, default='name of output directory to be created within feature_matrices')
#TODO: also add 6?
parser.add_argument('-r', '--rgi_version', type=str, default='None', choices=['None',"3", "5"], help='version of RGI, use None to not include RGI features in predictions') 
parser.add_argument('-rd','--rgi_directory',type=str,default='RGI_output',help='Directory containing RGI output files')
parser.add_argument('-t','--count_threshold',type=int,default=5,help='feature count threshold to be included as a feature')
parser.add_argument('-c', '--cluster_list_file',type=str,default='cluster_list.csv')
args = parser.parse_args()
input_dir = args.antismash_input_dir + "/"
output_dir = "feature_matrices/" + args.ouptut_dir +"/"
feature_count_threshold = args.count_threshold
RGI_type = args.rgi_version
RGI_input_dir = args.rgi_directory
cluster_list_file = args.cluster_list_file

#make output directories
if not os.path.isdir(output_dir):
    subprocess.run(["mkdir",output_dir])
if not os.path.isdir(output_dir + 'classifications'):
    subprocess.run(["mkdir",output_dir + 'classifications'])
if not os.path.isdir(output_dir + 'features'):
    subprocess.run(["mkdir",output_dir + 'features'])

def makeCountsandList(feature_dir):
    feature_list = []
    feature_total_counts = {}
    for c in feature_dir:
        for f in feature_dir[c]:
            if f not in feature_total_counts:
                feature_list.append(f)
                feature_total_counts[f] = 0
            feature_total_counts[f] += feature_dir[c][f]
    return (feature_list, feature_total_counts)

def filterLowCounts(feature_list, feature_counts, threshold):
    new_feature_list = []
    for f in feature_list:
        if feature_counts[f] >= threshold:
            new_feature_list.append(f)
    return new_feature_list

#write a matrix with counts of features to a file and a list of the features to a separate file
def writeFeatureMatrixFile(feature_dictionary, feature_list, output_dir, feature_name):
    matrix_outfile = open(output_dir + feature_name + ".csv", 'w')
    list_outfile = open(output_dir + feature_name + "_list.txt",'w')
    for cluster in feature_dictionary:
        for feature in feature_list:
            if feature in feature_dictionary[cluster]:
                matrix_outfile.write(str(feature_dictionary[cluster][feature]) + ",")
            else:
                matrix_outfile.write("0,")
        matrix_outfile.write("\n")
    matrix_outfile.close()
    
    for f in feature_list:
        list_outfile.write(f + "\n")
    list_outfile.close()
    
CDS_motifs = {}
CDS_motif_list = []
CDS_total_counts = {}
smCOGs = {}
smCOG_list = []
smCOG_total_counts = {}
pfam_list = []
pfam_id_list = []
pfam_counts = {}
NRPS_PKS = {}
NRPS_PKS_list = []
NRPS_PKS_total_counts = {}
tigrfam_counts = {}
tigrfam_list =[]
tigrfam_total_counts = {}
NRPS_PKS_substrate = {}
NRPS_PKS_substrate_list = []
#process files
for f in sorted(os.listdir(input_dir),key=str.lower):
    if ".gbk" not in f or ".gb" not in f: #skip files that do not have genbank extension
        continue
    record = SeqIO.read(open(input_dir + f, 'r'),"genbank")
    features = record.features
    cluster_name = f[0:f.index(".")]
    (cluster_pfam_counts, cluster_CDS_motifs, cluster_smCOGs, cluster_NRPS_PKS,cluster_tigrfam_counts,cluster_NRPS_PKS_substrate)  = readInputFiles.readAntismash6(features)
    pfam_counts[cluster_name] = cluster_pfam_counts
    smCOGs[cluster_name] = cluster_smCOGs
    CDS_motifs[cluster_name] = cluster_CDS_motifs
    NRPS_PKS[cluster_name] = cluster_NRPS_PKS
    tigrfam_counts[cluster_name]  =cluster_tigrfam_counts
    NRPS_PKS_substrate[cluster_name] = cluster_NRPS_PKS_substrate
    
#make feature and total count list
(pfam_list, pfam_total_counts) = makeCountsandList(pfam_counts)
(smCOG_list, smCOG_total_counts) = makeCountsandList(smCOGs)
(CDS_motif_list, CDS_total_counts) = makeCountsandList(CDS_motifs)
(NRPS_PKS_list, NRPS_PKS_total_counts) = makeCountsandList(NRPS_PKS)
(tigrfam_list, tigrfam_total_counts) = makeCountsandList(tigrfam_counts)
(NRPS_PKS_substrate_list, NRPS_PKS_substrate_total_counts) = makeCountsandList(NRPS_PKS_substrate)

#filter low frequency features
smCOG_list = filterLowCounts(smCOG_list, smCOG_total_counts, feature_count_threshold)
pfam_list = filterLowCounts(pfam_list, pfam_total_counts, feature_count_threshold)
CDS_motif_list = filterLowCounts(CDS_motif_list, CDS_total_counts, feature_count_threshold)
NRPS_PKS_list = filterLowCounts(NRPS_PKS_list, NRPS_PKS_total_counts, feature_count_threshold)
tigrfam_list = filterLowCounts(tigrfam_list, tigrfam_total_counts, feature_count_threshold)
NRPS_PKS_substrate_list = filterLowCounts(NRPS_PKS_substrate_list, NRPS_PKS_substrate_total_counts, feature_count_threshold)
#write data to files
cluster_list_out = open(output_dir + "cluster_list.txt", 'w')
for cluster in CDS_motifs:
    cluster_list_out.write(cluster + "\n")   
cluster_list_out.close()

writeFeatureMatrixFile(pfam_counts, pfam_list, output_dir + "features/","PFAM")
writeFeatureMatrixFile(CDS_motifs, CDS_motif_list, output_dir + "features/","CDS_motifs")
writeFeatureMatrixFile(smCOGs, smCOG_list, output_dir + "features/", "SMCOG")
writeFeatureMatrixFile(NRPS_PKS, NRPS_PKS_list, output_dir + "features/", "NRPS_PKS")
writeFeatureMatrixFile(tigrfam_counts, tigrfam_list, output_dir + "features/", "TIGR_FAM")
writeFeatureMatrixFile(NRPS_PKS_substrate, NRPS_PKS_substrate_list, output_dir + "features/", "NRPS_PKS_substrate")

#RGI parameters
e_value_threshold = 0.1 #e value must be less than this threshold for resistance marker to be considered a feature
threshold = 5 #number of times feature must occur in dataset to be included
bit_score_threshold = 40 #threshold for RGI5

if RGI_type == "3":
    resistance_genes = {}
    resistance_genes_list = []
    total_counts = {}
    for f in sorted(os.listdir(RGI_input_dir),key=str.lower):
        if ".txt" not in f:
            continue
        in_file = open(RGI_input_dir + f, 'r')
        cluster_name = f[0:f.index(".")]
        resistance_genes[cluster_name] = {}
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
                total_counts[best_hit] = 0
            total_counts[best_hit] += 1
            if best_hit not in resistance_genes[cluster_name]:
                resistance_genes[cluster_name][best_hit] = 0
            resistance_genes[cluster_name][best_hit] += 1
        in_file.close()
        
    new_resistance_genes_list = []
    for gene in resistance_genes_list:
        if total_counts[gene] >= threshold:
            new_resistance_genes_list.append(gene)
    resistance_genes_list = new_resistance_genes_list

    writeFeatureMatrixFile(resistance_genes, resistance_genes_list, output_dir + "features/", "CARD_gene")
    cluster_list = open(output_dir + "cluster_list_CARD.txt",'w')
    for cluster in resistance_genes:
        cluster_list.write(cluster + "\n")
    cluster_list.close()
    
elif RGI_type == "5":
    for f in sorted(os.listdir(RGI_input_dir),key=str.lower):
        if ".txt" not in f:
            continue
        in_file = open(RGI_input_dir + f, 'r')
        cluster_name = f[0:f.index(".")]
        resistance_genes[cluster_name] = {}
        for line in in_file:
            if "ORF_ID" in line:
                continue
            entries = line.split("\t")
            bit_score = float(entries[7])
            if bit_score < bit_score_threshold:
                continue
            best_hit = entries[8]
            if best_hit not in resistance_genes_list:
                resistance_genes_list.append(best_hit)
                total_counts[best_hit] = 0
                total_counts[best_hit] += 1
            if best_hit not in resistance_genes[cluster_name]:
                resistance_genes[cluster_name][best_hit] = 0
            resistance_genes[cluster_name][best_hit] += 1
        in_file.close()

    new_resistance_genes_list = []
    threshold = 5
    for gene in resistance_genes_list:
        if total_counts[gene] >= threshold:
            new_resistance_genes_list.append(gene)
    resistance_genes_list = new_resistance_genes_list

    resistance_output = open(output_dir + "features/" + "CARD_genes.csv",'w')
    cluster_list = open(output_dir + "cluster_list_CARD5.txt",'w')
    resistance_list_out = open(output_dir + "features/" + "CARD_gene_list.txt",'w')
    for cluster in resistance_genes:
        for gene in resistance_genes_list:
            if gene not in resistance_genes[cluster]:
                resistance_output.write("0,")
            else:
                #resistance_output.write("1,")
                resistance_output.write(str(resistance_genes[cluster][gene])+",")
        resistance_output.write("\n")
        cluster_list.write(cluster + "\n")
    cluster_list.close()
    resistance_output.close()
    for gene in resistance_genes_list:
        resistance_list_out.write(gene + "\n")
    resistance_list_out.close()

#read BGC classifications
cluster_info = open(cluster_list_file,'r',encoding='cp1252') #a spreadsheet with all clusters to be used and their classifications
is_antibacterial = {}
is_antifungal = {}
is_cytotoxic = {}
is_unknown = {}
targets_gram_pos = {}
targets_gram_neg = {}
for line in cluster_info:
    entries = line.split(",")
    cluster_name = entries[0]
    if entries[1] == "yes":
        is_antibacterial[cluster_name] = 1
    else:
        is_antibacterial[cluster_name] = 0
        
    if entries[2] == "yes":
        is_antifungal[cluster_name] = 1
    else:
        is_antifungal[cluster_name] = 0
    
    if "cytotoxic" in entries[3] or "antitumor" in entries[3] or "anti-tumor" in entries[3] or "proliferative" in entries[3]:
        is_cytotoxic[cluster_name] = 1
    else:
        is_cytotoxic[cluster_name] = 0
        
    if "unknown" in entries[3]:
        is_unknown[cluster_name] =1
    else:
        is_unknown[cluster_name]=0
        
    if "gram pos" in entries[7]:
        targets_gram_pos[cluster_name] = 1
    else:
        targets_gram_pos[cluster_name] = 0
    if "gram neg" in entries[7]:
        targets_gram_neg[cluster_name] = 1
    else:
        targets_gram_neg[cluster_name] = 0
        
#write classses
antibacterial_out = open(output_dir + "classifications/is_antibacterial.csv",'w')
antifungal_out = open(output_dir + "classifications/is_antifungal.csv",'w')
cytotoxic_out = open(output_dir + "classifications/is_cytotoxic.csv",'w')
unknown_out = open(output_dir + "classifications/is_unknown.csv",'w')
targets_gram_pos_out = open(output_dir + "classifications/targets_gram_pos.csv",'w')
targets_gram_neg_out = open(output_dir + "classifications/targets_gram_neg.csv",'w')

#read cluster list from features
cluster_list_input = open(output_dir + "/cluster_list.txt")
cluster_list = []
for line in cluster_list_input:
    c = line.replace("\n", "")
    cluster_list.append(c)
cluster_list_input.close()

#only write clusters that have feature data
for c in cluster_list:
    if c in is_antibacterial:
        antibacterial_out.write(str(is_antibacterial[c])+"\n")
        antifungal_out.write(str(is_antifungal[c]) + "\n")
        cytotoxic_out.write(str(is_cytotoxic[c]) + "\n")
        unknown_out.write(str(is_unknown[c]) + "\n")
        targets_gram_neg_out.write(str(targets_gram_neg[c])+"\n")
        targets_gram_pos_out.write(str(targets_gram_pos[c])+"\n")
    else:
        print("cluster " + c + " is not in class spreadsheet")
        antibacterial_out.write("0\n")
        antifungal_out.write("0\n")
        cytotoxic_out.write("0\n")
        unknown_out.write("1\n")
        targets_gram_neg_out.write("0\n")
        targets_gram_pos_out.write("0\n")
        
#check if there are clusters in spreadsheet missing in dataset
for c in is_antibacterial:
    if c not in cluster_list:
        print("cluster: " + c + " is in class speadsheet but not in antiSMASH files")
        
antibacterial_out.close()
antifungal_out.close()
cytotoxic_out.close()
unknown_out.close()
targets_gram_neg_out.close()
targets_gram_pos_out.close()