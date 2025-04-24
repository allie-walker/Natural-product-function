# -*- coding: utf-8 -*-
"""
Created on Thu Apr 24 12:27:14 2025

@author: Allison Walker
Scripts for getting feature matrices for antiSMASH6 output
"""

import os
from Bio import SeqIO
import readInputFiles

#parameters
feature_count_threshold = 5 #features must occur more then this number of times in dataset to be included
include_RGI = False
input_dir = "antiSMASH6_output/"
output_dir = "feature_matrices/antismash6_tigrfam/"
RGI_input_dir = "RGI_output/"

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

#process files
for f in sorted(os.listdir(input_dir),key=str.lower):
    if ".gbk" not in f or ".gb" not in f: #skip files that do not have genbank extension
        continue
    record = SeqIO.read(open(input_dir + f, 'r'),"genbank")
    features = record.features
    cluster_name = f[0:f.index(".")]
    (cluster_pfam_counts, cluster_CDS_motifs, cluster_smCOGs, cluster_NRPS_PKS,cluster_tigrfam_counts)  = readInputFiles.readAntismash6(features)
    pfam_counts[cluster_name] = cluster_pfam_counts
    smCOGs[cluster_name] = cluster_smCOGs
    CDS_motifs[cluster_name] = cluster_CDS_motifs
    NRPS_PKS[cluster_name] = cluster_NRPS_PKS
    tigrfam_counts[cluster_name]  =cluster_tigrfam_counts
    
#make feature and total count list
(pfam_list, pfam_total_counts) = makeCountsandList(pfam_counts)
(smCOG_list, smCOG_total_counts) = makeCountsandList(smCOGs)
(CDS_motif_list, CDS_total_counts) = makeCountsandList(CDS_motifs)
(NRPS_PKS_list, NRPS_PKS_total_counts) = makeCountsandList(NRPS_PKS)
(tigrfam_list, tigrfam_total_counts) = makeCountsandList(tigrfam_counts)

#filter low frequency features
smCOG_list = filterLowCounts(smCOG_list, smCOG_total_counts, feature_count_threshold)
pfam_list = filterLowCounts(pfam_list, pfam_total_counts, feature_count_threshold)
CDS_motif_list = filterLowCounts(CDS_motif_list, CDS_total_counts, feature_count_threshold)
NRPS_PKS_list = filterLowCounts(NRPS_PKS_list, NRPS_PKS_total_counts, feature_count_threshold)
tigrfam_list = filterLowCounts(tigrfam_list, tigrfam_total_counts, feature_count_threshold)

#write data to files
cluster_list_out = open(output_dir + "cluster_list_antiSMASH.txt", 'w')
for cluster in CDS_motifs:
    cluster_list_out.write(cluster + "\n")   
cluster_list_out.close()

writeFeatureMatrixFile(pfam_counts, pfam_list, output_dir + "features/","PFAM")
writeFeatureMatrixFile(CDS_motifs, CDS_motif_list, output_dir + "features/","CDS_motifs")
writeFeatureMatrixFile(smCOGs, smCOG_list, output_dir + "features/", "SMCOG")
writeFeatureMatrixFile(NRPS_PKS, NRPS_PKS_list, output_dir + "features/", "NRPS_PKS")
writeFeatureMatrixFile(tigrfam_counts, tigrfam_list, output_dir + "features/", "TIGR_FAM")

#RGI parameters
e_value_threshold = 0.1 #e value must be less than this threshold for resistance marker to be considered a feature
threshold = 5 #number of times feature must occur in dataset to be included

