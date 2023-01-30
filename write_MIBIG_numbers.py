# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 11:16:11 2022

@author: allison

This script extracts MIBIG identifiers from clusters and makes list
"""

#cluster_list_input = open("feature_matrices/antismash4rgi3/cluster_list.txt")
cluster_list_input = open("holdout_set/classifications/cluster_list.txt")
#sequence_file_dir = "antiSMASH_output/"
sequence_file_dir = "holdout_set/antiSMASH_output/"
#output_file = open("feature_matrices/antismash4rgi3/MIBIG_numbering.txt", 'w')
output_file = open("holdout_set/MIBIG_numbering.txt", 'w')


cluster_list = []
MIBIG_numbering = {}
for line in cluster_list_input:
    cluster_list.append(line.replace("\n", ""))
    
for c in cluster_list:
    try:
        bgc_in = open(sequence_file_dir + c + ".gbk")
    except:
        MIBIG_numbering[c] = "NONE"
        continue
    for line in bgc_in:
        MIBIG_numbering[c] = line.split()[1]
        break
    
for c in cluster_list:
    output_file.write(c + "," + MIBIG_numbering[c] + "\n")
output_file.close()