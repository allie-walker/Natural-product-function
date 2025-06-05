# -*- coding: utf-8 -*-
"""
Created on Tue Mar  4 21:08:52 2025

@author: Allison Walker
Rank BGCs by activity prediction scores
"""
from os import listdir
from os.path import isfile, join, exists
import math
import numpy as np
import math
import argparse

parser = argparse.ArgumentParser(description='+Ranks all BGCs in directory by prediction score')
parser.add_argument('indir', type=str, help='directory with predictions')
parser.add_argument('activity',type=str, choices=['antibacterial','antigrampos','antigramneg','antieuk','antitumor_cytotoxic','antifungal'],help='activity to sort by')
parser.add_argument('-o','--outfile', type=str, help='file to write results to',default="ranked_predictions.csv")
args = parser.parse_args()

predictions_dir = args.indir
activity = args.activity
outfile_base = args.outfile

def readActivityPredictionsFile(file, genome_name, avg_probabilities):
    infile = open(file)
    i = 0
    for line in infile:
        split_line = line.split()
        if i == 1:
            antibac_avg = (float(split_line[2]) + float(split_line[6])+float(split_line[9]))/3
        if i == 3:
            antigrampos_avg = (float(split_line[2]) + float(split_line[6])+float(split_line[9]))/3
        if i == 5:
            antigramneg_avg = (float(split_line[2]) + float(split_line[6])+float(split_line[9]))/3
        if i == 7:
            antieuk_avg = (float(split_line[2]) + float(split_line[6])+float(split_line[9]))/3
        if i == 9:
            antifungal_avg = (float(split_line[2]) + float(split_line[6])+float(split_line[9]))/3
        if i == 11:
            antitumor_avg = (float(split_line[2]) + float(split_line[6])+float(split_line[9]))/3
        i += 1
    return antibac_avg, antigrampos_avg, antigramneg_avg, antieuk_avg, antifungal_avg, antitumor_avg

def getPredictionFilesInDirectory(path, avg_probabilities):
    antibac_probs = {}
    antigrampos_probs = {}
    antigramneg_probs = {}
    antieuk_probs = {}
    antifungal_probs = {}
    antitumor_probs = {}

    for f in listdir(path):
        if isfile(join(path, f)):
            if ".txt" in f:
                try:
                 antibac_avg, antigrampos_avg, antigramneg_avg, antieuk_avg, antifungal_avg, antitumor_avg = readActivityPredictionsFile(join(path, f),  path[path.rfind("/")+1:len(path)], avg_probabilities)
                 avg_probabilities["antibacterial"][join(path, f)] = antibac_avg
                 avg_probabilities["antigrampos"][join(path, f)] = antigrampos_avg
                 avg_probabilities["antigramneg"][join(path, f)] = antigramneg_avg
                 avg_probabilities["antieuk"][join(path, f)] = antieuk_avg
                 avg_probabilities["antifungal"][join(path, f)] = antifungal_avg
                 avg_probabilities["antitumor_cytotoxic"][join(path, f)] = antitumor_avg
                except:
                    continue
        else:
            if "genometools" in f or "glimmer" in f:
                continue
            files = getPredictionFilesInDirectory(join(path, f), avg_probabilities)
    return avg_probabilities

avg_probabilities = {}
avg_probabilities["antibacterial"] = {}
avg_probabilities["antigrampos"] = {}
avg_probabilities["antigramneg"] = {}
avg_probabilities["antieuk"] = {}
avg_probabilities["antifungal"] = {}
avg_probabilities["antitumor_cytotoxic"] = {}
avg_probabilities = getPredictionFilesInDirectory(predictions_dir, avg_probabilities)

sorted_list = sorted(avg_probabilities[activity].items(), key=lambda item: item[1], reverse=True)
outfile = open(outfile_base,'w')

outfile.write("clustername,antibacterial,antigrampos,antigramneg,antieuk,antifungal,antitumor\n")
for cluster, score in sorted_list:
    outfile.write(cluster + ",")
    outfile.write(str(avg_probabilities["antibacterial"][cluster]) + ",")
    outfile.write(str(avg_probabilities["antigrampos"][cluster]) + ",")
    outfile.write(str(avg_probabilities["antigramneg"][cluster]) + ",")
    outfile.write(str(avg_probabilities["antieuk"][cluster]) + ",")
    outfile.write(str(avg_probabilities["antifungal"][cluster]) + ",")
    outfile.write(str(avg_probabilities["antitumor_cytotoxic"][cluster]) + "\n")
    
outfile.close()
