# -*- coding: utf-8 -*-
"""
Created on Fri Feb  4 12:13:47 2022

@author: Allison Walker

Make a graph showing the prediction results
"""
from os import listdir
from os.path import isfile, join, exists
import math
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
from Bio import SeqIO
import matplotlib.patches as mpatches

parser = argparse.ArgumentParser(description='Makes graph showing prediction results')
parser.add_argument('indir', type=str, help='directory with predictions')
parser.add_argument('-a', '--antismash_dir', type=str, help='directory with antismash results, required if you are using the color by BGC')
parser.add_argument('-o','--outdir',type=str, default="./",help="directory to write graph to")
parser.add_argument('-s','--sort',type=str,choices=["max_prob","max_actives"],default="max_prob",help="sort each genome graph, by max probabilty BGC or by total number of predicted active BGCs")
parser.add_argument('-c','--color',type=str,choices=["rainbow","type"],default="rainbow",help="color in rainbow or by BGC class")

args = parser.parse_args()
indir = args.indir + "/"
outdir = args.outdir + "/"
antismash_dir = args.antismash_dir
sort_mode = args.sort #max_prob to sort by maximum probability for entire genome
color = args.color
if color == "type" and antismash_dir == None:
    print("ERROR: need to specify antismash directory for color by type mode")
    exit()
elif antismash_dir != None:
    antismash_dir += "/"

def getPredictions(file):
    infile = open(file)
    avg_preds = {}
    indv_preds = {}
    activities = ["antibac", "antigrampos", "antigramneg", "antieuk", "antifungal", "antitumor"]
    count = 0
    for line in infile:
        if "probabilities" in line:
            continue
        split_line = line.split()
        tree = float(split_line[2])
        log = float(split_line[6])
        svm = float(split_line[9])
        indv_preds[activities[count]] = [tree, log, svm]
        avg_preds[activities[count]] = (tree + log + svm)/3
        count += 1

    return (avg_preds, indv_preds)
        
        

def getAllFilesInDirectory(path, files):
    for f in listdir(path):
        if isfile(join(path, f)):
            #print join(path, f)
            if ".txt" in f and ("cluster" in f or "region" in f) and ".log" not in f:
                #check if output exists
                fullpath = join(path,f)
                genome_name = path[path.rfind("/"):len(path)]
                files.append((genome_name, fullpath))
        else:
            if "genometools" in f or "glimmer" in f:
                continue
            files = getAllFilesInDirectory(join(path, f), files)
    return files

def getBGCType(gbk_filename):
    bgc_identites = {"RiPP":["bacteriocin","LAP","lanthipeptide", "TfuA-related","lassopeptide", "linaridin"],\
                     "NRPS":["NRPS", "NRPS-like"],\
                     "terpene":["terpene"], \
                      "PKS":["T1PKS","T2PKS","T3PKS","PKS-like","transAT-PKS"], \
                       "siderophore":["siderophore"],
                       "oligosaccharide":["oligosaccharide"]}
    record = SeqIO.read(open(gbk_filename),"genbank")
    for feature in record.features:
        if feature.type == "cand_cluster":
            if feature.qualifiers["kind"][0] != "single":
                bgc_type = "hybrid/other"
                break
            elif len(feature.qualifiers["product"]) > 1:
                bgc_type = "hybrid/other"
                break
            else:
                bgc_type = feature.qualifiers["product"][0]
    final_bgc_type = "hybrid/other"
    for i in bgc_identites:
        if bgc_type in bgc_identites[i]:
            final_bgc_type = i
    #if final_bgc_type == "hybrid/other":
     #   print(bgc_type)        
    return final_bgc_type


predictions_list = []
predictions_list = getAllFilesInDirectory(indir, predictions_list)
bgc_types = {}

if color == "type":
    for f in predictions_list:
        genome = f[0]
        if genome not in bgc_types:
            bgc_types[genome] = {}
        gbk_filename = antismash_dir + genome + f[1][f[1].rfind("/"):f[1].rfind(".")] + ".gbk"
        fname = f[1][f[1].rfind("/") + 1:f[1].rfind(".")]
        region = fname[fname.find("region")+6:len(fname)]
        contig = fname[0:fname.find(".")]
        fname = contig + "_" + region
        bgc_types[genome][fname] = getBGCType(gbk_filename)



predictions_avg = {}
predictions_idv = {}
predictions_max_genome = {}
predicted_actives = {}
order = {}
activities = ["antibac", "antigrampos", "antigramneg", "antieuk", "antifungal", "antitumor"]
possible_bgc_types = ["RiPP","NRPS","PKS","terpene","siderophore","oligosaccharide","hybrid/other"]

for a in activities:
    predictions_avg[a] = {}
    predictions_idv[a] = {}
    predicted_actives[a] = {}
    predictions_max_genome[a] = {}
    predictions_max_genome[a]= {}
    order[a] = []
for (genome, file) in predictions_list:
        fname = file[file.rfind("/") + 1:file.rfind(".")]
        region = fname[fname.find("region")+6:len(fname)]
        contig = fname[0:fname.find(".")]
        fname = contig + "_" + region
        for a in activities:
            if genome not in predictions_avg[a]:
                predictions_avg[a][genome] = {}
                predictions_idv[a][genome] = {}
        (avg, indv) = getPredictions(file)
        for a in avg:
            predictions_avg[a][genome][fname] = avg[a]
            predictions_idv[a][genome][fname] = indv[a]


for genome in predictions_avg[a]:
    for a in activities:
        predicted_actives[a][genome] = 0
        for f in predictions_avg[a][genome]:
            if predictions_avg[a][genome][f] > 0.5:
                predicted_actives[a][genome] += 1
        predictions_max_genome[a][genome] = sorted(predictions_avg[a][genome].items(), key=lambda item: item[1], reverse=True)[0][1]


for a in activities:    
    if sort_mode == "max_prob":
        order[a] = [g for g, m in sorted(predictions_max_genome[a].items(), key=lambda item: item[1], reverse = True)]
    else:
        order[a] = [g for g, m in sorted(predicted_actives[a].items(), key=lambda item: item[1], reverse = True)]


sns.set(style="whitegrid")     

for a in activities:
    fig, axs = plt.subplots(len(predictions_avg[a]),figsize=(15,5*len(predictions_avg[a])))   
    i = 0
    for g in order[a]:
        sorted_by_label = dict(sorted(predictions_avg[a][g].items(), key=lambda item: item[0]))
        
        x_values = []
        y_values = []
        if color == "rainbow":
            for l in sorted_by_label:
                for v in predictions_idv[a][g][l]:
                    x_values.append(l)
                    y_values.append(v)
            if len(predictions_avg[a]) > 1:
                sns.stripplot(ax=axs[i], x=x_values, y=y_values,  color="0", alpha=.35)
                sns.barplot(ax=axs[i],x=x_values, y=y_values, capsize=.1, errorbar="sd")
                axs[i].set_ylim(0, 1.0)
                axs[i].set_xticklabels(axs[i].get_xticklabels(),rotation = 90)
                axs[i].set_title(g[1:len(g)])
                axs[i].set_xlabel("BGC",fontsize=15)
                axs[i].set_ylabel("probability of activity",fontsize=15)
                i += 1        
            else:
                sns.stripplot(ax=axs, x=x_values, y=y_values,  color="0", alpha=.35)
                sns.barplot(ax=axs, x=x_values, y=y_values, capsize=.1, errorbar="sd")
                axs.set_ylim(0, 1.0)
                axs.set_xticklabels(axs.get_xticklabels(),rotation = 90)
                axs.set_title(g[1:len(g)])
                axs.set_xlabel("BGC",fontsize=15)
                axs.set_ylabel("probability of activity",fontsize=15)
        else:
            qualitative_colors = sns.color_palette("Set3", 7)
            bar_colors = []
            for l in sorted_by_label:
                bgc_type = bgc_types[g][l]
                type_index = possible_bgc_types.index(bgc_type)
                bar_colors.append(qualitative_colors[type_index])
            for l in sorted_by_label:
                for v in predictions_idv[a][g][l]:
                    x_values.append(l)
                    y_values.append(v)
            if len(predictions_avg[a]) > 1:
                sns.stripplot(ax=axs[i], x=x_values, y=y_values,  color="0", alpha=.35)
                sns.barplot(ax=axs[i],x=x_values, y=y_values, capsize=.1, errorbar="sd", palette=bar_colors)
                axs[i].set_ylim(0, 1.0)
                axs[i].set_xticklabels(axs[i].get_xticklabels(),rotation = 90)
                axs[i].set_title(g[1:len(g)])
                axs[i].set_xlabel("BGC",fontsize=15)
                axs[i].set_ylabel("probability of activity",fontsize=15)
                i += 1    
            else:
                sns.stripplot(ax=axs, x=x_values, y=y_values,  color="0", alpha=.35)
                sns.barplot(ax=axs,x=x_values, y=y_values, capsize=.1, errorbar="sd", palette=bar_colors)
                axs.set_ylim(0, 1.0)
                axs.set_xticklabels(axs.get_xticklabels(),rotation = 90)
                axs.set_title(g[1:len(g)])
                axs.set_xlabel("BGC",fontsize=15)
                axs.set_ylabel("probability of activity",fontsize=15)
    if color == "type":
        patches = []
        for i in range(0, len(possible_bgc_types)):
            patches.append(mpatches.Patch(color=qualitative_colors[i], label=possible_bgc_types[i]))
        if len(predictions_avg[a]) > 1:
            axs[-1].legend(handles=patches,fontsize=15)
            leg = axs[-1].get_legend()
        else:
            axs.legend(handles=patches,fontsize=15,loc='center right')
            leg = axs.get_legend()

    plt.tight_layout()
    plt.savefig(outdir + a + '_predictions.png')
    
