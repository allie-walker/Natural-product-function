# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 14:35:52 2019

@author: Allison Walker
Writes SLURM scripts to perform activity prediction
"""

from os.path import isfile, join, exists
from os import listdir
import math
import argparse

parser = argparse.ArgumentParser(description='Writes SLURM scripts for running activity prediction')
parser.add_argument('antismash_indir', type=str, help='path to directory with antiSMASH results')
parser.add_argument('antismash_version', type=int, choices=[4,5,6,7,8], help='version of antiSMASH inputs')
parser.add_argument('prediction_version', type=int, choices=[1,2], help='version of prediction method to use')
parser.add_argument('outdir', type=str, help='output directory')
parser.add_argument('files_per_job', type=int, help='the number of files to include per script')
parser.add_argument('-r','--rgi_indir', type=str, default=None)
parser.add_argument('-rv','--rgi_version', type=int,  choices=[None, 3,5, 6], default=None)
parser.add_argument('-t','--threads',type=int,default=1,help='Number of CPUs to request in jobs')
parser.add_argument('-m','--mem',type=int,default=8, help='Amount of RAM to request in jobs, in GB')
parser.add_argument('-hr','--hours',type=int,default=24,help='Hours to run job')

args = parser.parse_args()
antismash_indir = args.antismash_indir + "/"
rgi_indir = args.rgi_indir
if rgi_indir != None:
    rgi_indir += "/"
outdir = args.outdir + "/"
rgi_version = args.rgi_version
antismash_version = args.antismash_version
prediction_version = args.prediction_version

files_per_job = args.files_per_job
threads = args.threads
mem = args.mem
hours = args.hours

days = math.floor(hours/24)
extra_hours = hours%24

#check to make sure provided options meet requirements
if prediction_version == 1:
    if rgi_indir == "None":
        print("RGI indir is required for version 1 of the prediction software")
        exit()
    if rgi_version == "None":
        print("RGI version is required for version 1 of the prediction software")
        exit()
    if antismash_version > 5:
        print("Version 1 of the prediction software is only compatible with antiSMASH 4 and 5")
        exit()
    if rgi_version > 5:
        print("Version 1 of the prediction software is only compatible with RGI 3 and 5")
        exit()
if rgi_version != "None" and rgi_indir == "None":
    print("an input directory for RGI is required if you would like to use RGI features, if you do not want to use RGI features remove the argument for RGI version")
    exit()
if rgi_version == "None" and rgi_indir != "None":
    print("an iRGI version is required to use RGI features, if you do not want to use RGI features please delete the argument for the RGI directory")
    exit()
    
def getAllClusterListsInDirectory(path, files):
    for f in listdir(path):
        if isfile(join(path, f)):
            if ("cluster" in f or "region" in f) and (".gbk" in f or ".gb" in f or ".gbff" in f):
                fullpath = join(path, f)
                genome_name = fullpath[0:fullpath.rfind("/")]
                genome_name = genome_name[genome_name.rfind("/")+1:len(genome_name)]
                #check if output exists
                if not exists(outdir + genome_name + f[f.rfind("/")+1:f.rfind(".")]): 
                    files.append(join(path, f))
        else:
            if "genometools" in f or "glimmer" in f or "slurm" in f:
                continue
            files = getAllClusterListsInDirectory(join(path, f), files)
    return files


genome_names = []
cluster_list = []
cluster_list = getAllClusterListsInDirectory(antismash_indir, cluster_list)
total_files = len(cluster_list)
split = int(math.ceil(total_files/files_per_job))
j = 0
for i in range(0, split):
    outfile = open("run_all_predictions" + str(i) + ".sh",'w')
    outfile.write("#!/bin/bash \n")
    outfile.write("#SBATCH --nodes=1 \n")
    outfile.write("#SBATCH --ntasks=16 \n")
    outfile.write("#SBATCH -t" + str(int(days)) + "-" + str(extra_hours) + ":00:00" + "\n")
    outfile.write("#SBATCH --output=run_predictions_" + str(i) + ".out\n")
    outfile.write("#SBATCH --mem=" +str(mem) + "G                          # Memory total in MB (for all cores) \n")
    
    for k in range(0, files_per_job):
        if j > total_files -1:
            break
        cluster_name = cluster_list[j]
        
        genome_name = cluster_name[0:cluster_name.rfind("/")]
        genome_name = genome_name[genome_name.rfind("/")+1:len(genome_name)]
        cluster_name = cluster_name[cluster_name.rfind("/")+1:cluster_name.rfind('.')]
        antismash_input = cluster_list[j]
        if rgi_indir != None:
            rgi_input =  rgi_indir + genome_name + "/" + cluster_name + ".txt"
        if not exists(outdir + genome_name):
            outfile.write("mkdir " + outdir + genome_name + "\n")
        if prediction_version == 1:
            outfile.write("cluster_function_prediction.py '" + antismash_input + "' '" + rgi_input + "' --seed 0 --output " + outdir + genome_name + "/"+ " --antismash_version " + str(antismash_version) +" --rgi_version" + str(rgi_version) + " --no_SSN\n")
        elif prediction_version == 2:
            if rgi_indir != None:
                outfile.write("cluster_function_prediction.py '" + antismash_input + "' --rgi_results '" + rgi_input + "' --seed 0 --output " + outdir + genome_name + "/"+ " --antismash_version " + str(antismash_version) +" --rgi_version " + str(rgi_version) + "\n")
            else:
                outfile.write("cluster_function_prediction.py '" + antismash_input + "' --seed 0 --output " + outdir + genome_name + "/"+ " --antismash_version " + str(antismash_version) + "\n")
                
        j += 1
    outfile.close()
                        

