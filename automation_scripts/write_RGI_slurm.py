# -*- coding: utf-8 -*-
"""
Created on Sat May 31 22:08:32 2025

@author: Allison Walker

Writes a slurm script to run RGI on all fasta in a directory
"""

from os import listdir
from os.path import isfile, join, exists
import argparse
import math

parser = argparse.ArgumentParser(description='Writes SLURM scripts for running RGI on fastas in directory')
parser.add_argument('antiSMASH_dir', type=str, help='directory with fastas')
parser.add_argument('RGI_version', type=int,  choices=[3,5,6], help='RGI verison to use')
parser.add_argument('output_dir',type =str,help="output directory name")
parser.add_argument('-t','--threads',type=int,default=16,help='Number of CPUs to request in jobs')
parser.add_argument('-m','--mem',type=int,default=8,help='Amount of RAM to request in jobs, in GB')
parser.add_argument('-hr','--hours',type=int,default=24,help='Hours to run job')

args = parser.parse_args()
card_dir = args.output_dir + "/"
indir = args.antiSMASH_dir + "/"
RGI_version = args.RGI_version 
threads = args.threads
mem = args.mem
hours = args.hours

days = math.floor(hours/24)
extra_hours = hours%24


def getAllFilesInDirectory(path, files):
    for f in listdir(path):
        if isfile(join(path, f)):
            #print join(path, f)
            if ".fasta" in f and ("cluster" in f or "region" in f) and ".log" not in f:
                #check if output exists
                fullpath = join(path,f)
                genome_name = path[path.rfind("/"):len(path)]
                if not exists(indir + card_dir + genome_name + f[f.rfind("/")+1:f.rfind(".")]): #fix this to check if RGI has been run already!
                    files.append(join(path, f))
            #print files
        else:
            if "genometools" in f or "glimmer" in f:
                continue
            files = getAllFilesInDirectory(join(path, f), files)
    return files


outfile = open("run_all_rgi.sh",'w')
fasta_list = []
fasta_list = getAllFilesInDirectory(indir, fasta_list)



outfile.write("#!/bin/bash \n")
outfile.write("#SBATCH --nodes=1 \n")
outfile.write("#SBATCH --ntasks=16 \n")
outfile.write("#SBATCH -t "  + str(int(days)) + "-" + str(extra_hours) + "\n")
outfile.write("#SBATCH --output=run_rgi5.out\n")
outfile.write("#SBATCH --mem=8G                          # Memory total in MB (for all cores) \n")
for fasta in fasta_list:
    genome_name = fasta[0:fasta.rfind("/")]
    genome_name = genome_name[genome_name.rfind("/")+1:len(genome_name)]
    if not exists(indir + card_dir + genome_name):
        outfile.write("mkdir " + indir + card_dir + genome_name +"\n")
    outfile.write("rgi main -i " + fasta + " -o " + indir + card_dir + genome_name + "/" +fasta[fasta.rfind("/")+1:fasta.rfind(".")] + " --include_loose  --clean\n")
    
outfile.close()
    

