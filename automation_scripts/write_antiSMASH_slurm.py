# -*- coding: utf-8 -*-
"""
Created on Fri May 30 21:43:56 2025

@author: Allison Walker

Writes a slurm script to run antiSMASH on all genbank in a directory
"""

import argparse
import os
import math

parser = argparse.ArgumentParser(description='Writes SLURM scripts for running antiSMASH on genomes in directory')
parser.add_argument('genome_input_dir', type=str, default='directory with genomes')
parser.add_argument('antiSMASH_version', type=int,  choices=[4,5,6,7,8], default='directory with genomes')
parser.add_argument('output_dir',type =str,default="output directory name")
parser.add_argument('files_per_job', type=int, default='the number of files to include per script')
parser.add_argument('-t','--threads',type=int,default=16,help='Number of CPUs to request in jobs')
parser.add_argument('-m','--mem',type=int,default=8,help='Amount of RAM to request in jobs, in GB')
parser.add_argument('-hr','--hours',type=int,default=24,help='Hours to run job')

args = parser.parse_args()
genome_input_dir = args.genome_input_dir + "/"
antiSMASH_version = args.antiSMASH_version
files_per_job = args.files_per_job
threads = args.threads
output_dir = args.output_dir + "/"
mem = args.mem
hours = args.hours

days = math.floor(hours/24)
extra_hours = hours%24

input_list = []
already_run = []


for f in os.listdir(genome_input_dir): 
    if (".gbff" in f or ".gb" in f or ".gbk" in f):
        input_list.append(f)

for f in os.listdir(output_dir):
    already_run.append(f)  

new_input = []
for i in input_list:
	if i[0:i.rfind(".")] not in already_run:
		new_input.append(i)
input_list = new_input

total_files = len(input_list)

split = int(math.ceil(total_files/files_per_job))
j = 0
for i in range(0, split):
    outfile = open("run_all_antismash_" +str(i) + ".sh",'w')
    outfile.write("#!/bin/bash \n")
    outfile.write("#SBATCH --nodes=1 \n")
    outfile.write("#SBATCH --ntasks=" + str(threads) +" \n")
    outfile.write("#SBATCH -t" + str(int(days)) + "-" + str(extra_hours) + ":00:00" + "\n")
    outfile.write("#SBATCH --output=run_antismash_" + str(i) + ".out\n")
    outfile.write("#SBATCH --mem=8G                          # Memory total in MB (for all cores) \n")
    for k in range(0, files_per_job):
        if j > total_files -1:
            break
        infile_name = input_list[j]
        outfile.write("mkdir '" + output_dir  +infile_name[0:infile_name.rfind(".")] +"'\n")
        if antiSMASH_version == 4:
            outfile.write("antismash  --verbose --outputfolder '" + output_dir  +infile_name[0:infile_name.rfind(".")] + "' --statusfile status.txt --full-hmmer --borderpredict '" + genome_input_dir + infile_name + "'\n")
        elif antiSMASH_version == 5:
            outfile.write("antismash  --output-dir '" + output_dir  +infile_name[0:infile_name.rfind(".")] + "' --genefinding-tool prodigal --fullhmmer --cb-knownclusters '" + genome_input_dir + infile_name + "'\n")
        else:
            outfile.write("antismash  --output-dir '" + output_dir  +infile_name[0:infile_name.rfind(".")] + "' --genefinding-tool prodigal --fullhmmer --cb-knownclusters --cc-mibig --tigrfam --asf --pfam2go --smcog-trees '" + genome_input_dir + infile_name + "'\n")
        j += 1     
    outfile.close()