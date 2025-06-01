# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 14:01:09 2019

@author: Allison Walker
Converts antiSMASH output cluster genbank files to fasta files
"""

from Bio import SeqIO
from os import listdir
from os.path import isfile, join, exists
import argparse

parser = argparse.ArgumentParser(description='Converts antiSMASH cluster outputs to fasta file, will recursively search for cluster genbank files in the given input directory')
parser.add_argument('antismash_dir', type=str, default='directory with antismash results')
args = parser.parse_args()
directory_name = args.antismash_dir

def getAllFilesInDirectory(path, files):
    for f in listdir(path):
        if isfile(join(path, f)):
            if (".gbff" in f or ".gb" in f or ".gbk" in f) and ("region" in f or "cluster" in f):
                #check if output exists
                fullpath = join(path,f)
                if not exists(fullpath[0:fullpath.rfind(".")] + ".fasta"): 
                    files.append(join(path, f))
        else:
            if "genometools" in f or "glimmer" in f:
                continue
            files = getAllFilesInDirectory(join(path, f), files)
    return files

files = []
files = getAllFilesInDirectory(directory_name, files)

for f in files:
    record = SeqIO.read(open(f, 'rU'),"genbank")
    fasta_out = open(f[0:f.rfind(".")]+".fasta",'w')
    fasta_out.write(">")
    fasta_out.write(record.id + "|"+ record.description + "\n")
    for s in record.seq:
        fasta_out.write(s)
    fasta_out.write("\n")
    fasta_out.close()
