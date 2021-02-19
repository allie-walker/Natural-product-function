# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 11:14:09 2019

@author: Allison Walker
"""

from Bio import SeqIO
from Bio import Seq
from Bio import AlignIO
from Bio.Align.Applications import ClustalwCommandline
from Bio import Phylo
from Bio.Align.Applications import MuscleCommandline
from Bio import pairwise2
from xml.etree import ElementTree as ET
import numpy as np
import math
import subprocess

SSN_filenames = ["20799_thiolase_N_term_25_full_ssn.xgmml",
                 "21330_abc_35_full_ssn.xgmml",
                 "20870_acyl_transferase_pfam_100_full_ssn.xgmml",
                 "20923_aaa_30_full_ssn.xgmml",
                 "20928_abc2_30_full_ssn.xgmml",
                 "20929_acyl-coa_dehydrogenase_c_30_full_ssn.xgmml",
                 "20930_acyl-coa_dehydrogenase_n_20_full_ssn.xgmml",
                 "20931_alcohol_dehydrogenase_groes_20_full_ssn.xgmml",
                 "20932__alpha_beta_hydrolase_family_40_full_ssn.xgmml",
                 "20935_aminotransferase_40_full_ssn.xgmml",
                 "20944__beta-ketoacyl_synthase__c_55_full_ssn.xgmml",
                 "20941_beta-ketoacyl_synthase__n_100_full_ssn.xgmml",
                 "20942_cytochrome_p450_40_full_ssn.xgmml",
                 "20943_degt_dnrj_eryc1_strs_aminotransferase_60_full_ssn.xgmml",
                 "20946_enoyl-_acyl_carrier_protein_45_full_ssn.xgmml",
                 "20947_erythronolide_synthase_docking_9_full_ssn.xgmml",
                 "20948_fad_binding_60_full_ssn.xgmml",
                 "20949_glycosyl_transferase_family_2_20_full_ssn.xgmml",
                 "20953_glycosyltransferase_family_28_25_full_ssn.xgmml",
                 "20951_glycosyl_transferases_group_1_35_full_ssn.xgmml",
                 "20954_glycosyltransferase_like_family_2_25_full_ssn.xgmml",
                 "20955__glyoxalase_bleomycin_resistance_35_full_ssn.xgmml",
                 "20957_kr_domain_55_full_ssn.xgmml",
                 "20960_lanthionine_synthetase_c-like_30_full_ssn.xgmml",
                 "20963_major_facilitator_superfamily_45_full_ssn.xgmml",
                 "20977_methyltransferase_small_domain_20_full_ssn.xgmml",
                 "20976_methyltransferase_domain_45_full_ssn.xgmml",
                 "20979_nad_dependent_epimerase_45_full_ssn.xgmml",
                 "20982_ndp-hexose_2_3-dehydratase_70_full_ssn.xgmml",
                 "20986_o-methyltransferase_40_full_ssn.xgmml",
                 "20987_oxidoreductase_family__c_25_full_ssn.xgmml",
                 "20988_oxidoreductase_family_30_full_ssn.xgmml",
                 "20989__phosphopantetheine_attachment_25_full_ssn.xgmml",
                 "20991_polyketide_cyclase_50_full_ssn.xgmml",
                 "20992_polyketide_synthase_dehydratase_65_full_ssn.xgmml",
                 "20993_protein_of_unknown_function_20_full_ssn.xgmml",
                 "20998_short_chain_dehydrogenase_50_full_ssn.xgmml",
                 "21000_snoal-like_30_full_ssn.xgmml",
                 "21001_spab_c-terminal_25_full_ssn.xgmml",
                 "21003_sugar__and_other__transporter_20_full_ssn.xgmml",
                 "21005_transcriptional_regulatory_18_full_ssn.xgmml",
                 "21004_thioesterase_domain_40_full_ssn.xgmml",
                 "21006_ubie_coq5_methyltransferase_25_full_ssn.xgmml",
                 "21008_udp-glucoronosyl_25_full_ssn.xgmml",
                 "21009_ycao-like_40_full_ssn.xgmml",
                 "21010_zinc-binding_dehydrogenase_40_full_ssn.xgmml",
                 "23780_pyridine_nucleotide_oxidoreductase_full_ssn.xgmml"]

SSN_map_filenames = ["20800_20799_thiolase_N_term_25_full_ssn_coloredssn_mapping_table",
                     "21332_21330_abc_35_full_ssn_coloredssn_mapping_table",
                     "20878_20870_acyl_transferase_pfam_100_full_ssn_coloredssn_mapping_table",
                     "21013_20923_aaa_30_full_ssn_coloredssn_mapping_table",
                     "21015_20928_abc2_30_full_ssn_coloredssn_mapping_table",
                     "21016_20929_acyl-coa_dehydrogenase_c_30_full_ssn_coloredssn_mapping_table",
                     "21017_20930_acyl-coa_dehydrogenase_n_20_full_ssn_coloredssn_mapping_table",
                     "21018_20931_alcohol_dehydrogenase_groes_20_full_ssn_coloredssn_mapping_table",
                     "21019_20932__alpha_beta_hydrolase_family_40_full_ssn_coloredssn_mapping_table",
                     "21021_20935_aminotransferase_40_full_ssn_coloredssn_mapping_table",
                     "21028_20944__beta-ketoacyl_synthase__c_55_full_ssn_coloredssn_mapping_table",
                     "21025_20941_beta-ketoacyl_synthase__n_100_full_ssn_coloredssn_mapping_table",
                     "21026_20942_cytochrome_p450_40_full_ssn_coloredssn_mapping_table",
                     "21027_20943_degt_dnrj_eryc1_strs_aminotransferase_60_full_ssn_coloredssn_mapping_table",
                     "21029_20946_enoyl-_acyl_carrier_protein_45_full_ssn_coloredssn_mapping_table",
                     "21030_20947_erythronolide_synthase_docking_9_full_ssn_coloredssn_mapping_table",
                     "21031_20948_fad_binding_60_full_ssn_coloredssn_mapping_table",
                     "21032_20949_glycosyl_transferase_family_2_20_full_ssn_coloredssn_mapping_table",
                     "21035_20953_glycosyltransferase_family_28_25_full_ssn_coloredssn_mapping_table",
                     "21033_20951_glycosyl_transferases_group_1_35_full_ssn_coloredssn_mapping_table",
                     "21036_20954_glycosyltransferase_like_family_2_25_full_ssn_coloredssn_mapping_table",
                     "21037_20955__glyoxalase_bleomycin_resistance_35_full_ssn_coloredssn_mapping_table",
                     "21039_20957_kr_domain_55_full_ssn_coloredssn_mapping_table",
                     "21040_20960_lanthionine_synthetase_c-like_30_full_ssn_coloredssn_mapping_table",
                     "21043_20963_major_facilitator_superfamily_45_full_ssn_coloredssn_mapping_table",
                     "21045_20977_methyltransferase_small_domain_20_full_ssn_coloredssn_mapping_table",
                     "21044_20976_methyltransferase_domain_45_full_ssn_coloredssn_mapping_table",
                     "21047_20979_nad_dependent_epimerase_45_full_ssn_coloredssn_mapping_table",
                     "21049_20982_ndp-hexose_2_3-dehydratase_70_full_ssn_coloredssn_mapping_table",
                     "21052_20986_o-methyltransferase_40_full_ssn_coloredssn_mapping_table",
                     "21053_20987_oxidoreductase_family__c_25_full_ssn_coloredssn_mapping_table",
                     "21054_20988_oxidoreductase_family_30_full_ssn_coloredssn_mapping_table",
                     "21055_20989__phosphopantetheine_attachment_25_full_ssn_coloredssn_mapping_table",
                     "21056_20991_polyketide_cyclase_50_full_ssn_coloredssn_mapping_table",
                     "21057_20992_polyketide_synthase_dehydratase_65_full_ssn_coloredssn_mapping_table",
                     "21058_20993_protein_of_unknown_function_20_full_ssn_coloredssn_mapping_table",
                     "21063_20998_short_chain_dehydrogenase_50_full_ssn_coloredssn_mapping_table",
                     "21064_21000_snoal-like_30_full_ssn_coloredssn_mapping_table",
                     "21065_21001_spab_c-terminal_25_full_ssn_coloredssn_mapping_table",
                     "21067_21003_sugar__and_other__transporter_20_full_ssn_coloredssn_mapping_table",
                     "21069_21005_transcriptional_regulatory_18_full_ssn_coloredssn_mapping_table",
                     "21068_21004_thioesterase_domain_40_full_ssn_coloredssn_mapping_table",
                     "21070_21006_ubie_coq5_methyltransferase_25_full_ssn_coloredssn_mapping_table",
                     "21071_21008_udp-glucoronosyl_25_full_ssn_coloredssn_mapping_table",
                     "21072_21009_ycao-like_40_full_ssn_coloredssn_mapping_table",
                     "21073_21010_zinc-binding_dehydrogenase_40_full_ssn_coloredssn_mapping_table",
                     "23815_23780_pyridine_nucleotide_oxidoreductase_full_ssn_coloredssn_mapping_table"]

ssn_thresholds = [25,35,100,30,30,30,20,20,40,40,55,100,40,60,45,9,60,
                  20,25,35,25,35,55,30,45,20,45,45,70,40,25,30,
                  25,50,65,20,50,30,25,20,18,40,25,25,40,40,30]

def findAllPFAMSeqs(genbank_file, pfam_name):
    cluster_domains = []
    score_cutoff = 20

    record = SeqIO.read(open(genbank_file, 'rU'),"genbank")
    features = record.features
    for feature in features:
        if feature.type == "PFAM_domain":
            score = float(feature.qualifiers["score"][0])
            if score <score_cutoff:
                continue
            domain_description = feature.qualifiers["description"][0]
            #domain_description = domain_description.lower().replace("\n", " ")
            
            if pfam_name == domain_description:
                cluster_domains.append(feature.qualifiers["translation"])
    return cluster_domains

def makePFAMNameDictionary(SSN_pfam_names, SSN_feature_names):
    prev_feature_name = ""
    pfam_name_dictionary = {}
    i = 0
    for name in SSN_feature_names:
        if name[0:name.rfind("_")] == prev_feature_name:
            continue
        pfam_name_dictionary[SSN_pfam_names[i]] = name[0:name.rfind("_")]
        prev_feature_name = name[0:name.rfind("_")]
        i+=1
    return pfam_name_dictionary

def generateSSNFeatureMatrix(clusters, SSN_pfam_names, SSN_feature_names, included_clusters, blast_exe,genome_name, data_path):
    pfam_name_dictionary =makePFAMNameDictionary(SSN_pfam_names, SSN_feature_names)
    SSN_feature_matrix = np.zeros((len(clusters), len(SSN_feature_names)))
    j = 0
    for c in clusters:
        print(c)
        filename = c 

        i = 0
        for pfam in SSN_pfam_names:
            pfam_sequences = findAllPFAMSeqs(filename, pfam)
            for sequence in pfam_sequences:
                if len(sequence) == 0:
                    continue
                #print pfam
                #print sequence
                ssn_seq_filename = data_path + "SSN/gene_sequences/" + pfam_name_dictionary[pfam] +"_domains.fasta"
                ssn_membership =findSSNMembership(c, sequence, ssn_seq_filename, i, included_clusters[pfam_name_dictionary[pfam]], blast_exe,genome_name,data_path)
                print(ssn_membership)
                for cluster in ssn_membership:
                    #print pfam_name_dictionary[pfam] + "_" + cluster
                    if pfam_name_dictionary[pfam] + "_" + cluster not in SSN_feature_names:
                        print("no match")
                        continue
                    #print "match"
                    SSN_feature_matrix[j,SSN_feature_names.index(pfam_name_dictionary[pfam] + "_" + cluster)] +=1  
                    #print c 
                    #print cluster
                    #print ""
            i+= 1
        print(SSN_feature_matrix[j,:])
        j += 1
    #print SSN_feature_matrix
    return SSN_feature_matrix

def getGeneName(network_fn):
    tree = ET.parse(network_fn)
    root = tree.getroot()
    node_to_gene_ids = {}
    for child in root:
        if "node" in child.tag:
            node_id = child.attrib['id']
            #print node_id
            gene_id = child[1][0].attrib['value']
            node_to_gene_ids[node_id] = gene_id
    return node_to_gene_ids

def parseBLAST(blast_out):
    return

def findSSNMembership(cluster_name, sequence, ssn_seq_filename, ssn_index, included_clusters, blast_exe, genome_name,data_path):
    #read all the sequences from the sequence file, discard any with the same name as current cluster
    #print cluster_name
    #print ssn_seq_filename
    seq_file = open(ssn_seq_filename, 'r')
    sequences = {}
    for line in seq_file:
        if ">" in line:
            gene_name = line.replace("\n","").replace(">","").replace("\r","")
        else:
            sequences[gene_name] = line.replace("\n","")
    #group sequences by SSN cluster
    sequences_in_SSN = {}
    gene_ids_in_SSN = {}
    SSN_filename = SSN_filenames[ssn_index]
    node_to_gene = getGeneName(data_path+"SSN/" +SSN_filename)
    SSN_map_file =open(data_path+"SSN/" +SSN_map_filenames[ssn_index] + ".txt")
    for line in SSN_map_file:
        #print line
        if "UniProt" in line:
            continue
        node_id = line.split()[0]
        cluster_id = line.split()[1]
        if cluster_id not in included_clusters:
            continue
        gene_id = node_to_gene[node_id]
        if cluster_name in gene_id:
            continue
        if cluster_id not in sequences_in_SSN:
            sequences_in_SSN[cluster_id] = []
            gene_ids_in_SSN[cluster_id] = []
        sequences_in_SSN[cluster_id].append(sequences[gene_id])
        gene_ids_in_SSN[cluster_id].append(gene_id)
    SSN_map_file.close()
    
    #figure out which clusters sequences with e-value score lower than threshold
    clusters_with_match = []
    threshold = ssn_thresholds[ssn_index]
    waiting = True
    while waiting:
        try:
            temp_fasta = open(data_path+"temp_file/"+genome_name+"seq1.fasta", 'w')
            waiting = False
        except Exception as e:
            print(e)
            print("file error!!")
            
            
    temp_fasta.write(">"+cluster_name+"\n")
    temp_fasta.write(sequence[0])
    temp_fasta.close()
    for cluster in sequences_in_SSN:
        #print cluster
        for seq in sequences_in_SSN[cluster]:
            waiting = True
            while waiting:
                try:
                    temp_fasta = open(data_path+"temp_file/"+genome_name+"seq2.fasta", 'w')
                    waiting = False
                except:
                    print("file error!!")
                    
                                    
            temp_fasta.write(">"+cluster +"\n")
            temp_fasta.write(seq)
            temp_fasta.close()
            gap_open = '11' 
            gap_extend = '1'
            comp_based = 2
            cline = '{0} -query {1}temp_file/{2}seq1.fasta -subject {1}temp_file/{2}seq2.fasta -gapopen {3} -gapextend {4} -comp_based_stats {5} -use_sw_tback -outfmt \"6\" -max_hsps 1 -evalue {6}'.format(blast_exe,data_path,genome_name,gap_open,gap_extend,comp_based,5)
            #cline = 'blastp -query /home/asw23/antismash_on_full_genome/temp_file/seq1.fasta'
            child = subprocess.Popen(cline,stdout=subprocess.PIPE, shell=True)
            my_out,my_err = child.communicate()
            score = 0
            #print sequence[0]
            #print ""
            #print seq
            #print my_out
            if len(my_out.split()) >=12:
                score = float(my_out.split()[11])
                score = -1*math.log10((2**(-1*score))*len(sequence[0])*len(seq))
            
            if score > threshold:
                    #print threshold
                    #print "score: " + str(score)
                    clusters_with_match.append(cluster)
                    break
            
            cline = '{0} -query {1}temp_file/{2}seq2.fasta -subject {1}temp_file/{2}seq1.fasta -gapopen {3} -gapextend {4} -comp_based_stats {5} -use_sw_tback -outfmt \"6\" -max_hsps 1 -evalue {6}'.format(blast_exe,data_path,genome_name,gap_open,gap_extend,comp_based,5)
            #child = subprocess.Popen(cline,shell=True,stdout=subprocess.PIPE,close_fds=True)
            child = subprocess.Popen(cline,stdout=subprocess.PIPE,shell=True)
            my_out,my_err = child.communicate()
            if len(my_out.split()) >=12:
                score = float(my_out.split()[11])
                score = -1*math.log10((2**(-1*score))*len(sequence[0])*len(seq))
           
            if score > threshold:
                    clusters_with_match.append(cluster)
                    break
    #(cluster_counts, cluster_list) = readSSNFile(name)
    #print clusters_with_match
    return clusters_with_match