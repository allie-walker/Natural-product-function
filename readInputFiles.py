# -*- coding: utf-8 -*-
"""
Created on Fri Jan  8 10:45:34 2021

@author: allis
"""
import cluster_function_prediction_tools as tools
import numpy as np
import readFeatureFiles
import warnings

def checkForHits(feature_vector, pfam_indices, resistance_indices):
    pfam_vector = feature_vector[0,pfam_indices[0]:pfam_indices[1]]
    resistance_vector = feature_vector[0,resistance_indices[0]:resistance_indices[1]]
    if np.sum(pfam_vector) == 0:
        warnings.warn("no pfam features found, make sure to run antiSMASH with --fullhmmer option")
        warnings.warn("if antiSMASH was run with fullhmmer then this BGC does not have enough similarity to the training set for useful predictions")
    if np.sum(resistance_vector) == 0 and resistance_indices[1] != 0:
        warnings.warn("no resistance features found, double check RGI input file")
   

def processSecMetFeature(feature):
    subtype = ""
    pks_signature = ""
    minowa = ""
    consensus = ""
    stachelhaus = ""
    nrpspredictor = ""
    pHMM = ""
    predicat = ""
    sandpuma = ""
    for f in feature.qualifiers["sec_met"]:
        if "NRPS/PKS subtype" in f:
            subtype = f.split(":")[1]
        if "Substrate specificity predictions" in f:
            if "PKS_AT" in f:
                predictions = f.split(":")[4]
                pks_signature = predictions.split(",")[0].split()[0]
                minowa = predictions.split(",")[1].split()[0]
                consensus = predictions.split(",")[2].split()[0]
            if "AMP-binding" in f:
                predictions = f.split(":")[5]
                stachelhaus = predictions.split(",")[0].split()[0]
                nrpspredictor = predictions.split(",")[1].split()[0]
                pHMM = predictions.split(",")[2].split()[0]
                predicat = predictions.split(",")[3].split()[0]
                sandpuma = predictions.split(",")[4].split()[0]
    return (subtype, pks_signature, minowa, consensus, stachelhaus, nrpspredictor, pHMM, predicat, sandpuma)

def readAntismash4(as_features):
    score_cutoff = 20
    CDS_motif_list = []
    smCOG_list = []
    pfam_list = []
    
    CDS_motifs = {}
    smCOGs = {}
    pfam_counts = {}
    pks_nrp_subtypes = {}
    pk_monomers_signature = {}
    pk_monomers_minowa = {}
    pk_monomers_consensus = {}
    nrp_monomers_stachelhaus = {}
    nrp_monomers_nrps_predictor= {}
    nrp_monomers_pHMM = {}
    nrp_monomers_predicat = {}
    nrp_monomers_sandpuma= {}
    
    all_results = {}
    
    for feature in as_features:
        subtype = ""
        pks_signature = ""
        minowa = ""
        consensus = ""
        stachelhaus = ""
        nrpspredictor = ""
        pHMM = ""
        predicat = ""
        sandpuma = ""
        if feature.type == "CDS" and "sec_met" in feature.qualifiers:
           (subtype, pks_signature, minowa, consensus, stachelhaus, nrpspredictor, pHMM, predicat, sandpuma)= processSecMetFeature(feature)
        if subtype != "":
            if subtype not in pks_nrp_subtypes:
                pks_nrp_subtypes[subtype] = 0
            pks_nrp_subtypes[subtype] += 1
        if pks_signature != "" and "no_call" not in pks_signature and "N/A" not in pks_signature and "n/a" not in pks_signature:
            if pks_signature not in pk_monomers_signature:
                pk_monomers_signature[pks_signature] = 0
            pk_monomers_signature[pks_signature] += 1
        if minowa != "" and "no_call" not in minowa and "N/A" not in minowa and "n/a" not in minowa:
            if minowa not in pk_monomers_minowa:
                pk_monomers_minowa[minowa] = 0
            pk_monomers_minowa[minowa] += 1
        if consensus != "" and "no_call" not in consensus and "N/A" not in consensus and "n/a" not in consensus:
            if consensus not in pk_monomers_consensus:
                pk_monomers_consensus[consensus] = 0
            pk_monomers_consensus[consensus] += 1
        if stachelhaus != "" and "no_call" not in stachelhaus and "N/A" not in stachelhaus and "n/a" not in stachelhaus:
            if stachelhaus not in nrp_monomers_stachelhaus:
                nrp_monomers_stachelhaus[stachelhaus] = 0
            nrp_monomers_stachelhaus[stachelhaus] += 1
        if nrpspredictor != "" and "no_call" not in nrpspredictor and "N/A" not in nrpspredictor and "n/a" not in nrpspredictor:
            if nrpspredictor not in nrp_monomers_nrps_predictor:
                nrp_monomers_nrps_predictor[nrpspredictor] = 0
            nrp_monomers_nrps_predictor[nrpspredictor] += 1
        if pHMM != "" and "no_call" not in pHMM and "N/A" not in pHMM and "n/a" not in pHMM:
            if pHMM not in nrp_monomers_pHMM:
                nrp_monomers_pHMM[pHMM] = 0
            nrp_monomers_pHMM[pHMM] += 1
        if predicat != "" and "no_call" not in predicat and "N/A" not in predicat and "n/a" not in predicat:
            if predicat not in nrp_monomers_predicat:
                nrp_monomers_predicat[predicat] = 0
            nrp_monomers_predicat[predicat] += 1
        if sandpuma != "" and "no_call" not in sandpuma and "N/A" not in sandpuma and "n/a" not in sandpuma:
            if sandpuma not in nrp_monomers_sandpuma:
                nrp_monomers_sandpuma[sandpuma] = 0
            nrp_monomers_sandpuma[sandpuma] += 1
            
        if feature.type == "CDS_motif":
            note_text = feature.qualifiers['note'][0]
            if "(" not in note_text:
                continue
            motif_name = note_text[0:note_text.index("(")-1]
            if motif_name not in CDS_motif_list:
                CDS_motif_list.append(motif_name)
            if motif_name not in CDS_motifs:
                CDS_motifs[motif_name] = 0
            CDS_motifs[motif_name] += 1
        elif feature.type == "CDS":
            if "note" in feature.qualifiers:
                for note in feature.qualifiers["note"]:
                    if "smCOG" in note:
                        if ":" not in note or "(" not in note:
                            continue
                        smCOG_type = note[note.index(":")+2:note.index("(")-1]
                        if smCOG_type not in smCOG_list:
                            smCOG_list.append(smCOG_type)
                        if smCOG_type not in smCOGs:
                            smCOGs[smCOG_type] = 0
                        smCOGs[smCOG_type] += 1
        elif feature.type == "PFAM_domain":
            score = float(feature.qualifiers["score"][0])
            if score <score_cutoff:
                continue
            domain_description = feature.qualifiers["description"][0]
            if domain_description not in pfam_list:
                pfam_list.append(domain_description)
            if domain_description not in pfam_counts:
                pfam_counts[domain_description] = 0
            pfam_counts[domain_description] += 1
    
    all_results = {"CDS_motifs": CDS_motifs, "pk_consensus": pk_monomers_consensus, "nrp_stachelhaus": nrp_monomers_stachelhaus, \
                   "PFAM": pfam_counts, "nrp_pHMM": nrp_monomers_pHMM, "pks_nrps_type": pks_nrp_subtypes, "pk_minowa": pk_monomers_minowa, \
                   "pk_signature": pk_monomers_signature, "nrp_predicat": nrp_monomers_predicat, "nrp_nrpspredictor": nrp_monomers_nrps_predictor, \
                   "SMCOG": smCOGs, "nrp_sandpuma": nrp_monomers_sandpuma}
    return all_results

def readAntismash5(as_features):
    score_cutoff = 20
    CDS_motifs = {}
    smCOGs = {}
    pfam_counts = {}
    pk_monomers_consensus = {}
    for feature in as_features:
        consensus = ""
        if feature.type == "aSDomain":
            if "specificity" not in feature.qualifiers:
                continue
            for f in feature.qualifiers["specificity"]:
                if "consensus" in f:
                    consensus = f.split(":")[1].replace('"','')
        if consensus != "" and "no_call" not in consensus and "N/A" not in consensus and "n/a" not in consensus:
            if consensus not in pk_monomers_consensus:
                pk_monomers_consensus[consensus] = 0
            pk_monomers_consensus[consensus] += 1
            
        if feature.type == "CDS_motif":
            if 'label' not in feature.qualifiers: #this happens for ripp sequences
                continue
            motif_name = feature.qualifiers['label'][0]
            if motif_name not in CDS_motifs:
                CDS_motifs[motif_name] =0
            CDS_motifs[motif_name] += 1
        elif feature.type == "CDS":
            if "gene_functions" in feature.qualifiers:
                for note in feature.qualifiers["gene_functions"]:
                    if "SMCOG" in note:
                        if ":" not in note or "(" not in note:
                            continue
                        smCOG_type = note[note.index(":")+2:note.rfind("(")-1]
                        if smCOG_type not in smCOGs:
                            smCOGs[smCOG_type] = 0
                        smCOGs[smCOG_type] += 1
        elif feature.type == "PFAM_domain":
            score = float(feature.qualifiers["score"][0])
            if score <score_cutoff:
                continue
            domain_description = feature.qualifiers["description"][0]
            pfam_id = feature.qualifiers["db_xref"][0]
            pfam_id = pfam_id[pfam_id.find(" ")+1:len(pfam_id)]
            if domain_description not in pfam_counts:
                pfam_counts[domain_description] = 0
            pfam_counts[domain_description] += 1                
            
    return (pfam_counts, CDS_motifs, smCOGs, pk_monomers_consensus)



def readRGIFile3(rgi_infile):
    e_value_threshold = 0.1
    resistance_genes = {}
    resistance_genes_list = []
    for line in rgi_infile:
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
        if best_hit not in resistance_genes:
            resistance_genes[best_hit] = 0
        resistance_genes[best_hit] += 1
    rgi_infile.close()        
    return resistance_genes

def readRGIFile5(rgi_infile):
    bit_score_threshold = 40
    e_value_threshold = 0.1
    resistance_genes = {}
    use_bit_score = True
    for line in rgi_infile:
        if "ORF_ID" in line:
            if line.split("\t")[7] == "Best_Hit_evalue":
                warnings.warn('training set was generated using bit scores but RGI output has e-values, will use e-value threshold but this could increase error in predictions')
                use_bit_score = False
            continue
        entries = line.split("\t")
        bit_score = float(entries[7])
        if use_bit_score:
            if bit_score < bit_score_threshold:
                continue
        elif bit_score > e_value_threshold:
            continue
        best_hit = entries[8]
        if best_hit not in resistance_genes:
            resistance_genes[best_hit] = 0
        resistance_genes[best_hit] += 1
    rgi_infile.close()
    return resistance_genes

def readInputFiles(as_features, as_version, rgi_infile, rgi_version, training_features, training_features_types, feature_list_by_type, data_path):
    #TODO: need to read sandpuma features for antismash6
    #TODO: look into how to handle if there are multiple feature names reapeated in different sources?

    CDS_motifs = {}
    smCOGs = {}
    pfam_counts = {}
    pks_nrp_subtypes = {}
    pk_monomers_signature = {}
    pk_monomers_minowa = {}
    pk_monomers_consensus = {}
    nrp_monomers_stachelhaus = {}
    nrp_monomers_nrps_predictor= {}
    nrp_monomers_pHMM = {}
    nrp_monomers_predicat = {}
    nrp_monomers_sandpuma= {}
    try:
            
        """if as_version == 4:
            print("here1")
            print(data_path+"features/pks_nrps_type_list.txt")
            used_pks_nrps_type_list = readFeatureFiles.readFeatureList(data_path+"/features/pks_nrps_type_list.txt")
            used_pk_signature_list = readFeatureFiles.readFeatureList(data_path+"/features/pk_signature_list.txt")
            used_pk_minowa_list = readFeatureFiles.readFeatureList(data_path+"/features/pk_minowa_list.txt")
            used_pk_consensus_list = readFeatureFiles.readFeatureList(data_path+"/features/pk_consensus_list.txt")
            
            used_nrp_stachelhaus_list = readFeatureFiles.readFeatureList(data_path+"/features/nrp_stachelhaus_list.txt")
            used_nrp_nrps_predictor_list = readFeatureFiles.readFeatureList(data_path+"/features/nrp_nrpspredictor_list.txt")
            used_nrp_pHMM_list = readFeatureFiles.readFeatureList(data_path+"/features/nrp_pHMM_list.txt")
            used_nrp_predicat_list = readFeatureFiles.readFeatureList(data_path+"/features/nrp_predicat_list.txt")
            used_nrp_sandpuma_list = readFeatureFiles.readFeatureList(data_path+"/features/nrp_sandpuma_list.txt")
        
        
            used_pfam_list = readFeatureFiles.readFeatureList(data_path+"/features/PFAM_list.txt")
            used_CDS_list = readFeatureFiles.readFeatureList(data_path+"/features/CDS_motifs_list.txt")
            used_smCOG_list = readFeatureFiles.readFeatureList(data_path+"/features/SMCOG_list.txt")
            print("here2")
        if as_version == 5:
            used_pfam_list = readFeatureFiles.readFeatureList(data_path+"/features/pfam_list5.txt")
            used_smCOG_list = readFeatureFiles.readFeatureList(data_path+"/features/smCOG_list5.txt")
            used_CDS_list = readFeatureFiles.readFeatureList(data_path+"/features/CDS_motifs_list5.txt")    
            used_pk_consensus_list = readFeatureFiles.readFeatureList(data_path+"/features/pk_nrp_consensus_list5.txt")
        #TODO: as version 6
        if rgi_version == 3:
            used_resistance_genes_list = readFeatureFiles.readFeatureList(data_path+"/features/CARD_gene_list.txt")
        else:
            used_resistance_genes_list = readFeatureFiles.readFeatureList(data_path+"/features/CARD5_gene_list.txt")"""
    except:
        print("did not find file containing training data, please keep script located in directory downloaded from github")
        exit()

    if as_version == 4:
        as_results = readAntismash4(as_features)
    else:
        (pfam_counts, CDS_motifs, smCOGs, pk_monomers_consensus) = readAntismash5(as_features)
    if rgi_version == 3:
        resistance_genes = readRGIFile3(rgi_infile)
    elif rgi_version == 5:
         resistance_genes = readRGIFile5(rgi_infile)
         
    #TODO: make sure features get added in correct order!
    
    
    test_features = np.zeros((1, training_features.shape[1]))
    i = 0
    card_start_index = 0
    card_end_index =0
    pfam_start_index = 0
    pfam_end_index =0
    for t in training_features_types:
        print(t)
        if t in as_results:
            if t == "PFAM":
                pfam_start_index = i
            (test_features, i) = tools.addToFeatureMatrix(test_features, i, as_results[t], feature_list_by_type[t])
            if t == "PFAM":
                pfam_end_index = i
        elif "CARD" in t:
            print(resistance_genes)
            card_start_index = i
            (test_features, i) = tools.addToFeatureMatrix(test_features, i, resistance_genes, feature_list_by_type[t])
            card_end_index = i
        else:
            print(t)
            raise Exception("Unknown feature type")

    checkForHits(test_features, (pfam_start_index, pfam_end_index), (card_start_index, card_end_index))
    return test_features