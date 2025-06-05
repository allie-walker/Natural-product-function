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

def readAntismash6(as_features):
    #works for versions 6-8
    score_cutoff = 20 #PFAM HMM score to be counted
    CDS_motifs = {}
    smCOGs = {}
    pfam_counts = {}
    tigrfam_counts = {}
    NRPS_PKS= {}
    NRPS_PKS_substrate = {}

    for feature in as_features:
        consensus = ""
        if feature.type == "aSDomain":
            if feature.qualifiers["aSTool"][0] == "tigrfam":
                score = float(feature.qualifiers["score"][0])
                if score < score_cutoff:
                    continue
                domain_description = feature.qualifiers["description"][0]
                tigrfam_id = feature.qualifiers["identifier"][0]
                tigrfam_id = tigrfam_id[tigrfam_id.find(" ")+1:len(tigrfam_id)]
                if domain_description not in tigrfam_counts:
                    tigrfam_counts[domain_description] = 0
                tigrfam_counts[domain_description] += 1
            else:
                if 'label' not in feature.qualifiers: 
                    continue
                motif_name = feature.qualifiers['label'][0]
                if motif_name not in NRPS_PKS:
                    NRPS_PKS[motif_name] = 0
                NRPS_PKS[motif_name] += 1
                if "specificity" not in feature.qualifiers:
                    continue
                for substrate_pred in feature.qualifiers["specificity"]:
                    if "consensus" not in substrate_pred:
                        continue
                    if substrate_pred not in NRPS_PKS_substrate:
                        NRPS_PKS_substrate[substrate_pred] = 0
                    NRPS_PKS_substrate[substrate_pred] += 1
            
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
    return (pfam_counts, CDS_motifs, smCOGs, NRPS_PKS,tigrfam_counts, NRPS_PKS_substrate)


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

    if as_version == 4:
        as_results = readAntismash4(as_features)
    #TODO: add antismash 6, 7, 8
    elif as_version==5:
        (pfam_counts, CDS_motifs, smCOGs, pk_monomers_consensus) = readAntismash5(as_features)
        as_results = {}
        as_results["PFAM5"] = pfam_counts
        as_results["CDS_motifs"] = CDS_motifs
        as_results["SMCOG"] = smCOGs
        as_results["pk_nrp_consensus"] = pk_monomers_consensus
    else:
        (pfam_counts, CDS_motifs, smCOGs, NRPS_PKS,tigrfam_counts, NRPS_PKS_substrate) = readAntismash6(as_features)
        as_results = {}
        as_results["PFAM"] = pfam_counts
        as_results["CDS_motifs"] = CDS_motifs
        as_results["SMCOG"] = smCOGs
        as_results["NRPS_PKS"] = NRPS_PKS
        as_results["TIGR_FAM"] = tigrfam_counts
        as_results["NRPS_PKS_substrate"] = NRPS_PKS_substrate
    if rgi_version == 3:
        resistance_genes = readRGIFile3(rgi_infile)
    elif rgi_version == 5 or rgi_version==6:
         resistance_genes = readRGIFile5(rgi_infile)
    
    
    #TODO: make sure features get added in correct order!
    test_features = np.zeros((1, training_features.shape[1]))
    i = 0
    card_start_index = 0
    card_end_index =0
    pfam_start_index = 0
    pfam_end_index =0
    for t in training_features_types:
        if t in as_results:
            if t == "PFAM" or t == "PFAM5":
                pfam_start_index = i
            (test_features, i) = tools.addToFeatureMatrix(test_features, i, as_results[t], feature_list_by_type[t])
            if t == "PFAM" or t == "PFAM5":
                pfam_end_index = i
        elif "CARD" in t or "CARD5" in t:
            card_start_index = i
            (test_features, i) = tools.addToFeatureMatrix(test_features, i, resistance_genes, feature_list_by_type[t])
            card_end_index = i
        else:
            raise Exception("Unknown feature type " + t)
    
        

    checkForHits(test_features, (pfam_start_index, pfam_end_index), (card_start_index, card_end_index))
    return test_features

def readInputFilesLocs(as_features, as_version, rgi_infile, rgi_version, training_features, training_features_types, feature_list_by_type, data_path):
    features = {}
    if as_version == 5:
        (gene_data, pfam_data, cds_motif_data, smcog_data, pks_nrps_data) = readAntismash5Loc(as_features)
        features["PFAM5"] = pfam_data
        features["CDS_motifs"] = cds_motif_data
        features["SMCOG"] = smcog_data
        features["pk_nrp_consensus"] = pks_nrps_data
    else:
        (gene_data, pfam_data, cds_motif_data, smcog_data, pks_nrps_data, pks_nrps_substrate_data, tigrfam_data) = readAntismash6Loc(as_features)
        features["PFAM"] = pfam_data
        features["CDS_motifs"] = cds_motif_data
        features["SMCOG"] = smcog_data
        features["NRPS_PKS"] = pks_nrps_data
        features["NRPS_PKS_substrate"] = pks_nrps_substrate_data
        features["TIGR_FAM"] = tigrfam_data
    if rgi_version != None:
        rgi_data = readRGIFile5Loc(rgi_infile)
        features["CARD_genes"] = rgi_data
    return (gene_data, features)

def readAntismash5Loc(as_features):
    score_cutoff = 20
    CDS_motifs = {}
    smCOGs = {}
    pfam_counts = {}
    pk_monomers_consensus = {}
    pfam_data = []
    cds_motif_data = []
    smcog_data = []
    pks_nrps_data = []
    gene_data = []
    pfam_num = 1
    cds_num = 1
    smcog_num = 1
    nrps_num = 1
    gene_num = 1
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
            nrps_pks_dict = {}
            nrps_pks_dict['nrps+_pks_num'] = nrps_num 
            nrps_pks_dict['start'] = feature.location.start
            nrps_pks_dict['end'] = feature.location.end
            nrps_pks_dict['strand'] = feature.location.strand
            nrps_pks_dict["description"] = consensus
            pks_nrps_data.append(nrps_pks_dict)
            nrps_num +=1
        if feature.type == "CDS_motif":
            if 'label' not in feature.qualifiers: #this happens for ripp sequences
                continue
            motif_name = feature.qualifiers['label'][0]
            if motif_name not in CDS_motifs:
                CDS_motifs[motif_name] =0
            CDS_motifs[motif_name] += 1
            cds_motif_dict = {}
            cds_motif_dict['cds_num'] = cds_num 
            cds_motif_dict['start'] = feature.location.start
            cds_motif_dict['end'] = feature.location.end
            cds_motif_dict['strand'] = feature.location.strand
            cds_motif_dict["description"] = motif_name
            cds_motif_data.append(cds_motif_dict)
            cds_num += 1
        elif feature.type == "CDS":
            gene_dict = {}
            gene_dict["gene_num"] = gene_num
            gene_dict["start"] = feature.location.start
            gene_dict['end'] = feature.location.end
            gene_dict['strand'] = feature.location.strand
            gene_dict['product'] = feature.qualifiers["product"]
            if 'locus_tag' in feature.qualifiers:
                gene_dict['locus_tag'] = feature.qualifiers["locus_tag"]
            else:
                gene_dict['locus_tag'] = ""
            gene_data.append(gene_dict)
            gene_num += 1
            if "gene_functions" in feature.qualifiers:
                for note in feature.qualifiers["gene_functions"]:
                    if "SMCOG" in note:
                        if ":" not in note or "(" not in note:
                            continue
                        smCOG_type = note[note.index(":")+2:note.rfind("(")-1]
                        smCOG_description = note[note.index(":")+1:note.rfind("(")-1]
                        if smCOG_type not in smCOGs:
                            smCOGs[smCOG_type] = 0
                        smCOGs[smCOG_type] += 1
                        smcog_dict = {}
                        smcog_dict['smcog_num'] = smcog_num 
                        smcog_dict['start'] = feature.location.start
                        smcog_dict['end'] = feature.location.end
                        smcog_dict['strand'] = feature.location.strand
                        smcog_dict["description"] = smCOG_description
                        smcog_data.append(smcog_dict)
                        smcog_num = +1
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
            pfam_dict = {}
            pfam_dict['pfam_num'] = pfam_num 
            pfam_dict['start'] = feature.location.start
            pfam_dict['end'] = feature.location.end
            pfam_dict['strand'] = feature.location.strand
            pfam_dict["description"] = domain_description
            pfam_dict['pfam_id'] = pfam_id
            pfam_data.append(pfam_dict)
            pfam_num += 1                          
    return (gene_data, pfam_data, cds_motif_data, smcog_data, pks_nrps_data)

def readAntismash6Loc(as_features):
    #works for versions 6-8
    score_cutoff = 20 #PFAM HMM score to be counted
    CDS_motifs = {}
    smCOGs = {}
    pfam_counts = {}
    tigrfam_counts = {}
    NRPS_PKS= {}
    NRPS_PKS_substrate = {}
    
    pfam_data = []
    cds_motif_data = []
    smcog_data = []
    pks_nrps_data = []
    pks_nrps_substrate_data = []
    tigrfam_data = []
    gene_data = []
    pfam_num = 1
    cds_num = 1
    smcog_num = 1
    nrps_num = 1
    nrps_substrate_num = 1
    tigrfam_num = 1
    gene_num = 1
    for feature in as_features:
        consensus = ""
        if feature.type == "aSDomain":
            if feature.qualifiers["aSTool"][0] == "tigrfam":
                score = float(feature.qualifiers["score"][0])
                if score < score_cutoff:
                    continue
                domain_description = feature.qualifiers["description"][0]
                tigrfam_id = feature.qualifiers["identifier"][0]
                tigrfam_id = tigrfam_id[tigrfam_id.find(" ")+1:len(tigrfam_id)]
                if domain_description not in tigrfam_counts:
                    tigrfam_counts[domain_description] = 0
                tigrfam_counts[domain_description] += 1
                tigrfam_dict = {}
                tigrfam_dict['tigrfam_num'] = tigrfam_num 
                tigrfam_dict['start'] = feature.location.start
                tigrfam_dict['end'] = feature.location.end
                tigrfam_dict['strand'] = feature.location.strand
                tigrfam_dict["description"] = tigrfam_id
                tigrfam_data.append(tigrfam_dict)
                tigrfam_num +=1
            else:
                if 'label' not in feature.qualifiers: 
                    continue
                motif_name = feature.qualifiers['label'][0]
                if motif_name not in NRPS_PKS:
                    NRPS_PKS[motif_name] = 0
                NRPS_PKS[motif_name] += 1
                nrps_pks_dict = {}
                nrps_pks_dict['nrps+_pks_num'] = nrps_num 
                nrps_pks_dict['start'] = feature.location.start
                nrps_pks_dict['end'] = feature.location.end
                nrps_pks_dict['strand'] = feature.location.strand
                nrps_pks_dict["description"] = motif_name
                pks_nrps_data.append(nrps_pks_dict)
                nrps_num +=1
                if "specificity" not in feature.qualifiers:
                    continue
                for substrate_pred in feature.qualifiers["specificity"]:
                    if "consensus" not in substrate_pred:
                        continue
                    if substrate_pred not in NRPS_PKS_substrate:
                        NRPS_PKS_substrate[substrate_pred] = 0
                    NRPS_PKS_substrate[substrate_pred] += 1
                    nrps_pks_substrate_dict = {}
                    nrps_pks_substrate_dict['nrps+_pks_num'] = nrps_substrate_num 
                    nrps_pks_substrate_dict['start'] = feature.location.start
                    nrps_pks_substrate_dict['end'] = feature.location.end
                    nrps_pks_substrate_dict['strand'] = feature.location.strand
                    nrps_pks_substrate_dict["description"] = substrate_pred
                    pks_nrps_substrate_data.append(nrps_pks_substrate_dict)
                    nrps_substrate_num +=1
        if feature.type == "CDS_motif":
            if 'label' not in feature.qualifiers: #this happens for ripp sequences
                continue
            motif_name = feature.qualifiers['label'][0]
            if motif_name not in CDS_motifs:
                CDS_motifs[motif_name] =0
            CDS_motifs[motif_name] += 1
            cds_motif_dict = {}
            cds_motif_dict['cds_num'] = cds_num 
            cds_motif_dict['start'] = feature.location.start
            cds_motif_dict['end'] = feature.location.end
            cds_motif_dict['strand'] = feature.location.strand
            cds_motif_dict["description"] = motif_name
            cds_motif_data.append(cds_motif_dict)
            cds_num += 1
        elif feature.type == "CDS":
            gene_dict = {}
            gene_dict["gene_num"] = gene_num
            gene_dict["start"] = feature.location.start
            gene_dict['end'] = feature.location.end
            gene_dict['strand'] = feature.location.strand
            gene_dict['product'] = feature.qualifiers["product"]
            if 'locus_tag' in feature.qualifiers:
                gene_dict['locus_tag'] = feature.qualifiers["locus_tag"]
            else:
                gene_dict['locus_tag'] = ""
            gene_data.append(gene_dict)
            gene_num += 1
            if "gene_functions" in feature.qualifiers:
                for note in feature.qualifiers["gene_functions"]:
                    if "SMCOG" in note:
                        if ":" not in note or "(" not in note:
                            continue
                        smCOG_type = note[note.index(":")+2:note.rfind("(")-1]
                        if smCOG_type not in smCOGs:
                            smCOGs[smCOG_type] = 0
                        smCOGs[smCOG_type] += 1
                        smCOGs[smCOG_type] += 1
                        smcog_dict = {}
                        smcog_dict['smcog_num'] = smcog_num 
                        smcog_dict['start'] = feature.location.start
                        smcog_dict['end'] = feature.location.end
                        smcog_dict['strand'] = feature.location.strand
                        smcog_dict["description"] = smCOG_type
                        smcog_data.append(smcog_dict)
                        smcog_num = +1
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
            pfam_dict = {}
            pfam_dict['pfam_num'] = pfam_num 
            pfam_dict['start'] = feature.location.start
            pfam_dict['end'] = feature.location.end
            pfam_dict['strand'] = feature.location.strand
            pfam_dict["description"] = domain_description
            pfam_dict['pfam_id'] = pfam_id
            pfam_data.append(pfam_dict)
            pfam_num += 1                            
    return (gene_data, pfam_data, cds_motif_data, smcog_data, pks_nrps_data, pks_nrps_substrate_data, tigrfam_data)


def readRGIFile5Loc(rgi_infile):
    rgi_data = []
    rgi_num = 1
    bit_score_threshold = 40
    e_value_threshold = 0.1
    resistance_genes = {}
    use_bit_score = True
    for line in rgi_infile:
        if "ORF_ID" in line:
            if line.split("\t")[7] == "Best_Hit_evalue":
                #warnings.warn('training set was generated using bit scores but RGI output has e-values, will use e-value threshold but this could increase error in predictions')
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
        
        gene_start = float(entries[2])
        gene_end = float(entries[3])
        gene_orientation = entries[4]
        if gene_orientation == "-":
          gene_orientation = -1
        else:
          gene_orientation = 1
        rgi_dict = {}
        rgi_dict['rgi_num'] = rgi_num 
        rgi_dict['start'] = gene_start
        rgi_dict['end'] = gene_end
        rgi_dict['strand'] = gene_orientation
        rgi_dict["description"] = best_hit
        rgi_data.append(rgi_dict)
        rgi_num = +1
        #print(best_hit)
    rgi_infile.close()
    return rgi_data