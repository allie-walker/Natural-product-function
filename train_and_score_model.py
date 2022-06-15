# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 14:10:16 2022

@author: allison
Trains model parameters based on training set
Training set should be represented by a set of csv files with features and csv files with data labels
"""

import numpy as np
import random
import readFeatureFiles
from sklearn.linear_model import SGDClassifier
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.svm import SVC
from sklearn import preprocessing
import matplotlib.pyplot as plt
from matplotlib import cm

def readFeatureMatrix(filename):
    in_file = open(filename,'r')
    feature_matrix = []
    for line in in_file:
        if "," not in line:
            continue
        entries = line.split(",")
        feature_matrix.append([])
        for i in range(0,len(entries)-1):
            feature_matrix[len(feature_matrix)-1].append(int(entries[i]))
    in_file.close()
    return np.array(feature_matrix)

#reads file with a list of feature names
def readFeatureList(filename):
    in_file = open(filename,'r')
    feature_list = []
    for line in in_file:
        feature_list.append(line.replace("\n",""))
    return feature_list

#reads file with list of classifications for gene clusters
def readClassesMatrix(filename):
    class_matrix = []
    in_file = open(filename,'r')
    for line in in_file:
        if len(line)>0:
            class_matrix.append(int(line))
    in_file.close()
    return np.array(class_matrix)  

#gets indices of all natural products with known activities, based on a list classifying activities as known (0)
#or unknown (1)
def getNotUnknownIndices(is_unknown):
    is_not_unknown_indices = []
    for i in range(0, len(is_unknown)):
        if not is_unknown[i]:
            is_not_unknown_indices.append(i)
    return is_not_unknown_indices

def splitData(x_vars, y_vars, validation_index, num_splits):
    split_y = np.array_split(y_vars, num_splits)
    split_x = np.array_split(x_vars, num_splits)
    training_x = np.array([])
    training_y = np.array([])
    for i in range(0, 10):
        if i == validation_index:
            continue
        if training_x.shape[0] == 0:
            training_x = split_x[i]
            training_y = split_y[i]
        else:
            training_x = np.concatenate((training_x, split_x[i]), axis=0)
            training_y = np.concatenate((training_y, split_y[i]),axis=0)
    return (training_x, training_y, split_x[validation_index], split_y[validation_index])

def makeRandomOrder(seed, y_vars):
    np.random.seed(seed)
    y_vars_range = np.arange(y_vars.shape[0])
    np.random.shuffle(y_vars_range)
    return y_vars_range

#test classifier in ten-fold cross-validation
#returns avg accuracy
def validateClassifier(outfile, classifier, features, y_vars, regression):
    avg_accuracy = 0.0
    for i in range(0, 10):
        (training_x, training_y, val_x, val_y) = splitData(features, y_vars, i, 10)
        #scale features if regression classifier
        if regression:
            min_max_scaler = preprocessing.MinMaxScaler()
            training_x = min_max_scaler.fit_transform(training_x)
            val_x = min_max_scaler.transform(val_x)
        classifier.fit(training_x, training_y)
        score = classifier.score(val_x, val_y)
        avg_accuracy += score
        outfile.write(str(score) + ",")
    avg_accuracy /= 10.0
    return avg_accuracy

#parameters
#classification to train classifiers on 
#options: antibacterial, antieuk (defined as antifungal, antitumor, or cytotoxic), antifungal, cytotoxic_antitumor, antigramneg, antigrampos
classification = "antibacterial" 
#set random seed so results are consistent
random.seed(1)
#include features from SSN 
include_SSN = True

pfam_features = readFeatureMatrix("feature_matrices/PFAM.csv")
card_features = readFeatureMatrix("feature_matrices/CARD_gene.csv")
smCOG_features = readFeatureMatrix("feature_matrices/SMCOG.csv")
SSN_features = readFeatureMatrix("feature_matrices/SSN.csv")
CDS_features = readFeatureMatrix("feature_matrices/CDS_motifs.csv")

pks_nrps_type_features = readFeatureFiles.readFeatureMatrix("feature_matrices/pks_nrps_type.csv")
pk_signature_features = readFeatureFiles.readFeatureMatrix("feature_matrices/pk_signature.csv")
pk_minowa_features = readFeatureFiles.readFeatureMatrix("feature_matrices/pk_minowa.csv")
pk_consensus_features = readFeatureFiles.readFeatureMatrix("feature_matrices/pk_consensus.csv")

nrp_stachelhaus_features = readFeatureFiles.readFeatureMatrix("feature_matrices/nrp_stachelhaus.csv")
nrp_nrpspredictor_features = readFeatureFiles.readFeatureMatrix("feature_matrices/nrp_nrpspredictor.csv")
nrp_pHMM_features = readFeatureFiles.readFeatureMatrix("feature_matrices/nrp_pHMM.csv")
nrp_predicat_features = readFeatureFiles.readFeatureMatrix("feature_matrices/nrp_predicat.csv")
nrp_sandpuma_features = readFeatureFiles.readFeatureMatrix("feature_matrices/nrp_sandpuma.csv")

pfam_list = readFeatureFiles.readFeatureList("feature_matrices/PFAM_list.txt")
card_list = readFeatureFiles.readFeatureList("feature_matrices/CARD_gene_list.txt")
smCOG_list = readFeatureFiles.readFeatureList("feature_matrices/SMCOG_list.txt")
SSN_list = readFeatureFiles.readFeatureList("feature_matrices/SSN_list.txt")
CDS_list = readFeatureFiles.readFeatureList("feature_matrices/CDS_motifs_list.txt")

pks_nrps_type_list = readFeatureFiles.readFeatureList("feature_matrices/pks_nrps_type_list.txt")
pk_signature_list = readFeatureFiles.readFeatureList("feature_matrices/pk_signature_list.txt")
pk_minowa_list = readFeatureFiles.readFeatureList("feature_matrices/pk_minowa_list.txt")
pk_consensus_list = readFeatureFiles.readFeatureList("feature_matrices/pk_consensus_list.txt")

nrp_stachelhaus_list = readFeatureFiles.readFeatureList("feature_matrices/nrp_stachelhaus_list.txt")
nrp_nrpspredictor_list = readFeatureFiles.readFeatureList("feature_matrices/nrp_nrpspredictor_list.txt")
nrp_pHMM_list = readFeatureFiles.readFeatureList("feature_matrices/nrp_pHMM_list.txt")
nrp_predicat_list = readFeatureFiles.readFeatureList("feature_matrices/nrp_predicat_list.txt")
nrp_sandpuma_list = readFeatureFiles.readFeatureList("feature_matrices/nrp_sandpuma_list.txt")

#read classes
is_antibacterial = readFeatureFiles.readClassesMatrix("feature_matrices/is_antibacterial.csv")
is_antifungal = readFeatureFiles.readClassesMatrix("feature_matrices/is_antifungal.csv")
is_cytotoxic = readFeatureFiles.readClassesMatrix("feature_matrices/is_cytotoxic.csv")
is_unknown = readFeatureFiles.readClassesMatrix("feature_matrices/is_unknown.csv")
targets_gram_pos = readFeatureFiles.readClassesMatrix("feature_matrices/targets_gram_pos.csv")
targets_gram_neg = readFeatureFiles.readClassesMatrix("feature_matrices/targets_gram_neg.csv")
full_cluster_list = readFeatureFiles.readClusterList("feature_matrices/cluster_list_CARD.txt")
is_not_unknown_indices = readFeatureFiles.getNotUnknownIndices(is_unknown)
target_unannotated = is_antibacterial*((targets_gram_pos+targets_gram_neg)<1)
is_not_unknown_indices_gram =  readFeatureFiles.getNotUnknownIndices(is_unknown + target_unannotated)

is_antibacterial = (is_antibacterial >= 1).astype(int)
is_antieuk = ((is_antifungal + is_cytotoxic)>=1).astype(int)
is_gram_pos = (targets_gram_pos >= 1).astype(int)
is_gram_neg = (targets_gram_neg >= 1).astype(int)

#concatenate all features into one matrix
features = np.concatenate((pfam_features, card_features), axis=1)
features = np.concatenate((features,  smCOG_features), axis=1)
features = np.concatenate((features,  CDS_features), axis=1)
if include_SSN:
    features = np.concatenate((features,  SSN_features), axis=1)
features = np.concatenate((features,  pks_nrps_type_features), axis=1)
features = np.concatenate((features,  pk_signature_features), axis=1)
features = np.concatenate((features,  pk_minowa_features), axis=1)
features = np.concatenate((features,  pk_consensus_features), axis=1)
features = np.concatenate((features,  nrp_stachelhaus_features), axis=1)
features = np.concatenate((features,  nrp_nrpspredictor_features), axis=1)
features = np.concatenate((features,  nrp_pHMM_features), axis=1)
features = np.concatenate((features,  nrp_predicat_features), axis=1)
features = np.concatenate((features,  nrp_sandpuma_features), axis=1)

#process features for chosen classification
y_vars = []
if classification == "antibacterial":
    y_vars = is_antibacterial
    y_vars = y_vars[is_not_unknown_indices]
    features = features[is_not_unknown_indices,:]
if classification == "antieuk":
    y_vars = is_antieuk
    y_vars = y_vars[is_not_unknown_indices]
    features = features[is_not_unknown_indices,:]
if classification == "antifungal":
    y_vars = (is_antifungal >= 1).astype(int)
    y_vars = y_vars[is_not_unknown_indices]
    features = features[is_not_unknown_indices,:]
if classification == "cytotoxic_antitumor":
    y_vars = (is_cytotoxic >= 1).astype(int)
    y_vars = y_vars[is_not_unknown_indices]
    features = features[is_not_unknown_indices,:]
if classification == "antigramneg":
    y_vars = is_gram_neg
    y_vars = y_vars[is_not_unknown_indices_gram]
    features = features[is_not_unknown_indices_gram,:]
if classification == "antigrampos":
    y_vars = is_gram_pos
    y_vars = y_vars[is_not_unknown_indices_gram]
    features = features[is_not_unknown_indices_gram,:]
    
#randomize order of features
new_order = makeRandomOrder(0, y_vars)
y_vars = y_vars[new_order]
features = features[new_order,:]

#parameter values to test for logistic regression
l1_ratios = [.5, .2, .05, .1, .01, .001, .0001]
alphas = [.5, .3, .2, .1, .01,  .001,  .0001,.00001,.000001]
log_accuracies = {}
max_iter = 100
loss = "log"
reg= "elasticnet"

#test parameters for logistic regression
for l1_ratio in l1_ratios:
    log_accuracies[l1_ratio] =  {}
    accuracy_outfile = open("classifier_optimization/log_"+classification+loss+"_"+reg+str(l1_ratio)+"_alpha_vs_accuracy.txt",'w')
    for a in alphas:
        classifier = SGDClassifier(loss=loss,penalty=reg,max_iter=max_iter,alpha=a,l1_ratio=l1_ratio,tol=None)
        accuracy_outfile.write(str(a) + ",")
        log_accuracies[l1_ratio][a]= validateClassifier(accuracy_outfile, classifier, features, y_vars, True)
        accuracy_outfile.write("\n")
    accuracy_outfile.close()
    
    
#parameter values to test for SVM classifier
kernels = ['linear', 'rbf']
c_values = [100, 10, 1, .5, .1, .01]
gammas = [.01, .1, 1, 10] #for rbf only

svm_accuracies = {}
#test parameters
for k in kernels:
    if k == "linear":
        svm_accuracies[k] = {}
        accuracy_outfile = open("classifier_optimization/SVM_"+classification + "_" +k+"_c_vs_accuracy.txt",'w')
        for c in c_values:
            classifier = SVC(kernel=k,C=c)
            accuracy_outfile.write(str(c) + ",")
            svm_accuracies[k][c]=validateClassifier(accuracy_outfile, classifier, features, y_vars, True)  
            accuracy_outfile.write("\n")
        accuracy_outfile.close()
    if k == "rbf":
        svm_accuracies[k] = {}
        for g in gammas:
            accuracy_outfile = open("classifier_optimization/SVM_"+classification + "_" +k+"_"+str(g)+"_c_vs_accuracy.txt",'w')
            svm_accuracies[k][g] = {}
            for c in c_values:
                classifier = SVC(kernel=k,C=c, gamma=g)
                accuracy_outfile.write(str(c) + ",")
                svm_accuracies[k][g][c] = validateClassifier(accuracy_outfile, classifier, features, y_vars,True)
                accuracy_outfile.write("\n")
            accuracy_outfile.close()
            
random.seed(1)            
#parameters for random forest classifier
n_estimators = [1, 5, 10, 15, 25, 50, 100]
max_depths = [10, 20, 50, 100, 1000, None]
bootstrap = True
max_features = 'auto'
criterion = 'gini'

#optimize random forest classifier
forest_accuracies = {}
for d in max_depths:
    forest_accuracies[d] = {}
    accuracy_outfile = open("classifier_optimization/RT_"+classification + "_" +str(d)+"_estimators_vs_accuracy.txt",'w')
    for n in n_estimators:
        accuracy_outfile.write(str(n)+",")
        classifier = ExtraTreesClassifier(n_estimators=n, max_depth=d,bootstrap=bootstrap, max_features=max_features, criterion=criterion, random_state=0)
        forest_accuracies[d][n] = validateClassifier(accuracy_outfile, classifier, features, y_vars, False)
        accuracy_outfile.write("\n")
    accuracy_outfile.close()
    
#visualize results
fig, axs = plt.subplots(3, figsize=(20, 10))

cmap = cm.get_cmap('jet')
colors = []
for i in range(0, 7):
    colors.append(cmap(i*1.0/7.0))
    
i = 0

for l1_ratio in l1_ratios:
    y_values = []
    for a in alphas:
        y_values.append(log_accuracies[l1_ratio][a])
    axs[0].plot(alphas, y_values,  color=colors[i], label="l1: " + str(l1_ratio))
    i += 1

axs[0].legend()
axs[0].set_ylabel('accuracy')
axs[0].set_xlabel('alpha')
axs[0].set_xscale('log')
axs[0].set_title('Log Regression Elastic Net')

y_values = []
for c in c_values:
    y_values.append(svm_accuracies['linear'][c])
axs[1].plot(c_values, y_values, color=colors[0], label='linear')

i = 1    
for g in gammas:
    y_values = []
    for c in c_values:
        y_values.append(svm_accuracies['rbf'][g][c])
    axs[1].plot(c_values, y_values, color=colors[i], label='rbf gamma:' + str(g))
    i += 1

axs[1].legend()
axs[1].set_ylabel('accuracy')
axs[1].set_xlabel('C')
axs[1].set_xscale('log')
axs[1].set_title('SVM')

i = 0
for d in max_depths:
    y_values = []
    for n in n_estimators:
        y_values.append(forest_accuracies[d][n])
    axs[2].plot(n_estimators,y_values, color=colors[i], label='max depth :' + str(d))
    i += 1

axs[2].legend()
axs[2].set_ylabel('accuracy')
axs[2].set_xlabel('Num estimators')
axs[2].set_title('Extra Random Trees')

fig.tight_layout(pad=2.0)

#print optimal parameters
max_accuracy_log = 0
optimal_alpha = 0
optimal_l1_ratio = 0
for l1_ratio in l1_ratios:
    for a in alphas:
        if log_accuracies[l1_ratio][a] > max_accuracy_log:
            max_accuracy_log = log_accuracies[l1_ratio][a]
            optimal_alpha = a
            optimal_l1_ratio = l1_ratio
            
print("Log regression maximum accuracy: " + str(max_accuracy_log))
print("Optimal l1 ratio: " + str(optimal_l1_ratio) + " optimal alpha: " + str(optimal_alpha))
print("")

max_accuracy_svm = 0
optimal_kernel = ""
optimal_gamma = 0
optimal_c = 0
for c in c_values:
    if svm_accuracies['linear'][c] > max_accuracy_svm:
        max_accuracy_svm = svm_accuracies['linear'][c]
        optimal_kernel = 'linear'
        optimal_c = c
    for g in gammas:
        if svm_accuracies['rbf'][g][c] > max_accuracy_svm:
            max_accuracy_svm = svm_accuracies['rbf'][g][c]
            optimal_kernel = 'rbf'
            optimal_c = c
            optimal_gamma = g
            
print("SVM maximum accuracy: " + str(max_accuracy_svm))
if optimal_kernel == 'linear':
    print("Optimal kernel: linear optimal c: " + str(optimal_c))
else:
    print("Optimal kernel: rbf optimal c: " + str(optimal_c) + " optimal gamma: " + str(optimal_gamma))
print("")

max_accuracy_forest = 0
optimal_depth = 0
optimal_n = 0
for d in max_depths:
    for n in n_estimators:
        if forest_accuracies[d][n] > max_accuracy_forest:
            max_accuracy_forest = forest_accuracies[d][n]
            optimal_depth = d
            optimal_n = n
print("Extra random trees maximum accuracy: " + str(max_accuracy_forest))
print("Optimal depth: " + str(d) + " optimal n estimators: " + str(optimal_n))