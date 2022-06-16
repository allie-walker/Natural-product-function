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
from sklearn.metrics import balanced_accuracy_score
from sklearn.metrics import roc_curve
from sklearn.metrics import precision_recall_curve

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

def randomizeFeatures(seed, features):
    np.random.seed(seed)
    random_features = features
    for i in range(0, features.shape[1]):
        y_vars_range = np.arange(features.shape[0])
        np.random.shuffle(y_vars_range)
        random_features[y_vars_range,:]
    return random_features

def writeCurve(outfile, x_axis, y_axis):
    for i in range(0, min(x_axis.shape[0], y_axis.shape[0])):
        outfile.write(str(x_axis[i]) + "," + str(y_axis[i]) + "\n")
        
#parameters
#classification to train classifiers on 
#options: antibacterial, antieuk (defined as antifungal, antitumor, or cytotoxic), antifungal, cytotoxic_antitumor, antigramneg, antigrampos
#TODO: change to iterate through all classifications
classification = "antibacterial" 
#set random seed so results are consistent
random.seed(1)
#include features from SSN 
include_SSN = True

#TODO: add csv file list and label list as argument inputs, should be read from directory
#TODO: test antiSMASH6
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
#TODO: add argument with directory with class labels
#TODO: change is_unknown to give information about specific labels that migth be unknown
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
if classification == "antieuk":
    y_vars = is_antieuk
    y_vars = y_vars[is_not_unknown_indices]
if classification == "antifungal":
    y_vars = (is_antifungal >= 1).astype(int)
    y_vars = y_vars[is_not_unknown_indices]
if classification == "cytotoxic_antitumor":
    y_vars = (is_cytotoxic >= 1).astype(int)
    y_vars = y_vars[is_not_unknown_indices]
if classification == "antigramneg":
    y_vars = is_gram_neg
    y_vars = y_vars[is_not_unknown_indices_gram]
if classification == "antigrampos":
    y_vars = is_gram_pos
    y_vars = y_vars[is_not_unknown_indices_gram]
    
#randomize order of features
new_order = makeRandomOrder(0, y_vars)
y_vars = y_vars[new_order]
features = features[new_order,:]

#parameter values to test for logistic regression
#TODO: change variable search to a grid search
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


#scoring model starts here
#parameters
#TODO: iterate through all classifications
classification = "antibacterial" 
#set random seed so results are consistent
random.seed(1)
#include features from SSN 
include_SSN = True 
#parameters for classifiers
#TODO: get optimal parameters from search above
log_params = {}
log_params["antibacterial"] = {"l1_ratio":.05,"alpha":.01}
log_params["antigrampos"] = {"l1_ratio":.001,"alpha":.001}
log_params["antigramneg"] = {"l1_ratio":.05,"alpha":.001}
log_params["antieuk"] = {"l1_ratio":.001,"alpha":.001}
log_params["antifungal"] = {"l1_ratio":.0001,"alpha":.01}
log_params["cytotoxic_antitumor"] = {"l1_ratio":.001,"alpha":.001}

svm_params = {}
svm_params["antibacterial"] = {"kernel":'rbf',"C":10,"gamma":0.1}
svm_params["antigrampos"] = {"kernel":'rbf',"C":10,"gamma":.01}
svm_params["antigramneg"] = {"kernel":'rbf',"C":10,"gamma":0.01}
svm_params["antieuk"] = {"kernel":'linear',"C":.1}
svm_params["antifungal"] = {"kernel":'rbf',"C":10,"gamma":0.1}
svm_params["cytotoxic_antitumor"] = {"kernel":'rbf',"C":100,"gamma":0.1}

forest_params = {}
forest_params["antibacterial"] = {"depth":100,"n":50}
forest_params["antigrampos"] = {"depth":100,"n":50}
forest_params["antigramneg"] = {"depth":100,"n":25}
forest_params["antieuk"] = {"depth":None,"n":25}
forest_params["antifungal"] = {"depth":50,"n":50}
forest_params["cytotoxic_antitumor"] = {"depth":50,"n":100}


#make randomly shuffled version of features
features_rand = randomizeFeatures(0, features)

y_vars = []
if classification == "antibacterial":
    y_vars = is_antibacterial
    y_vars = y_vars[is_not_unknown_indices]
    
if classification == "antieuk":
    y_vars = is_antieuk
    y_vars = y_vars[is_not_unknown_indices]

    
if classification == "antifungal":
    y_vars = (is_antifungal >= 1).astype(int)
    y_vars = y_vars[is_not_unknown_indices]
    
if classification == "cytotoxic_antitumor":
    y_vars = (is_cytotoxic >= 1).astype(int)
    y_vars = y_vars[is_not_unknown_indices]
    
if classification == "antigramneg":
    y_vars = is_gram_neg
    y_vars = y_vars[is_not_unknown_indices_gram]
    
if classification == "antigrampos":
    y_vars = is_gram_pos
    y_vars = y_vars[is_not_unknown_indices_gram]
    
#reorder features
new_order = makeRandomOrder(0, y_vars)
y_vars = y_vars[new_order]
features = features[new_order,:]
features_rand = features_rand[new_order, :]

#initialize classifiers
opt_log_params = log_params[classification]
log_classifier = SGDClassifier(loss='log',penalty='elasticnet',max_iter=100,alpha=opt_log_params["alpha"],l1_ratio=opt_log_params["l1_ratio"],tol=None)

opt_svm_params = svm_params[classification]
if opt_svm_params['kernel'] == "linear":
    svm_classifier = SVC(kernel="linear",C=opt_svm_params["C"], probability=True)
else:
    svm_classifier = SVC(kernel="rbf",C=opt_svm_params["C"], gamma=opt_svm_params["gamma"], probability=True)
    
opt_forest_params = forest_params[classification]
tree_classifier = ExtraTreesClassifier(bootstrap=True,max_features="auto",n_estimators=opt_forest_params["n"], max_depth=opt_forest_params["depth"], random_state=0)

#do analysis, write to file, and visualize
output_fname_base = "classifier_metrics/" + classification + "_"

#output files
roc_out = open(output_fname_base + "roc.txt", 'w')
pr_out = open(output_fname_base + "precision_recall.txt", 'w')
accuracy_out = open(output_fname_base + "accuracy.txt", "w")
balanced_accuracy_out = open(output_fname_base + "balanced_accuracy.txt", "w")
#TODO: also output precision, recall, and aroc
#TODO: make separate output file for each classification type

#cross validation
#store results in dictionary for visualization
roc_curves = {}
roc_curves["log"] = {}
roc_curves["svm"] = {}
roc_curves["tree"] = {}
roc_curves["rnd_log"] = {}
roc_curves["rnd_svm"] = {}
roc_curves["rnd_tree"] = {}

pr_curves = {}
pr_curves["log"] = {}
pr_curves["svm"] = {}
pr_curves["tree"] = {}
pr_curves["rnd_log"] = {}
pr_curves["rnd_svm"] = {}
pr_curves["rnd_tree"] = {}

accuracies = {}
accuracies["log"] = []
accuracies["svm"] = []
accuracies["tree"] = []
accuracies["rnd_log"] = []
accuracies["rnd_svm"] = []
accuracies["rnd_tree"] = []

b_accuracies = {}
b_accuracies["log"] = []
b_accuracies["svm"] = []
b_accuracies["tree"] = []
b_accuracies["rnd_log"] = []
b_accuracies["rnd_svm"] = []
b_accuracies["rnd_tree"] = []

for i in range(0, 10):
    (training_x, training_y, val_x, val_y) = splitData(features, y_vars, i, 10)
    (rnd_training_x, rnd_training_y, rnd_val_x, rnd_val_y) = splitData(features_rand, y_vars, i, 10)
    min_max_scaler = preprocessing.MinMaxScaler()
    scaled_training_x = min_max_scaler.fit_transform(training_x)
    scaled_val_x = min_max_scaler.transform(val_x)
    scaled_rnd_training_x = min_max_scaler.fit_transform(rnd_training_x)
    scaled_rnd_val_x = min_max_scaler.transform(rnd_val_x)
    
    #metrics for classifiers fit to real data
    svm_classifier.fit(scaled_training_x, training_y)
    tree_classifier.fit(training_x, training_y)
    log_classifier.fit(scaled_training_x, training_y)
    
    accuracy_out.write(str(svm_classifier.score(scaled_val_x, val_y))+",")
    accuracy_out.write(str(log_classifier.score(scaled_val_x, val_y))+",")
    accuracy_out.write(str(tree_classifier.score(val_x, val_y))+",")
    accuracies["log"].append(log_classifier.score(scaled_val_x, val_y))
    accuracies["svm"].append(svm_classifier.score(scaled_val_x, val_y))
    accuracies["tree"].append(tree_classifier.score(val_x, val_y))
    
    balanced_accuracy_out.write(str(balanced_accuracy_score(val_y,svm_classifier.predict(scaled_val_x))) + ",")
    balanced_accuracy_out.write(str(balanced_accuracy_score(val_y,log_classifier.predict(scaled_val_x))) + ",")
    balanced_accuracy_out.write(str(balanced_accuracy_score(val_y,tree_classifier.predict(val_x))) + ",")
    b_accuracies["log"].append(balanced_accuracy_score(val_y,log_classifier.predict(scaled_val_x)))
    b_accuracies["svm"].append(balanced_accuracy_score(val_y,svm_classifier.predict(scaled_val_x)))
    b_accuracies["tree"].append(balanced_accuracy_score(val_y,tree_classifier.predict(val_x)))
    
    tree_probabilities = tree_classifier.predict_proba(val_x)
    log_probabilities = log_classifier.predict_proba(scaled_val_x)
    svm_probabilities = svm_classifier.predict_proba(scaled_val_x)
    
    fpr, tpr, thresholds = roc_curve(val_y, tree_probabilities[:,1])
    roc_curves["tree"][i] = (fpr, tpr)
    fpr, tpr, thresholds = roc_curve(val_y, svm_probabilities[:,1])
    roc_curves["svm"][i] = (fpr, tpr)
    fpr, tpr, thresholds = roc_curve(val_y, log_probabilities[:,1])
    roc_curves["log"][i] = (fpr, tpr)
    
    precision, recall, thresholds = precision_recall_curve(val_y, tree_probabilities[:,1])
    pr_curves["tree"][i] = (recall, precision)
    precision, recall, thresholds = precision_recall_curve(val_y, svm_probabilities[:,1])
    pr_curves["svm"][i] = (recall, precision)
    precision, recall, thresholds = precision_recall_curve(val_y, log_probabilities[:,1])
    pr_curves["log"][i] = (recall, precision)
    
    
    #metrics for classifiers fit to randomly scrambled data
    svm_classifier.fit(scaled_rnd_training_x, rnd_training_y)
    tree_classifier.fit(rnd_training_x, rnd_training_y)
    log_classifier.fit(scaled_rnd_training_x, rnd_training_y)
    
    accuracy_out.write(str(svm_classifier.score(scaled_rnd_val_x, rnd_val_y))+",")
    accuracy_out.write(str(log_classifier.score(scaled_rnd_val_x, rnd_val_y))+",")
    accuracy_out.write(str(tree_classifier.score(rnd_val_x, rnd_val_y))+",")
    accuracies["rnd_log"].append(log_classifier.score(scaled_rnd_val_x, rnd_val_y))
    accuracies["rnd_svm"].append(svm_classifier.score(scaled_rnd_val_x, rnd_val_y))
    accuracies["rnd_tree"].append(tree_classifier.score(rnd_val_x, rnd_val_y))
    
    balanced_accuracy_out.write(str(balanced_accuracy_score(rnd_val_y,svm_classifier.predict(scaled_rnd_val_x))) + ",")
    balanced_accuracy_out.write(str(balanced_accuracy_score(rnd_val_y,log_classifier.predict(scaled_rnd_val_x))) + ",")
    balanced_accuracy_out.write(str(balanced_accuracy_score(rnd_val_y,tree_classifier.predict(rnd_val_x))) + ",")
    b_accuracies["rnd_log"].append(balanced_accuracy_score(rnd_val_y,log_classifier.predict(scaled_rnd_val_x)))
    b_accuracies["rnd_svm"].append(balanced_accuracy_score(rnd_val_y,svm_classifier.predict(scaled_rnd_val_x)))
    b_accuracies["rnd_tree"].append(balanced_accuracy_score(rnd_val_y,tree_classifier.predict(rnd_val_x)))
    
    tree_probabilities = tree_classifier.predict_proba(rnd_val_x)
    log_probabilities = log_classifier.predict_proba(scaled_rnd_val_x)
    svm_probabilities = svm_classifier.predict_proba(scaled_rnd_val_x)
    
    fpr, tpr, thresholds = roc_curve(rnd_val_y, tree_probabilities[:,1])
    roc_curves["rnd_tree"][i] = (fpr, tpr)
    fpr, tpr, thresholds = roc_curve(rnd_val_y, svm_probabilities[:,1])
    roc_curves["rnd_svm"][i] = (fpr, tpr)
    fpr, tpr, thresholds = roc_curve(rnd_val_y, log_probabilities[:,1])
    roc_curves["rnd_log"][i] = (fpr, tpr)
    
    precision, recall, thresholds = precision_recall_curve(rnd_val_y, tree_probabilities[:,1])
    pr_curves["rnd_tree"][i] = (recall, precision)
    recision, recall, thresholds = precision_recall_curve(rnd_val_y, svm_probabilities[:,1])
    pr_curves["rnd_svm"][i] = (recall, precision)
    precision, recall, thresholds = precision_recall_curve(rnd_val_y, log_probabilities[:,1])
    pr_curves["rnd_log"][i] = (recall, precision)
    
    balanced_accuracy_out.write("\n")
    accuracy_out.write("\n")
    
accuracy_out.close()
balanced_accuracy_out.close()

#write ROC and pr curves for one trial
writeCurve(roc_out, roc_curves["svm"][0][0],roc_curves["svm"][0][1])
writeCurve(roc_out, roc_curves["log"][0][0],roc_curves["log"][0][1])
writeCurve(roc_out, roc_curves["tree"][0][0],roc_curves["tree"][0][1])
writeCurve(roc_out, roc_curves["rnd_svm"][0][0],roc_curves["rnd_svm"][0][1])
writeCurve(roc_out, roc_curves["rnd_log"][0][0],roc_curves["rnd_log"][0][1])
writeCurve(roc_out, roc_curves["rnd_tree"][0][0],roc_curves["rnd_tree"][0][1])
roc_out.close()

writeCurve(pr_out, pr_curves["svm"][0][0],pr_curves["svm"][0][1])
writeCurve(pr_out, pr_curves["log"][0][0],pr_curves["log"][0][1])
writeCurve(pr_out, pr_curves["tree"][0][0],pr_curves["tree"][0][1])
writeCurve(pr_out, pr_curves["rnd_svm"][0][0],pr_curves["rnd_svm"][0][1])
writeCurve(pr_out, pr_curves["rnd_log"][0][0],pr_curves["rnd_log"][0][1])
writeCurve(pr_out, pr_curves["rnd_tree"][0][0],pr_curves["rnd_tree"][0][1])
pr_out.close

#visualize results
fig, axs = plt.subplots(3, figsize=(15, 12))
cmap = cm.get_cmap('jet')
colors = []
for i in range(0, 6):
    colors.append(cmap(i*1.0/6.0))
    
mean_log_acc = np.mean(b_accuracies["log"])
mean_svm_acc = np.mean(b_accuracies["svm"])
mean_tree_acc = np.mean(b_accuracies["tree"])
mean_rnd_log_acc =np.mean(b_accuracies["rnd_log"])
mean_rnd_svm_acc = np.mean(b_accuracies["rnd_svm"])
mean_rnd_tree_acc = np.mean(b_accuracies["rnd_tree"])

mean_log_sd = np.std(b_accuracies["log"])
mean_svm_sd = np.std(b_accuracies["svm"])
mean_tree_sd = np.std(b_accuracies["tree"])
mean_rnd_log_sd =np.std(b_accuracies["rnd_log"])
mean_rnd_svm_sd = np.std(b_accuracies["rnd_svm"])
mean_rnd_tree_sd = np.std(b_accuracies["rnd_tree"])

bar_labels = ["log", "svm", "tree", "rnd log", "rnd svm", "rnd tree"]
x_pos = np.arange(len(bar_labels))
means = [mean_log_acc, mean_svm_acc, mean_tree_acc, mean_rnd_log_acc, mean_rnd_svm_acc, mean_rnd_tree_acc]
error = [mean_log_sd, mean_svm_sd, mean_tree_sd, mean_rnd_log_sd, mean_rnd_svm_sd, mean_rnd_tree_sd]
axs[0].bar(x_pos, means, yerr=error, align='center', capsize=10, color=(colors[0], colors[2], colors[5],colors[0], colors[2], colors[5]))
axs[0].set_xticks(x_pos)
axs[0].set_xticklabels(bar_labels)
axs[0].set_title("Balanced Accuracy of Classifiers")
axs[0].set_ylabel("Balanced Accuracy")

axs[1].plot(roc_curves["log"][0][0], roc_curves["log"][0][1], color=colors[0], label='log')
axs[1].plot(roc_curves["svm"][0][0], roc_curves["svm"][0][1], color=colors[2], label='svm')
axs[1].plot(roc_curves["tree"][0][0], roc_curves["tree"][0][1], color=colors[5], label='tree')
for i in range(1, 10):
    axs[1].plot(roc_curves["log"][i][0], roc_curves["log"][i][1], color=colors[0])
    axs[1].plot(roc_curves["svm"][i][0], roc_curves["svm"][i][1], color=colors[2])
    axs[1].plot(roc_curves["tree"][i][0], roc_curves["tree"][i][1], color=colors[5])
axs[1].set_title("ROC")
axs[1].set_ylabel("TPR")
axs[1].set_xlabel("FPR")
axs[1].legend()


axs[2].plot(pr_curves["log"][0][0], pr_curves["log"][0][1], color=colors[0], label='log')
axs[2].plot(pr_curves["svm"][0][0], pr_curves["svm"][0][1], color=colors[2], label='svm')
axs[2].plot(pr_curves["tree"][0][0], pr_curves["tree"][0][1], color=colors[5], label='tree')
for i in range(1, 10):
    axs[2].plot(pr_curves["log"][i][0], pr_curves["log"][i][1], color=colors[0])
    axs[2].plot(pr_curves["svm"][i][0], pr_curves["svm"][i][1], color=colors[2])
    axs[2].plot(pr_curves["tree"][i][0], pr_curves["tree"][i][1], color=colors[5])
axs[2].set_title("precision-recall curve")
axs[2].set_ylabel("precision")
axs[2].set_xlabel("recall")
axs[2].legend()
fig.tight_layout(pad=2.0)