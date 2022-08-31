# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 14:10:16 2022

@author: Allison Walker
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
from sklearn.preprocessing import MinMaxScaler
from sklearn.pipeline import Pipeline
from sklearn.model_selection import GridSearchCV, KFold
import matplotlib.pyplot as plt
from matplotlib import cm
from sklearn.metrics import balanced_accuracy_score
from sklearn.metrics import roc_curve
from sklearn.metrics import precision_recall_curve
import sys
import joblib
import train_model_tools

#writes cv results from GridSearchCV
def writeCVResults(filename, cv_results, var1, var2,num_cv=5):
    param_sets = cv_results["params"]
    var1_list = []
    var2_list = []
    var_combo_index = {}
    i= 0
    for p in param_sets:
        if p[var1] not in var1_list:
            var1_list.append(p[var1])
            var_combo_index[p[var1]] = {}
        if p[var2] not in var2_list:
            var2_list.append(p[var2])
        var_combo_index[p[var1]][p[var2]] = i
        i += 1
    for v in var1_list:
        outfile = open(filename + "_" + str(v) + ".txt", 'w')
        for v2 in var2_list:
            outfile.write(str(v2) + ",")
            for i in range(0, num_cv):
                score = cv_results["split"+str(i)+"_test_score"][var_combo_index[v][v2]]
                outfile.write(str(score) + ",")
            outfile.write("\n")
        outfile.close()
            
def writeCVResults3Vars(filename, cv_results, var1, var2, var3, num_cv=5):
    param_sets = cv_results["params"]
    var1_list = []
    var2_list = []
    var3_list = []
    var_combo_index = {}
    i = 0
    for p in param_sets:
        if p[var1] not in var1_list:
            var1_list.append(p[var1])
            var_combo_index[p[var1]] = {}
        if var2 not in p:
            if p[var3] not in var3_list:
                var3_list.append(p[var3])
            if "NA" not in var_combo_index[p[var1]]:
                var_combo_index[p[var1]]["NA"] = {}
            var_combo_index[p[var1]]["NA"][p[var3]] = i           
        else:
            if p[var2] not in var2_list:
                var2_list.append(p[var2])    
            if p[var2] not in var_combo_index[p[var1]]:
                var_combo_index[p[var1]][p[var2]] = {}
            if p[var3] not in var3_list:
                var3_list.append(p[var3])
            var_combo_index[p[var1]][p[var2]][p[var3]] = i  
        i += 1
    for v1 in var1_list:
        for v2 in var_combo_index[v1]:
            outfile = open(filename + "_" + str(v1) + "_" + str(v2) + ".txt", 'w')
            for v3 in var_combo_index[v1][v2]:
                outfile.write(str(v3) + ",")
                for i in range(0, num_cv):
                    score = cv_results["split"+str(i)+"_test_score"][var_combo_index[v1][v2][v3]]
                    outfile.write(str(score) + ",")
                outfile.write("\n")
            outfile.close()

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
#TODO: change to iterate through all classifications
classifications = ["antibacteria", "antieuk","antifungal", "cytotoxic_antitumor", "antigramneg","antigrampos"]
classification = "antibacterial" 
#set random seed so results are consistent
random.seed(1)



training_set_dir = "feature_matrices/antismash4rgi3"
training_set_name = training_set_dir[training_set_dir.find("/"):len(training_set_dir)]
#TODO: add feature directory as an argument
#TODO: test antiSMASH6
feature_dir = training_set_dir  + "/features/"
feature_type_list = readFeatureFiles.getFeatureFilesList(feature_dir)
features = readFeatureFiles.readFeatures(feature_dir, feature_type_list)
feature_list = readFeatureFiles.readFeatureNames(feature_dir, feature_type_list)

#read classes
#TODO: change so that these are read from training set directory
#TODO: change is_unknown to give information about specific labels that migth be unknown
is_antibacterial = readFeatureFiles.readClassesMatrix("feature_matrices/antismash4rgi3/classifications/is_antibacterial.csv")
is_antifungal = readFeatureFiles.readClassesMatrix("feature_matrices/antismash4rgi3/classifications/is_antifungal.csv")
is_cytotoxic = readFeatureFiles.readClassesMatrix("feature_matrices/antismash4rgi3/classifications/is_cytotoxic.csv")
is_unknown = readFeatureFiles.readClassesMatrix("feature_matrices/antismash4rgi3/classifications/is_unknown.csv")
targets_gram_pos = readFeatureFiles.readClassesMatrix("feature_matrices/antismash4rgi3/classifications/targets_gram_pos.csv")
targets_gram_neg = readFeatureFiles.readClassesMatrix("feature_matrices/antismash4rgi3/classifications/targets_gram_neg.csv")
full_cluster_list = readFeatureFiles.readClusterList(training_set_dir + "/cluster_list.txt")
is_not_unknown_indices = readFeatureFiles.getNotUnknownIndices(is_unknown)
target_unannotated = is_antibacterial*((targets_gram_pos+targets_gram_neg)<1)
is_not_unknown_indices_gram =  readFeatureFiles.getNotUnknownIndices(is_unknown + target_unannotated)

is_antibacterial = (is_antibacterial >= 1).astype(int)
is_antieuk = ((is_antifungal + is_cytotoxic)>=1).astype(int)
is_gram_pos = (targets_gram_pos >= 1).astype(int)
is_gram_neg = (targets_gram_neg >= 1).astype(int)


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
#TODO: change variable search to a grid search
log_params = [{'log__loss':['log'], 'log__penalty':['elasticnet'], 'log__max_iter':[100],\
               'log__alpha':[.5, .3, .2, .1, .01,  .001,  .0001,.00001,.000001],\
              'log__l1_ratio':[.5, .2, .05, .1, .01, .001, .0001], 'log__tol':[None]}]
log_pipeline = Pipeline([('mms',MinMaxScaler()),('log',SGDClassifier())])
gs_log = GridSearchCV(log_pipeline,param_grid=log_params,scoring='accuracy',cv=5)
gs_log.fit(features, y_vars)
writeCVResults("classifier_optimization/log_"+classification,gs_log.cv_results_,"log__l1_ratio", "log__alpha")


    
#parameter values to test for SVM classifier
c_values = [100, 10, 1, .5, .1, .01]
svm_params = [{'svm__kernel':['linear'], 'svm__C':c_values, 'svm__probability':[True]}, \
              {'svm__kernel':['rbf'], 'svm__C':c_values, 'svm__gamma':[.01, .1, 1, 10], 'svm__probability':[True]}]
svm_pipeline = Pipeline([('mms',MinMaxScaler()),('svm',SVC())])
gs_svm = GridSearchCV(svm_pipeline, param_grid=svm_params, scoring='accuracy',cv=5)
gs_svm.fit(features, y_vars)
writeCVResults3Vars("classifier_optimization/svm_"+classification,gs_svm.cv_results_,"svm__kernel","svm__gamma", "svm__C")

            
random.seed(1)            
#parameters for random forest classifier
tree_params = [{'max_features':['auto'], 'criterion':['gini'], 'bootstrap':[True],\
                'max_depth':[10, 20, 50, 100, 1000, None], \
                'n_estimators':[1, 5, 10, 15, 25, 50, 100]}]
gs_tree = GridSearchCV(ExtraTreesClassifier(), param_grid=tree_params, scoring='accuracy',cv=5)
gs_tree.fit(features, y_vars)
writeCVResults("classifier_optimization/RF_"+classification,gs_tree.cv_results_,"max_depth", "n_estimators")


#
#print optimal parameters
print("Logistic regression")
print("optimal accuracy: " + str(gs_log.best_score_))
print("best params: " + str(gs_log.best_params_))
print()
print("Support vector machine")
print("optimal accuracy: " + str(gs_svm.best_score_))
print("best params: " + str(gs_svm.best_params_))
print()
print("Random forest")
print("optimal accuracy: " + str(gs_tree.best_score_))
print("best params: " + str(gs_tree.best_params_))
print()


#asses model
best_log = gs_log.best_estimator_
best_svm = gs_svm.best_estimator_
best_tree = gs_tree.best_estimator_

#set random seed so results are consistent
random.seed(1)
#include features from SSN 
include_SSN = True 
#parameters for classifiers
#TODO: get optimal parameters from search above


#do analysis, write to file, and visualize
b_accuracies, roc_curves, pr_curves = train_model_tools.assessModel(classification, best_svm, best_log, best_tree, features, y_vars)


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
fig.savefig('classifier_optimization/' + classification + '.pdf')  

#train model on entire dataset and write model
best_svm.fit(features, y_vars)
best_tree.fit(features, y_vars)
best_log.fit(features, y_vars)

#print train accuracy
print("SVM train balanced accuracy: " +str(balanced_accuracy_score(y_vars,best_svm.predict(features))))
print("Log train balanced accuracy: " + str(balanced_accuracy_score(y_vars,best_log.predict(features))))
print("Log train balanced accuracy: " + str(balanced_accuracy_score(y_vars,best_tree.predict(features))))

outfilename = "trained_models/" + training_set_name +"_" + classification +".sav"
joblib.dump([best_svm, best_tree, best_log], outfilename)