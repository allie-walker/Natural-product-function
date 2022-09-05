# -*- coding: utf-8 -*-
"""
Created on Fri Aug 26 15:27:28 2022

@author: allison
"""

from sklearn.model_selection import KFold
import random
from sklearn.metrics import balanced_accuracy_score
from sklearn.metrics import roc_curve
from sklearn.metrics import precision_recall_curve
import numpy as np

def randomizeFeatures(seed, features):
    np.random.seed(seed)
    random_features = features
    y_vars_range = np.arange(features.shape[0])
    np.random.shuffle(y_vars_range)
    random_features = random_features[y_vars_range,:]        
    return random_features

def writeCurve(outfile, x_axis, y_axis):
    for i in range(0, min(x_axis.shape[0], y_axis.shape[0])):
        outfile.write(str(x_axis[i]) + "," + str(y_axis[i]) + "\n")

def assessModel(output_fname_base, svm_model, log_model, tree_model, features, y_vars):
    random.seed(1)
    features_rand = randomizeFeatures(0, features)   
    
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

    kf = KFold(n_splits=10, shuffle=True, random_state=0)
    data_split = kf.split(X=features, y=y_vars)
    i = 0
    for train_index, test_index in data_split:
        training_x, val_x = features[train_index], features[test_index]
        training_y, val_y = y_vars[train_index], y_vars[test_index]
        rnd_training_x, rnd_val_x = features_rand[train_index], features_rand[test_index]
        log_model.fit(training_x, training_y)
        svm_model.fit(training_x, training_y)
        tree_model.fit(training_x, training_y)
        
        accuracy_out.write(str(svm_model.score(val_x, val_y))+",")
        accuracy_out.write(str(log_model.score(val_x, val_y))+",")
        accuracy_out.write(str(tree_model.score(val_x, val_y))+",")
        accuracies["log"].append(log_model.score(val_x, val_y))
        accuracies["svm"].append(svm_model.score(val_x, val_y))
        accuracies["tree"].append(tree_model.score(val_x, val_y))
        
        balanced_accuracy_out.write(str(balanced_accuracy_score(val_y,svm_model.predict(val_x))) + ",")
        balanced_accuracy_out.write(str(balanced_accuracy_score(val_y,log_model.predict(val_x))) + ",")
        balanced_accuracy_out.write(str(balanced_accuracy_score(val_y,tree_model.predict(val_x))) + ",")
        b_accuracies["log"].append(balanced_accuracy_score(val_y,log_model.predict(val_x)))
        b_accuracies["svm"].append(balanced_accuracy_score(val_y,svm_model.predict(val_x)))
        b_accuracies["tree"].append(balanced_accuracy_score(val_y,tree_model.predict(val_x)))
        
        tree_probabilities = tree_model.predict_proba(val_x)
        log_probabilities = log_model.predict_proba(val_x)
        svm_probabilities = svm_model.predict_proba(val_x)
        
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
        svm_model.fit(rnd_training_x, training_y)
        tree_model.fit(rnd_training_x, training_y)
        log_model.fit(rnd_training_x, training_y)
        
        accuracy_out.write(str(svm_model.score(rnd_val_x, val_y))+",")
        accuracy_out.write(str(log_model.score(rnd_val_x, val_y))+",")
        accuracy_out.write(str(tree_model.score(rnd_val_x, val_y))+",")
        accuracies["rnd_log"].append(log_model.score(rnd_val_x, val_y))
        accuracies["rnd_svm"].append(svm_model.score(rnd_val_x, val_y))
        accuracies["rnd_tree"].append(tree_model.score(rnd_val_x, val_y))
        
        balanced_accuracy_out.write(str(balanced_accuracy_score(val_y,svm_model.predict(rnd_val_x))) + ",")
        balanced_accuracy_out.write(str(balanced_accuracy_score(val_y,log_model.predict(rnd_val_x))) + ",")
        balanced_accuracy_out.write(str(balanced_accuracy_score(val_y,tree_model.predict(rnd_val_x))) + ",")
        b_accuracies["rnd_log"].append(balanced_accuracy_score(val_y,log_model.predict(rnd_val_x)))
        b_accuracies["rnd_svm"].append(balanced_accuracy_score(val_y,svm_model.predict(rnd_val_x)))
        b_accuracies["rnd_tree"].append(balanced_accuracy_score(val_y,tree_model.predict(rnd_val_x)))
        
        tree_probabilities = tree_model.predict_proba(rnd_val_x)
        log_probabilities = log_model.predict_proba(rnd_val_x)
        svm_probabilities = svm_model.predict_proba(rnd_val_x)
        
        fpr, tpr, thresholds = roc_curve(val_y, tree_probabilities[:,1])
        roc_curves["rnd_tree"][i] = (fpr, tpr)
        fpr, tpr, thresholds = roc_curve(val_y, svm_probabilities[:,1])
        roc_curves["rnd_svm"][i] = (fpr, tpr)
        fpr, tpr, thresholds = roc_curve(val_y, log_probabilities[:,1])
        roc_curves["rnd_log"][i] = (fpr, tpr)
        
        precision, recall, thresholds = precision_recall_curve(val_y, tree_probabilities[:,1])
        pr_curves["rnd_tree"][i] = (recall, precision)
        recision, recall, thresholds = precision_recall_curve(val_y, svm_probabilities[:,1])
        pr_curves["rnd_svm"][i] = (recall, precision)
        precision, recall, thresholds = precision_recall_curve(val_y, log_probabilities[:,1])
        pr_curves["rnd_log"][i] = (recall, precision)
        
        balanced_accuracy_out.write("\n")
        accuracy_out.write("\n")
        i += 1
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
    return b_accuracies, roc_curves, pr_curves