# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 14:43:23 2019

@author: Allison Walker
"""
import os
from sklearn.svm import SVC
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.linear_model import SGDClassifier
from sklearn import preprocessing
import numpy as np

def checkIfFileExists(filename, file_description):
    #check if file exists
    if not os.path.isfile(filename):
        print(file_description + " file does not exist, please enter a valid file")
        return False
    return True
    
    
def removeReturnChars(feature_list):
    for i in range(0, len(feature_list)):
        feature_list[i] = feature_list[i].replace("\n","").replace("\r","")
    return feature_list

def addToFeatureMatrix(feature_matrix, i, feature_counts, feature_list):
    for f in feature_list:
        if f in feature_counts:
            feature_matrix[0,i] = feature_counts[f]
        else:
            feature_matrix[0,i] = 0

        i += 1
    return (feature_matrix, i)

def treePrediction(training_features, training_y, test_features, params, seed):
    np.random.seed(seed)
    tree_classifier = ExtraTreesClassifier(bootstrap=True,max_features="auto",n_estimators=params["n"], max_depth=params["depth"])
    tree_classifier.fit(training_features, training_y)
    tree_probabilities = tree_classifier.predict_proba(test_features)
    return tree_probabilities

def logPrediction(training_features, training_y, test_features, params, seed):
    np.random.seed(seed)
    log_classifier = SGDClassifier(loss='log',penalty='elasticnet',max_iter=100,tol=None,alpha=params["alpha"],l1_ratio=params["l1_ratio"])
    min_max_scaler = preprocessing.MinMaxScaler()
    scaled_training_x = min_max_scaler.fit_transform(training_features)
    scaled_test_x = min_max_scaler.transform(test_features)
    log_classifier.fit(scaled_training_x, training_y)
    log_probabilities = log_classifier.predict_proba(scaled_test_x)
    return log_probabilities

def svmPrediction(training_features, training_y, test_features, params, seed):
    np.random.seed(seed)
    if params['kernel'] == "linear":
        svm_classifier = SVC(kernel="linear",C=params["C"], probability=True)
    else:
        svm_classifier = SVC(kernel="rbf",C=params["C"], gamma=params["gamma"], probability=True)
    min_max_scaler = preprocessing.MinMaxScaler()
    scaled_training_x = min_max_scaler.fit_transform(training_features)
    scaled_test_x = min_max_scaler.transform(test_features)
    svm_classifier.fit(scaled_training_x, training_y)
    svm_probabilities = svm_classifier.predict_proba(scaled_test_x)
    return svm_probabilities

def writeProbabilitiesToFile(outfile, classification_name, tree_prob, log_prob, svm_prob):
    outfile.write("probabilities of " + classification_name +" activity:\n")
    outfile.write("tree classifier: " + str(tree_prob[0,1])  + " logistic regression classifier: " + str(log_prob[0,1]) + " svm classifier: " + str(svm_prob[0,1]) + "\n")
    return