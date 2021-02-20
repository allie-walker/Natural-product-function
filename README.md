# Natural-product-function

Natural_product_BGC_activity_prediction.ipnyb is an python notebook containing code to optimize, train, and use machine learning models to predict natural product
activity from biosynthetic gene clusters. We have also provided a command line version of the script in the cluster_function_prediction.py. The activity prediction program requires ouptut files from antiSMASH and Resistance Gene Identifier (RGI) (the txt file not the JSON). We used antiSMASH version 4.0 and RGI version 3.2.1 with database version 1.1.9, to produce the results described in our manuscript. Since our initial development of the tool newer versions of antiSMASH and RGI were released, we have added support for antiSMASH 5 and RGI 5 to the command line version of our tool and are planning to add support to the python notebook version. Initial analysis showed no significant change in classifier accuracy when updating to the newer versions of antiSMASH and RGI. 

Code is developed to run in python 3 (3.8.3). Packages required, version used for publication in parentheses:<br/>
Biopython (1.77)<br/>
Numpy (1.18.5)<br/>
Scikit learn (0.23.1)<br/>
Matplotlib (3.2.2)<br/>

To generate the proper files, antiSMASH 4 should be run with the command:

antismash --verbose --outputfolder OUTPUT_FOLDER --statusfile status.txt --full-hmmer --borderpredict INPUT_FILE

and RGI 3 shoud be run with the command:

rgi -i INPUT_FASTA_FILE -o OUTPUT_FILE --loose_criteria YES --verbose ON


Instructions for python notebook:
Throughout the cells the classification variable determines which binary classification problem predictions are made for. Options are "antibacterial", "antieuk" (which predicts antifungal, antitumor, or cytotoxic activity), "antifungal", "cytotoxic_antitumor", "antigramneg", or "antigrampos.

Cells 1-4 contain methods for processing the training data. These only need to be run if the training set data have been updated.

Cell 5 contains methods used by cells 6-7 and should be run before running those cells.

Cell 6 determines optimal parameters for the different classifiers.

Cell 7 performs cross validation and evaluates classifers by various metrics.

Cell 8 uses the logisitic regression coefficients to analyze which features are most important for the prediction.


Cells 9 and 10 are the only cells that need to be run to make predictions for novel BGCs.

Cell 9 extracts the features from the antiSMASH and rgi output files for the BGC. antiSMASH output genbank files should be saved in the BGCs_for_prediction/antiSMASH folder and rgi output text files should be saved in the BGCs_for_prediction/CARD folder. antiSMASH and RGI files must have matching names so the program knows that they are for the same cluster (e.g. if the antiSMASH genbank file is called cluster1.gbk, the corresponding RGI file must be callsed cluster1.txt).

Predictions can be run either with or without sequence similarity network (sub-PFAM domain) features. Running cell 9 with SSN domains makes the program significantly slower as genes from the cluster need to be blasted against all training set SSN proteins to determine which, if any, sub-PFAM domain the gene belongs to. Running it without SSN features does not significantly decrease the accuracy of the predictions. To run with out SSN features set include_SSN = False in both cell 9 and 10. To run with SSN features set include_SSN = True in both cell 9 and 10 and set the blast_exe variable in cell 9 equal to the path to the blastp executable on your computer.

Cell 10 runs the predictions on the features extracted by cell 9. The classification variable determines which binary prediction will be made, options are "antibacterial", "antieuk" (which predicts antifungal, antitumor, or cytotoxic activity), "antifungal", "cytotoxic_antitumor", "antigramneg", or "antigrampos.


Instructions for command line tool:
requirements.txt contains all required packages for the command line tool.

To run the command line tool enter:

python cluster_function_prediction.py "filename of antismash cluster gbk file" "filename of rgi txt file"
  
optional arguments include:

--output - directory where predictions should be saved, current directory by default

--seed - random seed for predictions

--no_SSN - use this flag to run without SSN features, makes predictions faster with minimal accuracy cost

--blastp_path - path to blastp executable, by default assumes blastp is already in PATH variable

--write_features - directory to write features to, if not defined features will not be saved

--antismash_version antiSMASH version used to generate antiSMASH input file, support for 4 or 5, default 5

--rgi_version RGI version to generate RGI input file, support for 3 or 5, default 5

We have also supplied a Dockerfile that can be used to create a Docker image to run the command line tool. To create the docker image use the command:
docker build -t function_prediction_docker_img .

To run the docker use the command:
docker run --volume EXACT_PATH_TO_DIRECTORY_WITH_INPUTS_AND_OUTPUTS:/src/input function_prediction_docker_image ./input/ANTISMASH_FILENAME ./input/RGI_FILENAME --output ./input/OUTPUT_FOLDERNAME
additional arguments from the command line tool can also be added to this command
