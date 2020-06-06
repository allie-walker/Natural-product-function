# Natural-product-function

Natural_product_BGC_activity_prediction.ipnyb is an python notebook containing code to optimize, train, and use machine learning models to predict natural product
activity from biosynthetic gene clusters. It requires ouptut files from antiSMASH and Resistance Gene Identifier (RGI).

To generate the proper files, antiSMASH should be run with the command:

antismash --verbose --outputfolder OUTPUT_FOLDER --statusfile status.txt --full-hmmer --borderpredict INPUT_FILE

and RGI shoud be run with the command:

rgi -i INPUT_FASTA_FILE -o OUTPUT_FILE --loose_criteria YES --verbose ON