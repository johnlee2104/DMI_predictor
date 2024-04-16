# DMI_predictor
This repository contains codes to develop, test, and implement a DMI (domain-motif interface) predictor. This program automates the search of domains and motifs in the sequences of interacting proteins. Following that, the program uses a reference list of DMI types to pair potential DMI among the domain and motif matches. These potential DMIs are subsequently scored with a random forest model to allow for the ranking of potential DMI based on their likelihood to be functional.

The predictor requires some files to run:

20220311_elm_classes.tsv: List of motif types that the program uses to find motif matches.

20220311_elm_interaction_domains_complete.tsv: List of DMI types for the pairing of potential DMIs.

all_pfam/smart_domains_with_frequency.txt: Frequencies of different Pfam domains as features for the scoring of potential DMIs.

final_PRS_RRSv4_3_used_to_fit_model.tsv: The positive reference set (PRS) and the third replicate of random reference set (RRSv4) used to fit the random forest and median imputer for deployment.

final_RF_model_with_RRSv4_3.joblib: Fitted random forest model.

final_median_imputer_with_RRSv4_3.joblib: Fitted median imputer.

interpro_9606_pfam/smart_matches_20210122.json: JSON file downloaded from InterPro for all domain matching in human protein sequences in SwissProt.

The scripts to generate different versions of RRS, run the program, analyse features and trained models, as well as precompute different features, are stored in scripts folder.

DMI_prediction: Codes to run the DMI predictor

RRS_formation: Codes to generate different versions of RRS

features_analysis: Codes and notebooks to analyse features among the PRS and RRSs

features_precomputation_annotation: Codes to precompute features needed by the random forest model
model_fitting_evaluation: Codes to automate the training and testing of random forest models using different RRS versions