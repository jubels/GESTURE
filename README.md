# GESTURE
GESTURE (Gene-Expression based Simulated Treatment Using similaRity between patiEnts) implements the concept of Simulated Treatment Learning for gene expression datasets in R. This software may be used for non-commercial purposes; the software license is provided in the license.md file. 

# Input files
GESTURE needs a gene expression dataset and survival dataset as input. It also needs you to define gene sets you want to test.  
The gene expression dataset should be a numeric matrix with the genes as rows and the samples as columns. The survival dataset is preferably a dataframe, with the samples in the rows. It should have a column with a numeric survival variable, a binary censoring variable and a treatment variable, in this order. A potential 4th column can contain a stratification variable. 
Colnames of the gene expression dataset need to match with the rownames of the survival dataset.
The gene sets need to be in a named list. Obviously the names for the genes in the gene set list need to match with the rownames of the gene expression data. 

# Running GESTURE 
Four functions are needed to run GESTURE and predict the classification of independent samples: 
- gestureTrain.r
- trainFinalClassifiers.r
- gestureFinalClassifier.r
- gesturePredict.r

The script runGesture.r provides the code to run the entire algorithm, which relies on the 4 functions above.

Note that the gestureTrain function is relatively slow, with 12 repeats per run, training classifiers for 30 gene sets takes ~1 hour (although this obviously also depends on the size of your dataset). I made this a seperate function to make it easier to run the algorithm parallel (the training of the final classifier and predicting new patients is much faster). If you want to run GESTURE on many gene sets (>100), I recommend seperate runs of gestureTrain with subsets of your gene sets. Then merge the output from the seperate runs and use that as input for gestureFinalClassifier.r. 

# Additional files in the repository 
The files survivalBortezomib, FoldsBortezomib and FoldsLenalidomide are associated with the manuscript about GESTURE. The survivaldata contained within the survivalBortezomib file is associated with GEO datasets GSE2658 and GSE19784. 30 patients from the Total Therapy 3 study used in the manuscript are not included in the GSE2658 dataset, these can be found in ArrayExpress dataset E-TABM-1138. The GEO and ArrayExpress IDs are in the "sampleID" column in the survival data. The files FoldsBortezomib and FoldsLenalidomide contain the assignment of samples to different folds in the crossvalidation contained within the manuscript, to ensure reproducibility of the results. 
