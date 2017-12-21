# GESTURE
GESTURE (Gene-Expression based Simulated Treatment Using similaRity between patiEnts) implements the concept of Simulated Treatment Learning for gene expression datasets in R. 

# Input files
GESTURE needs a gene expression dataset and survival dataset as input. 
The gene expression dataset should be a numeric matrix with the genes as rows and the samples as columns. The survival dataset is preferably a dataframe, with the samples in the rows. It should have a column with a numeric survival variable, a censoring variable and a treatment variable, in this order. A potential 4th column can contain a stratification variable. 
Colnames of the gene expression dataset need to match with the rownames of the survival dataset.

# Running GESTURE 
Four functions are needed to run GESTURE and predict the classification of independent samples: 
- stlTrain.r
- trainFinalClassifiers.r
- stlFinalClassifier.r
- stlPredict.r

The script runGesture.r provides the code to run the entire algorithm, which relies on the 4 functions above.
