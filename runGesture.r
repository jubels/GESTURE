#### source the functions
source("gestureTrain.r") 
source("trainFinalClassifiers.r")
source("gestureFinalClassifier.r") 
source("gesturePredict.r") 

#### load the required libraries
library(caret)
library(survival)
library(lsa)
library(fields)

#### load whichever survivaldata, gene expression and gene sets you want to use

#### gep = 			matrix with samples in columns and observations in rows 
#### surv = 		dataframe or matrix with in columns time, status, treatment and optionally stratification variable (in that order)
#### sets = 		named list of gene sets to be tested with entries matching rownames of gep
#### toi = 			name of treatment of interest, which should have superior survival within class 1 of the classification
#### z = 			maximum of high scoring patients used as prototype for classification
#### g = 			maximum number of neighbours used to determine threshold for gamma 
#### nbrs = 		number of neighbours taken into account when calculating zPFS, default is 10
#### min.size = 	minimum fraction of patiens classified as 1, optional
#### min.HR =		minimum HR in class 1, to avoid extreme HRs, optional
#### folds = 		list of length repeats, containing lists with data partitioned into three folds, does not need to be provided  
#### stratify =		set to  "YES" when stratification variable needs to be used in Cox model 
#### repeats =		number of repeats, STL will use this number of repeats in both Run 1 and Run 2 
#### alpha = 		threshold for p-value associated with HR in class 1 with optimal combination k/gamma

trainingResults = gestureTrain(surv, gep, sets, toi, z, g, nbrs, min.size, min.HR, folds, stratify, repeats,  alpha, seed)

#### In addition to survival data, gene epxression data and the gene sets, the function FinalClassifier takes the output of gestureTrain as input (results). 
#### Make sure to use the same survivaldata, gene expression data and gene sets as used for the gestureTrain function as input. 
#### The maxUse parameter controls how many gene sets can be maximally used in the final ensemble classifier. Default is 1000. 
final = FinalClassifier(surv, gep, sets, results, maxUse) 

#### gesturePredict predicts labels for independent patients. gep.train is the gene expression the model was trained on, gep.validation should contain the gene expression
#### for the samples you want to predict a label for. The input sets should be the same as used in training. 
#### gesturePredict takes the output of FinalClassifier as input (finalResults). It outputs a vector that assigns each sample in the columns of gep.validation the label 'benefit' or 'no benefit'
prediction = gesturePredict(gep.train, gep.validation, sets, finalResults) 


