FinalClassifier = function(surv, gep, sets, results, maxUse = 1000){ 
	class.c = results[[1]]
	results.all = results[[2]]
	settings = results[[3]]
	repeats = settings$repeats
	toi = settings$toi
	stratify = settings$stratify
	### base meanHR on repeats from Run 1
	mean.HR <- apply(results.all[1,1:repeats,], 2, mean)
	
	if(any(is.na(apply(results.all[1,,],2,mean)))){
		mean.HR = mean.HR[!is.na(apply(results.all[1,,],2,mean))]
		class.c = class.c[,,!is.na(apply(results.all[1,,],2,mean))] 
		warning(paste("NAs found in results; ", length(which(is.na(apply(results.all[1,,],2,mean)))), " geneset(s) excluded from analysis", sep = ""))
	}	
	surv[which(surv[,3] != toi),3] <- "a"
	surv[which(surv[,3] == toi),3] <- "b"
	
	if(length(sets) < maxUse){
		maxUse = length(sets)
	}	
	
	training <- array(NA, c(repeats, maxUse,4))
	#### determine optimal number of gene sets
	for(x in 1:repeats){
		### select classification of the repeats from Run 2
		class.ab <- class.c[c(settings$folds[[(x+repeats)]][[1]],settings$folds[[(x+repeats)]][[2]]),(x+repeats),]
		class.d <- class.c[settings$folds[[(x+repeats)]][[3]],(x+repeats),]
		
		### select best geneset
		top <- which(mean.HR <= sort(mean.HR, decreasing = FALSE)[1])[1]
		surv.t <- surv[c(settings$folds[[(x+repeats)]][[1]], settings$folds[[(x+repeats)]][[2]]),]
		surv.v <- surv[settings$folds[[(x+repeats)]][[3]],]
		if(stratify == "YES"){
			training[x,1,1] <- summary(coxph(Surv(surv.t[,1],surv.t[,2])~surv.t[,3] + strata(surv.t[,4]) , data = surv.t, subset = class.ab[,top] == 1))$coef[2]
			training[x,1,2] <- summary(coxph(Surv(surv.t[,1],surv.t[,2])~surv.t[,3] + strata(surv.t[,4]) , data = surv.t, subset = class.ab[,top] == 1))$coef[5]
			training[x,1,3] <- summary(coxph(Surv(surv.v[,1],surv.v[,2])~surv.v[,3] + strata(surv.v[,4]) , data = surv.v, subset = class.d[,top] == 1))$coef[2]
			training[x,1,4] <- summary(coxph(Surv(surv.v[,1],surv.v[,2])~surv.v[,3] + strata(surv.v[,4]) , data = surv.v, subset = class.d[,top] == 1))$coef[5]
		}
		if(stratify == "NO"){
			training[x,1,1] <- summary(coxph(Surv(surv.t[,1],surv.t[,2])~surv.t[,3], data = surv.t, subset = class.ab[,top] == 1))$coef[2]
			training[x,1,2] <- summary(coxph(Surv(surv.t[,1],surv.t[,2])~surv.t[,3], data = surv.t, subset = class.ab[,top] == 1))$coef[5]
			training[x,1,3] <- summary(coxph(Surv(surv.v[,1],surv.v[,2])~surv.v[,3], data = surv.v, subset = class.d[,top] == 1))$coef[2]
			training[x,1,4] <- summary(coxph(Surv(surv.v[,1],surv.v[,2])~surv.v[,3], data = surv.v, subset = class.d[,top] == 1))$coef[5]
		}	
		
		### add gene sets sequentially, based on their mean HR
		for(g in 2:maxUse){
			use <- which(mean.HR <= sort(mean.HR, decreasing = FALSE)[g])
    
			use.c <- class.ab[,use]
			use.d <- class.d[,use]
    
			surv.t$score <- apply(use.c, 1, sum)
			surv.v$score <- apply(use.d, 1, sum)
    
			results.e <- matrix(NA, nrow(surv.t), 4)
			count <- 0
			for(i in unique(sort(surv.t$score, decreasing = FALSE))){
				count <- count + 1
				surv.t$solution <- rep(0, nrow(surv.t))
				surv.t$solution[which(surv.t$score >= i)] <- 1
				if(length(unique(surv.t[which(surv.t$solution == 1),3])) == 2 & length(unique(surv.t[which(surv.t$solution == 0),3])) == 2){
					if(stratify == "YES"){
						results.e[count,1] <- summary(coxph(Surv(surv.t[,1],surv.t[,2]) ~surv.t[,3] + strata(surv.t[,4]), data = surv.t, subset = solution == 1))$coef[2]	
						results.e[count,2] <- summary(coxph(Surv(surv.t[,1],surv.t[,2]) ~surv.t[,3] + strata(surv.t[,4]), data = surv.t, subset = solution == 1))$coef[5]
					}
					if(stratify == "NO"){
						results.e[count,1] <- summary(coxph(Surv(surv.t[,1],surv.t[,2]) ~surv.t[,3], data = surv.t, subset = solution == 1))$coef[2]	
						results.e[count,2] <- summary(coxph(Surv(surv.t[,1],surv.t[,2]) ~surv.t[,3], data = surv.t, subset = solution == 1))$coef[5]
					}	
					results.e[count,3] <- mean(surv.t$solution)
					results.e[count,4] <- i
				}	 
			} 
    
    
			cutoff <- results.e[which(results.e[,1] == min(results.e[which(results.e[,3] > settings$min.size),1]))[1], 4 ]
    
    
			surv.t$solution <- rep(0, nrow(surv.t))
			surv.t$solution[which(surv.t$score >= cutoff)] <- 1
			surv.v$solution <- rep(0, nrow(surv.v))
			surv.v$solution[which(surv.v$score >= cutoff)] <- 1
			if(stratify == "YES"){
				training[x,g,1] <- summary(coxph(Surv(surv.t[,1],surv.t[,2])~surv.t[,3] + strata(surv.t[,4]), data = surv.t, subset = solution == 1))$coef[2]
				training[x,g,2] <- summary(coxph(Surv(surv.t[,1],surv.t[,2])~surv.t[,3] + strata(surv.t[,4]), data = surv.t, subset = solution == 1))$coef[5]
				training[x,g,3] <- summary(coxph(Surv(surv.v[,1],surv.v[,2])~surv.v[,3] + strata(surv.v[,4]), data = surv.v, subset = solution == 1))$coef[2]
				training[x,g,4] <- summary(coxph(Surv(surv.v[,1],surv.v[,2])~surv.v[,3] + strata(surv.v[,4]), data = surv.v, subset = solution == 1))$coef[5]
			}
			if(stratify == "NO"){
				training[x,g,1] <- summary(coxph(Surv(surv.t[,1],surv.t[,2])~surv.t[,3], data = surv.t, subset = solution == 1))$coef[2]
				training[x,g,2] <- summary(coxph(Surv(surv.t[,1],surv.t[,2])~surv.t[,3], data = surv.t, subset = solution == 1))$coef[5]
				training[x,g,3] <- summary(coxph(Surv(surv.v[,1],surv.v[,2])~surv.v[,3], data = surv.v, subset = solution == 1))$coef[2]
				training[x,g,4] <- summary(coxph(Surv(surv.v[,1],surv.v[,2])~surv.v[,3], data = surv.v, subset = solution == 1))$coef[5]
			}
		
		}
  
  
	}	

	x = 1:maxUse
	#### calculate median HR over all repeats
	y = apply(training[,,3], 2, median)
	
	
	#### fit regression model
	model <- loess(y ~x) 
	delta = NULL
	for(i in 2:maxUse){ 
		delta[i] <- fitted(model)[(i -1)] - fitted(model)[i]
	}
	
	#### select optimal number of gene sets
	g = (which(delta < 0.0001)[1]) 
	
	if(is.na(g)){
		g = maxUse
		warning("No optimal performance found; maximum available number of genesets used")
	}
	#### base final ranking on the meanHR from Run2
	mean.HR <- apply(results.all[1,(repeats+1):(repeats*2),], 2, mean)
	
	use <- which(mean.HR <= sort(mean.HR, decreasing = FALSE)[g])
	#### train Final Classifier
	finalClassifier = trainFinalClassifiers(gep = gep, surv = surv, sets = sets[use], toi = settings$toi, z = settings$z, g = settings$g, min.size = settings$min.size, stratify = settings$stratify, nbrs = settings$nbrs, seed = settings$seed)
	
	
	use.c <- finalClassifier[[1]]
	
	#### treshold the classification score
	surv$score <- apply(use.c, 1, sum)
	results.e <- matrix(NA, nrow(surv), 4)
	count <- 0
	for(i in unique(sort(surv$score, decreasing = FALSE))){
		count <- count + 1
		surv$solution <- rep(0, nrow(surv))
		surv$solution[which(surv$score >= i)] <- 1
		if(length(unique(surv[which(surv$solution == 1),3])) == 2 & length(unique(surv[which(surv$solution == 0),3])) == 2){
			if(stratify == "YES"){
				results.e[count,1] <- summary(coxph(Surv(surv[,1],surv[,2])~surv[,3] + strata(surv[,4]), data = surv, subset = solution == 1))$coef[2]	
				results.e[count,2] <- summary(coxph(Surv(surv[,1],surv[,2])~surv[,3] + strata(surv[,4]), data = surv, subset = solution == 1))$coef[5]	
			}
			if(stratify == "NO"){
				results.e[count,1] <- summary(coxph(Surv(surv[,1],surv[,2])~surv[,3], data = surv, subset = solution == 1))$coef[2]	
				results.e[count,2] <- summary(coxph(Surv(surv[,1],surv[,2])~surv[,3], data = surv, subset = solution == 1))$coef[5]	
			}
			results.e[count,3] <- mean(surv$solution)
			results.e[count,4] <- i
		} 
	} 
	threshold <- results.e[which(results.e[,1] == min(results.e[which(results.e[,3] > settings$min.size & results.e[,2] < settings$alpha),1]))[1], 4 ]
	expectedPerformance = fitted(model)[g] 
	selectedSets = use 
	finalResults = list(finalClassifier, threshold, selectedSets, expectedPerformance, training) 
	names(finalResults) = c("finalClassifier", "threshold", "selectedSets", "expectedPerformance", "training")
	return(finalResults) 
}


