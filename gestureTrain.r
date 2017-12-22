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

#### note: function relies on R packages 'lsa', 'fields' and 'survival'. If folds are not provided R package 'caret' is required.

gestureTrain <- function(gep, surv, sets, toi, z = 50, g = 30, nbrs = 10, min.size = 0, min.HR = 0, folds = NA, stratify = "NO", repeats = 12, alpha = 0.05, seed)
{
	set.seed(seed)
	### recode treatment variable to ensure right comparison in Cox model 
	surv[which(surv[,3] != toi),3] <- "a"
	surv[which(surv[,3] == toi),3] <- "b"
	centroids.all <- vector("list", length(sets))
	repeats.use = repeats*2
	div <- c(1:2, 1:2)
	if(is.na(folds)){
		#### making folds for cross validation
		folds <- list()
		for (i in 1:repeats.use){
			HR.1 <- 1               
			HR.2 <- 2
			HR.3 <- 3
			while(abs(HR.1 - HR.2) > 0.05 | abs(HR.2 - HR.3) > 0.05 | abs(HR.3 - HR.1) > 0.05){
				folds.cv <- createFolds(surv[,3], k = 3)
				survival.1 <- surv[folds.cv$Fold1,]
				survival.2 <- surv[folds.cv$Fold2,]
				survival.3 <- surv[folds.cv$Fold3,]
				if(stratify == "YES"){
					HR.1 <- summary(coxph(Surv(survival.1[,1], survival.1[,2])~survival.1[,3] + strata(survival.1[,4])))$coef[2]
					HR.2 <- summary(coxph(Surv(survival.2[,1], survival.2[,2])~survival.2[,3] + strata(survival.2[,4])))$coef[2]
					HR.3 <- summary(coxph(Surv(survival.3[,1], survival.3[,2])~survival.3[,3] + strata(survival.3[,4])))$coef[2]
				}
				if(stratify == "NO"){
					HR.1 <- summary(coxph(Surv(survival.1[,1], survival.1[,2])~survival.1[,3]))$coef[2]
					HR.2 <- summary(coxph(Surv(survival.2[,1], survival.2[,2])~survival.2[,3]))$coef[2]
					HR.3 <- summary(coxph(Surv(survival.3[,1], survival.3[,2])~survival.3[,3]))$coef[2]
				}			
			}
		folds[[i]] <- folds.cv
		}
	}
	for(l in 1:length(sets)){
		centroids.path <- list()
		if(length(sets) > 0){	
			for(f in 1:repeats.use){
				#### dividing data  for cross validation 
				dset.z <- gep[sets[[l]],folds[[f]][[1]]]	
				dset.g <- gep[sets[[l]],folds[[f]][[2]]]
				dset.v <- gep[sets[[l]],folds[[f]][[3]]]
				surv.z <- surv[folds[[f]][[1]],]	
				surv.g <- surv[folds[[f]][[2]],]	
				surv.v <- surv[folds[[f]][[3]],]
				
				### calculating euclidean distance matrix for fold A
				model <- rdist(t(dset.z))
				
				rownames(model) <- rownames(surv.z)
				colnames(model) <- rownames(surv.z) 
				
				
				treatment <- rownames(surv.z)[which(surv.z[,3] == "b")]
				other <- rownames(surv.z)[which(surv.z[,3] == "a")]

				distanceA <- model[other,treatment] 
				
				### calculating deltaPFS
				surv.diff <- matrix(NA, length(other), length(treatment))
				for (i in 1:length(treatment)){
					for (j in 1:length(other)){
						if(surv.z[treatment[i],2] == 1 & surv.z[other[j],1] > surv.z[treatment[i],1]){
							surv.diff[j,i] <- surv.z[treatment[i],1] - surv.z[other[j],1] 
						} 
						if(surv.z[other[j],2] == 1 & surv.z[other[j],1] < surv.z[treatment[i],1]){
							surv.diff[j,i] <- surv.z[treatment[i],1] - surv.z[other[j],1] 
						}
					}
				}

				rownames(surv.diff) <- other
				colnames(surv.diff) <- treatment

				true.nbr <- matrix(NA, nbrs, length(treatment)) 
				true.dpfs <- rep(NA, length(treatment))
				rand.nbr <- array(NA, c(nbrs, length(treatment), 1000)) 
				rand.dpfs <- matrix(NA, 1000, length(treatment))
				zpfs <- rep(NA, length(treatment))
		
				colnames(true.nbr) <- treatment
				colnames(rand.dpfs) <- treatment
				names(true.dpfs) <- treatment
				names(zpfs) <- treatment
				
				for(i in 1:length(treatment)){
					true.nbr[,i] <- rownames(distanceA)[sort(distanceA[,i], decreasing = FALSE, index = TRUE)$ix[1:nbrs]]
					true.dpfs[i] <- mean(surv.diff[true.nbr[,i],i], na.rm = TRUE)
				} 
					
				#### calculating RPFS
				for (r in 1:1000){
					for (i in 1:length(treatment)){
						rand.nbr[,i,r] <- other[sample(1:length(other), nbrs, replace = FALSE)] 
						rand.dpfs[r,i] <- mean(surv.diff[rand.nbr[,i,r],i], na.rm = TRUE)
					}
				} 
					
				#### calculating zPFS
				for (i in 1:length(treatment)){
					zpfs[i] <- (true.dpfs[i] - mean(rand.dpfs[,i], na.rm = TRUE)) / sd(rand.dpfs[,i], na.rm = TRUE) 
				}
				
				p.values <- matrix(NA, z, g)
				HRs <- matrix(NA, z, g)
				size <- matrix(NA, z, g)
				
				### calculating distance matrix for fold B
				distanceB <- rdist(t(cbind(dset.g, dset.z)))
				
				
				
				rownames(distanceB) <- c(colnames(dset.g), colnames(dset.z))
				colnames(distanceB) <- c(colnames(dset.g), colnames(dset.z))
				
				### optimizing k and gamma 
				for(b in 1:z){
					high <- treatment[which(true.dpfs > 0 & zpfs >= sort(zpfs[which(true.dpfs > 0)], decreasing = TRUE)[b])]
					if(length(high) > 0){
						use.dist <- distanceB[high, colnames(dset.g)]
						if(length(high) == 1){
							use.dist <- t(as.matrix(use.dist))
						}
						classification <- matrix(0, z, ncol(dset.g)) 
						high.g <- NULL
						for (o in 1:g){
							if(length(high.g)/ncol(dset.g) < 0.9){
								limits <- apply(use.dist,1, sort)[o,]
								for(d in 1:length(high)){
									classification[d,which(use.dist[d,] <= limits[d])] <- 1
								} 					
								high.g <- which(apply(classification, 2, sum) > 0)
								if(length(unique(surv.g[high.g,3])) == 2){
									size[b,o] <- length(high.g)/nrow(surv.g)
									surv.high.g <- surv.g[high.g,]
									if(stratify == "YES"){
										cox.model <- summary(coxph(Surv(surv.high.g[,1],surv.high.g[,2])~surv.high.g[,3] + strata(surv.high.g[,4])))
										p.values[b,o] <- cox.model$coef[5] 
										HRs[b,o] <- cox.model$coef[2] 
									}
									if(stratify == "NO"){
										cox.model <- summary(coxph(Surv(surv.high.g[,1],surv.high.g[,2])~surv.high.g[,3]))
										p.values[b,o] <- cox.model$coef[5] 
										HRs[b,o] <- cox.model$coef[2] 
									}	
								}	
							}	
						}
					}	
				}
				
				### select combinations k and gamma that do not violate size/HR constraints 
				optimal.options <- HRs[which(p.values < alpha & size > min.size & HRs > min.HR)]
				if(length(optimal.options) == 0){
					optimal.options <- HRs[which(size > min.size & HRs > min.HR)]
				}	
					
				### validate optimal k and gamma in fold C		
				if(length(optimal.options) > 0){
					optimal <- which(HRs == min(optimal.options, na.rm = T), arr.ind = T)
					high <- treatment[which(true.dpfs > 0 & zpfs >= sort(zpfs[which(true.dpfs > 0)], decreasing = TRUE)[optimal[1,1]])]
					distanceC <- rdist(t(cbind(dset.v, dset.z[,high])))	
					rownames(distanceC) <- c(colnames(dset.v), high)
					colnames(distanceC) <- c(colnames(dset.v), high)
					use.dist <- distanceB[high, colnames(dset.g)]
					use.dist.v <- distanceC[high, colnames(dset.v)]
					if(length(high) == 1){
						use.dist <- t(as.matrix(use.dist))
						use.dist.v <- t(as.matrix(use.dist.v))
					}	
					optimal.limits <- apply(use.dist,1, sort)[optimal[1,2],]
					classification.v <- matrix(0,length(high), nrow(surv.v))
					for(b in 1:nrow(use.dist.v)){
						classification.v[b,which(use.dist.v[b,] <= optimal.limits[b])] <- 1
					} 	
					names(optimal.limits) <- high 
					centroids.path[[f]] <- optimal.limits
					surv.v$solution <- rep(0, nrow(surv.v))
					surv.v$solution[which(apply(classification.v, 2, sum) > 0)] <- 1 
					if(length(unique(surv.v[which(surv.v$solution == 1),3])) == 2 & length(unique(surv.v[which(surv.v$solution == 0),3])) == 2){ 
						if(stratify == "YES"){
							cox.1 <- summary(coxph(Surv(surv.v[,1], surv.v[,2])~surv.v[,3] + strata(surv.v[,4]), subset = surv.v$solution == 1))
							cox.0 <- summary(coxph(Surv(surv.v[,1], surv.v[,2])~surv.v[,3] + strata(surv.v[,4]), subset = surv.v$solution == 0))
							attr(centroids.path[[f]], "perf") <- c(cox.1$coef[2], cox.1$coef[5], mean(surv.v$solution))
						}
				
						if(stratify == "NO"){
							cox.1 <- summary(coxph(Surv(surv.v[,1], surv.v[,2])~surv.v[,3], subset = surv.v$solution == 1))
							cox.0 <- summary(coxph(Surv(surv.v[,1], surv.v[,2])~surv.v[,3], subset = surv.v$solution == 0))
							attr(centroids.path[[f]], "perf") <- c(cox.1$coef[2], cox.1$coef[5], mean(surv.v$solution))
						}
					}
			
				}	
				
			}
				 
		}	
			
	
		centroids.all[[l]] <- centroids.path
	}
	
	### gather results of all classifiers in an array
	results.all = array(NA, c(3,repeats.use, length(sets)))
	for(l in 1:length(centroids.all)){
		for(f in 1:repeats.use){	
			if(length(attr(centroids.all[[l]][[f]], "perf")) == 3){
				results.all[1,f,l] <- attr(centroids.all[[l]][[f]], "perf")[1]
				results.all[2,f,l] <- attr(centroids.all[[l]][[f]], "perf")[2]
				results.all[3,f,l] <- attr(centroids.all[[l]][[f]], "perf")[3]
 
			}	
				
		}	
	}
	
	settings <- list(z = z, g = g, nbrs = nbrs, min.size = min.size, min.HR = min.HR, folds = folds, toi = toi, stratify = stratify, repeats = repeats, alpha = alpha, seed = seed)
	
	### classify all patients with all classifiers
	class.c <- array(NA, c(ncol(gep), dim(results.all)[2],length(sets)))
	rownames(class.c) <- colnames(gep)
	for(f in 1:repeats.use){
        for(l in 1:length(sets)){
                dset.val <- gep[sets[[l]], ]
                if(!is.na(results.all[1,f,l])){
                        high <- names(centroids.all[[l]][[f]])
                        dset.z <- gep[sets[[l]],high]
                        if(length(high) == 1){
                                dset.z <- as.matrix(dset.z)
                        }
                        optimal.limits <- as.vector(centroids.all[[l]][[f]])
                        cos.dist.v <- rdist(t(cbind(dset.val, dset.z)))
                        rownames(cos.dist.v) <- c(colnames(dset.val), high)
                        colnames(cos.dist.v) <- c(colnames(dset.val), high)
                        use.dist.v <- cos.dist.v[high, colnames(dset.val)]
                        if(length(high) == 1){
                                use.dist.v <- t(as.matrix(use.dist.v))
                        }

                        classification.v <- matrix(0,length(high), ncol(dset.val))
                        for(b in 1:nrow(use.dist.v)){
                                classification.v[b,which(use.dist.v[b,] <= optimal.limits[b])] <- 1
                        }
                        class.c[which(apply(classification.v, 2, sum) > 0),f,l] <- 1
                        class.c[which(apply(classification.v, 2, sum) == 0),f,l] <- 0
                }

        }
	}

	results = list(class.c, results.all, settings)
	names(results) = c("classification", "performance", "settings")
	return(results)	
}

	
