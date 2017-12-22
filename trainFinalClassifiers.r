trainFinalClassifiers <- function(gep, surv, sets, toi, z = 50, g = 30, nbrs = 10, min.size = 0, min.HR = 0, folds = NA, stratify = "NO",  alpha = 0.05, seed)
{
	set.seed(seed)
	
	centroids.all <- vector("list", length(sets))
	repeats = 1
	div <- c(1:2, 1:2)
	
	if(is.na(folds)){
		#### making folds for cross validation
		folds <- list()
		for (i in 1:repeats){
			HR.1 <- 1               
			HR.2 <- 2
		
			while(abs(HR.1 - HR.2) > 0.05 ){
				folds.cv <- createFolds(surv[,3], k = 2)
				survival.1 <- surv[folds.cv$Fold1,]
				survival.2 <- surv[folds.cv$Fold2,]
			
				if(stratify == "YES"){
					HR.1 <- summary(coxph(Surv(survival.1[,1], survival.1[,2])~survival.1[,3] + strata(survival.1[,4])))$coef[2]
					HR.2 <- summary(coxph(Surv(survival.2[,1], survival.2[,2])~survival.2[,3] + strata(survival.2[,4])))$coef[2]
				}
				if(stratify == "NO"){
					HR.1 <- summary(coxph(Surv(survival.1[,1], survival.1[,2])~survival.1[,3]))$coef[2]
					HR.2 <- summary(coxph(Surv(survival.2[,1], survival.2[,2])~survival.2[,3]))$coef[2]
				}			
			}
		folds[[i]] <- folds.cv
		}
	}
	for(l in 1:length(sets)){
		centroids.path <- list()
		if(length(sets) > 0){	
			for(f in 1:repeats){
				#### dividing data  for cross validation 
				dset.z <- gep[sets[[l]],folds[[f]][[1]]]	
				dset.g <- gep[sets[[l]],folds[[f]][[2]]]
			
				surv.z <- surv[folds[[f]][[1]],]	
				surv.g <- surv[folds[[f]][[2]],]	
		
				
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
				
				
				distanceB <- rdist(t(cbind(dset.g, dset.z)))
				
				rownames(distanceB) <- c(colnames(dset.g), colnames(dset.z))
				colnames(distanceB) <- c(colnames(dset.g), colnames(dset.z))
				
				### optimzing z and gamma 
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
				
				optimal.options <- HRs[which(p.values < alpha & size > min.size & HRs > min.HR)]
				if(length(optimal.options) == 0){
					optimal.options <- HRs[which(size > min.size & HRs > min.HR)]
				}	
					
				### validate k and gamma 		
				if(length(optimal.options) > 0){
					optimal <- which(HRs == min(optimal.options, na.rm = T), arr.ind = T)
				
					high <- treatment[which(true.dpfs > 0 & zpfs >= sort(zpfs[which(true.dpfs > 0)], decreasing = TRUE)[optimal[1,1]])]
			
					use.dist <- distanceB[high, colnames(dset.g)]
				
					if(length(high) == 1){
						use.dist <- t(as.matrix(use.dist))
					
					}	
					optimal.limits <- apply(use.dist,1, sort)[optimal[1,2],]
					
					names(optimal.limits) <- high 
					centroids.path[[f]] <- optimal.limits
				
			
				}	
				
			}
				 
		}	
			
	
		centroids.all[[l]] <- centroids.path
	}
	
	### classify all patients according to the trained classifier
	class.c <- matrix(NA, ncol(gep), length(sets))
	rownames(class.c) <- colnames(gep)
    for(l in 1:length(sets)){
		dset.val <- gep[sets[[l]], ]
        high <- names(centroids.all[[l]][[1]])
        dset.z <- gep[sets[[l]],high]
        if(length(high) == 1){
            dset.z <- as.matrix(dset.z)
        }
        optimal.limits <- as.vector(centroids.all[[l]][[1]])
        dist.v <- rdist(t(cbind(dset.val, dset.z)))
        rownames(dist.v) <- c(colnames(dset.val), high)
        colnames(dist.v) <- c(colnames(dset.val), high)
        use.dist.v <- dist.v[high, colnames(dset.val)]
        if(length(high) == 1){
			use.dist.v <- t(as.matrix(use.dist.v))
        }

        classification.v <- matrix(0,length(high), ncol(dset.val))
        for(b in 1:nrow(use.dist.v)){
            classification.v[b,which(use.dist.v[b,] <= optimal.limits[b])] <- 1
        }
        class.c[which(apply(classification.v, 2, sum) > 0),l] <- 1
        class.c[which(apply(classification.v, 2, sum) == 0),l] <- 0
    }

        
	
	settings <- list(z = z, g = g, nbrs = nbrs, min.size = min.size, min.HR = min.HR, folds = folds, stratify = stratify, repeats = repeats,  alpha = alpha, seed = seed)
	results = list(class.c, centroids.all, settings)
	return(results)	
}

	
