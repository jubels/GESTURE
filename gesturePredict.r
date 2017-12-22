gesturePredict = function(gep.train, gep.validation, sets, finalResults){
	centroids.all = finalResults$finalClassifier[[2]]
	threshold = finalResults$threshold 
	use = finalResults$selectedSets
	
    class.c <- matrix(NA, ncol(gep.validation), length(use))
	rownames(class.c) <- colnames(gep.validation)
	for(l in 1:length(use)){
		dset.val <- gep.validation[sets[[use[l]]], ]
		high <- names(centroids.all[[l]][[1]])
        dset.z <- gep.train[sets[[use[l]]],high]
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
	score = apply(class.c, 1, sum)
	label = rep("no benefit", ncol(gep.validation)) 
	label[which(score >= threshold)] = "benefit" 
	names(label) = colnames(gep.validation)
	return(label)
} 	