
#classifierTraining
#performanceEstimation
#network...

geNetClassifier <- function(eset, sampleLabels, plotsName=NULL, buildClassifier=TRUE, estimateGError=FALSE, calculateNetwork=TRUE,   labelsOrder=NULL, geneLabels=NULL, numGenesNetworkPlot = 100, minGenesTrain=1, maxGenesTrain  = 100, continueZeroError=FALSE, numIters = 6, lpThreshold = 0.95, numDecimals=3, removeCorrelations=FALSE, correlationsThreshold=0.8, removeInteractions=FALSE, interactionsThreshold=0.5, minProbAssignCoeff=1, minDiffAssignCoeff=0.8,    IQRfilterPercentage = 0, skipInteractions=TRUE,    precalcGenesNetwork = NULL, precalcGenesRanking=NULL, returnAllGenesRanking=TRUE, verbose=TRUE)
# numDecimals=3 (solo necesario en GE)
# Precalculated genesNetwork / genesRanking -> Use at your own risk!
{    

    #### Fixed variables  (could be added as arguments in the future)
    geneAdd <- "V1.b"                                                                # "V0" "V1" o "V1.b"
    geneSelection<-"maxSoutliers"                                        # "max" "mean" "maxSoutliers"
    outlThreshold <- 1                                                                # "SD" q consideramos outliers. (1=64%) x defecto era 1.5 (86%). 2 (94%) es casi igual q "max"            (0=median?)
    returnTopGenesNetwork <-calculateNetwork                # Si se calcula, se devuelve

    ###################################
    ###     Checking arguments
    ###################################
    # sampleLabels format: 
    if(is.character(sampleLabels) && (length(sampleLabels) ==1))
    {
        if( is(eset, "ExpressionSet")  && (sampleLabels %in% colnames(pData(eset)))) {
        sampleLabels <- pData(eset)[,sampleLabels, drop=FALSE]
        }    else{
            stop("The sampleLabels should be either a factor, or contain the name of the phenoData column containing the labels.") 
        }
    }
    if(is.data.frame(sampleLabels) || is.matrix(sampleLabels))
    {
        if(dim(sampleLabels)[2] == 1)  
        {
            tempSamplesLabels <- sampleLabels
            sampleLabels <- as.factor(sampleLabels[,1])
            names(sampleLabels) <- rownames(tempSamplesLabels)
        }
    }
    if(class(sampleLabels) != "factor") { warning("The argument 'classification sampleLabels' had to be converted into a factor.")}    
    sampleLabels <- factor(sampleLabels) #Just in case there are not samples of all the original labels    
    
    # Eset should be a matrix. If it is a whole ExpressionSet, we extract the expressions matrix with exprs()
    if(is(eset, "ExpressionSet")) eset <- exprs(eset) else if (!is.matrix(eset)) stop("The first argument should be an expression matrix or an ExpressionSet.")
    if(length(unique(rownames(eset))) > length(rownames(eset))) stop("The row names of the expression matrix (gene ID) should be unique.")
    
    # Check wether the eset and the sampleLabels match
    numSamples <- ncol(eset)
    if(numSamples != length(sampleLabels)) stop("The number of labels does not match the number of samples.")
    if(!is.null(names(sampleLabels)))
    {
        if(sum(!names(sampleLabels) %in% colnames(eset))>0 ) stop("The names of the labels do not match the samples.")
    }else 
    {
        names(sampleLabels)<-colnames(eset)
        warning("The data labels vector is not named, it will be assumed the labels are in order: the first label applies to the first sample... ", immediate. = TRUE)
    }    
    if(!is.null(labelsOrder))
    {
        if(!is.vector(labelsOrder) && !is.factor(labelsOrder)) stop("The labels order should be either a vector or a factor")
        if(any(!labelsOrder %in% levels(sampleLabels)) || any(!levels(sampleLabels) %in% labelsOrder)) 
        {
            warning("The labelsOrder doesn't match the samples labels. It will be ignored.")
            labelsOrder <- NULL
        }
    }
    
    # Calculate basic variables
    classes <- levels(sampleLabels)
    if(!is.null(labelsOrder)) classes <- labelsOrder
    numClasses <- length(classes)
       if(numClasses < 2) stop("The classification can not be done with less than two classes.(There is nothing to classify)")
       if(maxGenesTrain  < numClasses) stop(paste("The classification can not be performed with less genes than classes. Increase the argument 'maxGenesTrain'.",sep=""))
    numElemClass <- table(sampleLabels)
    

    # Check the rest of the arguments
    # ShowNames checked in function (matrix, dataframe or vector)
    if(!is.numeric(numDecimals)){ numDecimals <- 3
    }else {
        if (numDecimals<0 || numDecimals>7) 
        {
        numDecimals <- 3
        warning("The argument 'numDecimals' should be a number between 0 and 7. The default value (3) will be used.", immediate. = TRUE)
        }
    }
    if(!is.numeric(IQRfilterPercentage) || (IQRfilterPercentage>=1 || IQRfilterPercentage <0)) stop("The filter percentage should be a probability (a number between 0 and 1).")
    if(!is.numeric(minGenesTrain) || (minGenesTrain<1 || minGenesTrain>=maxGenesTrain)) stop("The minimum number of genes per class to train the classifier with should be a number higher than zero and lower than the maxGenesTrain.")
    if(!is.numeric(maxGenesTrain) || maxGenesTrain<3) stop("The maximum number of genes per class to train the classifier with should be a number higher than two.")
    if(!is.numeric(numIters) || numIters<3) stop("The process for calculating the number of genes should be repeated at least three times (recommended: 10).")
    if(!is.numeric(lpThreshold) || (lpThreshold>=1 || lpThreshold <0)) stop("Lp threshold should be a probability (a number between 0 and 1).")
    if(!is.null(precalcGenesNetwork) && !is.list(precalcGenesNetwork)){ stop("The argument genesNetwork should be a list.")
                                                                    }else{ genesNetwork <- precalcGenesNetwork}
    if(!is.logical(removeInteractions) || !is.logical(removeCorrelations)) stop("The arguments removeInteractions and removeCorrelations should be either TRUE or FALSE.")
    if(!is.logical(skipInteractions) || (skipInteractions && removeInteractions)) skipInteractions <- FALSE
    if(!is.logical(buildClassifier))  stop("The argument buildClassifier should be either TRUE or FALSE.")
    if(!is.logical(estimateGError))  stop("The argument estimateGError should be either TRUE or FALSE.")
    if(!is.logical(calculateNetwork))  stop("The argument calculateNetwork should be either TRUE or FALSE.")
    if(!is.logical(returnTopGenesNetwork)) returnTopGenesNetwork <- FALSE
    if(skipInteractions && calculateNetwork) warning("The network will be calculated, but without interactions", immediate. = TRUE)
    if((removeInteractions || removeCorrelations) && (!calculateNetwork && is.null(genesNetwork))) 
    {
        warning("In order to remove correlations or interactions, the genes network needs to be calculated.")
        calculateNetwork<-TRUE
    }
    if(calculateNetwork && !is.null(genesNetwork)) 
    {
        warning("The genes network was provided. It will not be recalculated.")
        calculateNetwork <- FALSE
    }
    if((!buildClassifier && !estimateGError) && (calculateNetwork && !returnTopGenesNetwork)) 
    {
        warning("The genes network will be calculated but not used (buildClassifier and estimateGError are False) so it will be returned although it was not requested.", immediate. = TRUE)
        returnTopGenesNetwork <- TRUE # Otherwise only the genes ranking will be returned (the network will be lost & not used)
    }
    if ((!calculateNetwork&& is.null(genesNetwork)) && returnTopGenesNetwork) stop("The genes network cannot be returned without providing or calculating it.")
    if(!is.numeric(numGenesNetworkPlot) || numGenesNetworkPlot< 2) stop("The number of genes to plot in the network should be a number higher than two.")
    if(!is.numeric(interactionsThreshold) || (interactionsThreshold>=1 || interactionsThreshold <=0)) stop("The threshold for the interactions should be a probability (a number between 0 and 1).")
    if(!is.numeric(correlationsThreshold) || (correlationsThreshold>=1 || correlationsThreshold <=0)) stop("The threshold for the correlations should be a probability (a number between 0 and 1).")
    if(!is.numeric(minProbAssignCoeff)) stop("'minProbAssignCoeff' should be a coefficient to modify the minimum probability required to assign the sample to a class.")
    if(!is.numeric(minDiffAssignCoeff)) stop("'minDiffAssignCoeff' should be a coefficient to modify the required difference between probabilites to assign the sample to a class.")
        if(minProbAssignCoeff<0 || ((numClasses != 2) &&(minProbAssignCoeff>(numClasses/2)))) stop("'minProbAssignCoeff' should be between 0 and half of the number of classes.")
        if(minDiffAssignCoeff<0 || minDiffAssignCoeff>numClasses) stop("'minDiffAssignCoeff' should be between 0 and the number of classes.")
    if(!is.null(plotsName) && !is.character(plotsName))  stop("The plots file name is not a valid name.")
    if(!is.logical(verbose)) verbose <- TRUE
    if(!is.logical(returnAllGenesRanking)) returnAllGenesRanking <- FALSE
    
    # Can the genes ranking provided belong to this dataset?
    genesRankingGlobal <- precalcGenesRanking
    if(!is(genesRankingGlobal, "GenesRanking"))
    {
        if (!is.null(genesRankingGlobal)) warning("The genesRanking provided is not valid. It will be re-calculated.")
        genesRankingGlobal <- NULL        
    }else
    {
        if (!any(gClasses(genesRankingGlobal) %in% classes)) # The genes ranking and the classes given sould be in the same order?     
        {
            if  (numClasses!=2 && length(gClasses(genesRankingGlobal)) != 1)        
            {
                genesRankingGlobal <- NULL
                warning("The genesRanking given as argument does not match the expression matrix. It will be re-calculated.")
            }
        }
    }
    
    #################################################
    ###     Calculate variables, ranking and network
    #################################################
    
    # Check if there are the same number of samples for all the classes
    sameNumSamplesClass <- length(table(numElemClass))==0
     if (!sameNumSamplesClass) warning("It is recommended to have the same number of samples in each class in order to obtain balanced external validation stats.")  
    if(estimateGError && sameNumSamplesClass) if((unique(numElemClass)%%5)!=0) warning("Since the number of samples is not multiple of 5, some samples might be used as test in several cross-validation loops when estimating the generalization error of the classifier.")
    
    
    # Loops number for crossvalidation. Intern: classifier training, extern: stats (1 = only training, 5 = stats, 5+1=both)
    numCV <- cvNumbers(numElemClass)
    if (estimateGError && (numCV[2] ==1)) stop("There are not enough samples of each class to calculate the generalization error of classifiers trained with this data.\n The error of each specific classifier built may be very different, so you should check the errors for each of them individually.")
    numCV.intern <- numCV[1]
    numCV.extern <- ifelse(estimateGError,numCV[2],0)
    ifelse (buildClassifier, numCV.total <- numCV.extern + 1, numCV.total <- numCV.extern) # If buildClassifier and Calculate statistics : +1 loop
    
    # Filter data, calculate Posterior Probabilities and rank genes
    if (verbose){ message(paste(format(Sys.time(), "%H:%M:%S"),"- Filtering data and calculating the genes ranking...")); flush.console()}
    esetFiltered <- eset[iqr.filter(eset,IQRfilterPercentage),]
    if(dim(esetFiltered)[1]< numClasses) stop(paste("Applying a filter percentage of ",IQRfilterPercentage," there are not enough genes left to perform the classiffication. Try with a lower filter percentage.",sep=""))
    if(!is.null(geneLabels)) geneLabels <- extractGeneLabels(geneLabels, rownames(esetFiltered))

    # Make sure the columns in filtered eset are ordered by sample type
    indexes<-NULL
    for(i in 1:numClasses) 
    { 
        indexes <- c(indexes, which(sampleLabels==classes[i]))
    }
    esetFiltered <- esetFiltered[, indexes, drop=FALSE]
    sampleLabels <- sampleLabels[indexes] 
    
    
    # Check whether the genesRanking provided matches the filtered eset
    if(!is.null(genesRankingGlobal) && (sum(!rownames(genesRankingGlobal@postProb) %in% rownames(eset)) > 0))     # There are genes in the genesRanking which arent in the eset -> Recalculate
    {
        genesRankingGlobal <- NULL
        warning("The genesRanking given as argument doesn't match the dataset. Recaculculating the genesRanking...", immediate. = TRUE)
    }
    if (!is.null(genesRankingGlobal) &&    any(!rownames(esetFiltered) %in% rownames(genesRankingGlobal@postProb)) )    # There are missing genes from the eset  in the genesranking -> Warning (Most likely ok, just filtered)
    {
        warning("The genesRanking does not contain all the genes in the expression matrix.", immediate. = TRUE)
    }
    # Calculate genesRanking (global) if it was not provided
    if(is.null(genesRankingGlobal))  genesRankingGlobal <- PEB(esetFiltered, sampleLabels, labelsOrder= labelsOrder)
       if(!is.null(geneLabels)) genesRankingGlobal <- setProperties(genesRankingGlobal, geneLabels=geneLabels)
    
    # Default lpThreshold= 0.95. 
    lp <-     numSignificantGenes(genesRankingGlobal, lpThreshold=lpThreshold)
    lpMaxGenes <- lp
    if (estimateGError)
    {
        lpMaxGenes <- sapply(lpMaxGenes, function(x) max(x,4*maxGenesTrain)) # To asure there are enough in the internal CV loops
    } 
    maxGenesTrain <- rep(maxGenesTrain, length(lp))
    names(maxGenesTrain) <- names(lp)
    
    # To assure that there is at least the requested number of genes in the network
    lpMaxGenes <- sapply(lpMaxGenes, function(x) max(x, numGenesNetworkPlot))        
    #To asure there are enough genes in the global genes ranking (4xmaxGenesTrain) 
    lpMaxGenes <- apply(rbind(lpMaxGenes, genesRankingGlobal@numGenesClass), 2, function(x) min(x))
    
    # Add expression difference to the genes ranking
    topGenes <- as.vector(getRanking(getTopRanking(genesRankingGlobal, lpMaxGenes), showGeneID=TRUE, showGeneLabels=FALSE)$geneID)
    topGenes <- topGenes[which(!is.na(topGenes))]
    meanExprDiff <- difMean(esetFiltered[topGenes,], sampleLabels)
    genesRankingGlobal <- setProperties(genesRankingGlobal, meanDif=meanExprDiff)
    
    # Calculate genes with Correlations & Interactions
    rankENSG <- getRanking(genesRankingGlobal, showGeneLabels=FALSE, showGeneID=TRUE)$geneID
    if (is.null(genesNetwork) && calculateNetwork)                                
    {    
        # Calculate Correlations
        if (verbose) { message(paste(format(Sys.time(), "%H:%M:%S"),"- Calculating correlations between genes...")) ; flush.console()}
        for( cl in gClasses(genesRankingGlobal))
        {        
            if(lpMaxGenes[cl] > 0)    
            {
               nodes <- rankENSG[1:lpMaxGenes[cl],cl]
               edges <- correlation.net(esetFiltered, nodes, lpMaxGenes[cl], method="pearson", threshold=correlationsThreshold)
               genesNetwork <- c(genesNetwork, genesNetwork=list(new("GenesNetwork", nodes = nodes,edges = edges)))
            }else {genesNetwork <- c(genesNetwork, genesNetwork=list(NULL))}
            names(genesNetwork)[length(genesNetwork)] <- cl
        }
    
        # Calculate Interations (mutual information)
        if (!skipInteractions || removeInteractions)
        {
            if (verbose) { message(paste(format(Sys.time(), "%H:%M:%S"),"- Calculating interactions between genes...")) ; flush.console()}
            for( cl in 1:length(gClasses(genesRankingGlobal)))
            {
                if(lpMaxGenes[cl] >0)   genesNetwork[[cl]]@edges <- rbind(genesNetwork[[cl]]@edges, interaction.net(esetFiltered, rankENSG[1:lpMaxGenes[cl],cl], lpMaxGenes[cl], method="clr", estimator="mi.empirical", threshold=interactionsThreshold))
            }
        }
    }

    # If they have to be removed...
    genesRedundancy <- NULL
    if(removeCorrelations && removeInteractions)              
    {
        genesRedundancy <- remove.redundancy(esetFiltered, rankENSG, lpMaxGenes, genesNetwork) #Ambos        
    }
    if(removeCorrelations && !removeInteractions)
    {
        genesRedundancy <- remove.redundancy(esetFiltered, rankENSG, lpMaxGenes, genesNetwork, relation="Correlation - pearson") # Correlation
    }
    if(!removeCorrelations && removeInteractions)
    {
        genesRedundancy <- remove.redundancy(esetFiltered, rankENSG, lpMaxGenes, genesNetwork, relation="Interaction - clr")     #Interaction
    }
    numNonRedundantGenes <- sapply(genesRedundancy$nonRedundantGenes, function(x) length(x))

    if(!is.null(genesRedundancy)  && (buildClassifier||estimateGError))  
    {
        availableGenes <-  sapply(genesRedundancy$nonRedundantGenes, function(x) length(x))
        
        lessMaxGenes <- (maxGenesTrain > availableGenes)
        if(any(lessMaxGenes))
        {
            warning(paste("There are not ",maxGenesTrain[1]," non-redundant genes for some of the classes. These available genes will be used instead. If they are not enough, try lowering 'lpThreshold':",sep=""))
            maxGenesTrain[which(lessMaxGenes)] <- availableGenes[which(lessMaxGenes)]
            print(maxGenesTrain[which(lessMaxGenes)])
        }
    }
    
    if(estimateGError)
    {    
        # Generalization Error stats
         cvSamples <- cvSplitSamples(numCV.extern, sampleLabels)                        
        mxcf <- matrix(data=0,  nrow=numClasses, ncol=numClasses+1) 
        globalError <- vector(length=numCV.extern) # Almacenara el error para cada iteracion del bucle externo de crossValidation
        globalResults <- NULL        
        genesUsed <- NULL 
    }
    
    ###################################################
    ###     Cross validation and classifier training
    ###################################################
    trainGenes      <- NULL #Needed in case (!buildClassifier && !estimateGError)
    numBestGenes    <- NULL
    maxGenesTrainGlobal <- maxGenesTrain
    if(buildClassifier || estimateGError)
    {
        # Inizialization
        numTrainGenes <- matrix(nrow=numCV.total,ncol=ncol(genesRankingGlobal@ord))                
                            if(estimateGError) rownames(numTrainGenes)<-  c(paste("CV ", 1:numCV.total,":", sep=""))
                            if(buildClassifier) rownames(numTrainGenes)[nrow(numTrainGenes)] <- "Classifier"
                            colnames(numTrainGenes)<- colnames(genesRankingGlobal@ord)
        esetFilteredDataFrame <- as.data.frame(esetFiltered)
        showWarningMaxGenes <- NULL
        numGenesTPglobal <- NULL
        numGenesClassGE <- NULL
        genesRankingLoops <- NULL
        
        if (verbose) { message(paste(format(Sys.time(), "%H:%M:%S"),"- All required parameters have been calculated and checked. Building classifier...")); flush.console()}
        
        # External Cross Validation / Train Classifier    
        for( i in 1:numCV.total)         
        {
            maxGenesTrain <- maxGenesTrainGlobal
            if (buildClassifier && (i== numCV.total))    # Building classifier.  Train samples =  all samples
            {
                genesRankingLoop <- genesRankingGlobal
                trainSamples <- 1:numSamples
            
            }else # CV loop
            {
                # Divide samples in Train and Test
                testSamples <- cvSamples$test[[i]]
                trainSamples <- cvSamples$train[[i]]
        
                if(sameNumSamplesClass)
                {
                        # Make sure all the classes have the same number of train samples (remove some if any has more than others)
                        minNumSamplesClass <- min(summary(sampleLabels[trainSamples]))
                        removedSamples <- NULL         # To avoid removing the same ones all the time
                        for (clNum in classes)
                        {
                            numSamplesClass <- summary(sampleLabels[trainSamples])[clNum]
                            if(numSamplesClass > minNumSamplesClass) 
                            {
                                classRemovedSamples<- removedSamples[which(sampleLabels[removedSamples] == clNum)]
                                classSamples <- trainSamples[which(sampleLabels[trainSamples] == clNum)]
                                
                                tableRemovedSamples<-table(classRemovedSamples)
                                if(sum(!classSamples%in%classRemovedSamples)==0)
                                {
                                    minRemoved <-min(tableRemovedSamples,numElemClass) 
                                }else{
                                    minRemoved <- 0
                                }
                                
                                sample2Remove <- classSamples[sample(1:numSamplesClass,1)]
                                while( (sample2Remove%in%names(tableRemovedSamples)) && (tableRemovedSamples[which(names(tableRemovedSamples)==sample2Remove)]>minRemoved))
                                {
                                    sample2Remove <- classSamples[sample(1:numSamplesClass,1)]
                                }
                                
                                trainSamples <- trainSamples[-which(trainSamples == sample2Remove)]
                                testSamples <- c(testSamples, sample2Remove)
                                
                                removedSamples<- c(removedSamples, sample2Remove)
                            }
                        }
                }        
                # Genes Ranking for these samples    
                genesRankingLoop <- PEB(esetFiltered[,trainSamples], sampleLabels[trainSamples])
            }
            genesRankingLoops<- c(genesRankingLoops, list(genesRankingLoop))
                            
            #Comprobar maxGenesTrain
            lessMaxGenes <- (maxGenesTrain > numGenes(genesRankingLoop))
            if(any(lessMaxGenes))
            {
                maxGenesTrain[which(lessMaxGenes)] <- numGenes(genesRankingLoop)[which(lessMaxGenes)]
                showWarningMaxGenes<- rbind(showWarningMaxGenes, maxGenesTrain)
                if(i== numCV.total) rownames(showWarningMaxGenes)[dim(showWarningMaxGenes)[1]] <- "Classifier"
                else rownames(showWarningMaxGenes)[dim(showWarningMaxGenes)[1]] <- paste("iter", i, sep="")
            }
                
            # Select the best number of genes from the Ranking to train the classifier
            underMin <- which(numGenes(genesRankingLoop)<minGenesTrain)
            minGenesTrainLoop <- rep(minGenesTrain, length(numGenes(genesRankingLoop)))
            minGenesTrainLoop[underMin] <- numGenes(genesRankingLoop)[underMin]

            numGenesClass <- matrix(data=minGenesTrainLoop, nrow=numIters, ncol=length(numGenes(genesRankingLoop)), byrow=TRUE)
            colnames(numGenesClass) <- gClasses(genesRankingLoop)
            if(ncol(genesRankingLoop@ord)==1) colnames(numGenesClass) <- "BothClasses"
            
            for(k in 1:numIters)
            {
                # Select initial number of genes (minGenesTrain) from all the classes
                if (is.null(genesRedundancy))
                {
                     genes <- getRanking(genesRankingLoop, showGeneID=TRUE)$geneID[1:minGenesTrain,, drop=FALSE] # En vez de minGenesTrainLoop, se quitan los NA despues
                } else 
                {
                    genes <- sapply(genesRedundancy$nonRedundantGenes, function(x) x[1:minGenesTrain]) 
                }
                if(any(is.na(genes))) genes <- genes[-which(is.na(genes))]
                
                # Initialize variables with the result from the initial genes
                numGenesTP <- matrix(nrow=0, ncol=ncol(genesRankingLoop@ord)+1)
                colnames(numGenesTP) <- c(colnames(genesRankingLoop@ord), "TruePositives")   # NOT inverse of error: % of True Positives out of the total (including NotA).
                cvInternSamples <- cvSplitSamples(numCV.intern, sampleLabels[trainSamples])
                geneSelectClassifier <- linear.SVM(esetFiltered[genes,trainSamples, drop=FALSE], sampleLabels[trainSamples], cvInternSamples, minProbAssignCoeff=minProbAssignCoeff, minDiffAssignCoeff=minDiffAssignCoeff)
                numGenesTP <- rbind(numGenesTP, c(numGenesClass[k,], geneSelectClassifier$sensitGlobal))
                
                # While error is > 0 and none of the classes have reached maxGenesTrain
                numGenes2add <- numGenesTP[nrow(numGenesTP),-ncol(numGenesTP)]  
                worstClasses <- 1:numClasses
                while((any(numGenes2add< maxGenesTrain) && any(which(numGenes2add < maxGenesTrain) %in% worstClasses)) &&  (geneSelectClassifier$sensitGlobal<1 || continueZeroError)) 
                {
                    # Add one gene to the classes with errors   
                    if (numClasses > 2 ) 
                    {
                        if(geneAdd=="V0")
                        {
                            worstClasses <- c(1:length(classes))    #V0: Adds to all classes
                        }else
                        {
                            worstClasses <- which(geneSelectClassifier$sensitClass<1) #V1: Classes which don't have all in the diagonal
                            if(geneAdd=="V1.b") worstClasses <- unique(c(worstClasses, which(apply(geneSelectClassifier$confMatrix,2,sum)>diag(geneSelectClassifier$confMatrix)))) #V1.b: Classes with which they were confussed
                            if(length(worstClasses)==0 && continueZeroError)  worstClasses <- 1:numClasses
                        }                                                    
                    }else worstClasses <- 1    #If only two classes: both
                    # Remove those which reached maxGenesTrain (or without more available genes)
                    if (is.null(genesRedundancy))
                    {
                        worstClasses <- worstClasses[ numGenes2add[worstClasses] < maxGenesTrain[worstClasses]] 
                    }else
                    {
                        worstClasses <- worstClasses[ numGenes2add[worstClasses] < sapply(genesRedundancy$nonRedundantGenes, length)[worstClasses]] 
                        worstClasses <- worstClasses[ numGenes2add[worstClasses] < maxGenesTrain[worstClasses]] 
                    }
                            
                    if(length(worstClasses)>0) # initialized: worstClasses <- 1:numClasses
                    {
                        # Add line & genes to numGenesTP
                        numGenesTP <- rbind(numGenesTP, c(numGenes2add,1))
                        numGenesTP[nrow(numGenesTP),-ncol(numGenesTP)][worstClasses] <-  numGenes2add[worstClasses] +1
                    
                        # Numero de gene por clase: Ultima fila numGenesTP, sin ceros
                        numGenes2add <- numGenesTP[nrow(numGenesTP),-ncol(numGenesTP)]         
                        numGenes2addNonZero <- numGenes2add[numGenes2add>0]                            # 1:numGenes2add , si es 0...
                            
                        # Select the genes
                        if (is.null(genesRedundancy))
                        {
                            genes <- NULL
                            for (cl in names(numGenes2addNonZero))
                            {
                                genes<- c(genes, getRanking(genesRankingLoop, showGeneID=TRUE)$geneID[1:numGenes2addNonZero[cl],cl])
                            }
                        } else
                        {
                            genes<-NULL
                            notEnoughNonRedundantGenes <- sapply( genesRedundancy$nonRedundantGenes, length)[which(numGenes2addNonZero > sapply( genesRedundancy$nonRedundantGenes, length)[names(numGenes2addNonZero)])]
                            numGenes2addNonZero[names(notEnoughNonRedundantGenes)] <- notEnoughNonRedundantGenes
                            numGenes2addNonZero <- numGenes2addNonZero[numGenes2addNonZero>0]
                            for (cl in names(numGenes2addNonZero))
                            {
                                genes <- c(genes, genesRedundancy$nonRedundantGenes[[cl]] [1:numGenes2addNonZero[cl]])
                            }
                        } 

                        # Train the classifier
                        geneSelectClassifier <- linear.SVM(esetFiltered[genes, trainSamples, drop=FALSE], sampleLabels[trainSamples], cvInternSamples, minProbAssignCoeff=minProbAssignCoeff, minDiffAssignCoeff=minDiffAssignCoeff) # V4: Usando NA tb internamente
                        # Take note of its accuracy
                        numGenesTP[nrow(numGenesTP),"TruePositives"] <- geneSelectClassifier$sensitGlobal
                    }
                }
                
                # Find the number of genes with best performance
                numGenesClass[k,] <-  rbind(NULL, numGenesTP[which(numGenesTP[,"TruePositives"] == max(numGenesTP[,"TruePositives"])),])[1,1:(ncol(numGenesTP)-1)]    # Best True Positives
                if (buildClassifier && ( i== numCV.total)) numGenesTPglobal <- c(numGenesTPglobal, list(numGenesTP))
            }
            if (estimateGError && !(buildClassifier && (i== numCV.total)))
            {
                numGenesClassGE <- c(numGenesClassGE, list(numGenesClass))
                names(numGenesClassGE)[length(numGenesClassGE)] <- paste("CV", eval(i), sep="")
            }
            
            #Select the best combination of genes
            # Max number of genes
            if(geneSelection=="max")    numTrainGenes[i,] <- apply(numGenesClass,2,function(x) max(x))    
            if(geneSelection=="mean")    numTrainGenes[i,] <- apply(numGenesClass,2,function(x) round(mean(x)))
            if(geneSelection=="maxSoutliers") numTrainGenes[i,]<-apply(numGenesClass, 2, function(x) max(x[which(x<=(mean(x)+(outlThreshold*sd(x))))]))

            # Obtain the matrix of best genes 
            numBestGenes <- max(max(lp), max(numTrainGenes[i,]))  # using lp for global Ranking
            
            if (dim(genesRankingLoop@ord)[1] < numBestGenes) numBestGenes <- dim(genesRankingLoop@ord)[1]

            bestGenes <-  getRanking(getTopRanking(genesRankingLoop, numBestGenes), showGeneID=TRUE)$geneID

            if(!is.null(genesRedundancy))             
            {                
                nonRedundantBestGenes <-matrix(ncol=length(genesRedundancy$nonRedundantGenes), nrow=max(numNonRedundantGenes))
                    colnames(nonRedundantBestGenes) <- names(genesRedundancy$nonRedundantGenes)

                # All non redundant genes are over lpThreshold (requirement in calculation)
                for(cl in 1:ncol(bestGenes)) # Just in case the nonRedundantGenes are not in order 
                {
                    nonRedundantBestTemp  <- bestGenes[which(bestGenes[,cl] %in% genesRedundancy$nonRedundantGenes[[cl]]), cl]
                    if(length(nonRedundantBestTemp)>0) nonRedundantBestGenes[1: length(nonRedundantBestTemp),cl] <- nonRedundantBestTemp
                }        
                numNonRedundantBest <- apply(nonRedundantBestGenes, 2, function(x) length(na.omit(x)))
                bestGenes <-  nonRedundantBestGenes 
                    if (!is.matrix(bestGenes)) bestGenes <- cbind(NULL,bestGenes) 
                    if(is.null(colnames(bestGenes))&& numClasses ==2) colnames(bestGenes) <- colnames(numTrainGenes)
            }    

            # Select genes for training classifier
            if (max(numTrainGenes[i,]) > dim(bestGenes)[1])
            {
                #Nunca se deberia entrar aqui...
                warning ("There arent enough genes to build the classifier with the minimum error.")
                trainGenes <- bestGenes        
            }else
            {
                trainGenes <- matrix(ncol=length(numTrainGenes[i,]), nrow= max(numTrainGenes[i,]))
                for (cl in 1:length(numTrainGenes[i,]))
                {
                    if(numTrainGenes[i,cl] > 0) trainGenes[1:numTrainGenes[i,cl],cl] <- bestGenes[1:numTrainGenes[i,cl],cl]
                }    
            }
            colnames(trainGenes) <- colnames(numTrainGenes)        
            buildGenesVector <- trainGenes[which(trainGenes!="NA")]
            #if (!is.matrix(trainGenes)) trainGenes <- cbind(NULL,trainGenes) 
            

            # Train classifier
            finalClassifier  <-  svm(x=t(esetFilteredDataFrame[buildGenesVector, trainSamples, drop=FALSE]), y=sampleLabels[trainSamples], C=1, kernel="linear", probability=TRUE )           
            
                            # If there is only 1 gene, $SV is not labeled -> SV.dif & query.predictor do not work
                            if (is.null(dimnames(finalClassifier$SV)) && numClasses ==2)
                            {
                                    dimnames(finalClassifier$SV)<- vector("list",2)
                                    colnames(finalClassifier$SV) <- buildGenesVector
                            }
                                                                
            if (estimateGError && !(buildClassifier && ( i== numCV.total)))    # Loop for estimateGError
            {        
                    genesUsed<- c(genesUsed, trainGenes=list(trainGenes))
                    
                    #Evaluate the classifier built
                    queryResult <- queryGeNetClassifier(finalClassifier, esetFiltered[,testSamples], minProbAssignCoeff=minProbAssignCoeff, minDiffAssignCoeff=minDiffAssignCoeff, verbose=FALSE)
                    
                    prediction             <- factor()
                    levels(prediction)     <- c(levels(sampleLabels), "NotAssigned")
                    for (l in 1:length(queryResult$class)) prediction[l] <- queryResult$class[l]
                    
                    testLabels<-sampleLabels[testSamples]
                    levels(testLabels)     <- c(levels(sampleLabels), "NotAssigned")    
                    
                    assignedSamples <- which(prediction!="NotAssigned")
                    globalError[i] <- sum(prediction[assignedSamples] != testLabels[assignedSamples]) 
                    globalResults <- c(globalResults, queryResult)
                    
                    # Matriz de confusion - Real x prediction    
                    mxcfi <- table(testLabels, prediction)[1:numClasses,]    
                    mxcf <- mxcf + mxcfi
                    
                    if (verbose) { message(paste(format(Sys.time(), "%H:%M:%S")," - ",i," out of " ,ifelse(buildClassifier,numCV.total-1,numCV.total)," cross-validation loops finished.", sep="")); flush.console()}
            }    
        }
        if(!is.null(showWarningMaxGenes)) 
        {
            warning("The maximum numbers of genes passed as argument is bigger than the number of available genes in some of the internal loops. These were used instead:")
            print(showWarningMaxGenes)
        }
    }

    #########################################################
    ###     Calculate: CV statistics, netClassificationGenes, classificationGenesDetails
    #########################################################
        
    # CV statistics 
    if(estimateGError)
    {    
        # General Stats
        numTestedSamples <- sum(mxcf)

        confusionMatrixStats <- externalValidation.stats (mxcf, numDecimals=numDecimals)
        predictionStats <- querySummary(globalResults, numDecimals=numDecimals, showNotAssignedSamples=TRUE, verbose=FALSE)
        probMatrix <- externalValidation.probMatrix(globalResults, sampleLabels, numDecimals=numDecimals)
      
        # Gene Stats
        genesStats<-NULL
        if(numClasses > 2)
        {
            classList<-classes
        }else
        {
            classList <- "BothClasses"
        }

        allGenesDetails <- lapply(genesRankingLoops[1:numCV.extern], function(x) genesDetails(x))   # $iter$cl[gen, "ranking")
        for(cl in classList)
        {            
            allGenesClass <- NULL
            for(iter in 1:numCV.extern)
            {
                genesClass<- genesUsed[iter]$trainGenes[,cl][which(genesUsed[iter]$trainGenes[,cl]!="NA")]   #GenesUsed: TrainGenes from GE loops
                if(length(genesClass)>0) allGenesClass <- rbind(allGenesClass,cbind(genesClass, 1:length(genesClass)))
            }
            # If buildClassifier, add the final train genes to the genes class in case there is any missing one in the GE loops
            if(buildClassifier)                         
            {
                genesClass <- trainGenes[!is.na(trainGenes[, cl]), cl] 
                genesClass <- genesClass[!genesClass %in% allGenesClass[,1]]
                if(length(genesClass)>0) allGenesClass <- rbind(allGenesClass, cbind(genesClass, rep(NA,length(genesClass))))
            }

            colnames(allGenesClass)<- c("gen","rank")    
            timesChosen <- table(allGenesClass[,1, drop=FALSE], dnn="timesChosen")
            # Substract from those with NA
            if(any(is.na(allGenesClass[,2]))) 
            {
                naGenes <- allGenesClass[which(is.na(allGenesClass[,2])),1]
                timesChosen[names(timesChosen)%in% naGenes] <- timesChosen[names(timesChosen)%in% naGenes, drop=FALSE]-1
            }
            
            classTable <- data.frame(cbind(timesChosen), chosenRankMean=NA, chosenRankSD=NA, gERankMean=NA, gERankSD=NA)
            if(!is.null(geneLabels)) classTable <- cbind(geneSymbols=extractGeneLabels(geneLabels, rownames(classTable)), classTable)
            for(gene in unique(allGenesClass[,1]))  
            {                                
                geneRank <- as.integer(allGenesClass[which(allGenesClass[,1]==gene),2])
                classTable[gene,"chosenRankMean"] <- round(mean(geneRank),2)
                if(classTable[gene,"timesChosen"]>1) classTable[gene,"chosenRankSD"] <- round(sd(geneRank),2)
                else classTable[gene,"chosenRankSD"] <- 0
                
                geneRanks <- sapply(allGenesDetails[1:numCV.extern], function(x) x[[cl]][gene,"ranking"])
                classTable[gene,"gERankMean"] <- round(mean(geneRanks, na.rm=TRUE),2)
                classTable[gene,"gERankSD"] <- round(sd(geneRanks, na.rm=TRUE),2)
            }
            classTable <- classTable[order(-classTable[,"timesChosen"], classTable[,"gERankMean"]),, drop=FALSE]

            genesStats<- c(genesStats, list(classTable))
        }

        names(genesStats) <- classList
        generalizationError  <- new("GeneralizationError", 
                                                                accuracy=confusionMatrixStats$global, 
                                                                sensitivitySpecificity=confusionMatrixStats$byClass, 
                                                                confMatrix=mxcf, 
                                                                probMatrix=probMatrix, 
                                                                querySummary=predictionStats, 
                                                                classificationGenes.stats=genesStats, 
                                                                classificationGenes.num=numTrainGenes[1:numCV.extern,,drop=FALSE] )
    }
        
    # If genesNetwork was calculated
    # add Discriminant Power and Redundancy information to the genes ranking
    # Redundancy calculation
    isRedundant <- NULL
    if(!is.null(genesNetwork) && (buildClassifier || estimateGError))
    {
        if((removeCorrelations && removeInteractions)  || (skipInteractions && removeCorrelations))
        {
            nonRedundant <- genesRedundancy$nonRedundantGenes
        }else # If the correlations + interactions (both) hasnt been removed, they has to be calculated
        {
            nonRedundant <- remove.redundancy(esetFiltered, rankENSG, lpMaxGenes, genesNetwork)$nonRedundantGenes        
        }
        
        genesLpMax <- getRanking(getTopRanking(genesRankingGlobal, numGenesClass=lpMaxGenes), showGeneID=TRUE)$geneID
        genesLpMax <- as.vector(genesLpMax)
        genesLpMax <- genesLpMax[!is.na(genesLpMax)]
        
        isRedundant[1:length(genesLpMax)] <- TRUE
        names(isRedundant) <- genesLpMax
        
        for(cl in 1:length(nonRedundant)) isRedundant[names(isRedundant) %in% nonRedundant[[cl]]] <- FALSE
        
        if(!all(is.na(isRedundant))) genesRankingGlobal <- setProperties(genesRankingGlobal, isRedundant=isRedundant)
    }
    
    # Discriminant Power
    discriminantPower <- NULL
    if(buildClassifier) 
    {
        for(cl in 1:dim(trainGenes)[2])
        {
            clGenes <- which(trainGenes[,cl]!="NA")
            if(length(clGenes)>0)
            {
                tempDP <- t(sapply(trainGenes[clGenes,cl, drop=FALSE], function(x) SV.dif(finalClassifier, x, originalGeneNames=rownames(esetFiltered))))

                tempDP <- data.frame(discriminantPower= tempDP[,"discriminantPower", drop=FALSE], discrPwClass= tempDP[,"discrPwClass",drop=FALSE])
                tempDP[,1]<-as.numeric(tempDP[,1])
                tempDP[,2]<-as.character(tempDP[,2])

                discriminantPower <- rbind(discriminantPower, tempDP)
            }
        }
    }
    
    # Add properties to classificationGenes
    classificationGenesRanking <- NULL
    if(buildClassifier) 
    {
        # GE Rank Mean
        gERankMean <- NULL
        if(estimateGError) 
        {
            for(cl in 1:length(genesStats))
            {
                classGenes <- trainGenes[which(trainGenes[,cl] %in% rownames(genesStats[[cl]])),cl]
                names(classGenes) <- classGenes
                classGenes <- genesStats[[cl]][classGenes,"gERankMean"] 
                gERankMean <- c(gERankMean, classGenes)
            }
        }
        
        classificationGenesRanking <- extractGenes (genesRankingGlobal, trainGenes)
        classificationGenesRanking <- setProperties(classificationGenesRanking, discriminantPower = discriminantPower, gERankMean = gERankMean) # isRedundant = isRedundant,
    }

    

    #########################################################
    ###     Verbose, plots and returns
    #########################################################
    
    ####### Verbose (Print results)
    if(buildClassifier && verbose)
    {
        # Print basic results
        message(paste(format(Sys.time(), "%H:%M:%S"),"Classifier built for the following classes:"))
        for (i in 1:numClasses) message(paste("\tC",i,": ",classes[i],sep=""))
        message(paste("Total number of genes included in the classifier: ",sum(numTrainGenes[numCV.total,]),". Number per class: ",sep=""));
        print(numTrainGenes[numCV.total,])            
        message(paste("\nClassifier trained with ",numSamples," samples.",sep=""))
        if(IQRfilterPercentage>0) message(paste("Gene expression matrix filtered to remove ",(IQRfilterPercentage*100),"% of the genes (IQR filter). Number of genes left: ",dim(esetFiltered)[1],sep=""))
        message(paste("Gene associations calculated between genes with posterior probability over  ",lpThreshold,sep=""))        # ," (top ",100-(lpThreshold*100),"% of the genes ranking)"
        message(paste("Classifier built removing gene associations/redundancy: ",(removeInteractions || removeCorrelations),sep=""))
        message(paste("\tCorrelations removed: ",ifelse(removeCorrelations, paste("Pearsons Correlation Coeff. >", correlationsThreshold, sep="") , "FALSE"),sep=""))
        message(paste("\tInteractions removed: ",ifelse(removeInteractions, paste("Mutual Information Coeff. >", interactionsThreshold, sep="") , "FALSE"),sep=""))        
        message(paste("Maximum number of genes explored to find the optimum classifier:",sep=""))
        print( apply(t(sapply(numGenesTPglobal, function(x) x[nrow(x),-ncol(x)])), 2, max))

        #if utilizado genesNetwork/genesing... avisar...
    }
    if (estimateGError && verbose)
    {
        message("\nGeneralization error:\n")
        print(generalizationError@accuracy)
        print(generalizationError@sensitivitySpecificity)
        message("\nFor more details use overview(Object@generalizationError)\n")
    }
    
    # Avisar si se han utilizado genes por debajo del lpThreshold 
    # Default lpThreshold: 0.95
    # lp: Number of genes per class with PostProb over lpThreshold    
    if(buildClassifier || estimateGError)
    {
        firstWarning <- TRUE
        for(i in 1: length(lp))
        {
            if(any(numTrainGenes[,i]> lp[i]))
            {
                    if(firstWarning)
                    {
                        warning (paste("For these classes, it was needed to take some genes with posterior probability under the ", lpThreshold," threshold:",sep=""))
                        firstWarning<-FALSE
                    }
                    if(estimateGError) message(paste(names(lp[i]), ": ", lp[i]," genes over threshold. Between ",min(numTrainGenes[,i])," and ",max(numTrainGenes[,i])," genes needed to train the classifier.", sep=""))
                    else message(paste(names(lp[i]), ": ", lp[i]," genes over threshold. ", numTrainGenes[,i]," genes needed to train the classifier.", sep=""))
            }
        }
    }
    
    #### RETURNS
        # testCall <- function(){return (match.call()) }
        # geNetClassifierReturn <- new("GeNetClassifierReturn", call=testCall())
    geNetClassifierReturn <- new("GeNetClassifierReturn", call=match.call())

    if(buildClassifier)
    { 
        geNetClassifierReturn@classifier <- list(SVMclassifier=finalClassifier)
        geNetClassifierReturn@classificationGenes <- classificationGenesRanking
    }
    
    if(estimateGError)
    {
        geNetClassifierReturn@generalizationError <- generalizationError                                                                                                                                 
    }
 
    # Return ranking of all genes or only significant 
    # (If significant remove correlations or interactions if required)
        
    if (returnAllGenesRanking) 
    {    
        geNetClassifierReturn@genesRanking <- genesRankingGlobal
        geNetClassifierReturn@genesRankingType <- "all"
    }else 
    {    # Ranking of the significant genes
        if(is.null(genesRedundancy)) 
        {                                                                                                                                    
            # Select genes over lpThreshold
            plotGenesRankingRanking <- getTopRanking(genesRankingGlobal, numSignificantGenes(genesRankingGlobal))
        }else            
        {
            if(! exists("nonRedundantBestGenes"))
            {
                nonRedundantBestGenes <- matrix(ncol=length(names(genesRedundancy$nonRedundantGenes)), nrow=max(sapply(genesRedundancy$nonRedundantGenes, length)))
                for( nr in 1:length(genesRedundancy$nonRedundantGenes))nonRedundantBestGenes[1:length(genesRedundancy$nonRedundantGenes[[nr]]),nr] <- genesRedundancy$nonRedundantGenes[[nr]]
                colnames(nonRedundantBestGenes) <- names(genesRedundancy$nonRedundantGenes)
                numNonRedundantBest <- apply(nonRedundantBestGenes, 2, function(x) length(na.omit(x)))
            }
            # All non redundant genes are over lpThreshold (requirement in calculation)
            newOrd <- matrix(nrow=max(numNonRedundantBest) , ncol=ncol(genesRankingGlobal@ord) )
            colnames(newOrd)<-colnames(genesRankingGlobal@ord) 
            index<-0
            for(cl in colnames(genesRankingGlobal@ord)) 
            {
                    newOrd[1:numNonRedundantBest[cl],cl]<- (index+1):(index+numNonRedundantBest[cl])
                    index <- index+numNonRedundantBest[cl]
            }
            nonRedundantBestGenesVector <-  as.vector(nonRedundantBestGenes[which(nonRedundantBestGenes!="NA")])
            plotGenesRankingRanking <- new("GenesRanking",  postProb=genesRankingGlobal@postProb[nonRedundantBestGenesVector,], meanDif=genesRankingGlobal@meanDif[nonRedundantBestGenesVector,], numGenesClass=numNonRedundantBest , ord=newOrd)
        }    
        # (Ranking de la ultima vuelta: La que entrena el clasificador final)
        geNetClassifierReturn@genesRanking <- plotGenesRankingRanking
        if(is.null(genesRedundancy)) {
            geNetClassifierReturn@genesRankingType <- "significant" 
        } else {
            geNetClassifierReturn@genesRankingType <- "significantNonRedundant" 
        }
    }
    
    # Network and Relations
    if(!is.null(genesNetwork))
    {
        # Add full genesRanking & networks to return if requested
        if (returnTopGenesNetwork)
        { 
            geNetClassifierReturn@genesNetwork <- genesNetwork
            geNetClassifierReturn@genesNetworkType <- "topGenes"
        }
    }
    
    ####### Plots
    if(!is.null(plotsName))
    {        
        plotGeNetClassifierReturn(geNetClassifierReturn, fileName=plotsName, lpThreshold=lpThreshold, numGenesLpPlot=numBestGenes, numGenesNetworkPlot=numGenesNetworkPlot, geneLabels=geneLabels, verbose=FALSE)
        if(!is.null(classificationGenesRanking)) plotExpressionProfiles(eset, classificationGenesRanking, fileName=paste(plotsName,"_expressionProfiles.pdf",sep=""), sampleLabels=sampleLabels, showMean=TRUE, labelsOrder=labelsOrder,verbose=FALSE, type=c("lines","boxplot"))
        plotCVgE <- FALSE
        if(buildClassifier || (estimateGError && plotCVgE)) plotErrorNumGenes(numGenesTPglobal, numGenesClass, numTrainGenes, numGenesClassGE, paste(plotsName,"_errorNumGenes.pdf",sep=""), plotCVgE=plotCVgE) 
        if(estimateGError)
        {
            pdf(paste(plotsName,"_assingments.pdf",sep=""))
            plotAssignments(queryResult=globalResults, realLabels=sampleLabels, minProbAssignCoeff=minProbAssignCoeff, minDiffAssignCoeff=minDiffAssignCoeff)
            dev.off()
        }
        if (verbose) { message(paste("The plots were saved in ",getwd()," with the prefix '",plotsName,"'.",sep="")); flush.console()}
    }
        
    return(geNetClassifierReturn)
}

    
