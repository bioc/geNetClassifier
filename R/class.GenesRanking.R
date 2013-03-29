

########################################
###
### CLASS:      GenesRanking
###
########################################
### METHODS:
### - initialize
### - show
### - plot.GenesRanking
### - getRanking 	
### - overview		# Renombrar a summary/overview/preview?
### - gClasses		    # getClasses exists in findClass {methods}...
### - numGenes
### - numSignificantGenes
### - getTopRanking
### - extractGenes		
### - setProperties
### - genesDetails 	
###    -- Funcion summary?
########################################

setClass(Class="GenesRanking",
		representation=representation(
							postProb="matrix",
							meanDif="matrix",
							numGenesClass="numeric",
							ord="matrix",
							geneLabels="character",
							discriminantPower="data.frame",
							isRedundant="logical",
							gERankMean="numeric")
		)
	
	# To initialize with content, all slots should be provided.
	# To initialize an empty object, the classes should be provided.
	# testGR<-new("GenesRanking", postProb=classificationGenesRanking@postProb,meanDif=classificationGenesRanking@meanDif, numGenesClass=classificationGenesRanking@numGenesClass,ord=classificationGenesRanking@ord,  discriminantPower=discriminantPower)
	setMethod("initialize", "GenesRanking", function(.Object, classes=NULL, postProb=NULL,meanDif=NULL, numGenesClass=NULL,ord=NULL,  geneLabels=NULL, discriminantPower=NULL, isRedundant=NULL, gERankMean=NULL) 
	{
		if(!any(!c(is.matrix(postProb),is.numeric(numGenesClass),is.matrix(ord)))) 	# All are provided with right data type
		{
			if(is.matrix(numGenesClass)) numGenesClass <- numGenesClass[1,,drop=T]
			# Check dimensions 
			# Dim[1] (gene number)
			if(max(numGenesClass)>dim(ord)[1]) stop("There cant be more genes by class than availables in ord.")
			if(dim(ord)[1]>max(numGenesClass)) ord <- ord[1:max(numGenesClass),] # They should be equal
			# Dim[2] (classes)
			if(dim(postProb)[2] != length(numGenesClass)+1) stop("The number of classes in postProb does not match numGenesClass.")

			if(dim(ord)[2] != length(numGenesClass)) stop("The number of classes in ord does not match numGenesClass.")
			
			# Check names
			if(any(!colnames(ord) %in% names(numGenesClass)) || any(!names(numGenesClass) %in% colnames(ord))) stop("The classes in ord do not match the ones in numGenesClass.")
			if(any(!colnames(postProb) %in% c("Null Hipothesis", names(numGenesClass))) || any(!names(numGenesClass) %in% colnames(postProb))) stop("The classes in postProb do not match the ones in numGenesClass.")
			
			if(is.matrix(meanDif))
			{
				if(length(numGenesClass)>1 &&(dim(meanDif)[2] != length(numGenesClass))) stop("The number of classes in meanDif does not match numGenesClass.")
				if(dim(postProb)[1] != dim(meanDif)[1]) stop("The number of genes in postProb and in meanDif does not match")
				.Object@meanDif  <- meanDif
			}
			.Object@postProb <- postProb
			.Object@numGenesClass<- numGenesClass
			.Object@ord		 <- ord
		}	
		else 
		{
			if(!is.null(classes)) # Not all are provided, but the classes are. An empty genesRanking will be created.
			{
				.Object@postProb <- matrix(data=0, nrow=0, ncol=length(classes)+1);
				colnames(.Object@postProb)<- c("Null Hipothesis", classes)
				.Object@meanDif  <- matrix(data=0, nrow=0, ncol=length(classes));
				colnames(.Object@meanDif)<-classes
				.Object@numGenesClass<-rep(as.numeric(0),length(classes));
				names(.Object@numGenesClass)<-classes
				.Object@ord		 <-  matrix(data=0, nrow=0, ncol=length(classes));
				colnames(.Object@ord)<-classes
			}
			else 
			{
				if(any(c(is.matrix(postProb),is.matrix(meanDif),is.numeric(numGenesClass),is.matrix(ord)))) stop("In order to initialize the object, all slots should be provided (postProb, meanDif, numGenesClass,ord). Otherwise create an empty object and insert slot by slot.")
				else stop('To create an empty GenesRanking, please provide the classes. new("GenesRanking", classes=c("a","b",...))')
			};
		}
		
		## Add optional arguments
		rankingGenes <- rownames(.Object@postProb)
		
		if(!is.null(geneLabels))
		{
			
			.Object@geneLabels <- extractGeneLabels(geneLabels, rankingGenes)
		}
		
		if(!is.null(discriminantPower))
		{
			if(!is.data.frame(discriminantPower)) stop("discriminantPower should be a data.frame")
			if (any(!rownames(discriminantPower) %in% rankingGenes) ||  any(!rankingGenes %in% rownames(discriminantPower))) stop("The genes in discriminantPower do not match the ranking genes.")
			.Object@discriminantPower <- discriminantPower
		}
		
		if(!is.null(isRedundant))
		{
			if(!is.logical(isRedundant)) stop("isRedundant should be logical.")
			if (any(!names(isRedundant) %in% rankingGenes) ||  any(!rankingGenes %in% names(isRedundant))) stop("The genes in isRedundant do not match the ranking genes.")
			.Object@isRedundant  <- isRedundant
		}
		
		if(!is.null(gERankMean))
		{
			if(!is.numeric(gERankMean)) stop("gERankMean should be numeric.")
			if (any(!names(gERankMean) %in% rankingGenes)) stop("The genes in gERankMean do not match the ranking genes.")
			
			# Add NA for missing genes
			notInCV <- rankingGenes[which(!rankingGenes %in% names(gERankMean))]
			gERankMean <- c(gERankMean, rep(NA, length(notInCV)))
			names(gERankMean)[((length(gERankMean)-length(notInCV))+1):length(gERankMean)]<- notInCV
			
			.Object@gERankMean<- gERankMean
		}		
		.Object
	})
	
	# When typing the GenesRanking object, doesn't show the whole object content, but the ranking of the top10 genes of each class
	setMethod("show", "GenesRanking", function(object)
	{
		nRows <- ifelse(nrow(object@ord)>100, 10, nrow(object@ord))
		if(nRows > 0) 
		{
			cat("Top ranked genes for the classes: " ,colnames(object@ord),"\n", sep=" ")
			print(getRanking(object, showGeneLabels=TRUE)[[1]][1:nRows,,drop=F])		# By default shows the geneName if available

			if(nrow(object@ord)>100){
				cat("...\n\nNumber of ranked significant genes (posterior probability over threshold):\n\t",colnames(object@ord),"\n\t",numSignificantGenes(object, lpThreshold=0.75),sep=" ")
				cat("\nTo see the whole ranking (",nrow(object@ord)," rows) use: getRanking(...)",sep="")
			}
			cat("\nDetails of the top X ranked genes of each class: genesDetails(..., nGenes=X)\n")
			#cat("\nTo see an overview of the whole object (top 5 rows of all the slots) use: gSlots(object)\n")		
		}else print("The genes ranking contains no genes.")
	})
	
	plot.GenesRanking <- function(x, y="missing", numGenesPlot=1000, plotTitle="Significant genes", lpThreshold = 0.75, ...)
	{ 
		calculateGenesRanking(precalcGenesRanking=x, numGenesPlot=numGenesPlot, plotTitle=plotTitle, lpThreshold = lpThreshold)
	}

	# Note: showGeneID may also refer to probes or proteins... or whatever the ID in eset is.
	if(!exists("getRanking")) setGeneric(name="getRanking", def=function(object, showGeneLabels=TRUE, showGeneID=FALSE) standardGeneric("getRanking"))
    setMethod(f="getRanking", signature="GenesRanking", definition=function(object, showGeneLabels=TRUE, showGeneID=FALSE) 
	{
		if(!is.logical(showGeneLabels) || !is.logical(showGeneID)) stop("showGeneLabels and showGeneID should be either TRUE or FALSE.")
		if(!showGeneLabels && !showGeneID) stop("Either showGeneLabels or showGeneID should be TRUE, otherwise there is nothing to show.")
		if(length(object@geneLabels)==0)	# Check if it is possible to show the names
		{
			showGeneLabels <- FALSE
			showGeneID <- TRUE
		}
			
		ret <- NULL
		if(showGeneLabels)
		{	
			geneLabels <- apply(object@ord,2,function(x) rownames(object@postProb[x, ,drop=F]))
			if(is.null(geneLabels))			geneLabels <- matrix(nrow=0, ncol=length(gClasses(object)), byrow=F)
			if(!is.matrix(geneLabels)) 	geneLabels <- matrix(geneLabels, ncol=length(gClasses(object)), byrow=F)
			
			geneLabels  <- apply(geneLabels, 2, function(x){object@geneLabels[x]})
			if(is.character(geneLabels))	geneLabels <- matrix(geneLabels, ncol=length(gClasses(object)), byrow=F)
			
			if(is.null(colnames(geneLabels))) colnames(geneLabels) <- gClasses(object)
			ret[["geneLabels"]] <- geneLabels
		}
		if(showGeneID) 
		{
			geneIDs <- apply(object@ord,2,function(x) rownames(object@postProb[x, ,drop=F]))
			if (is.null(geneIDs))   		geneIDs <- matrix(nrow=0, ncol=length(gClasses(object)), byrow=F)
			if (!is.matrix(geneIDs)) geneIDs <- matrix(geneIDs, ncol=length(gClasses(object)), byrow=F)
			if(is.null(colnames(geneIDs))) colnames(geneIDs) <- gClasses(object)  # (Only for 1 row or empty)
			ret[["geneID"]]  <- geneIDs			
		}
		
        return (ret)
    })
	
	if(!exists("overview")) setGeneric(name="overview", def=function(object) standardGeneric("overview")) # Renombrar a summary/overview?
    setMethod(f="overview", signature="GenesRanking", definition=function(object) 
	{
		nRowsShow <- 5
		nRowsOrd <- 5
		if(dim(object@ord)[1] < nRowsOrd) nRowsOrd <- dim(object@ord)[1] 
		if(dim(object@postProb)[1] < nRowsOrd) nRowsOrd <- dim(object@postProb)[1] 

		cat("GenesRanking slots:", sep="")
		cat("\n@numGenesClass:\n", sep="")
		print(object@numGenesClass)
		cat("\n@ord[1:",nRowsOrd,",,drop=F]:\n", sep="")
		print(object@ord[1:nRowsOrd,,drop=F])
		cat("... (",nrow(object@ord)," rows)\n\n", sep="")
		cat("\n@postProb[1:",nRowsShow,",]:\n", sep="")
		print(object@postProb[1:nRowsShow,])
		cat("... (",nrow(object@postProb)," rows)\n\nUsage sample: \nPosterior probability of the top10 ranked gene of the second class (",colnames(object@ord)[2],"): object@postProb[object@ord[1:10,2],]\n",sep="")
		

		# Optional properties
		if(nrow(object@meanDif) > 0) 
		{
			cat("\n\n@meanDif[1:",nRowsShow,",]:\n", sep="")
			print(object@meanDif[1:nRowsShow,])	
			cat("... (",nrow(object@meanDif)," rows)\n\n", sep="")
		}
		if(nrow(object@discriminantPower) > 0)
		{
			cat("\n@discriminantPower[1:",nRowsShow,",]:\n", sep="")
			print(object@discriminantPower[1:nRowsShow,])
			cat("... (",nrow(object@discriminantPower)," rows)\n\n", sep="")
		}
		if(length(object@isRedundant) > 0)
		{
			cat("@isRedundant[1:",nRowsShow,"]:\n", sep="")
			print(object@isRedundant[1:nRowsShow])
			cat("... (length: ",length(object@isRedundant),")\n\n", sep="")
		}
		if(length(object@gERankMean) > 0)
		{
			cat("\n@gERankMean[1:",nRowsShow,"]:\n", sep="")
			print(object@gERankMean[1:nRowsShow])
			cat("... (length: ",length(object@gERankMean),")\n\n", sep="")
		
		}
    })
	
	if(!exists("gClasses")) setGeneric(name="gClasses", def=function(object) standardGeneric("gClasses"))
    setMethod(f="gClasses", signature="GenesRanking", definition=function(object) 
	{
        return (colnames(object@ord))
    })
	
	# Returns the number of available ranked genes per class
	if(!exists("numGenes")) setGeneric(name="numGenes", def=function(object) standardGeneric("numGenes"))
    setMethod(f="numGenes", signature="GenesRanking", definition=function(object) 
	{
		if(class(object)!="GenesRanking") stop("The first argument should be a 'GenesRanking' object.")
		nGenes<- apply(object@ord,2,function(x){sum(!is.na(x))})
		if(!length(nGenes)) nGenes <- 0 
		return(nGenes)
    })
	
	# Returns the number of ranked genes over the threshold (not exactly the same as the lp returned by the main classifier)
	# lp from classifier: number of genes per class over the threshold
	# numSignificantGenes: Number of RANKED genes per class over the threshold (some classes with many genes over the threshold migth have "given" some genes to classes with really few genes
	if(!exists("numSignificantGenes")) setGeneric(name="numSignificantGenes", def=function(object, lpThreshold=0.75, numSignificantGenesType="ranked") standardGeneric("numSignificantGenes"))
    setMethod(f="numSignificantGenes", signature="GenesRanking", definition=function(object, lpThreshold=0.75, numSignificantGenesType="ranked") 
	{
		if(class(object)!="GenesRanking") stop("The first argument should be a 'GenesRanking' object.")
        if(!is.numeric(lpThreshold) || (lpThreshold>=1 || lpThreshold <0)) stop("Lp threshold should be a probability (a number between 0 and 1).")
		if(!is.character(numSignificantGenesType)) numSignificantGenesType <- "ranked"
		if(!numSignificantGenesType %in% c("ranked", "global"))  numSignificantGenesType <- "ranked"
		
		# Calculate genes over lpThreshold
		if(numSignificantGenesType=="ranked")
		{
			nGenes <- rep(0, length(gClasses(object)))
			names(nGenes) <- gClasses(object)
			for(cl in 1:length(nGenes)) nGenes[cl] <- sum(object@postProb[object@ord[1:object@numGenesClass[cl],cl],cl+1] > lpThreshold)	
			nGenes[which(is.na(nGenes))] <- 0
		}
		else
		{
			if (length(gClasses(object)) > 2) 
			{
				nGenes <- apply(object@postProb[,-1], 2, function(x) length(which(x>lpThreshold)))   
				names(nGenes) <- gClasses(object)	
			} else 
			{
			nGenes<- length(which(object@postProb[,-1]>lpThreshold))   
			names(nGenes) <- "BothClasses"
			} 
		}			
		nGenes[which(is.na(nGenes))] <- 0
		return(nGenes)
    })
	

	# Returns a new Genes ranking with the selected top ranked genes of each class
	if(!exists("getTopRanking")) setGeneric(name="getTopRanking", def=function(object, numGenesClass) standardGeneric("getTopRanking"))
    setMethod(f="getTopRanking", signature="GenesRanking", definition=function(object, numGenesClass) 
	{
		if(class(object)!="GenesRanking") stop("The first argument should be a 'GenesRanking' object.")
		if(!is.numeric(numGenesClass)) stop("The second argument should be the number of genes per class to extract.")
		if(is.matrix(numGenesClass)) numGenesClass <- numGenesClass[1,,drop=T]
		if(length(numGenesClass)==1){
			numGenesClass <- rep(numGenesClass, dim(object@ord)[2])
			names(numGenesClass)<-colnames(object@ord)
		}
		numAvailableGenes <- numGenes(object)
		tooMany <- which(numGenesClass>numAvailableGenes)
		numGenesClass[tooMany] <- numAvailableGenes[tooMany]
		
		topGenes <- getRanking(object, showGeneLabels=FALSE, showGeneID=TRUE)[[1]][1:max(numGenesClass),,drop=F] # Needs the ID
		
								## TOP GENES, no hay que quitar los q no se usan?!?!
								# revisar con leucemias
		
		if(max(numGenesClass) == 0)  topGenes <- topGenes[0,, drop=FALSE]

		newOrd <- matrix(nrow=max(numGenesClass) , ncol=ncol(object@ord) )
		colnames(newOrd)<-colnames(object@ord) 
		index<-0
		for(cl in colnames(object@ord)) 
		{
				if(numGenesClass[cl] > 0) 
				{
					topGenes[-c(1:numGenesClass[cl]),cl] <- NA
					newOrd[1:numGenesClass[cl],cl] <- (index+1):(index+numGenesClass[cl])
				}else{
						topGenes[,cl] <- NA
				}
				index <- index+numGenesClass[cl]
		}

		topGenes <- topGenes[which(topGenes!="NA")]
		ret <- new("GenesRanking",  postProb=object@postProb[topGenes,], numGenesClass=numGenesClass , ord=newOrd)
		if(nrow(object@meanDif) > 0) meanDif <- object@meanDif[topGenes,]  else meanDif <- object@meanDif
		if(length(object@geneLabels) > 0 &&  any(!is.na(object@geneLabels))) geneLabels <- object@geneLabels[topGenes] else geneLabels <-  object@geneLabels
		#if(ncol(object@discriminantPower) > 0 &&  any(!is.na(object@discriminantPower))) discriminantPower <- object@discriminantPower[topGenes,] else discriminantPower <-  object@discriminantPower
		if(ncol(object@discriminantPower) > 0) discriminantPower <- object@discriminantPower[topGenes,] else discriminantPower <-  object@discriminantPower
		#if(length(object@isRedundant) > 0 && any(!is.na(object@isRedundant))) isRedundant <- object@isRedundant[topGenes] else isRedundant <- object@isRedundant
		if(length(object@isRedundant) > 0) isRedundant <- object@isRedundant[topGenes] else isRedundant <- object@isRedundant
		#if(length(object@gERankMean) > 0 && any(!is.na(object@gERankMean))) gERankMean <- object@gERankMean[topGenes] else gERankMean <- object@gERankMean
		if(length(object@gERankMean) > 0) gERankMean <- object@gERankMean[topGenes] else gERankMean <- object@gERankMean
		ret <- setProperties(ret, geneLabels=geneLabels, discriminantPower=discriminantPower, meanDif=meanDif, isRedundant=isRedundant, gERankMean=gERankMean)
		
		return( ret)
	})
	
	# Returns a new Genes ranking with the selected genes of each class
	if(!exists("extractGenes")) setGeneric(name="extractGenes", def=function(object, genes) standardGeneric("extractGenes"))
    setMethod(f="extractGenes", signature="GenesRanking", definition=function(object, genes) 
	{
		if(class(object)!="GenesRanking") stop("The first argument should be a 'GenesRanking' object.")
		if(!is.matrix(genes)) stop("The genes should be in a matrix, ordered by class.")
		if(any(!colnames(genes) %in% gClasses(object))) stop("The genes matrix classes, do not match the genesRanking.")
		genes <- genes [, gClasses(object), drop=F]
		geneList <- as.vector(na.omit(as.vector(genes)))
		if(any(!geneList %in% rownames(object@postProb))) stop("The requested genes are not in the genesRanking.")
		
		newPostProb <- object@postProb[geneList, ]
		#if(!is.vector(geneList)) genesOrder <- calculateOrder(newPostProb, untie="bestPostProb")		# WARNING: The gene class is not known. It may differ, if it was initially assigned to a class with low post prob.
		newOrd	<- matrix(ncol=ncol(genes), nrow=nrow(genes))
		colnames(newOrd) <- colnames(genes)
		ordList <- apply(genes, 2, function(x) which(geneList %in% x))
		if(is.matrix(ordList)) { newOrd <- ordList
		}else
		{
			for(o in 1:length(ordList))
			{
				if (length(ordList[[o]])>0) newOrd[1:length(ordList[[o]]), names(ordList)[o]] <- ordList[[o]]
			}
		}

		newNGenesClass <- apply(newOrd, 2, function(x) length(na.omit(x)))
		
		ret <- new("GenesRanking",  postProb=newPostProb, numGenesClass=newNGenesClass, ord=newOrd)
		if(nrow(object@meanDif) > 0) meanDif <- object@meanDif[geneList,]  else meanDif <- object@meanDif
		if(length(object@geneLabels) > 0 &&  any(!is.na(object@geneLabels))) geneLabels <- object@geneLabels[geneList] else geneLabels <-  object@geneLabels
		if(ncol(object@discriminantPower) > 0) discriminantPower <- object@discriminantPower[geneList,] else discriminantPower <-  object@discriminantPower
		if(length(object@isRedundant) > 0) isRedundant <- object@isRedundant[geneList] else isRedundant <- object@isRedundant
		if(length(object@gERankMean) > 0) gERankMean <- object@gERankMean[geneList] else gERankMean <- object@gERankMean
		ret <- setProperties(ret, geneLabels=geneLabels, discriminantPower=discriminantPower, meanDif=meanDif, isRedundant=isRedundant, gERankMean=gERankMean)
		
		return( ret)
	})
		
	# Method to assign values to the properties discriminantPower, isRedundant and gERankMean
	#	testGR <- setProperties(testGR, discriminantPower=discriminantPower)
	if(!exists("setProperties")) setGeneric(name="setProperties", def=function(object, geneLabels=NULL, discriminantPower=NULL, meanDif=NULL, isRedundant=NULL, gERankMean=NULL) standardGeneric("setProperties"))
	setMethod(f="setProperties", signature="GenesRanking", definition=function(object, geneLabels=NULL, discriminantPower=NULL, meanDif=NULL, isRedundant=NULL, gERankMean=NULL) 
	{
		rankingGenes <- rownames(object@postProb)
		
		if(!is.null(geneLabels))
		{
			object@geneLabels <- extractGeneLabels(geneLabels, rankingGenes)
		}
		
		if(!is.null(discriminantPower))
		{
			if(!is.data.frame(discriminantPower)) stop("discriminantPower should be a data.frame")
			if (any(!rownames(discriminantPower) %in% rankingGenes)) stop("The genes in discriminantPower do not match the ranking genes.")
			
			# Add NA for missing genes
			notInProp <- rankingGenes[which(!rankingGenes %in% rownames(discriminantPower))]
			
			addedDP <- matrix(ncol=2, nrow= length(notInProp))
			rownames(addedDP) <- notInProp
			colnames(addedDP) <- colnames(discriminantPower)
			discriminantPower <- rbind(discriminantPower, addedDP)
			
			object@discriminantPower <- discriminantPower
		}
		
		if(!is.null(meanDif))
		{
			if(!is.matrix(meanDif)) stop("meanDif should be a matrix")
			if (any(!rownames(meanDif) %in% rankingGenes)) stop("The genes in meanDif do not match the ranking genes.")
			
			if(nrow(meanDif) > 0)
			{
				# Add NA for missing genes
				notInProp <- rankingGenes[which(!rankingGenes %in% rownames(meanDif))]
				
				addedMD <- matrix(ncol= ncol(meanDif), nrow= length(notInProp))
				rownames(addedMD) <- notInProp
				colnames(addedMD) <- colnames(meanDif)
				meanDif <- rbind(meanDif, addedMD)
			} 
			object@meanDif <- meanDif
		}
		
		if(!is.null(isRedundant))
		{
			if(!is.logical(isRedundant)) stop("isRedundant should be logical.")
			if (any(!names(isRedundant) %in% rankingGenes)) stop("The genes in isRedundant do not match the ranking genes.")
			
			# Add NA for missing genes
			notInProp <- rankingGenes[which(!rankingGenes %in% names(isRedundant))]
			isRedundant <- c(isRedundant, rep(NA, length(notInProp)))
			if(length(notInProp) > 0) names(isRedundant)[((length(isRedundant)-length(notInProp))+1):length(isRedundant)]<- notInProp
			
			object@isRedundant  <- isRedundant
		}
		
		if(!is.null(gERankMean))
		{
			if(!is.numeric(gERankMean)) stop("gERankMean should be numeric.")
			#if (any(!names(gERankMean) %in% rankingGenes)) warning("Some genes in gERankMean are not in the ranking.")
			gERankMean <- gERankMean[which(names(gERankMean) %in% rankingGenes)]
			
			# Add NA for missing genes
			notInProp <- rankingGenes[which(!rankingGenes %in% names(gERankMean))]
			gERankMean <- c(gERankMean, rep(NA, length(notInProp)))
			if(length(notInProp) > 0) names(gERankMean)[((length(gERankMean)-length(notInProp))+1):length(gERankMean)] <- notInProp
			
			object@gERankMean<- gERankMean
		}		
		object
	})
	
	# Metodo q devuelve la tabla de genes:	
	if(!exists("genesDetails")) setGeneric(name="genesDetails", def=function(object, nGenes=NULL, numDecimals=4, classes=NULL, genes=NULL) standardGeneric("genesDetails"))
	setMethod(f="genesDetails", signature="GenesRanking", definition=function(object, nGenes=NULL, numDecimals=4, classes=NULL, genes=NULL) #, showRawPostProb=TRUE) 
	{
		####################
		# Checking arguments
		####################
		#classes: c("AML",...)
		if (is.null(classes)){ classes<-gClasses(object)
		}else if(any(!classes %in% gClasses(object))) stop("The requested classes are not in the ranking object.")
		
		#genes c("geneid1","geneid2",...)
		if(!is.null(genes) && !is.character(genes)) stop("The genes list should be of type character.")
		
		#nGenes: number or c(1,3,5,6,8)
		if(is.null(nGenes))
		{ 
			nGenes <- object@numGenesClass
		}else
		{
			if(is.numeric(nGenes))
			{
				if(length(nGenes)==1 || length(nGenes)==length(classes))
				{
					if (length(nGenes)==1) nGenes <- rep(nGenes, length(classes))

					if (is.null(names(nGenes))) names(nGenes) <- classes
					else if(any(names(nGenes)!=classes)) stop("nGenes class names do not match the ranking/provided class names.")
				}
				else {stop("nGenes should contain the number of genes per class to show.")}
			}
			else stop("nGenes should be numeric.")
		}
		
		####################
		# Creating table
		####################
		dfList <- NULL
		ranking <- getRanking(object, showGeneLabels=FALSE, showGeneID=TRUE)$geneID[0:max(nGenes),,drop=F]    
		for(cl in classes)
		{	
			# Gene Name or ID
			classGenes <- ranking[0:nGenes[cl],cl, drop=F]	
			# Gene Ranking
			classRanking <- c(1:nGenes[cl])
			if(nGenes[cl]==0) classRanking<-classRanking[0]
			genesDF<- data.frame(ranking=classRanking, class=rep(cl,nGenes[cl]))
			rownames(genesDF)<-classGenes
			
			# PostProb
			# if(showRawPostProb) 
			# {
				posteriorProb <- signif(object@postProb[classGenes, cl], numDecimals+1)
				genesDF <- cbind(genesDF, postProb=posteriorProb)
			# } else 
			# {
				# lPosteriorProb <- abs(signif(1-object@postProb[classGenes, cl], numDecimals+1))	# 1- very very very small post prob... results negative!!
				# genesDF <- cbind(genesDF, l_postProb=lPosteriorProb)
			# }
		
			# Add optional properties
			if(nrow(object@meanDif) > 0) 
			{
				exprsMeanDiff <- round(object@meanDif[classGenes, ifelse(length(classes)>1,cl,1)], numDecimals)
				exprsUpDw <- ifelse(as.numeric(exprsMeanDiff)>=0, "UP", "DOWN")
				genesDF <- cbind(genesDF, exprsMeanDiff=exprsMeanDiff, exprsUpDw=exprsUpDw)
			}
			if(any(dim(object@discriminantPower) > 0) && any(!is.na(object@discriminantPower))) genesDF <- cbind(genesDF, object@discriminantPower[rownames(genesDF),])   #discriminantPower=
			if(length(object@gERankMean) > 0 &&  any(!is.na(object@gERankMean))) genesDF <- cbind(ranking=genesDF[,1], gERankMean=object@gERankMean[rownames(genesDF)], genesDF[,-1])
			if(length(object@isRedundant) > 0 && any(!is.na(object@isRedundant))) genesDF <- cbind(genesDF, isRedundant=object@isRedundant[rownames(genesDF)])
			
			if(length(object@geneLabels) > 0 &&  any(!is.na(object@geneLabels))) genesDF <- cbind(GeneName=object@geneLabels[rownames(genesDF)], genesDF)
			
			dfList <- c(dfList, list(genesDF))
		}
		#if(length(classes)==1){   	names(dfList)<- classes #dfList <- dfList[[1]] (Hay q hacer caso especial al usarlo para 2clases)
		names(dfList)<- classes
		
		####################
		# Filter genes 			(keep only the requested ones)
		####################
		if(!is.null(genes))
		{
			for(cl in names(dfList))
			{
				genesClass <- genes[which(genes%in%rownames(dfList[[cl]]))]
				dfList[[cl]] <-  rbind(dfList[[cl]][genesClass,, drop=F])
				if(length(genesClass)==1) rownames(dfList[[cl]] )<-genesClass
			}
		}
		
		return(dfList)
    })
	
	
