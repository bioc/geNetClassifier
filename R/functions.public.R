
## PUBLIC FUNCTIONS ##  
#
# clasificador						(in file: classifier.main.r)
# queryGeNetClassifier
#    externalValidation.stats 
#    externalValidation.probMatrix
#    querySummary 
#
# PLOTS
# plotAssignments
# plotExpressionProfiles
# plotDiscriminantPower			
#	 (private)  errorNumGenes.plot  
#
# calculateGenesRanking
# plotNetwork 			
#
######################
  	
 # queryGeNetClassifier: 
  #
  # Args:
  #   classifier: 
  #   eset: 
  #   verbose: If TRUE, prints _____; if not, not. Default is TRUE.
  #
  # Returns:
  #   
  # WARNING!: The arrays should have been normalized with the samples used for the classifier training.
queryGeNetClassifier <- function(classifier, eset,	minProbAssignCoeff=1, minDiffAssignCoeff=0.8, verbose=TRUE)
{
	# Comprobacion de parametros
	if(is(eset, "ExpressionSet")) eset <- exprs(eset) else if (!is.matrix(eset)) stop("The last argument should be either an expression matrix or an ExpressionSet.")
	esetTdf <- data.frame(t(eset))
	if(is(classifier, "GeNetClassifierReturn")){
		if("classifier" %in% names(classifier)) { classifier <- classifier@classifier$SVMclassifier
		}else stop("'classifier' doesn't contain a trained classifier.")
	}
	if(!is(classifier, "svm")) classifier <- classifier$SVMclassifier
	if(!is(classifier, "svm")) stop("The first argument should be the classifier returned by geNetClassifier.")
	numClasses <- length(classifier$levels)
	if(!is.numeric(minProbAssignCoeff)) stop("'minProbAssignCoeff' should be a coefficient to modify the probability required to assign the sample to a class.")
	if(!is.numeric(minDiffAssignCoeff)) stop("'minDiffAssignCoeff' should be a coefficient to modify the required difference between probabilites to assign the sample to a class.")
		if(minProbAssignCoeff<0 || ((numClasses != 2) &&(minProbAssignCoeff>(numClasses/2)))) stop("'minProbAssignCoeff' should be between 0 and half of the number of classes.")
		if(minDiffAssignCoeff<0 || minDiffAssignCoeff>numClasses) stop("'minDiffAssignCoeff' should be between 0 and the number of classes.")
	genes <- colnames(classifier$SV)
	if(sum(!genes %in% colnames(esetTdf))>0) 
	{	
		m <- rbind(rep(0,length(genes)))	#Transformed into dataframe colnames to make sure they match the $SV
		colnames(m) <- genes
		m <- data.frame(m)
		genes<- colnames(m)
								
		if(sum(!genes %in% colnames(esetTdf))>0) stop("The expression matrix provided does not have the required genes.")
	}
					
	# Calculate minimum probability and difference between probabilities
	rand <- 1/length(classifier$levels)
	minProb <- 2*rand * minProbAssignCoeff
	minDiff <- rand * minDiffAssignCoeff
	if(verbose) { message(paste("Coefficients for assignment: Minimum Probability to be assigned = ",round(minProb,2), ifelse(minProbAssignCoeff==1, "(default)",""),".\n Minimum difference between the probabilities of first and second most likely classes  = ", round(minDiff,2),ifelse(minDiffAssignCoeff==1, "(default)",""), sep=""))  ; flush.console()}
					
	# Calculo de la clasificacion
	esetSelection <- esetTdf[,genes, drop=F]
								#								if (!is.data.frame(esetSelection)) esetSelection <- t(cbind(NULL,esetSelection)) #To avoid error when there is only 1gen/class
	prob <- t(attributes(predict(classifier, esetSelection, probability=TRUE ))$probabilities)
	if(is.null(names(prob))) colnames(prob)<-rownames(esetTdf) # if 2 classes... not labeled. Needed for mxcf...
	
	classes <- factor(apply(prob, 2, function(x) {assignment.conditions(x, minProb, minDiff)}))
	ret <- list(call=match.call(), class= classes, probabilities=prob)
	return(ret)
}


# Calculates stats from the confussion matrix. i.e. sensitivity and specificity of the predictor for each class, global accuracy and call rate (rate of assigned samples)
# Receives a confussion matrix (actual(rows) x prediction(cols)) --> Columns and rows in same order. "NotAssigned" last column
# 100% Sensitivity = Recognizes all positives for the class
# 100% Specificity = Recognizes all negatives for the class

#Renombrado de class.accuracy
externalValidation.stats <- function(confussionMatrix, numDecimals=2) #Confussion matrix
{
	mxcf <- confussionMatrix
	if(!is.matrix(mxcf))stop("The argument should be a confussion matrix.")
	if(!is.numeric(numDecimals)){ numDecimals <- 2
	}else { if (numDecimals<0 || numDecimals>7) numDecimals <- 2 }
	
	if ("NotAssigned" %in% rownames(mxcf))
	{
		mxcf<- t(mxcf)
		warning("The confussion matrix should have the real class in rows and the assigned class in cols. The matrix provided didn't seem to be in the right order so it was transposed:")
	}
	if (!"NotAssigned" %in% colnames(mxcf))  mxcf <- cbind(mxcf, NotAssigned=rep(0, dim(mxcf)[1]))
	
	classes<-unique(c(rownames(mxcf),colnames(mxcf)))
	if(any(classes=="NotAssigned")) classes <- classes[-which(classes=="NotAssigned")]
	nclasses <- length(classes)
	
	#Comprobar si falta alguna clase en filas o columnas
	showWarning<-FALSE
	if (any(!classes %in% colnames(mxcf))) 
	{
		missingClasses <- classes[which(!classes %in% colnames(mxcf))]
		mxcf <- cbind(mxcf, matrix(ncol=length(missingClasses), nrow=dim(mxcf)[1], data=0))
		colnames(mxcf)[which(colnames(mxcf)=="")]<-missingClasses

		showWarning <- TRUE
	}
	if (any(!classes %in% rownames(mxcf))) 
	{
		missingClasses <- classes[which(!classes %in% rownames(mxcf))]
		mxcf <- rbind(mxcf, matrix(nrow=length(missingClasses), ncol=dim(mxcf)[2], data=0))
		rownames(mxcf)[which(rownames(mxcf)=="")]<-missingClasses

		showWarning <- TRUE
	}
	
	#Just in case they are not in order (Diagonal=hits). Will use only the real classes.
	mxcf<- mxcf[,c(classes,"NotAssigned")]																										
	mxcf<- mxcf[c(classes),]		
	
	if(showWarning) 
	{
		warning("There were missing columns or rows in the confussion matrix, empty ones were added.", immediate.=TRUE)
	}
	
	#nclasses <- dim(mxcf)[1]
	numSamples <- sum(mxcf)
	numNA <- sum(mxcf[,nclasses+1])
	
	byClass <- matrix(nrow=nclasses, ncol=4)
	rownames(byClass)<- classes
	colnames(byClass)<-c("Sensitivity","Specificity", "MCC", "CallRate")

	falseNegatives <- array(0,dim=nclasses)
	falsePositives <- array(0,dim=nclasses)
	trueNegatives <- array(0,dim=nclasses)
	truePositives <- array(0,dim=nclasses)
	
	for (i in 1:dim(mxcf)[1])
	{	
		for(j in 1:nclasses)  #dim(mxcf)[2])  Para incluir NA
		{
			if(i!=j) 
			{
				falseNegatives[i] <- falseNegatives[i] + mxcf[i,j]
				falsePositives[j] <- falsePositives[j] + mxcf[i,j]
			}
			trueNegatives[-c(i,j)] <- trueNegatives[-c(i,j)] + mxcf[i,j]
		}
		truePositives[i] <- mxcf[i,i]
	}

	for (i in 1:nclasses) #We need another loop in order to have the whole trueNegatives ready
	{
		#Sensitivity
		byClass[i,1] <- round(100*(truePositives[i]/(truePositives[i]+falseNegatives[i])) ,numDecimals)
		#Specificity
		byClass[i,2] <- round(100*(trueNegatives[i]/(trueNegatives[i]+falsePositives[i])) ,numDecimals)
		#Matthews Correlation Coefficient
		byClass[i,3] <- round(100*( ((truePositives[i]*trueNegatives[i])-(falsePositives[i]*falseNegatives[i])) /  sqrt( (truePositives[i]+falsePositives[i])*(truePositives[i]+falseNegatives[i])*(trueNegatives[i]+falsePositives[i])*(trueNegatives[i]+falseNegatives[i])  ))  ,numDecimals)   
		#Call Rate
		if(i <= dim(mxcf)[1]) byClass[i,4] <- round(100*( (sum(mxcf[i,])- mxcf[i,dim(mxcf)[2]])/sum(mxcf[i,])) ,numDecimals)
	}
	
	global = cbind("Accuracy"=round(100*sum(truePositives)/(numSamples-numNA),numDecimals) , "CallRate"=round(100*((numSamples-numNA)/numSamples),numDecimals) )
	rownames(global) <- "Global"

	return( list(byClass=byClass, global=global, confMatrix=mxcf) )
}

# Returns the matrix with the average probabilities of assigning a sample to each class (only of assigned samples)
# Can receive the result from executing queryGeNetClassifier, or a list of several: queryResult<-c(assignment1, assignment2)
externalValidation.probMatrix<- function(queryResult, realLabels, numDecimals=2)
{
	# Comprobar y concatenar si hay varias queries
	 globalQueryResult <- queryResultCheck(queryResult)
	 
	 #Comprobar el resto de los argumentos
	if(!is.factor(realLabels)) { warning("The second argument (real labels) had to be converted into a factor.", immediate. = TRUE)}	
	realLabels <- factor (realLabels)
	if(is.null(names(realLabels))) 
	{
		if (length(realLabels) == length( names(queryResult$class))) 
		{	
			names(realLabels) <- names(queryResult$class)
		} else stop("The sample's labels vector is not named.")
		warning("The real data labels vector is not named, it will be assumed the labels are in order: the first label applies to the first predicted sample... ", immediate. = TRUE)
	}
	
	if(!is.numeric(numDecimals)) {numDecimals <- 2
	}else if (numDecimals<0 || numDecimals>7) 
	{
		numDecimals <- 2
		warning("The argument 'numDecimals' should be a number between 0 and 7. The default value (2) will be used.", immediate. = TRUE)
	}
	
	#Check if there are "real labels" for all the prediction samples 
	if(sum(!names(globalQueryResult$class) %in% names(realLabels)) >0) stop("There are samples for which the real label was not provided.") 
	#El numero de samples no encaja
	if((length(globalQueryResult$class)!=dim(globalQueryResult$probabilities)[2]) || (sum(!names(globalQueryResult$class) %in% colnames(globalQueryResult$probabilities))>0 )) {stop("The samples in $class and in $probabilities do not match.")}
	
		
	#Calcular estadisticas
	predClasses <- c(levels(realLabels), rownames(globalQueryResult$probabilities)[which(!rownames(globalQueryResult$probabilities) %in%  levels(realLabels))])
	
	probMatrix <- matrix(0,nrow=length(levels(realLabels)), ncol=length(predClasses))
	rownames(probMatrix)<- levels(realLabels)
	colnames(probMatrix)<- predClasses

	for (label in levels(realLabels))
	{
		classAssignments <- globalQueryResult$class[names(realLabels)[which(realLabels==label)]] #Prob for the class samples, even if the prediction was wrong
		assignedSamples <- names(classAssignments)[	which( classAssignments!= "NotAssigned")]
		if (!is.null(assignedSamples))
		{
			if(length(assignedSamples)>1) probMatrix[label,] <- apply(globalQueryResult$probabilities[,assignedSamples], 1, function(x) {mean(x)})[colnames(probMatrix)]
			else probMatrix[label,] <- globalQueryResult$probabilities[,assignedSamples][colnames(probMatrix)]
		}
	}   

	ret<- round(probMatrix,numDecimals)
	return(ret)
}	
	
# Gives basic stats of the probabilities with wich the samples were assigned to the class
# Can receive the result from executing queryGeNetClassifier, or a list of several: queryResult<-c(prediction1, prediction2)
querySummary <- function(queryResult, showNotAssignedSamples=TRUE, numDecimals=2, verbose=TRUE)
{
	# Comprobar y concatenar si hay varias queries
	 globalQueryResult <- queryResultCheck(queryResult)
	 
	 #Comprobar el resto de los argumentos
	if(!is.logical(showNotAssignedSamples)) showNotAssignedSamples <- TRUE
	if(!is.numeric(numDecimals)) numDecimals <- 2
	if(!is.logical(verbose)) verbose <- FALSE
	else if (numDecimals<0 || numDecimals>7) 
	{
		numDecimals <- 2
		warning("The argument 'numDecimals' should be a number between 0 and 7. The default value (2) will be used.", immediate. = TRUE)
	}
	
	#Calcular estadisticas
	if(length(globalQueryResult$class)!=dim(globalQueryResult$probabilities)[2]) {}#El numero de samples no encaja

	numSamples <- length(globalQueryResult$class)
	classes <- rownames(globalQueryResult$probabilities)
	numClasses <- length(classes)

	stats <- cbind(c(rep(0,numClasses)), c(rep(1,numClasses)), c(rep(0,numClasses)),c(rep(NA,numClasses)),c(rep(NA,numClasses)))
	rownames(stats)<-c(classes)
	colnames(stats)<-c("Count", "MinProb", "MaxProb", "Mean", "SD")
	notAssigned<-0
	
	for (c in 1:numClasses)
	{
		temp<-NULL
		for (i in 1:numSamples)
		{
			if (globalQueryResult$class[i] == classes[c])
			{
				prob <- globalQueryResult$probabilities[c,i]
				temp<-c(temp,prob)
				if(prob < stats[c,2]) stats[c,2] <- round(prob, numDecimals)
				if(prob > stats[c,3]) stats[c,3] <- round(prob, numDecimals)
				stats[c,1] <- stats[c,1] + 1
			}
			else if (c==1 && (globalQueryResult$class[i] == "NotAssigned")) notAssigned <- notAssigned+1
		}
		if (stats[c,1]!=0) 
		{
			stats[c,4] <- round(mean(temp), numDecimals)
			stats[c,5] <- round(sd(temp), numDecimals)
		}
		else
		{
			stats[c,2] <- NA
			stats[c,3] <- NA
		}
	}
	
	# Info about NotAssigned samples (most likely class & probs)
	notAssignedSamples="All samples have been assigned."
	if (notAssigned > 0 && showNotAssignedSamples)
	{
		highestProb		<- apply(globalQueryResult$probabilities[,which(globalQueryResult$class == "NotAssigned"), drop=F], 2, function(x) round(max(x), numDecimals))    #--> mayor probabilidad por sample
		highestProbClass <- rownames(globalQueryResult$probabilities)[apply(globalQueryResult$probabilities[,which(globalQueryResult$class == "NotAssigned"), drop=F], 2, function(x) which(order(x, decreasing=TRUE)==1))] #Clase con la mayor probabilidad

		nextProbIndex 	<- cbind(apply(globalQueryResult$probabilities[,which(globalQueryResult$class == "NotAssigned"), drop=F], 2, function(x) which(order(x, decreasing=TRUE)==2)), which(globalQueryResult$class == "NotAssigned"))
		nextProb 		<- apply(nextProbIndex, 1, function(x) round(globalQueryResult$probabilities[x[1],x[2]], numDecimals))
		nextClass 		<- rownames(globalQueryResult$probabilities)[apply(globalQueryResult$probabilities[,which(globalQueryResult$class == "NotAssigned"), drop=F], 2, function(x) which(order(x, decreasing=TRUE)==2))]

		notAssignedSamples <- cbind(highestProbClass=highestProbClass, as.data.frame(highestProb), nextProb=nextProb, nextClass=nextClass)
	}
	callRate<-round(100*(1-(notAssigned/sum(stats[,1],notAssigned))),numDecimals)
	
	# Verbose
	if(verbose)
	{
		message(paste("The query contains ", samplesQueried=numSamples, " samples. ",sum(stats[,1])," were assigned to a class resulting on a call rate of ", callRate,"%. \n", sep=""))
		flush.console()
	}
	
	# Returns
	ret<- list(callRate=callRate, assigned=stats)
	if(showNotAssignedSamples) ret <- c(ret, notAssigned=list(notAssignedSamples))
	
	return(ret)
}

############################
#  Plots the assignment probabilities
############################
# Changes in version 1.2:
# function plotPoints no longer required. Optimized speed.
# New argument: pointSize 
# Background color: Blue
# Cairo
plotAssignments <- function(queryResult, realLabels, minProbAssignCoeff=1, minDiffAssignCoeff=0.8, totalNumberOfClasses=NULL, pointSize=0.8)							
{
	# Comprobar y concatenar si hay varias queries
	queryResult <- queryResultCheck(queryResult)
	
	# Comprobar  los argumentos
	if(!is.factor(realLabels)) { warning("The second argument (real labels) had to be converted into a factor.", immediate. = TRUE)}	
	realLabels <- factor (realLabels)
	if(is.null(names(realLabels))) 
	{
		if (length(realLabels) == length(names(queryResult$class))) 
		{	
			names(realLabels) <- names(queryResult$class)
		} else stop("The sample's labels vector is not named.")
		warning("The real sample's labels  vector is not named, it will be assumed the labels are in order: the first label applies to the first predicted sample... ", immediate. = TRUE)
	}
	if(is.null(totalNumberOfClasses)) {numClasses <- length(levels(realLabels))
	}else{
		if(!is.numeric(totalNumberOfClasses) || (totalNumberOfClasses<length(levels(realLabels)))) stop ("totalNumberOfClasses should be the number of classes for which the classifier was originaly trained.")
		numClasses <- totalNumberOfClasses
	}
	if(!is.numeric(minDiffAssignCoeff)) stop("'minDiffAssignCoeff' should be a coefficient to modify the required difference between probabilites to assign the sample to a class.")
	if(minProbAssignCoeff<0 || ((numClasses != 2) &&(minProbAssignCoeff>(numClasses/2)))) stop("'minProbAssignCoeff' should be between 0 and half of the number of classes.")
	if(minDiffAssignCoeff<0 || minDiffAssignCoeff>numClasses) stop("'minDiffAssignCoeff' should be between 0 and the number of classes.")
	if(!is.numeric(pointSize)) stop("'pointSize' should be numeric.")
	
	# Settings:
	abLineColor <- "midnightblue"
	bgColor <- "aliceblue"
	correctColor <- "#3F9728"
	incorrectColor <- "red"
	minX <- 0.2
	
	# Calculate asignments
	rand <- 1/numClasses
	minProb <- 2*rand * minProbAssignCoeff
	minDiff <- rand * minDiffAssignCoeff
	
	# Prepare canvas
	plot(c(minX,1), c(0,1), type="n", xlab="Probability of the most likely class", ylab="Difference with next class", frame=FALSE, main="Thresholds to assign query samples")
	if(numClasses>2) 
	{ 
		rect(minProb,minDiff,1,1,  col=bgColor, border=bgColor) 
		abline(v=minProb, h=minDiff, col=abLineColor, lty="dashed") 
	} else
	{	
		rect(0, minDiff,1,1,  col=bgColor, border=bgColor) 
		abline(h=minDiff, col=abLineColor, lty="dashed") 
		axis(1)
		axis(2)
	}
	
	text(0.7, 0.95, labels="Assigned", col=abLineColor, cex=0.8)
	if(numClasses>2)  text(0.3, 0.95, labels="Not Assigned", col="#606362", cex=0.8)
	text(minProb-0.03, minDiff + 0.01, labels="minProb", col=abLineColor, srt = 90, pos=4, cex=0.8)
	text(minX+0.09, minDiff+0.01, labels="minDiff", col=abLineColor, pos=1, cex=0.8)
	legend("bottomright", "(x,y)", legend=c("Correct", "Incorrect"), title = "Most likely class", text.width = strwidth("1,000,000"),  xjust = 1, yjust = 1, lty = 0, pch=16, col=c(correctColor, incorrectColor), cex=0.8)
	
	# Plot probabilities
	realLabels <- as.character(realLabels[colnames(queryResult$probabilities)])	
	prob <- apply(queryResult$probabilities, 2, function(x)
					{														
						largest <- which(x == max(x))
						nextProb <- max(x[-largest])
						class <- names(x)[largest]
				
				
						return(c(biggestProb=x[largest], nextProb=nextProb, assignedClass=class))
					})
	rownames(prob) <- c("biggestProb", "nextProb", "assignedClass")
	prob <- rbind(prob, realLabels=realLabels)	
				
	correct <- which(prob["assignedClass",] == prob["realLabels",])
	incorrect <- which(prob["assignedClass",] != prob["realLabels",])
	
	biggestProb <- as.numeric(prob["biggestProb",])
	nextProb <- biggestProb - as.numeric(prob["nextProb",])
	points(biggestProb[correct], nextProb[correct], col=correctColor, pch=16, cex=pointSize)
	points(biggestProb[incorrect], nextProb[incorrect], col=incorrectColor, pch=16, cex=pointSize)								

	if(!names(dev.cur()) %in% c("pdf", "Cairo")) 
	{
		print("To identify a sample on the plot click on it. Press ESC or  right-click on the plot screen to finish.")
		id <- identify(coordinates, labels=paste(rownames(coordinates)," (",prob["realLabels", rownames(coordinates)], ")", sep=""))
	}
}


## Calculates the genes ranking and plots the significant genes  
#####
# Options: 	- Arguments:	 a genesRanking
#          					 an eset+sampleLabels (to calculate genesRanking)
#		   	- plotLp=FALSE don't plot
#		  	- returnRanking: (FALSE/NULL, "full", "lp"/"significant"/"lpThreshold"/TRUE)
calculateGenesRanking <- function(eset=NULL, sampleLabels=NULL, numGenesPlot=1000, plotTitle="Significant genes", plotLp=TRUE, lpThreshold = 0.75, numSignificantGenesType="ranked", returnRanking="full", nullHiphothesisFilter=0.95,  nGenesExprDiff=1000, geneLabels=NULL,  precalcGenesRanking=NULL, IQRfilterPercentage= 0, verbose=TRUE)
{ 
	###################################
	# Checking arguments
	###################################
	genesRanking <- precalcGenesRanking
	if(!is.null(precalcGenesRanking)) returnRanking <- FALSE
	if(is.null(genesRanking) && (is.null(eset) && is.null(sampleLabels))) stop("Please, provide either a precalculated genesRanking, or an expressionSet and its sample labels to calculate it")
	if(!is.null(genesRanking) && (!is.null(eset) || !is.null(sampleLabels))) warning("Since a genesRanking was provided, eset and sampleLabels are not needed. They will be ignored.", immedite.=TRUE)
	
	# Arguments always required:
	if(!is.numeric(numGenesPlot)) stop("The argument numGenesPlot should be a number.")
	if(!is.numeric(lpThreshold) || (lpThreshold>=1 || lpThreshold <0)) stop("The threshold should be a percentage (a number between 0 and 1).")
	if(!is.character(numSignificantGenesType)) numSignificantGenesType <- "ranked"
	if(!numSignificantGenesType %in% c("ranked", "global"))  numSignificantGenesType <- "ranked"
	if(!is.numeric(nullHiphothesisFilter) || (nullHiphothesisFilter>1 || nullHiphothesisFilter <0)) stop("The nullHiphothesisFilter threshold should be a percentage (a number between 0 and 1).")
	if(!is.numeric(nGenesExprDiff)) stop("nGenesExprDiff should be a number.")
	if(!is.character(plotTitle)) stop("The title is not of type 'character'.")
	if(!is.logical(plotLp)) plotLp <- TRUE
	if(!is.logical(verbose)) verbose <- TRUE
	
	# sampleLabels format: 
	if(is.character(sampleLabels) && (length(sampleLabels) ==1))
	{
		if( is(eset, "ExpressionSet")  && (sampleLabels %in% colnames(pData(eset)))) {
		sampleLabels <- pData(eset)[,sampleLabels, drop=F]
		}	else{
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
	if(class(sampleLabels) != "factor") { 
		#warning("The argument 'classification sampleLabels' had to be converted into a factor.", immediate. = TRUE)
	}	
	sampleLabels <- factor (sampleLabels) #Just in case there are not samples of all the original labels	
	
	# If genesRanking was not provided, calculate:
	if(is.null(genesRanking)) 
	{
		# Check eset
		if(is(eset, "ExpressionSet")) eset <- exprs(eset) else if (!is.matrix(eset)) stop("The argument 'eset' should be an expression matrix or an ExpressionSet.")
		# Check sampleLabels
		if(class(sampleLabels) != "factor") { warning("The argument 'sampleLabels' had to be converted into a factor.", immediate. = TRUE)}	
		sampleLabels <- factor (sampleLabels) #Just in case there aren't samples of all the original labels	
		if(dim(eset)[2] != length(sampleLabels)) stop("The number of labels does not match the number of samples.")
		if(!is.null(names(sampleLabels)))
		{
			if(sum(!names(sampleLabels) %in% colnames(eset))>0 ) stop("The names of the labels do not match the samples.")
		}else 
		{
			names(sampleLabels)<-colnames(eset)
			warning("The data labels vector is not named, it will be assumed the labels are in order: the first label applies to the first sample... ", immediate. = TRUE)
			#print(sampleLabels)
		}	
		# Check filter pecentage
		if(!is.numeric(IQRfilterPercentage) || (IQRfilterPercentage>=1 || IQRfilterPercentage <0)) stop("The filter percentage should be a probability (a number between 0 and 1).")
							
		# IQR filter
		giqrs <- iqr.filter(eset, percentage=IQRfilterPercentage)														
		esetFiltered <- eset[giqrs,]
	
		# Calculate the genes ranking:
		genesRanking <- PEB(esetFiltered, sampleLabels, nullHiphothesisFilter=nullHiphothesisFilter)							
		
		# Add geneLabels
		if(!is.null(geneLabels)) geneLabels <- extractGeneLabels(geneLabels, rownames(esetFiltered))
		if(!is.null(geneLabels)) genesRanking <- setProperties(genesRanking, geneLabels=geneLabels)
		
		# Add expression
		if((!is.null(returnRanking))&&(nGenesExprDiff>0))
		{
				numClasses <- length(gClasses(genesRanking))
				topGenes <- as.vector(getRanking(getTopRanking(genesRanking, nGenesExprDiff), showGeneID=T, showGeneLabels=F)$geneID)
				topGenes <- topGenes[which(!is.na(topGenes))]
				meanExprDiff<- difMean(esetFiltered[topGenes,], numClasses, (dim(esetFiltered)[2]/numClasses))
				colnames(meanExprDiff) <-  gClasses(genesRanking)
				genesRanking <- setProperties(genesRanking, meanDif=meanExprDiff)
		}
	}								
	
	# Calculate genes over lpThreshold
	lp <- numSignificantGenes(genesRanking, lpThreshold=lpThreshold, numSignificantGenesType=numSignificantGenesType)
	

	if(plotLp)
	{
		# Extract from the genesRanking the required postProb and ord (only the required numGenesPlot...)
		if (numGenesPlot > dim(genesRanking@ord)[1]) numGenesPlot <- dim(genesRanking@ord)[1]
		if (length(gClasses(genesRanking)) > 2){ ord <- genesRanking@ord[1:numGenesPlot,]  
		}else ord <- cbind(NULL,genesRanking@ord[1:numGenesPlot] ) #Para que sea una matriz de una columna
		postProb <- matrix(nrow=numGenesPlot, ncol=ncol(ord))   
		for(i in 1:ncol(ord))  {
		postProb[,i] <- genesRanking@postProb[ord[,i],i+1]  
		}

		numClasses <- ncol(postProb) # If there are only 2 classes, postProb only has 1 column
		plot(postProb[,1], type="n", ylab="Posterior Probability", xlab="Gene Rank", xlim=c(1,numGenesPlot), ylim=c(0,1))
		if((numClasses>3 && numClasses<10) && library(RColorBrewer,logical.return=TRUE))
		{
			cols<- brewer.pal(numClasses,"Set1")	
		}else cols <- rainbow(numClasses)
		pchs <- c(15:(15+ncol(postProb)))

		if(numGenesPlot < 200) {interval <- 1:numGenesPlot
		} else {
			interval <- seq(0, numGenesPlot, numGenesPlot/20)
			interval[1]<-1		
		}
		
		for(i in 1:numClasses)
		{
			lines(interval,postProb[interval,i], type="b", col=cols[i], pch=pchs[i]) 
		}

		title(plotTitle)
		abline(h=lpThreshold)
		text(numGenesPlot-(numGenesPlot/10), lpThreshold-0.02, paste("Threshold=",lpThreshold, sep=""), col="grey50", cex=0.8)
		legend("bottomleft", paste( gClasses(genesRanking)," (",lp," genes)",sep=""), lty=1, col=cols,  pch=pchs)
	}
	
	if(!is.null(returnRanking) && returnRanking!=FALSE)
	{
		if (tolower(returnRanking) == "full")	return (genesRanking)
		
		# Leave in ranking only significant genes (genes over lpThreshold)
		if ((returnRanking==TRUE || tolower(returnRanking)=="lpthreshold") || (tolower(returnRanking)=="lp" || tolower(returnRanking)=="significant"))
		{
			lpRanking <- getTopRanking(genesRanking, lp)
			return(lpRanking)
		}
	}
}


## Perfiles de expresion de los genes del clasificador en todas las muestras
# Saves the plots in a PDF file of name fileName
# Example: plotExpressionProfiles (basicLeukemias[,trainSamples], getRanking(getTopRanking(genesRankingGlobal, 5))$geneID, fileName="expresion.pdf", sampleLabels=basicLeukemias$LeukemiaType[trainSamples])
plotExpressionProfiles <- function(eset, genes=NULL, fileName=NULL, geneLabels=NULL, sampleLabels=NULL, sampleColors=NULL, labelsOrder=NULL, showSampleNames=FALSE, showMean= FALSE, sameScale=TRUE, verbose=TRUE)
{
	printText<-""
	
	#########################
	#Comprobacion de parametros
	#########################
		# sampleLabels format: 
	if(is.character(sampleLabels) && (length(sampleLabels) ==1))
	{
		if( is(eset, "ExpressionSet")  && (sampleLabels %in% colnames(pData(eset)))) {
		sampleLabels <- pData(eset)[,sampleLabels, drop=F]
		}	else{
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
	if(!is.null(sampleLabels)) {
		if(class(sampleLabels) != "factor") { warning("The argument 'sampleLabels' had to be converted into a factor.", immediate. = TRUE)}	
		sampleLabels <- factor (sampleLabels) #Just in case there are not samples of all the original labels	
	}
	
	# eset
	exprMatrix <- eset
	if(is(exprMatrix, "ExpressionSet")) exprMatrix <- exprs(exprMatrix) else if (!is.matrix(exprMatrix)) stop("The first argument should be an expression matrix or an ExpressionSet.")
	if(is.null(printText)) printText<-""
	
	# Other
	genesRanking <- NULL
	if(is.null(genes)) {
		genes <- rownames(exprMatrix)
		if(length(genes)>1000 && !is.null(fileName))  warning(paste("Plotting the expression profiles of the ",length(genes)," genes in the expression set.", sep=""), immediate. = TRUE)
	}
	if(is(genes, "GeNetClassifierReturn") && "classificationGenes" %in% names(genes)) {
		genes <- genes@classificationGenes
		warning("Plotting expression profiles of the classification genes. To plot other genes, set i.e. genes=...@genesRanking")
	}
	if(is(genes, "GenesRanking"))
	{
		if(length(genes@geneLabels) > 0 &&  any(!is.na(genes@geneLabels))) geneLabels <- genes@geneLabels
		genesRanking <- genes
		genes <- getRanking(genes, showGeneLabels=F, showGeneID=T)$geneID		
	}
	if(!is.matrix(genes) && !is.vector(genes)) stop ("The genes list should be either a vector or a matrix.")
	if(sum(!genes[which(genes!="NA")] %in% rownames(exprMatrix))!=0) stop ("The expression matrix doesn't contain all the genes.")
	genesVector <- unique(as.vector(genes))
	genesVector <- genesVector[!is.na(genesVector)]
	if(!is.null(geneLabels)) geneLabels<-extractGeneLabels(geneLabels, rownames(exprMatrix[genesVector,]))
    if(!is.null(fileName) && !is.character(fileName)) stop("The file name is not valid.")
	if(!is.null(fileName) && regexpr(".pdf", fileName) == -1) fileName <- paste(fileName, ".pdf", sep="") 
	numSamples <- dim(exprMatrix)[2] 
	if(!is.null(sampleLabels))
	{
		if(class(sampleLabels) != "factor") { warning("The argument 'sampleLabels' had to be converted into a factor.")}
		sampleLabels <- factor (sampleLabels) #Just in case there aren't samples of all the original labels	
		if(numSamples != length(sampleLabels)) stop("The number of labels doesn't match the number of samples.")  
		if(!is.null(names(sampleLabels)))
		{
			if(sum(!names(sampleLabels) %in% colnames(exprMatrix))>0 ) stop("The names of the labels do not match the samples.")
		}else 
		{
			names(sampleLabels)<-colnames(exprMatrix)
			warning("The data labels vector is not named, it will be assumed the labels are in order: the first label applies to the first sample... ", immediate. = TRUE)
			#print(sampleLabels)
		}	
	}
	if(!is.null(labelsOrder))
	{
		if(any(!labelsOrder %in% levels(sampleLabels)) || any(!levels(sampleLabels) %in% labelsOrder)) 
		{
			warning("The labelsOrder doesn't match the samples labels. It will be ignored.")
			labelsOrder <- NULL
		}
	}
	if(is.null(sampleColors)) sampleColors <- "red"
	if((length(sampleColors) > 1) && (numSamples != length(sampleColors))) warning("The number of sampleColors doesn't match the number of samples.")  
	

	if(!is.logical(showMean)) showMean <- TRUE
	if(!is.logical(sameScale)) showMean <- TRUE
	if(!is.logical(verbose)) verbose <- TRUE

	if(!is.null(sampleLabels))
	{
		classes <- levels(sampleLabels)
		if(!is.null(labelsOrder)) classes <- labelsOrder
		numClasses <- length(classes)
		matriz <- NULL
		indexes<- NULL
		for(i in 1:numClasses) #Por si no estan agrupados
		{
			indexes <- c(indexes,  which(sampleLabels==classes[i]))
		}
		matriz <- exprMatrix[genesVector, indexes, drop=F]
		sampleLabels <- sampleLabels[indexes]
		if(length(sampleColors) > 1) sampleColors <- sampleColors[indexes]
	}else
	{
		classes <- colnames(genes)			
		numClasses <- length(classes)		
		matriz <- exprMatrix[genesVector,, drop=F]
		if(is.null(rownames(matriz))) rownames(matriz)<- genesVector   # When there is only one gene: no names
	}

	#########################
	# Plot
	#########################	
	# Configure output (pdf or split window?)
	numGenesPlot <- configurePlotOutput(length(genesVector), fileName)
	if(numGenesPlot != length(genesVector))
	{
		warning(paste("Up to ",numGenesPlot," genes will be shown. To plot more genes specify a PDF output file name.",sep=""))
		if( numClasses == 0 || !is.matrix(genes) ) { 
			genesVector <-  genesVector[1:numGenesPlot]
		}
		else {
			numRows <- numGenesPlot/dim(genes)[2]
			while ((numRows < dim(genes)[1]) && (length(which(!is.na(genes[1:numRows,, drop=F]))) < numGenesPlot))
			{
					numRows <- numRows + 1
			}
			genesVector <- as.vector(genes[1:numRows,, drop=F])
			genesVector <-  genesVector[which(!is.na(genesVector))[1:numGenesPlot]]
		}
	}

	# Plot
	if(sameScale) 
	{
		minEset <- min(0, matriz)
		maxEset <- max(matriz)
		ylim <- c(minEset - (minEset*0.05), maxEset +(maxEset*0.05))
	}
	for(i in 1:numGenesPlot) 	# No length(genesVector): Puede haber mas genes q espacio
	{
		if(!is.na(genesVector[i]))
		# {
			# plot.new()
		# } else 
		{
			if(!sameScale) 
			{
				minEset <- min(0, matriz[genesVector[i],])
				maxEset <- max(0, matriz[genesVector[i],])
				ylim <- c(minEset - (minEset*0.05), maxEset +(maxEset*0.05))
			}
			plot(matriz[genesVector[i],], type="h", col=sampleColors, lwd=2, ylim=ylim, xlab="Sample index", ylab="Expression values")
			# if(dim(matriz)[2] < 20)
			# {	
				# points(matriz[genesVector[i],], type="l", col="grey", lwd=1)
			# }
			
			if(!is.null(colnames(genes))) {  geneClass <- unique(colnames(genes)[which(genes == genesVector[i],arr.ind=TRUE)[,2]])
			} else	{	geneClass<-"" }					#classes[which(genes == genesVector[i], arr.ind=TRUE)[2]]
			if(!is.null(geneLabels[genesVector[i]]) && !is.na(geneLabels[genesVector[i]])) { geneName<- paste(genesVector[i]," (" , geneLabels[genesVector[i]],")" ,sep="")
			}else geneName<-  genesVector[i]
			tit <- paste(geneClass[1], geneName, sep="\n" )
			title(tit, sub=eval(parse(text=printText)))
			
			if(showSampleNames) for (nSample in 1:length(colnames(matriz))) {text(nSample,ylim[2]-(ylim[2]*0.1),colnames(matriz)[nSample], pos=1, srt = 90, cex=0.5, col="grey")} 
			
			if(!is.null(sampleLabels))
			{
				prevLim <- 0
				for(j in 1:(numSamples))
				{
					if (is.na((sampleLabels[j] != sampleLabels[j+1]) ) || sampleLabels[j] != sampleLabels[j+1]) 
					{
						if(!is.na((sampleLabels[j] != sampleLabels[j+1]) ) ) abline(v=j+0.5, col="black")
						if (any(nchar(classes)>10) ) {text(prevLim+((j-prevLim)/2), ylim[2]-(ylim[2]*0.04), labels=paste("C",which(classes==sampleLabels[j]),sep=""))
						}else  text(prevLim+((j-prevLim)/2)+0.5, y=ylim[2]-(ylim[2]*0.04), labels=paste(sampleLabels[j],sep=""), pos=3) 
						if(showMean)
						{
							classMean <- mean(matriz[genesVector[i], (prevLim+1):j])
							lines(c(prevLim+1, j), c(classMean, classMean) , col="grey")
						}
						prevLim <- j
					}
				}
			}
		}
	}	

	if (!is.null(fileName))
	{
		dev.off()
		if (verbose){ message(paste("The plot was saved as ",getwd(),"/",fileName," (PDF file)",sep="")); flush.console()} 
	}	else 
	{
	
		if(i==1)
		{
			geneExprs <-t(matriz[genesVector[i],, drop=F])
			if(names(dev.cur())!="pdf") print("To identify a sample on the plot click on it. Press ESC or  right-click on the plot screen to finish.")
			identify(geneExprs, labels=rownames(geneExprs))
		}
	}	
}


## Representamos los valores escalados de los genes en los vectores soporte que sirven de limite entre las clases
## tenemos en cuenta los valores de los coeficientes de lagrange para cada una de las SV.
# El parametro correctedAlpha para corregir o no los valores teniendo en cuenta los alphai
# discriminant.power.plot(classifier, classificationGenes, classNames=c("ALL","AML","CLL","CML","NoLeu"), fileName="test.pdf",correctedAlpha=TRUE)
# discriminant.power.plot(classifier, classificationGenes, fileName="test.pdf")
# Classification genes: Genes por columnas (nombrecolumna= clase)
# discriminant.power.plot(classifier, colnames(classif$SV), fileName="test.pdf")
#> geneLabels
# ENSG00000164398 ENSG00000169575 ENSG00000143153 
#         "gen1"          "gen2"          "gen3" 
# games(geneLabels)       >> [1] "ENSG00000164398" "ENSG00000169575" "ENSG00000143153"
# 
# classifier: puede ser un svm o el objeto devuelto por la funcion principal
# classificationGenes: puede ser un c(), una matriz o un GenesRanking

plotDiscriminantPower <- function(classifier, classificationGenes=NULL , geneLabels=NULL, classNames=NULL, plotDP = TRUE, fileName= NULL, returnTable=FALSE, verbose=TRUE)
{	
	correctedAlpha<-FALSE  #Eliminado temporalmente
		
	################################
	# Check arguments
	################################
	if(!is.logical(plotDP)) plotDP <- TRUE
	if(!is.logical(returnTable)) stop ("returnTable should be either TRUE or FALSE.")
	if(!is.logical(verbose)) verbose <- TRUE
	if(!is.null(fileName) && !is.character(fileName)) stop("The file name is not valid.")
	if(!is.null(fileName) &&regexpr(".pdf", fileName) == -1) fileName <- paste(fileName, ".pdf", sep="") 
		
	# Classifier
	if(is(classifier, "GeNetClassifierReturn")){
		if("classificationGenes" %in% names(classifier))
		{
			if(is.null(classificationGenes))
			{
				classificationGenes <- classifier@classificationGenes
			}else {
				if(is.null(geneLabels) && is.character(classificationGenes))
				{
					if(length(classifier@classificationGenes@geneLabels) > 0 &&  any(!is.na(classifier@classificationGenes@geneLabels[classificationGenes])))				
						geneLabels <- classifier@classificationGenes@geneLabels[classificationGenes]
				}
			}
		}
		if("classifier" %in% names(classifier)) {classifier <- classifier@classifier$SVMclassifier
		}else stop("'classifier' doesn't contain a trained classifier.")
	}
	if(is.list(classifier) && ("SVMclassifier" %in% names(classifier))) classifier <- classifier$SVMclassifier
	if(!is(classifier,"svm")) stop("The first argument should be a svm classifier or the object returned by geNetClassifier.")
	
	# ClassificationGenes (GenesRanking)
	if(class(classificationGenes) == "GenesRanking") 
	{	
		if(sum(numGenes(classificationGenes)))
		{
			if(is.null(geneLabels) && (length(classificationGenes@geneLabels) > 0 &&  any(!is.na(classificationGenes@geneLabels)))) geneLabels <- classificationGenes@geneLabels
			classificationGenesRanking <- classificationGenes
			classificationGenes <- getRanking(classificationGenes, showGeneLabels=FALSE, showGeneID=TRUE)$geneID 
		}else {
			classificationGenes <- NULL
			classificationGenesRanking<-NULL
		}
	}else{
		classificationGenesRanking<-NULL
	}
	
	# If classificationGenes is not provided/valid, use the classifier's SV
	missingGenes <- !as.vector(classificationGenes[!is.na(classificationGenes)]) %in% colnames(classifier$SV)
	if(any(missingGenes)) 
	{
		missingGenes <- as.vector(classificationGenes[!is.na(classificationGenes)])[which(missingGenes)]
		# Transformed into dataframe colnames to make sure they match the $SV
			m <- matrix(ncol=length(missingGenes))	
			colnames(m) <- missingGenes
			m <- data.frame(m)
			missingGenes <- colnames(m)
		missingGenes <- missingGenes[which(!missingGenes %in% colnames(classifier$SV))]
		
		classificationGenes[which(classificationGenes %in% missingGenes)] <- NA
		
		 if(all(is.na(classificationGenes))) stop("The given 'classificationGenes' are not used by the classifier. Their Discriminant Power cannot be calculated.")
		if(length(missingGenes)>0) warning(paste("The following classificationGenes are not used by the classifier. Their Discriminant Power cannot be calculated: ", missingGenes, sep=""))
	}
	if(is.null(classificationGenes)) classificationGenes <- colnames(classifier$SV)

	# geneLabels
	if(is.matrix(geneLabels) || is.data.frame(geneLabels))
	{	
		if(dim(geneLabels)[1] != 1 && dim(geneLabels)[2] != 1)
		{ 
			stop("geneLabels should be a named vector or a one dimensional matrix.")
		}else 
		{
				if(dim(geneLabels)[2] == 1) 
				{ 
					geneLabels <- geneLabels[,1]
				}	else if(dim(geneLabels)[1] == 1) geneLabels <- geneLabels[1,]		
		}
	}
	if(is.factor(geneLabels))  
	{
		tmp <- names(geneLabels)
		geneLabels <- as.character(geneLabels)
		names(geneLabels) <- tmp
	}
	if(!is.null(geneLabels) && !is.character(geneLabels)) stop("geneLabels should be a vector containing the gene symbol.")
	if(!is.null(geneLabels) && is.null(names(geneLabels))) stop("names(geneLabels) can't be empty. It should contain the names used in classification genes.")
	#if(!is.null(geneLabels) && sum(!names(geneLabels) %in% classificationGenes[which(classificationGenes!="NA")] )>0) warning("Some geneLabels will not  be used.")
	if(!is.null(geneLabels) && sum(!classificationGenes[which(classificationGenes!="NA")] %in% names(geneLabels))>0) warning("geneLabels doesn't contain the symbol for all the classification genes.")
	
	# classNames
	if(is.null(classNames) || (length(classifier$levels) != length(classNames))){
		if (!is.null(classNames) && length(classifier$levels) != length(classNames))  {
			warning(paste("The number of classes provided don't match the classifier's. The default class names will be used instead.",sep=""), immediate. = TRUE)	}
		classNames <-classifier$levels 		
	}
	longClassNames <- any(nchar(classNames)>6)
	if(longClassNames)
	{
		for( i in 1:length(classNames)) #Add "C1:..."
		{
			if (nchar(classNames[i])>10 ) classNames[i] <- paste(substr(classNames[i] ,1,10), "...",sep="")	
			classNames[i] <- paste("C", i, ": ", classNames[i], sep="")
		}
	}
	numClasses <- ifelse(is.matrix(classificationGenes), length(classNames), length(classifier$levels))

	if(!is.matrix(classificationGenes)) { #if(verbose) warning("The 'classification genes' are not sorted by colums and classes, the gene class will not be shown .")
	}else
	{
			if(length(classifier$levels) == 2) {
				if (dim(classificationGenes)[2] != 1)  stop("The classes of the classifier and the classification genes provided don't match.")
			}else {	
				if(sum(!colnames(classificationGenes) %in% classifier$levels)>0) stop("The classes of the classifier provided and the classification genes don't match.")
			}			
			classificationGenes <- classificationGenes[,apply(classificationGenes, 2, function(x) !all(is.na(x)))] # Is there any class without genes?
	}
	
	nGenes <- length(classificationGenes[which(!is.na(classificationGenes))]) #sum(numGenes(classifier$classificationGenes))
	if(nGenes>dim(classifier$SV)[2]){ warning(paste("The given number of genes is bigger than the classifier's.",sep=""), immediate. = TRUE)}
			
		
	################################
	# Calculate discriminant power 
	################################
	discrPwList<-NULL
	if(is.matrix(classificationGenes))	# If it contains the genes by classes (columns)
	{
		for(cl in 1:dim(classificationGenes)[2])
		{
			discrPwList <- c(discrPwList, discrPwList=list(sapply(as.character(classificationGenes[which(classificationGenes[,cl]!="NA"),cl]), function(x) SV.dif(classifier, x, correctedAlpha=correctedAlpha))))
			names(discrPwList)[cl] <- colnames(classificationGenes)[cl]
		}
	}else
	{
		classificationGenes <- classificationGenes[which(!is.na(classificationGenes))]
		discrPwList <- list(sapply(as.character(classificationGenes), function(x) SV.dif(classifier, x, correctedAlpha=correctedAlpha)))
		classificationGenes <- as.matrix(classificationGenes)
	}


	###### Configure output
	numGenesPlot <- configurePlotOutput(nGenes, fileName)
	if(numGenesPlot != nGenes)
	{		
		warning(paste("Up to ",numGenesPlot ," genes will be shown. To plot more genes specify a PDF output file name.",sep=""))	

		numRows <- numGenesPlot/dim(classificationGenes)[2]
		while ((numRows < dim(classificationGenes)[1]) && (length(classificationGenes[which(!is.na(classificationGenes[1:numRows,, drop=F]))]) < numGenesPlot))
		{
				numRows <- numRows + 1
		}
				
		if(length(classificationGenes[which(!is.na(classificationGenes[1:numRows,, drop=F]))]) <= numGenesPlot)
		{
			classificationGenes <- classificationGenes[1:numRows,, drop=F]
		} else {
			classificationGenes <- classificationGenes[1:(numRows-1),, drop=F]
		}
	}

	################################
	# Plot
	################################
	if(plotDP)
	{
		supLim <- max(sapply(discrPwList,function(y) max(sapply(y["positive",], function(x) max(apply(x,2,sum))))))
		infLim <- min(sapply(discrPwList,function(y) min(sapply(y["negative",], function(x) min(apply(x,2,sum))))))
		lims <- round(max(supLim ,abs(infLim)),1)
		if (lims > 2.5) {lims <- ceiling(lims/3)*3
		}else 			lims <- max(lims, round((lims+0.1)/3,1)*3)
		lims <- c(-lims,lims)

		mycols <- colorRampPalette(c("blue","white"))(max(classifier$nSV+2))	
		
		for(c in 1:dim(classificationGenes)[2]) #numClasses
		{
			for(g in 1:dim(classificationGenes)[1]) 
			{
				# Get Data
				gene <- classificationGenes[g,c]
				if(!is.na(gene))
				 # {
					# if(is.null(fileName)) plot.new()
				 # } else 
				{
					if(!is.null(geneLabels[gene]) && !is.na(geneLabels[gene])) geneName<- paste(gene," (" ,geneLabels[gene],")" ,sep="")
					else geneName<-gene
					
					geneClass<-NULL
					if(!is.null( names(gene)))
					{
						geneClass <- names(gene)
						if(nchar(geneClass)>70) geneClass<- substr(geneClass,1,70)
					}
					
					pos <- discrPwList[[c]]["positive",][[gene]]
					neg <- discrPwList[[c]]["negative",][[gene]]
					
					if(is.null(pos) || is.null(neg))
					{
							if(length(discrPwList[[c]]["positive",]) ==1) pos<-discrPwList[[c]]["positive",][[1]]
							if(length(discrPwList[[c]]["negative",]) ==1) neg<-discrPwList[[c]]["negative",][[1]]
					}
				
					# Seg graph parameters
					tit<- paste(geneClass,"\n", geneName, "\n", sep="")
					if(longClassNames){ 	tit<- paste(tit,"DP: ", sep="")
					} else {							tit<- paste(tit,"Discriminant power: ", sep="") }			
					tit<- paste(tit, abs(round(discrPwList[[c]][,gene]$discriminantPower,2)), " (", discrPwList[[c]][,gene]$discrPwClass, ")", sep="")
					barplot(rep(0,numClasses),add=FALSE, ylim=lims, main=tit, col=mycols, width=0.9, space=0.1, cex.main=1)
					
					# Draw ab lines
					sum_pos <- apply(pos,2,sum)
					sum_neg <- apply(neg,2,sum)
					if (discrPwList[[c]][,gene]$discriminantPower > 0)
					{
						maxim <- which(sum_pos == max(sum_pos), arr.ind=TRUE)[1]
						if (max(sum_pos[-maxim]) != 0) 		sig <- max(sum_pos[-maxim]) 	#Si el siguiente valor no es cero
						else 								sig <- max(sum_neg[-maxim])
						abline(h=c(sum_pos[maxim], sig), col="red", lty="dashed")
					}else
					{
						minim <- which(sum_neg == min(sum_neg), arr.ind=TRUE)[1]
						if (min(sum_neg[-minim]) !=0) 		sig <- min(sum_neg[-minim])
						else 								sig <- min(sum_pos[-minim]) 
						abline(h=c(sum_neg[minim], sig), col="red", lty="dashed")
					}
					
					# Draw bars
					barplot(pos,add=TRUE, col=mycols, width=0.9, space=0.1, names.arg=rep("",length(classNames)))
					barplot(neg,add=TRUE, col=mycols, width=0.9, space=0.1, names.arg=rep("",length(classNames)))
					abline(h=0, lwd=2)
					if(!correctedAlpha) text(seq(1, length(classNames), by=1)-0.5, par("usr")[3] - 0.2, labels = classNames, srt = 90, pos = 4, xpd = TRUE) 
					if(correctedAlpha) text(seq(1, length(classNames), by=1)-0.5, par("usr")[3], labels = classNames, srt = 90, pos = 4, xpd = TRUE) 
				}
			}	
		}
	}
	if (!is.null(fileName))
	{
		dev.off()
		if (verbose){ message(paste("The SV plot was saved as ",getwd(),"/",fileName," (PDF file)",sep="")); flush.console()} 
	}		
		
	################################
	# Prepare return table
	################################	
	if(returnTable)
	{
		discrPwDF <- NULL
		for(cl in 1:length(discrPwList))
		{
			discrPwDF<- rbind(discrPwDF, cbind(t(discrPwList[[cl]][c("discriminantPower","discrPwClass"),]), originalClass=rep(names(discrPwList)[cl],dim(discrPwList[[cl]])[2])))
		}
		# Transform into Data Frame
		discrPwDF <- data.frame(discrPwDF)
		discrPwDF[,"discriminantPower"] <- as.numeric(discrPwDF[,"discriminantPower"])
		discrPwDF[,"discrPwClass"] <- as.character(discrPwDF[,"discrPwClass"])
		
		# Order by Discriminant power
		tempDpMatrix <- discrPwDF
		tempDpMatrix[,"discriminantPower"] <- abs (tempDpMatrix[,"discriminantPower"] )
		discrPwDF <- NULL
		for(cl in classifier$levels)
		{
			clGenes <- which(tempDpMatrix[,"discrPwClass"]==cl)
			discrPwDF <- rbind(discrPwDF, tempDpMatrix[clGenes[order(as.numeric(tempDpMatrix[clGenes,"discriminantPower"]),decreasing=T)],])
		}
		
		# Return
		if(is.null(classificationGenesRanking)) 
		{
			return(discrPwDF)
		}else{ 
			# Merge the geneDetails(genesRanking) with the discriminant power into one matrix
			gDetails<-genesDetails(classificationGenesRanking)
			genesDetailsDF<-NULL
			for(cl in 1:length(discrPwList))
			{
				genesDetailsDF <- rbind(genesDetailsDF, gDetails[[cl]])
			}
			genesDetailsDF <- cbind(discrPwClass=rep(NA,dim(discrPwDF)[1]), discriminantPower=rep(NA,dim(discrPwDF)[1]), genesDetailsDF[rownames(discrPwDF),]) # "normal" cbind doesnt work
			genesDetailsDF[,"discrPwClass"] <- as.character(discrPwDF[,"discrPwClass"])
			genesDetailsDF[,"discriminantPower"] <- as.numeric(discrPwDF[,"discriminantPower"])
			
			return(genesDetailsDF)	
		}
	}
}


# plotType="dynamic" (each can be modified), plotType="static" (1 image divided into classes), plotType="pdf" 
# if geneLabels exists, it will use these labels, not the ones in the Ranking object
# genesInfo: Data.frame containing info about the genes. Can be replaced by classificationGenes or genesRanking (recommended).
# If classificationGenes + genesRanking:  
# classificationGenes: Tiene q ser un genesRanking

plotNetwork  <- function(genesNetwork, classificationGenes=NULL, genesRanking=NULL, genesInfo=NULL,geneLabels=NULL, returniGraphs=FALSE, plotType="dynamic", fileName=NULL, plotAllNodesNetwork=TRUE, plotOnlyConnectedNodesNetwork=FALSE,  plotClassifcationGenesNetwork=FALSE, labelSize=0.5, vertexSize=NULL, width=NULL, height=NULL, verbose=TRUE)
{
	layoutList <- NULL
	if(!library(igraph, logical.return=TRUE)) 
	{	
		warning("The function plotNetwork()  requires the packge igraph but it could not be loaded.")
		graphList<-NULL
	}
	else
	{
		showWarning <- FALSE
		genesInfoList <- NULL
		
		#####################
		# Check arguments & prepare variables
		# - Arguments: check
		# - Network: check
		# - Genes Info: check and merge
		# - Gene Labels: check or extract
		#####################
		
		if(!is.logical(returniGraphs)) stop("returniGraphs should be either TRUE or FALSE.")
		if(!is.logical(plotAllNodesNetwork)) stop("plotAllNodesNetwork should be either TRUE or FALSE.")
		if(!is.logical(plotOnlyConnectedNodesNetwork)) stop("plotOnlyConnectedNodesNetwork should be either TRUE or FALSE.")
		if(!is.logical(plotClassifcationGenesNetwork)) stop("plotClassifcationGenesNetwork should be either TRUE or FALSE.")
		if(!is.character(plotType)) { stop("plotType is not valid.")
		}else{
			if(!plotType %in% c("dynamic", "static", "pdf")) stop("plotType should be either 'dynamic', 'static' or 'pdf'.")
			if ((plotType == "dynamic") && !("tcltk" %in% rownames(installed.packages())))
			{
				warning("tcltk package is required for dynamic plots. A static network plot will be drawn instead.")
				plotType <- "static"
			}
		}
		if(!is.null(fileName) && !is.character(fileName)) stop("The file name is not valid.")
		if(!is.null(fileName)) 
		{
			plotType <- "pdf"
			if(regexpr(".pdf", fileName) == -1) fileName <- paste(fileName, ".pdf", sep="") 
		}
		if(is.null(fileName) && plotType=="pdf") fileName <- "genesNetwork.pdf"
		
		if(is.null(width))
		{
			if(plotType =="dynamic") width <- 800
			if(plotType =="pdf") 		 width <- 7
		}
		if(is.null(height))
		{
			if(plotType =="dynamic") height <- 500
			if(plotType =="pdf") 		 height <- 7
		}
		if(!is.null(vertexSize) && !is.numeric(vertexSize)) vertexSize <- NULL
		if(!is.logical(verbose)) verbose <- TRUE

		if(!(returniGraphs || plotAllNodesNetwork || plotOnlyConnectedNodesNetwork || plotClassifcationGenesNetwork)) stop("No network plots have been requested.")
			if(!(plotAllNodesNetwork || plotOnlyConnectedNodesNetwork || plotClassifcationGenesNetwork)) warning("No network plots have been requested, only the iGraph will be returned.") #(else)
			
		# Check NETWORK format
		if(is(genesNetwork, "GeNetClassifierReturn"))
		{	
			if(is.null(classificationGenes) && ("classificationGenes" %in% names(genesNetwork))) 	classificationGenes <- genesNetwork@classificationGenes
			if((is.null(genesRanking) && is.null(genesInfo)) && ("genesRanking" %in% names(genesNetwork))) 
			{
				nGenes <- max( 100, numGenes(genesNetwork@classificationGenes))
				genesRanking <- getTopRanking(genesNetwork@genesRanking, nGenes)
				warning(paste("Plotting up to ", max(numGenes(genesRanking)), " genes of each class.", sep=""))
			}
			if("genesNetwork" %in% names(genesNetwork)) {	
				if(!is.null(genesRanking)) genesNetwork <- getSubNetwork(genesNetwork@genesNetwork, genesRanking)
				else		genesNetwork <- genesNetwork@genesNetwork
			}else stop("'genesNetwork' is the return of geNetClassifier, but doesn't contain a genesNetwork.")
		}
		
		if(is.list(genesNetwork)) {
			if(any(sapply(genesNetwork, is.null))) {genesNetwork<-genesNetwork[-which(sapply(genesNetwork, is.null))] }
			nwClasses <- names(genesNetwork)
		} else
		{
			if(!class(genesNetwork) == "GenesNetwork") stop("genesNetwork should be either a list or a GenesNetwork.")
			if((sum(c("class1", "class2") %in% colnames(genesNetwork@edges)) == 2 ) && nrow(genesNetwork@edges)>0)
			{
				nwClasses <- unique(as.vector(genesNetwork@edges[,c("class1", "class2")]))
				# if (nwClasses[1] == nwClasses[2]) nwClasses <- nwClasses[1]
			}else 	nwClasses <- "geneClass"
			
			genesNetwork <- list(genesNetwork)
			names(genesNetwork) <- nwClasses[1]
		}
		ntwColnames<-c("gene1","gene2","relation","value")
		for(cl in names(genesNetwork))
		{
			if(!is.matrix(genesNetwork[[cl]]@edges)) stop("genesNetwork should be either a GenesNetwork or a list of GenesNetwork.")
			if(any(!ntwColnames %in% colnames(genesNetwork[[cl]]@edges))) stop("genesNetwork column names are not what expected.")
		}

		# Check classificationGenes and Genes ranking format and EXTRACT its genes INFO.
		if(is.matrix(classificationGenes) && nrow(classificationGenes)==0) classificationGenes<-NULL
		if(class(classificationGenes) == "GenesRanking" && numGenes(classificationGenes) == 0)  classificationGenes <- NULL
		if(plotClassifcationGenesNetwork && is.null(classificationGenes)) warning("The classifcation genes network can only be plotted if the classification genes are provided.")
		if((!is.null(classificationGenes) && !is.null(genesRanking)) && !is.null(genesInfo)) stop("Please, provide either 'genesInfo' OR a genesRanking and classificationGenes.")
		if(!is.null(genesRanking) || !is.null(classificationGenes))
		{
			if(!is.null(genesRanking))
			{
				if(class(genesRanking) != "GenesRanking") stop("genesRanking should be an object of type GenesRanking.")
				#genesInfo <- genesDetails(getTopRanking(genesRanking, nRankedGenesOverThreshold(genesRanking)))
				genesInfo <- genesDetails(genesRanking)
			}
			
			if(!is.null(classificationGenes))
			{
				if(class(classificationGenes) != "GenesRanking") stop("classificationGenes should be an object of type GenesRanking (the classificationGenes object returned by the classifier).")

				classificationGenesInfo <- genesDetails(classificationGenes)[nwClasses]
				clGenes <- lapply(classificationGenesInfo, rownames)			
				for( cl in names(genesNetwork))	if(any(!clGenes[[cl]] %in% getNodes(genesNetwork[[cl]]))) showWarning <- TRUE    #Cambiado 0.99
				if(showWarning) warning("Not all the classificationGenes are available in the genesNetwork. They will be represented, but there may be missing relationships.") # Or error?
				
				if(is.null(genesInfo))  # Create genes info
				{
					genesInfo <- classificationGenesInfo
				}else# Add to the existing genes info
				{ 	
					missingColumnsInGlobal <- colnames(classificationGenesInfo[[1]])[which(!colnames(classificationGenesInfo[[1]]) %in% colnames(genesInfo[[1]]))]
					if(length(missingColumnsInGlobal)>0)
					{
						for( cl in nwClasses)
						{
							# Add missing columns
							temp <- cbind(genesInfo[[cl]], matrix(nrow=nrow(genesInfo[[cl]]),ncol=length(missingColumnsInGlobal)))
							colnames (temp) <- c(colnames(genesInfo[[cl]]), missingColumnsInGlobal)
							
							for(tempCol in colnames(temp))
							{
								if((is.factor(temp[,tempCol]))) levels(temp[,tempCol]) <- unique(c(levels(temp[,tempCol]), levels(classificationGenesInfo[[cl]][,tempCol])))
							}

							temp[rownames(classificationGenesInfo[[cl]]), ] <- classificationGenesInfo[[cl]][,colnames(temp)]
							genesInfo[[cl]] <- temp	
						}
					}
				}
			}
		}
		
		# Check genesInfo format & split into list
		if(!is.null(genesInfo))
		{
			if(!is.data.frame(genesInfo) && !is.list(genesInfo)) stop("genesInfo should be a list or a data.frame.")
			if(is(genesInfo, "list"))
			{
				if(any(! names(genesNetwork) %in% names(genesInfo))) { stop("The class names in genesInfo and genesNetwork do not match.")
				} else genesInfo <- genesInfo[names(genesNetwork)]
				genesInfoList <- genesInfo
			}

			# Or split into if it was only one table
			if(is.null(genesInfoList))
			{
				for(cl in names(genesNetwork))
				{
					genesInfoList[[cl]] <- genesInfo[which(rownames(genesInfo) %in% getNodes(genesNetwork[[cl]])), ]
				}
			}
					
			# Check wether all the requested genes (in genesInfo) are available in the network
			for(cl in names(genesNetwork))
			{
				if(any(!rownames(genesInfoList[[cl]]) %in% getNodes(genesNetwork[[cl]]))) 
				{
					 if(!showWarning) warning("Not all the genes given in the Ranking or Info are available in the genesNetwork. They will be plot, but there may be missing relationships.") # showWarning: The warning was already shown for the classification genes.
				}
			}
		}

		# Check GENELABELS or get from GenesInfo 
		if(!is.null(geneLabels)) 
		{
			geneLabels <- extractGeneLabels(geneLabels) 
		}else
		{
			if (!is.null(genesInfoList))
			{
				for(cl in names(genesInfoList))
				{
					if ("GeneName" %in% colnames(genesInfoList[[cl]])) 
					{
						availableNames <-  which(!is.na(as.vector(genesInfoList[[cl]][,"GeneName"])))
						classGeneLabels	<- as.vector(genesInfoList[[cl]][,"GeneName"])[availableNames]
						names(classGeneLabels) <- rownames(genesInfoList[[cl]])[availableNames]
						geneLabels <- c(geneLabels, classGeneLabels)
					}
				}
			}
		}
		
		# Add possible missing genes
		if(!is.null(genesInfoList))
		{	
			for( cl in names(genesNetwork))
			{
				classGenes <- getNodes(genesNetwork[[cl]])
				missingGenes <- classGenes[which(!classGenes %in% rownames(genesInfoList[[cl]]))]

				temp <- matrix(NA, nrow=length(missingGenes), ncol=length(colnames(genesInfoList[[cl]])))
				rownames(temp) <- missingGenes
				colnames(temp) <- colnames(genesInfoList[[cl]])
				
				genesInfoList[[cl]] <- rbind(genesInfoList[[cl]], temp)
			}
		}	
		
		############################
		# Prepare NETWORK list
		# - Add classification nodes Network 
		#	- Add connected nodes Network
		############################
		# Extract CLASSIFICATIONgenesNetwork if available/needed & add to list
		classificationGenesNetwork <- NULL
		classificationGenesID <- NULL
		if(!is.null(classificationGenes) && plotClassifcationGenesNetwork)
		{
			classificationGenesID <- getRanking(classificationGenes, showGeneLabels=F, showGeneID=T)$geneID[, nwClasses, drop=F]
			classificationGenesNetwork <- getSubNetwork(genesNetwork, classificationGenesID)
			names(classificationGenesNetwork) <- paste(names(classificationGenesNetwork), " - Classification Genes",sep="")
								
			clToAdd <- which(sapply(genesNetwork, function(x){length(getNodes(x))}) - sapply(classificationGenesNetwork, function(x){length(getNodes(x))}) != 0)
			if(length(clToAdd) == 0)
			{
				warning("Only the classification genes network was provided. Only 'AllNodesNetwork' will be plotted.")
				plotAllNodesNetwork <- TRUE
			} else
			{
				for( i in names(clToAdd))
				{
					pos <- clToAdd[i]+ length(which(clToAdd<clToAdd[i]))
					if(pos == 0) pos <- 1
					genesNetwork <- c(genesNetwork[1:pos], classificationGenesNetwork[clToAdd[i]], genesNetwork[-(1:pos)])

				   genesInfoList <- c(genesInfoList[1:pos], list(genesInfoList[[i]][classificationGenesID[,clToAdd[i]][!is.na(classificationGenesID[,i])],] ), genesInfoList[-(1:pos)])
				   names(genesInfoList)[pos+1] <- names(classificationGenesNetwork[clToAdd[i]])
				}
			}
		}

		# Extract onlyCONNECTEDNodesNetwork if required & add to list
		# (Needs to be added after classific. in order to add it right after the "full" network)
		if(plotOnlyConnectedNodesNetwork)
		{
			for(cl in names(genesNetwork))
			{
				connectedGenes <- unique(as.vector(getEdges(genesNetwork[[cl]])[,c("gene1","gene2")]))
				if(length(connectedGenes)>0)
				{
					pos <- which(names(genesNetwork)==cl)
					genesNetwork <- c(genesNetwork[1:pos], getSubNetwork(genesNetwork[[cl]], connectedGenes), genesNetwork[-(1:pos)])
					names(genesNetwork)[pos+1] <- paste(cl, "\n(Only connected nodes)", sep="")
					
					if(!is.null(genesInfoList)) genesInfoList <- c(genesInfoList[1:pos], list(genesInfoList[[pos]][connectedGenes,]), genesInfoList[-(1:pos)])
				}
			}
			if(!is.null(genesInfoList)) names(genesInfoList) <- names(genesNetwork)
		}
		
		if (!plotAllNodesNetwork)
		{
			genesNetwork <- genesNetwork[-which(names(genesNetwork) %in% nwClasses)]
		}

		############################
		# PLOT
		############################
		# Check plotType
		if (!tolower(plotType) %in% c("dynamic", "static", "pdf")) stop("plotType is not valid. Please specify either 'dynamic', 'static' or 'pdf'.")
		
		if (tolower(plotType) == "pdf") pdf(fileName, height=height, width=width)	# Open device
		
		if (tolower(plotType) == "static") # Divide window
		{
			numClasses<-length(genesNetwork)
			if( numClasses>25 )  stop("Too many classes to draw in a single plot. Use 'pdf' instead.")						
			cols <- ceiling(sqrt(numClasses))
			rows <- ifelse(sqrt(numClasses)<round(sqrt(numClasses)), ceiling(sqrt(numClasses)),round(sqrt(numClasses)))
			par(mfrow=c(rows,cols))
		}


		# For each class...
		graphList<-NULL
		for(nw in names(genesNetwork))
		{
			classGenes <- unique(c(genesNetwork[[nw]]@edges[,"gene1"],genesNetwork[[nw]]@edges[,"gene2"]))
			#### Create graph object ####
			if((is.null(genesInfoList) || nrow(genesInfoList[[nw]])==0 ) || any(!classGenes %in% rownames(genesInfoList[[nw]]))) 
			{
				if(length(genesNetwork[[nw]]@nodes)>0)	{		classGraph <- graph.data.frame(as.data.frame(genesNetwork[[nw]]@edges[,ntwColnames,drop=F]), vertices=data.frame(nodes=genesNetwork[[nw]]@nodes), directed=FALSE)
				} else																	classGraph <- graph.data.frame(as.data.frame(genesNetwork[[nw]]@edges[,ntwColnames,drop=F]), directed=FALSE)
			}else
			{
				classGraph <- graph.data.frame(as.data.frame(genesNetwork[[nw]]@edges[,ntwColnames,drop=F]), directed=FALSE, vertices=data.frame(nodes=rownames(genesInfoList[[nw]]),genesInfoList[[nw]]))		
			}
			if (vcount(classGraph) != 0) 
			{
				#### Set graph parameters #####
				# Layout
				if(is.null(layoutList)) 
				{
					graphLayout <- layout.fruchterman.reingold(classGraph, niterNumeric=500) # .grid is faster, but the result looks far worse.
				}else
				{
					graphLayout <- layoutList[[nw]]
				}
				

				# Vertex labels
				vertexLabels <- get.vertex.attribute(classGraph,"name")
				if(!is.null(geneLabels)) vertexLabels[which(vertexLabels %in% names(geneLabels))] <- geneLabels[vertexLabels[which(vertexLabels %in% names(geneLabels))]]
				nVertex<-length(vertexLabels)
				
				# Vertex colors: Expression
				vertexColors<-rep("#7094FF", nVertex) # Ligth blue
				if(!is.null(get.vertex.attribute(classGraph,"exprsMeanDiff")))
				{
					exprsDiff <- as.numeric(get.vertex.attribute(classGraph,"exprsMeanDiff"))
					reds <- colorRampPalette(c("white","red"))(8)[c(7,6,5,4,4,3,3,2,2,2)] 	# Colors: 5, length:10. Index=0 very overexpressed, 10=almost 0
					greens <- colorRampPalette(c("white","darkgreen"))(8)[c(7,6,5,4,4,3,3,2,2,2)]				# Index=0 very rexpressed, 10 = almost 0 ?

					#vertexColors[which(exprsDiff>0)] <- reds[apply(sapply((max(exprsDiff,na.rm=T)/length(reds)) * 1:length(reds), function(x){ exprsDiff<=x }), 1,sum)][which(exprsDiff>0)] 
					vertexColors[which(exprsDiff>0)] <- reds[apply(	matrix(data=sapply((max(exprsDiff+0.001,na.rm=T)/length(reds)) * 1:length(reds), function(x){ exprsDiff<=x }), nrow=length(exprsDiff), ncol=length(reds))		, 1,sum)]			[which(exprsDiff>0)] 					
					vertexColors[which(exprsDiff<0)] <- greens[apply(matrix(data=sapply((min(exprsDiff-0.001,na.rm=T)/length(greens))*1:length(greens), function(x){exprsDiff>=x}), nrow=length(exprsDiff), ncol=length(greens))	, 1,sum)][which(exprsDiff<0)]				
				}
			#	vertexColors[which(is.na(vertexColors))]

				# Vertex Size: Discriminant power				
				if(!is.null(vertexSize)) 
				{
					vSize <- vertexSize
				}else
				{
					if(nVertex <= 50) {vSize <- 15}  # vertexSize = Minimum size
					if(nVertex > 50)  {vSize <- 10}
					if(nVertex >= 100){vSize <- 5 }
				}
								
				vertSizeArr <- rep(vSize,nVertex)  # By default sightly smaller than the minimum DP
				if(!is.null(get.vertex.attribute(classGraph,"discriminantPower")))
				{
					discPower<-round(as.numeric(get.vertex.attribute(classGraph,"discriminantPower")))
					if(any(!is.na(discPower)))
					{
						incr = ifelse( max(discPower, na.rm=T)>min(discPower, na.rm=T), (vSize*0.6/(max(discPower, na.rm=T)-min(discPower, na.rm=T))), 1) 		# vSize + incr = max size
							#min(discPower)	[1] 8 =15
							#max(discPower)  [1] 21 =25
						vertSizeArr <- rep((vSize*0.8), nVertex) # By default sightly smaller than the minimum DP
						dpNotNA <- which(!is.na(discPower))
						vertSizeArr[dpNotNA] <- vSize + ((discPower[dpNotNA] - min(discPower, na.rm=T))*incr)	
					}
				}
				# Shape: Classification gene
				vertexShape<- rep("circle",nVertex)
				labelColor <- rep("black",nVertex)
				if(!is.null(get.vertex.attribute(classGraph,"discriminantPower")))
				{
					dpNotNA <- which(!is.na(discPower))
					vertexShape[dpNotNA] <- rep("square", length(dpNotNA)) # Only works with plot() (not with tkplot) 
				}	
				if(!is.null(classificationGenes)) # alguna comprobacion mas?
				{
					vertexShape[which(get.vertex.attribute(classGraph,"name")%in% as.vector(getRanking(classificationGenes, showGeneID=T)$geneID))] <- "square"
				}

				# Edge color (Relation type)
				relColors<- c("blue", "orange")    # It is assumed there are only two types of relations
				relColors <- ifelse( get.edge.attribute(classGraph,"relation")==levels(factor(get.edge.attribute(classGraph,"relation")))[1], relColors[1],relColors[2])

				#### Output plot ####
				if (tolower(plotType)=="dynamic") {
					if(ecount(classGraph) > 0)
					{
						tkplot(classGraph, layout=graphLayout, vertex.label=vertexLabels, vertex.label.family="sans",  vertex.color=vertexColors, vertex.frame.color=vertexColors, vertex.label.color=labelColor, vertex.size=vertSizeArr, edge.color=relColors, edge.width=2, canvas.width=width, canvas.height=height)  
						#vertex.label.font=2, #Error in BioConductor Check?
						# vertex.label.cex=labelSize,  : En valores menores de 1 (0.3, 0.8...) a veces da error.
					} else {
						plotType<-"static"
						warning("Dynamic plot cannot be drawn for a network without edges.")
					}
				}
				if (tolower(plotType)=="static" || tolower(plotType)=="pdf"  )	
				{
					plot(classGraph, layout=graphLayout, vertex.label=vertexLabels, vertex.label.family="sans",  vertex.label.cex=labelSize, vertex.color=vertexColors, vertex.frame.color=vertexColors, vertex.label.color=labelColor, vertex.size=vertSizeArr, edge.color=relColors, edge.width=2, vertex.shape=vertexShape, main=nw) #vertex.label.font=2,
				}
			}
			graphList <- c(graphList, graph=list(classGraph))
		}
		
		#### Add legend ####
		if (tolower(plotType) == "pdf" ) 
		{
			plot.new()
			title("Network legend")

			text(0,0.9,"Node color: Expression", pos=4, font=2) #font=2 (bold)
			   text(0.15,0.8,"Repressed", pos=4)
			  points(0.37,0.8, pch=16, col="#31A354", cex=5)#238B45
			  points(0.47,0.8, pch=16, col="#EDF8E9", cex=5)#C7E9C0
			  points(0.57,0.8, pch=16, col="#FEE5D9", cex=5)#FCBBA1
			  points(0.67,0.8, pch=16, col="#FC9272", cex=5)#CB181D
			   text(0.7,0.8,"Overexpressed", pos=4)
			  #points(0.3,0.8, pch=16, col="#FFD0D0", cex=5)
			   #text(0.35,0.8,"Slightly overexpr.", pos=4)
			  points(0.25,0.7, pch=16, col="#7094FF", cex=5)
			   text(0.3,0.7,"Unknown", pos=4)
			   
			text(0,0.6,"Node shape: Chosen/Not chosen for classification", pos=4, font=2)
			  points(0.25,0.5, pch=22, col="black", cex=5)
			   text(0.30,0.5,"Chosen", pos=4)
			  points(0.6,0.5, pch=1, col="black", cex=5)
			   text(0.65,0.5,"Not chosen/Unknown", pos=4)
			   
			text(0,0.4,"Node size: Discriminant power (if available)", pos=4, font=2)
			  points(0.25,0.3, pch=22, col="black", cex=7)
			   text(0.30,0.3,"High DP", pos=4)
			  points(0.6,0.3, pch=22, col="black", cex=3)
			   text(0.65,0.3,"Low DP", pos=4)
			   
			text(0,0.2,"Line color: Relation type", pos=4, font=2)
			  lines(c(0.25,0.5),c(0.05,0.05),col="blue", lty="solid", lwd=2)
			   text(0.30,0.1,"Correlation", pos=4)
			  lines(c(0.6,0.9),c(0.05,0.05),col="orange", lty="solid", lwd=2)
			   text(0.6,0.1,"Mutual information", pos=4)
			   
			lines(c(0,1),c(1,1),col="black", lty="solid", lwd=1)  #  __
			lines(c(0,0),c(0,1),col="black", lty="solid", lwd=1)  # |
			lines(c(1,1),c(0,1),col="black", lty="solid", lwd=1)  #    |
			lines(c(1,0),c(0,0),col="black", lty="solid", lwd=1)  #  __
			
			#### Close dev ####
			dev.off()
			if (verbose){ message(paste("The plot was saved as ",getwd(),"/",fileName," (PDF file)",sep="")); flush.console()} 
		}
		
		############################
		# Return graph
		############################
		if(returniGraphs) 
		{
			names(graphList) <- names(genesNetwork)	
		}
	}
	if(returniGraphs) return (graphList)
}
