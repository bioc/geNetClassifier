## I N D E X ##
# cvNumbers
# create.hypothesis
# iqr.filter
# PEB
#   calculateOrder
#   minimumDiff
#   gene.select
# difMean
# sampeling
# linear.SVM
# interaction.net	
# correlation.net   		
# remove.redundancy
# plotGeNetClassifierReturn
# configurePlotOutput
# queryResultCheck
#
# assignment.conditions
# SV.dif
# extractGeneLabels
# plotAssignmentPoint
# plotErrorNumGenes

############################### PRIVATE FUNCTIONS ##############################


## Estima los parametros de cross validation (numero de vueltas del bucle interno y externo) en base al numero de elemtos de cada clase  
## Si el numero de elementos por clase es menor que 7, no hay suficientes elementos para realizar el bucle externo: numCV.extern=1
cvNumbers <- function(numElemClass)
{
  if (numElemClass <=6)
  {
    numCV.intern <- numElemClass
    numCV.extern <- 1                # NO hacer la externa (no hay suficientes)
  }else{
    numCV.extern <- 5
    numCV.intern <- numElemClass-2
    if (numCV.intern > 8)
      numCV.intern <- 8
  }
  return (c(numCV.intern, numCV.extern))
}

#Funcion Crear hipotesis
# Ejemplo de uso: hipotesis <- create.hypothesis (datosCell$TIPO)
# datosCell es un ExpressionSet, $TIPO la variable utilizada para clasificar
create.hypothesis <- function(classLabels)
{
	classLabels <- factor(classLabels) #Transformed in case it wasn't a factor or the levels weren't updated
    numSamples <- length(classLabels)
    classNames <- levels(classLabels)
    numClasses <- length(classNames)

    # Hipotesis nula (todos iguales)
    hypothesis <- paste(rep(1,numSamples), sep=",", collapse=",")

	if (numClasses == 1){
		warning("There is only one class. Only the null hypothesis was created.", immediate. = TRUE)
		names(hypothesis)<- c("Null Hipothesis")
		return (hypothesis);
	}
	
    # Add hypothesis for each class
    for (i in classNames)
    {
         newHypothesis <- c()
         for (j in 1:numSamples)
         {
            if (classLabels[j] == i) newHypothesis <- c(newHypothesis, 2)
            else  newHypothesis <- c(newHypothesis, 1)
         }
         newHypothesis <- paste(newHypothesis, sep=",", collapse=",")
    	 hypothesis <- c(hypothesis, newHypothesis)
		 if (numClasses ==2 ) # Only one hypothesis needed. The another one is the complementary.
		 {
			names(hypothesis)<- c("Null Hipothesis", "BothClasses")
			return (hypothesis)
		 }
    }
	names(hypothesis)<- c("Null Hipothesis", classNames)
    return (hypothesis)
}

## filtrado IRQ
iqr.filter <- function(eset, percentage = 0.50)
{
	if(is.integer(eset[1,1]))		# Transforms in case of integer. !is.numeric doesnt work (it is true even if it is numeric)
	{
		newEset<-cbind(apply(eset,2,as.numeric))
		colnames(newEset)<-colnames(eset)
		rownames(newEset)<-rownames(eset)
		eset<-newEset
	}

  num <- dim(eset)[[2]]
  lowQ <- rowQ(eset, floor(0.25*num))
  upQ <- rowQ(eset, ceiling(0.75*num))
  iqrs<- upQ - lowQ
  giqrs <- iqrs > quantile(iqrs, probs= percentage)
  return(giqrs)
}

## Calculo de las probabilidades posteriores de los genes para cada hipotesis
## Recibe como argumento una matriz con los datos de expresion, el conjunto de hipotesis, el numero de clases, 
## el numero de elementos de cada clase y un booleano indicando si se van a devolver solo los genes sobreexpresados (TRUE) o tanto los genes soberexpresados como reprimidos (FALSE)

# En nuestro programa la opcion onlyOverExpressed no se utiliza
PEB <- function(eset, sampleLabels, labelsOrder=NULL,  nullHiphothesisFilter=0.95, untie= "bestRank")
{
	hypothesis 	<- create.hypothesis (sampleLabels)
	#if(!is.null(labelsOrder)) hypothesis <- hypothesis[c("Null Hipothesis",labelsOrder)]

	numClasses <- length(hypothesis) -1
	if(numClasses == 1) numClasses <- 2
	
	numElemClass 	<- table(table(sampleLabels))
	if(numElemClass[1]==numClasses){ 		numElemClass <- as.integer(names(numElemClass[1]))
						}else stop("There should be the same number of samples per class.")
	
	postProb <- gene.select(eset, hypothesis)
	colnames(postProb) <- names(hypothesis)
	postProb <- postProb[postProb[,"Null Hipothesis"] < nullHiphothesisFilter,]
	if(dim(postProb)[1] == 0) stop("The current genes don't differentiate the classes.")
	
	if(!is.null(labelsOrder))  postProb <-postProb[,c("Null Hipothesis", labelsOrder)]
	
	genesOrder <- calculateOrder(eset, sampleLabels, postProb, untie=untie)

	# Returns the order and the posterior probabilities
	ret <- new("GenesRanking",  postProb=postProb, numGenesClass=genesOrder$nGenes , ord=genesOrder$ord)
	return(ret)
}


# calculate order 
calculateOrder <- function (eset, sampleLabels, postProb,untie)		#untie="bestRank" / "bestPostProb"
{
	numClasses <- ncol(postProb) -1
	ranking <- apply(-postProb[,-1,drop=F], 2, rank, ties.method="min")
	
	# Desempatar:
	empates <- NULL
	if(numClasses > 2)
	{   
		tab<- apply(ranking, 2, table)
		if(is.matrix(tab)) 
		{
			temp<-NULL
			for(i in 1:dim(tab)[2]) temp[[i]]<-tab[,i]
			names(temp) <- colnames(tab)
			tab <- temp
		}
		empates <- lapply(tab, function(x)   names(which(x>1))) #if (any(x>1))
	}else 							empates <- list(BothClasses=names(which(table(ranking)>1)))
	
	if(!is.null(names(empates)))
	{
		for ( cl in names(empates))
		{
			#genes <- rownames(ranking)[which(ranking[,cl] %in% empates[[cl]])]
			#meanDiff <- minimumDiff (genes, cl, eset, sampleLabels)
			meanDiffs <- lapply( empates[[cl]], function(x) { minimumDiff (rownames(ranking)[which(ranking[,cl] == x)], cl, eset, sampleLabels)})
			names(meanDiffs)  <- empates[[cl]]
			
			meanDiffRanks <- lapply(meanDiffs, function(x) {rank(-x)-1}) # rank(-meanDiff)-1
			for(rnk in names(meanDiffRanks))
			{
				ranking[names(meanDiffRanks[[rnk]]),cl] <- as.numeric(rnk) + meanDiffRanks[[rnk]]
			}
		}	
	}
	
	# New order:
	ord <- apply(ranking, 2, order)	
	
	# Remove repeated genes from ord: 	Keep genes only in the class with best rank
	if(numClasses > 2 )
	{
		if(untie == "bestRank")	# The gene is kept in the class with best rank. (if there is a tie, maybe the postProbs should be compared, but not likely and would take longer)
		{
			for(gen in 1:dim(ord)[1]) 
			{
				geneRanks <- which(ord == ord[gen,1], arr.ind=TRUE)
				rmv <- geneRanks[-which(geneRanks[,1] == min(geneRanks[,1]), arr.ind=TRUE)[1],,drop=F]

				ord[rmv] <- 0
				# if (dim(rmv)[1] >0)
				# {
					# for(i in 1:dim(rmv)[1]) ord[ rmv[i,1],rmv[i,2]]<-0
				# }
			}
		}
		if(untie== "bestPostProb")  # The gene is kept in the class with best postProb. 
		{
			for(gen in 1:dim(ord)[1]) 
			{
				geneRanks <- which(ord == ord[gen,1], arr.ind=TRUE)
				genePostProb <- postProb[ord[gen,1],-1]
				rmv <- geneRanks[-which(geneRanks[,2]==which(genePostProb == max(genePostProb))),,drop=F]
				ord[rmv] <- 0
			}
		}

		nGenes<-apply(ord,2,function(x) { sum(x!=0)})
		# Reorder the matrix
		tmp<- matrix(ncol=ncol(ord), nrow=max(nGenes))
		colnames(tmp) <- colnames(ord)
		for(cl in names(nGenes))
		{	
			if(nGenes[cl]>0) tmp[1:nGenes[cl],cl]<-   ord[which(ord[,cl]!=0),cl]
		}
		ord <- tmp
	} else 
	{
		nGenes<-dim(ord)[1] 
		names(nGenes)<-colnames(ord)
	}
	return(list(ord=ord, nGenesClass=nGenes))
}
## Funcion para la seleccion de genes
## Recibe como argumentos un objeto exprSet y la lista de hipotesis que se quieren contrastar
## devuelve un objeto ''ebarraysPostProb'' con la probabilidad posterior para cada una de las hipotesis para cada gen
## Un objeto de tipo ebarraysPostProb is essentially a matrix of probabilities with number of rows given by the number of
##genes in 'data' and as many columns as the number of hypotheses for the fit. It additionally contains a slot `hypotheses'
##containing these hypotheses.
gene.select <- function (eset, hypothesis){
pattern <- ebPatterns(hypothesis)
ggOVR.em.out <- emfit(eset, family = "GG", hypotheses = pattern, num.iter = 10)
return(round(postprob(ggOVR.em.out, eset)$pattern,12))
}

# Returns the value of the difference between the mean expression of the given class, with the mean expression of the closest class
# This is used for untie the ranking
# minimumDiff (c("ENSG00000178445", "ENSG00000175097"), "ALL", eset, sampleLabels)
minimumDiff <- function (genesDiff, genesClass, eset, sampleLabels)
{
	if(any(!colnames(eset) %in% names(sampleLabels))) stop("The expressionSet doesn't contain the requested samples.")
	if(any(!genesDiff %in% rownames(eset))) stop("The expressionSet doesn't contain the requested genes.")
	if(genesClass %in% "BothClasses") genesClass <- levels(sampleLabels)[1]
	if(any(!genesClass %in% levels(sampleLabels))) stop("The genes' class is not one of the provided sample labels.")
	minDiff <- rep(0, length(genesDiff))
	names(minDiff)  <- genesDiff
	
	classes <-  levels(sampleLabels)
	numClasses <- length(classes)
	genesClass <- which(classes == genesClass)
	if(length(classes)==2) genesClass <- 1
	for( gen in genesDiff)
	{
		classExprsMean <- rep(0, numClasses)
		names(classExprsMean) <- classes

		for(cl in classes)
		{ 
			classExprsMean[cl] <- mean(eset[gen, names(sampleLabels) [which(sampleLabels == cl)]])
		}
		#minDiff[gen] <- min(abs(classExprsMean[genesClass] - classExprsMean[which(classExprsMean == max(classExprsMean[-genesClass]))]), abs(classExprsMean[genesClass] - classExprsMean[which(classExprsMean == min(classExprsMean[-genesClass]))]))
		minDiff[gen] <- min(abs(classExprsMean[genesClass] - classExprsMean[-genesClass]))
	}
	return(minDiff)
}


## Mean difference for each gene between a class and the other classes (calculated for all classes for each gene, returns a matrix)
# Used for knowing if the gene is up- or down-regulated
difMean <- function(eset, numClasses, numElemClass)
{
  dif <- matrix(data=0, nrow=dim(eset)[1], ncol=numClasses)
  rownames(dif) <- rownames(eset)
  samplesClass <- 1 : numElemClass
  for(cl in 1:numClasses)
  {
    # Para cada gen: Tomamos la media de todas las muestras de una clase
    # y le restamos media para ese gen en el resto de las clases quedandonos con el valor absoluto
	if (dim(eset)[1] > 1)
	{
		dif[,cl] <- apply(eset[,samplesClass],1,mean) - apply(eset[,-samplesClass],1, mean)
	}
	else
	{
		dif[,cl] <- mean(eset[,samplesClass]) - mean(eset[,-samplesClass])
		rownames(dif) <- rownames(eset)
	}

    # Podriamos quedarnos con el ratio entre las medias o medianas tambien en lugar de con la diferencia
    #   ratio[,cl] <- apply(exprs(rmaFiltrado)[,indices],1,median)/apply(exprs(rmaFiltrado)[,-indices],1,median)
    samplesClass <- ((cl*numElemClass)+1):((cl+1)*numElemClass)
  }
  return(dif)
}

## Returns a vector with a random order of the indexes of the matrix (treated as vector: only one dimension)
sampeling <- function(numElemClass, numClasses)
{
  mat <- matrix(data=0, nrow=numElemClass, ncol=numClasses)
  col <- sample(1:numClasses) ## randomizo el orden de las columnas
  for(i in 1:numClasses)
  {
    mat[,i] <- sample(((numElemClass*i)-(numElemClass-1)):(numElemClass*i)) #randomizo los elementos dentro de cada columna(filas)
  }
  sampleOVR <- as.vector(t(mat[,col]))
  return(sampleOVR)
}

## Funcion para construir el clasificador. LLeva a cabo la Cross validation
## Recibe como argumentos un objeto exprSet, una lista con la clase a la que pertenece cada muestra, el numero de clases,
## una lista con las muestras ordenadas de forma aleatoria para hacer el cross validation y el numero de grupos para la
## cross validation
## Devuelve una lista con el clasificador construido, la matriz de confusion y el error
## resultado$clasificador, resultado$matriz.confusion, resultado$error

linear.SVM <- function(trainEset, sampleLabels, numClasses, sampleOVR, numCrossValidation, minProbAssignCoeff=1, minDiffAssignCoeff=0.8)
{
  if(!is.matrix(trainEset)) trainEset <- rbind(NULL,trainEset)
  numSamples <- dim(trainEset)[[2]] #n es el numero de muestras totales
  numSamplesCV <- kfoldcv( numCrossValidation, numSamples ) #numSamplesCV array que nos indica el numero de muestras en cada uno de los grupos de cross validation
  rmaDataFrame <- as.data.frame(trainEset)
  mxcf <- matrix(data=0,  nrow=numClasses, ncol=numClasses) 
   colnames(mxcf) <- levels(sampleLabels)
   rownames(mxcf) <- levels(sampleLabels)
  for( i in 1:numCrossValidation)
  {
    if( i == 1)
      index.i <- sampleOVR[1:numSamplesCV[1]] #toma los primeros indices para escoger el conjunto de test
    else
      index.i <- sampleOVR[(sum(numSamplesCV[1:(i-1)]) + 1):sum(numSamplesCV[1:i])] # en las siguientes iteraciones toma un grupo para conjunto de test
	# SVM - TRAINING  (phase 1)
    svmClassif  <-  svm(x=t(rmaDataFrame[,-index.i]), y=sampleLabels[-index.i], C=1, kernel="linear", probability=T )
	
	# In case there are only two classes, svmClassif$SV doesn't contain the gene names
	if (is.null(dimnames(svmClassif$SV)) && numClasses ==2)
	{
		dimnames(svmClassif$SV)<- vector("list",2)
		colnames(svmClassif$SV) <- rownames(rmaDataFrame)
	}

    #predict.index.i  <-  predict(svmClassif, t(rmaDataFrame[,index.i]), probability=FALSE )
	predict.index.i  <- queryGeNetClassifier(svmClassif, as.matrix(rmaDataFrame[,index.i]), verbose=FALSE, minProbAssignCoeff=minProbAssignCoeff, minDiffAssignCoeff=minDiffAssignCoeff)$class
	predict.index.i  <- predict.index.i[predict.index.i!="NotAssigned"] # Ignorar los NA
	levels(predict.index.i)<- unique(c(levels(predict.index.i), levels(sampleLabels)))
	
    mxcf  <-  mxcf + table(predict.index.i, sampleLabels[names(predict.index.i)])[levels(sampleLabels),]
  }
  error.mxcf <- (1-(sum(diag( mxcf ))/sum( mxcf )))
  result <- list(classifier=svmClassif, confMatrix=mxcf, error=error.mxcf)	
  return(result)
}


# si lo voy a hacer para todas las clases juntas pero luego separando por las interacciones dentro de la misma clase hay que pasarle los genes en columnas para cada clase en lugar de una lista con todos los genes...
# eset <- as.data.frame(exprs(rmaFiltrado))
# si devolvemos el geneSymbol en lugar del ensg hay que utilizar gnes.human.Annotation
# genes deberia se el equivalente a rankENSG de manera que si es una matriz asignemos a cada gen su clase y si es un vector pues todos directamente

## Devuelve la red de interaccion de los genes pasados como parametro
interaction.net <- function(eset, genes, lp, method="clr", estimator="mi.empirical", threshold=0.5)
{
	net <- cbind(gene1=NULL, class1=character(0), gene2=NULL, class2=character(0), relation=character(0), value=numeric(0))
	
	#Comprobacion de parametros
	if(is(eset, "ExpressionSet")) eset <- exprs(eset) else if (!is.matrix(eset)) {warning("The first argument (eset) should be an expression matrix or an ExpressionSet.", immediate. = TRUE)
								return(net)	}
	if(is.vector(genes)) genes <- as.matrix(genes)
	if(!is.matrix(genes)) 		{warning("The argument 'genes' should be a matrix or a vector.", immediate. = TRUE)
								return(net)	}
	if(dim(genes)[1]>dim(eset)[1] || sum(which(dim(genes)==0))) {warning("The argument 'genes' can't be neither empty nor have more genes than the available.") 
												return(net)}
	if(!is.numeric(lp)) 		{warning("The argument 'lp' should be numeric.", immediate. = TRUE)
								return(net)	}
		if (sum(lp!=0)==0) return(net)
		if(dim(genes)[2]!=length(lp)){warning("The dimensions of the given genes and lp don't match.", immediate. = TRUE)
								return(net)	}
		if(sum(lp>dim(genes)[1])>0) {warning("The number of genes in lp is bigger than the available genes.", immediate. = TRUE)
								lp[1:length(lp)]<-dim(genes)[1]
								}
	#Calculos
    ensg<-c()
    classes <- c()
    for(i in 1:length(lp))
	{
		if(lp[i]!=0) 
		{
			ensg <- c(ensg, genes[1:lp[i],i])
			classes <- c(classes, rep(names(lp[i]), times=lp[i]))
		}
    }
    datos <- eset[as.vector(ensg),]
	if(!is.matrix(datos)) 
	{	
		datos<- rbind(NULL, datos)
		rownames(datos)<- as.vector(ensg)
	}	
	if(any(datos == -Inf)) stop("Expression matrix contains invalid values (-Inf).")
	names(classes) <- as.character(ensg)

	mi <- minet(t(datos), method=method, estimator=estimator, disc="equalwidth", nbins=sqrt(nrow(datos)))
	tmp <- which(mi > threshold, arr.ind=TRUE)
	tmp <- tmp[which(tmp[,1] < tmp[,2]),,drop=F]
	if(!is.matrix(class(tmp))) tmp<- rbind(NULL,tmp)

	net <- cbind(gene1=rownames(tmp), class1=classes[rownames(tmp)], gene2=colnames(mi)[tmp[,2]], class2=classes[colnames(mi)[tmp[,2]]], relation=rep(paste("Interaction -",method), times=nrow(tmp)), value=as.numeric(apply(tmp,1,function(x)mi[x[1],x[2]])))
	#rownames(net) <- apply(net[,c("gene1","gene2")], 1, function(x) paste(x[1],x[2],sep="_"))
	rownames(net) <- NULL
	return(net)
}

## Devuelve la red de correlaciones de los genes pasados como parametro
correlation.net <- function(eset, genes, lp, method="pearson", threshold=0.8)
##(correlation.net(esetFiltered, rankENSG[1:lp[i],i], lp[i], method="pearson", threshold=correlationsThreshold)
{
	net <- cbind(gene1=NULL, class1=character(0), gene2=NULL, class2=character(0), relation=character(0), value=numeric(0))
	
	#Comprobacion de parametros
	if(is(eset, "ExpressionSet")) eset <- exprs(eset) else if (!is.matrix(eset)) {stop("The first argument should be an expression matrix.", immediate. = TRUE)
								return(net)	}
	if(is.vector(genes)) genes <- as.matrix(genes)
	if(!is.matrix(genes) ) 		{stop("The argument 'genes' should be a matrix or a vector.", immediate. = TRUE)
								return(net)	}
		if(dim(genes)[1]>dim(eset)[1] || sum(which(dim(genes)==0))) {warning("The argument 'genes' can't be neither empty nor have more genes than the available.") 
												return(net)}
	if(!is.numeric(lp)) 		{stop("The argument 'lp' should be numeric.", immediate. = TRUE)
								return(net)	}
		if (sum(lp!=0)==0) return(net)
		if(dim(genes)[2]!=length(lp)){stop("The dimensions of the given genes and lp don't match.", immediate. = TRUE)
								return(net)	}
		if(sum(lp>dim(genes)[1])>0) {warning("The number of genes in lp is bigger than the available genes.", immediate. = TRUE)
								lp[1:length(lp)]<-dim(genes)[1] 
								}
	
	#Calculo de la correlacion
    ensg<-c()
    classes <- c()
    for(i in 1:length(lp)){
	  if(lp[i]!=0) 
	  {
		ensg <- c(ensg, genes[1:lp[i],i])
		classes <- c(classes, rep(names(lp[i]), times=lp[i]))
	  }
    }
    datos <- eset[as.vector(ensg),]
    names(classes) <- as.character(ensg)
  
	pearsoncor <- cor(t(datos), method=method)
	tmp <- which(pearsoncor > threshold, arr.ind=TRUE)
	tmp <- tmp[which(tmp[,1] < tmp[,2]),,drop=F] #Eliminamos los redundantes porque la matriz es simetrica
	if(!is.matrix(class(tmp))) tmp<- rbind(NULL,tmp)
	net <- cbind(gene1=colnames(pearsoncor)[tmp[,1]], class1=classes[colnames(pearsoncor)[tmp[,1]]], gene2=colnames(pearsoncor)[tmp[,2]], class2=classes[colnames(pearsoncor)[tmp[,2]]], relation=rep(paste("Correlation -",method), times=nrow(tmp)), value=as.numeric(apply(tmp,1,function(x)pearsoncor[x[1],x[2]])))
	#rownames(net) <- apply(net[,c("gene1","gene2")], 1, function(x) paste(x[1],x[2],sep="_"))
	rownames(net) <- NULL
	return(net)
}

## Elimina los genes correlacionados o que interactuan de la lista total de genes
remove.redundancy <- function(eset, genes, lp, net, relation=NULL)
{
	#Comprobacion de parametros
	if(is(eset, "ExpressionSet")) eset <- exprs(eset) else if (!is.matrix(eset))  	stop("The first argument should be an expression matrix or an expressionset.")	
	if(is.vector(genes)) genes <- as.matrix(genes)	
	if(!is.matrix(genes)) 	stop("The argument 'genes' should be a matrix.")
	if(!is.numeric(lp)) 	stop("The argument 'lp' should be numeric.")
			if(dim(genes)[2]!=length(lp)){warning("The dimensions of the given genes and lp don't match.", immediate. = TRUE)
								return(net)	}
	if(!class(net) == "GenesNetwork" && !is.list(net)) stop("The argument 'net' should be the GenesNetwork or a list of GenesNetwork returned by correlation.net or interaction.net.")	
	#if (is.matrix(net) && (sum(colnames(net) %in% c("gene1", "class1", "gene2", "class2", "relation", "value")) !=6 )) 	stop("The argument 'net' should be the matrix or a list of matrixes returned by correlation.net or interaction.net.")	
	if (is.list(net)&!is.null(net))
	{
		netList<- net
		net <- NULL
		i <- 1
		for( cl in 1:length(netList))
		{
			if(!is.null(netList[[cl]]))
			{
				subnet <- netList[[cl]]@edges
				if(!is.matrix(subnet) || (sum(!colnames(subnet) %in% c("gene1", "class1", "gene2", "class2", "relation", "value")) >0 ) ) 	stop("The argument 'net' should be the GenesNetwork or a list of GenesNetwork returned by correlation.net or interaction.net.")		
				if(sum(colnames(subnet) %in% c("class1", "class2")) ==0)  
				{
					classNameRep <- rep(names(netList)[i],dim(subnet)[1])
					subnet<- cbind(subnet, class1=classNameRep, class2=classNameRep)
				}
				net<- rbind(net,subnet)
			}
			i <- i+1
		}
	}else{
		if(is.null(net)) { return(net) }
		net <- net@edges
	}
	
	#Filter matrix to remove only the relation requested
	if (!is.null(relation))
	{
		if (relation != "Correlation - pearson" && relation != "Interaction - clr") warning ("The relation was not recognized. Check if the results were what you expected.", immediate=TRUE)
		net <- net[which(net[,"relation"]==relation),]			
	}
													
	#Calculos
	nonRedundantGenes <- list()
	removedGenes <- list()
		
	for( j in 1:length(lp))
	{
		tmp <- genes[1:lp[j],j]   
		tmp <- tmp[which(tmp!="NA")]
		i <- 1
		while(i<length(tmp)) 	# No cambiar por for, length(tmp) cambia en el bucle
		{						
		  if (tmp[i] %in% net[,c("gene1","gene2")]) 
		  {
			ind <- which(net[,c("gene1","gene2"), drop=F]==tmp[i] & net[,"class1"]==names(lp[j]) & net[,"class2"]==names(lp[j]), arr.ind=TRUE) # drop=F added. check it works!
			redundant <- unique(c(tmp[i], as.vector(net[ind[,1],c("gene1","gene2")])))[-1]
			if( length(which(tmp %in% redundant)) != 0)
			{
			  tmp <- tmp[-which(tmp %in% redundant)]
			  removedGenes[[names(lp[j])]] <- unique(c(removedGenes[[names(lp[j])]], redundant)) 
			}
		  }
		  i <- i+1 
		}
		nonRedundantGenes[[names(lp[j])]] <- tmp
	}
	ret <- list(nonRedundantGenes=nonRedundantGenes, removedGenes=removedGenes)
	return(ret)
}

# Used in main. Also exported as plot() for class Global Return
plotGeNetClassifierReturn <- function( geNetClassifierReturn, fileName=NULL, lpThreshold=0.75, numGenesLpPlot=1000, numGenesNetworkPlot=100, geneLabels=NULL, verbose=TRUE)
	{
		if(class(geNetClassifierReturn) != "GeNetClassifierReturn") stop("")
		if(!is.null(fileName) && !is.character(fileName))  stop("The plots name is not a valid name.")
		if(!is.numeric(lpThreshold) || (lpThreshold>=1 || lpThreshold <0)) stop("Lp threshold should be a probability (a number between 0 and 1).")
		if(!is.numeric(numGenesLpPlot)) stop("The number of genes to plot in the posterior probability plot should be a number .")
		if(!is.numeric(numGenesNetworkPlot) || numGenesNetworkPlot< 2) stop("The number of genes to plot in the network should be a number higher than two.")
		if(!is.logical(verbose)) verbose <- TRUE


		####### LP plot
		if(("genesRanking" %in% names(geNetClassifierReturn)) && numGenes(geNetClassifierReturn@genesRanking))
		{
			if(!is.null(fileName))
			{
				pdf(paste(fileName,"_postProb.pdf",sep=""))
			}
			calculateGenesRanking(precalcGenesRanking = geNetClassifierReturn@genesRanking, lpThreshold=lpThreshold, numGenesPlot=numGenesLpPlot, plotLp=TRUE, returnRanking=FALSE, verbose=FALSE)
			if(!is.null(fileName))
			{
				dev.off()
			}
		}		

		####### NETWORK plot
		if(("genesNetwork" %in% names(geNetClassifierReturn)) && length(geNetClassifierReturn@genesNetwork)) if ("igraph" %in% rownames(installed.packages()))  
		{
				# Extract requested network
				topGenes <- getRanking(getTopRanking(geNetClassifierReturn@genesRanking, numGenesClass=numGenesNetworkPlot), showGeneID=TRUE, showGeneLabels=FALSE)$geneID
				netTopGenes <- getSubNetwork(geNetClassifierReturn@genesNetwork, topGenes, showWarnings=FALSE)
				if(is.null(fileName))
				{
					plotType <-"dynamic"
					nwFileName <- ""
					plotOnlyConnectedNodesNetwork <- FALSE
				}else{
					plotType <-"pdf"
					nwFileName <- paste(fileName,"_network_top",numGenesNetworkPlot,"Genes.pdf",sep="")
					plotOnlyConnectedNodesNetwork <- TRUE
				}
				plotNetwork(genesNetwork=netTopGenes,  classificationGenes=geNetClassifierReturn@classificationGenes, genesRanking=getTopRanking(geNetClassifierReturn@genesRanking, numGenesNetworkPlot), plotAllNodesNetwork=TRUE, plotOnlyConnectedNodesNetwork=plotOnlyConnectedNodesNetwork,  plotClassifcationGenesNetwork=FALSE, genesInfo=NULL, geneLabels=geneLabels, returniGraphs=FALSE, plotType=plotType, labelSize=0.3, fileName=nwFileName, verbose=FALSE)
		}
			
		####### DISCRIMINANT POWER plot
		if (is.null(fileName)) 
		{
			dpFileName<- NULL
			x11()
		}
		else dpFileName <- paste(fileName,"_discriminantPower.pdf",sep="")
		
		if(("classifier" %in% names(geNetClassifierReturn)) && length(geNetClassifierReturn@classifier)) plotDiscriminantPower(geNetClassifierReturn@classifier, classificationGenes=geNetClassifierReturn@classificationGenes, fileName=dpFileName, returnTable=FALSE, verbose=FALSE)
				
		if (verbose && !is.null(fileName)) { message(paste("The plots were saved in ",getwd()," with the prefix '",fileName,"'.",sep="")); flush.console()}
}

# Opens the pdf file name, or, if null, divides the plot window for the appropiate number of plots
# Returns the number of genes for which the screen was prepared. Check if it is lower than provided to show warning.
# After the plot is drawn, the device the devide should be closed:
	# if (!is.null(fileName))
	# {
		# dev.off()
		# if (verbose){ message(paste("The SV plot was saved as ",getwd(),"/",fileName," (PDF file)",sep="")); flush.console()} 
	# }		
configurePlotOutput <- function(nGenes, fileName=NULL)
{
	
	if (!is.null(fileName))	
	{
		if(!is.character(fileName)) stop("The file name is not valid.")
		pdf(fileName)
	}else
	{
		if(nGenes>20)
		{ 
			nGenes <- 20
			#	warning(paste("Up to ",nGenes," genes will be shown. To plot more genes specify a PDF output file name.",sep=""))
		}
		if(nGenes<=20)
		{
			#Columns and rows of the plot
			if(nGenes == 4) rows <-2 		# 2x2
			else if (nGenes <= 5) rows <- 1
			else if (nGenes <= 10) rows <- 2
			else if (nGenes <= 15 || nGenes == 18) rows <- 3
			else rows <- 4
			
			cols <-ceiling(nGenes/rows)
			
			# Split plot window:
			par(mfcol=c(rows,cols))
		}
	}
	return(nGenes)
}



# Comprobar y concatenar si hay varias queries
# globalQueryResult <- queryResultCheck(queryResult)
queryResultCheck <- function(queryResult)
{
	if((!is.list(queryResult)) || (length(queryResult)<2)) stop('The argument should be the result from executing queryGeNetClassifier.')
	else if (!is.factor(queryResult$class) || !is.matrix(queryResult$probabilities))  stop('The argument should be the result from executing queryGeNetClassifier.')

	# Eliminar los $call de la prediccion
	temp<- queryResult
	queryResult<- NULL
	for(i in 1:length(temp))
	{
		if(is.null(temp[i]$call)) queryResult<- c(queryResult, temp[i])
	}

	# Concatenar
	if(length(queryResult)==2) globalQueryResult <- queryResult
	if(length(queryResult)>2)
	{
		globalQueryResult <- c(queryResult[1], queryResult[2])
		
		numQueryResults <- length(queryResult)/2
		for (i in 2:numQueryResults)
		{
			#Comprobar q coinciden etiquetas de probabilities
			# Probabilities: Tienen q ser exactamente las mismas (en numero y orden)
			if(dim(globalQueryResult$probabilities)[1] !=  dim(queryResult[i*2]$probabilities)[1])
			{
				warning("The classes of the predictors don't match.", immediate. = TRUE)
				return (list(stats=NULL, countNotAssigned=NULL,statsLegend=NULL))
			}
			else 
			{
				if (sum(rownames(globalQueryResult$probabilities) != rownames(queryResult[i*2]$probabilities))!=0) #Si las dimensiones son distintas da error
				{
					warning("The classes of the predictors don't match.", immediate. = TRUE)
					return (list(stats=NULL, countNotAssigned=NULL,statsLegend=NULL))
				}
			}
			
				# Unir las clases (puede que los niveles no coincidan si solo se hayan asignado algunas clases en cada resultado concreto)
			globalNames<- names(globalQueryResult$class)
			globalQueryResult$class <- factor(c(as.character(globalQueryResult$class),as.character(queryResult[(i*2)-1]$class)))
			names(globalQueryResult$class) <- c(globalNames, names(queryResult[(i*2)-1]$class))
			
			# Unir las probabilidades
			globalQueryResult$probabilities <- cbind(globalQueryResult$probabilities,queryResult[i*2]$probabilities)
		}		
	}								

	return(globalQueryResult)
}

## Establece los valores minimos de probabilidad de una muestra para que sea asignada a una determinada clase
# Para 2 clases, requiere al menos un 75%
assignment.conditions <- function(prob, minProb, minDiff)
{
	ord <- order(prob,decreasing=TRUE)
	ret <- names(prob)[ord[1]]

	if(prob[ord[1]]-prob[ord[2]]<=minDiff)
	{
			ret <- "NotAssigned"								
	}   
	if(length(prob)!=2 && prob[ord[1]]<=minProb)
	{
		ret <- "NotAssigned"
	}
	
	return(ret)
}

## Calcula la diferencia para un gen entre la suma de los valores SV para la clase que diferencian y la siguiente clase mas cercana a ese valor (valor discriminante...)
SV.dif <- function(classifier, gene, originalGeneNames=NULL, correctedAlpha=FALSE)
{
	numClasses <- length(classifier$levels)
	
	if(!gene %in% colnames(classifier$SV))   
	{
		m <- matrix(0)	#Transformed into dataframe colnames to make sure they match the $SV
		colnames(m) <- gene
		m <- data.frame(m)
		#gene <- gsub("^X","",colnames(data.frame(m)))   
		
		colnames(classifier$SV)[which(colnames(classifier$SV)%in% colnames(m))] <- gene
		#colnames(classifier$SV)<-  gsub("^X","",colnames(classifier$SV))	#If the genenames start by number -->DataFrame adds an X to the name
	}
	
	# same number of elements in plot as the number of support vectors, nrow=max(classifier$nSV)
	tmp <- matrix(ncol=numClasses, nrow=max(classifier$nSV))
	colnames(tmp) <- classifier$levels

	# Si se corrige alfa: Sumar los valores absolutos de los coeficientes y multiplicar por el del gen
	if(correctedAlpha)
	{
		alpha <- apply(classifier$coefs,1,function(x)sum(abs(x)))
		mult <- classifier$SV[,gene]*alpha
	}else 
	{
		mult <- classifier$SV
	}

	ind <- 1
	for(i in 1:numClasses)
	{
		# mult para los SV de la clase. Rellenado con la media en las clases que tienen menos.
		if(correctedAlpha)
		{
		  tmp[,i] <- mean(mult[ind:(ind+classifier$nSV[i]-1)]) 		 
		  tmp[1:classifier$nSV[i],i] <- mult[ind:(ind+classifier$nSV[i]-1)]
		} else 
		{
		  tmp[,i] <- mean(mult[ind:(ind+classifier$nSV[i]-1),gene]) 
		  tmp[1:classifier$nSV[i],i] <- mult[ind:(ind+classifier$nSV[i]-1),gene]
		}
		ind <- ind+classifier$nSV[i]
	}
	
	# Separate negative/positive
	neg <- tmp
	pos <- tmp
	for(i in 1:numClasses)
	{
		neg[which(neg[,i]>0),i] <- 0
		pos[which(pos[,i]<0),i] <- 0
	}

	# Max and min
	sum_pos <- apply(pos,2,sum)
	maxim <- which(sum_pos == max(sum_pos), arr.ind=TRUE)[1]
	sum_neg <- apply(neg,2,sum)
	minim <- which(sum_neg == min(sum_neg), arr.ind=TRUE)[1]
	
	# up regulated
	if (max(sum_pos[-maxim]) != 0){ 		sig <- max(sum_pos[-maxim]) 	#Si el siguiente valor no es cero
	}else 									sig <- max(sum_neg[-maxim])		#Si no hay mas "positivas" se busca la negativa
	up_reg <- max(sum_pos)-sig
	
	# down regulated
	if (min(sum_neg[-minim]) !=0){ 			sig <- min(sum_neg[-minim])
	}else 									sig <- min(sum_pos[-minim]) 
	dw_reg <- min(sum_neg)-sig

	if(numClasses>2)
	{
		# Choose biggest absolute value (The most diferent class' DP) 
		if(up_reg < abs(dw_reg)){ 	discriminantPower <- dw_reg
		}else					 					discriminantPower <- up_reg
	}else
	{
		# Choose smallest absolute value (Taking into account the most likely to be mistaken)
		if(up_reg > abs(dw_reg)){ 	discriminantPower <- dw_reg
		}else					 					discriminantPower <- up_reg
	}
	
	# Find save the class for wich it was calculated
	dpClass <- classifier$levels[ifelse(discriminantPower>0,which(sum_pos == max(sum_pos)),which(sum_neg == min(sum_neg)))]
	
	ret<- list(discriminantPower=discriminantPower, discrPwClass=dpClass ,positive=pos, negative=neg)
	return(ret)
}


# Checks the given annotation format and extracts the desired gene names from it
# Returns: Vector (name=geneID, content=geneShowName), or NULL
# genesAnnotation: Data frame or matrix containing the expressionSet "geneID" as first column or rowname, and the desired "geneName" (symbol) to show in plots and tables 
#
extractGeneLabels <- function(genesAnnotation, geneList=NULL)
{
	if(is.null(genesAnnotation)) 
	{
		ret <- NULL
	}else
	{
		if (is.matrix(genesAnnotation) || is.data.frame(genesAnnotation))
		{
			if(dim(genesAnnotation)[1]<dim(genesAnnotation)[2]) genesAnnotation <- t(genesAnnotation)  # Ordered as column
			if(dim(genesAnnotation)[2]>1)
			{
				if(!any(rownames(genesAnnotation)!=genesAnnotation[,1]))  genesAnnotation <- as.matrix(genesAnnotation[,2])		# The rows are named and the first column contain the geneID
				if(is.null(rownames(genesAnnotation))) 
				{
					rownames(genesAnnotation) <- genesAnnotation[,1]	# First column assumed to contain the geneID
					genesAnnotation <- as.matrix(genesAnnotation[,2])	# Second column assumed to contain the gene show name (symbol)
				} 
			}
		}else 
		{
			if (is.factor(genesAnnotation) || is.character(genesAnnotation))
			{	
				if(length(genesAnnotation) > 0) genesAnnotation <- as.matrix(genesAnnotation) # Transform into matrix
			}else stop("The genesAnnotation should be either a matrix, a dataframe or a vector.")
		}
		

		if(!is.null(geneList))
		{
			if(any(!geneList %in% rownames(genesAnnotation)))
			{
				#warning("Some of the given gene IDs are not available in the genes Annotation. Their alias will not be shown.", immediate.=TRUE)
				# return(c("a")[0])
			}
		}else { geneList <- rownames(genesAnnotation) }
			
		if(!is.null(rownames(genesAnnotation))) 
		{
			ret <- as.character(genesAnnotation[geneList,]) # Returns vector
			names(ret) <- geneList
		}else ret <- genesAnnotation
	}
	return(ret)
}

# Aux function for plotAssignments
plotAssignmentPoint<-function(prob, correctColor, incorrectColor)
{
	realLabels <- prob["realLabels"]
	prob <- prob[-length(prob)]
	
	ord <- order(prob,decreasing=TRUE)
	color <- ifelse(realLabels ==  names(prob)[ord[1]], correctColor, incorrectColor)
	size <- ifelse(realLabels ==  names(prob)[ord[1]], 0.8, 1.2)
	
	prob<-as.numeric(prob)
	cords <- c(prob[ord[1]], prob[ord[1]]-prob[ord[2]])
	points(cords[1], cords[2], col=color, pch=16, cex=size)
	return(cords)
}


	# ARREGLAR O ELIMINAR
		plotErrorNumGenes <- function(numGenesTPglobal, numGenesClass, numTrainGenes, numGenesClassGE, plotsName, plotCVgE=FALSE)
		{
			if(!is.null(numGenesTPglobal) && (!is.vector(numGenesTPglobal) && !is.matrix(numGenesTPglobal)))  stop("numGenesTPglobal should be a numeric vector.")
			numClasses <- ncol(numGenesClass)

			pdf(plotsName)
	
			# buildClassifier=TRUE
			if(!is.null(numGenesTPglobal))
			{
				######
				# Plot of errorNumGenes (all in one)
				if((length(numGenesTPglobal)>3 && length(numGenesTPglobal)<10) && library(RColorBrewer,logical.return=TRUE))
				{
					colores<- brewer.pal(length(numGenesTPglobal),"Dark2")	
				}else colores <- rainbow(length(numGenesTPglobal))
				maxError <- round(1-min(sapply(numGenesTPglobal, function(x){min(x[,numClasses+1])})), 1)
				maxGenes <- ceiling(max(sapply(numGenesTPglobal, function(x){max(apply(x[,1:numClasses, drop=F],1,sum))}))/20)*20
				numGenesTPglobal <- lapply(numGenesTPglobal, function(x)  { cbind(x, totalGenes=apply(x[,1:numClasses,drop=F],1,sum), error=1-x[,numClasses+1])})
			
				plot.new( )
				plot.window(xlim=c(0, maxGenes), ylim=c(0,maxError))
				title( main="Gene-selection iterations",  xlab="Total number of genes", ylab="Error rate") #sub="(Error  vs Total gene number)",
				axis(1)
				axis(2)
				for(iter in 1:length(numGenesTPglobal))
				{
					errorNumGenes <- numGenesTPglobal[[iter]][,c("totalGenes", "error"), drop=F]
					#maxError <- max(errorNumGenes[,"error"])
					minError <- which(errorNumGenes == min(errorNumGenes[,"error"]), arr.ind=T)
					lines(errorNumGenes, type="l", col=colores[iter], lwd=2)
					points(errorNumGenes[minError[1],1], errorNumGenes[minError[1],2], type="b", pch=16, col=colores[iter]) 
					text(maxGenes*0.60,maxError-(iter*(maxError*0.05)), paste("Error ",round(errorNumGenes[minError[1],2],2), ": ",errorNumGenes[minError[1],1]," genes",sep=""), col=colores[iter], pos=4)
				}	
				
				######
				# Gene selection (barplot) # Classifier
				barplot(numGenesClass, beside=TRUE, main="Number of genes selected in each iteration", density=50, col=colores, sub=paste("Total number of selected genes: ", sum(numTrainGenes["Classifier",])), ylim=c(0, max(numGenesClass)), border=colores, xlab="Classes", ylab="Number of genes")

				numTrainGenesPlot <-numGenesClass
				for(cl in 1:dim(numTrainGenes)[2])
				{			
					numTrainGenesPlot[-which(numTrainGenesPlot[,cl] == numTrainGenes["Classifier",cl])[1],cl]<- 0
				}
				barplot(numTrainGenesPlot, beside=TRUE,add=TRUE, col=colores, border=colores)
				
				######
				# Plot max genes recorridos
				maxGenesChecked <- NULL
				for(i in 1:length(numGenesTPglobal)) {
					maxGenesChecked <- rbind(maxGenesChecked, numGenesTPglobal[[i]][dim(numGenesTPglobal[[i]])[1],1:numClasses])
				}
				barplot(maxGenesChecked, beside=TRUE, main="Total number of genes explored in each iteration", col=colores, density=50, xlab="Classes", ylab="Number of genes")
				barplot(numGenesClass, beside=TRUE,add=TRUE, col=colores)
			}
			
			# GeneralizationError
			if(!is.null(numGenesClassGE) && plotCVgE)
			{
				######
				# Gene selection (barplot) # GE
				nIters <- dim(numGenesClassGE[[1]])[1]
				if((length(numGenesClassGE)==5) && library(RColorBrewer,logical.return=TRUE))
				{
					colores<- brewer.pal(9,"Purples")[-c(1,2,3,9)]	
				}else colores <- rainbow(length(numGenesClassGE))
				
				numGenesClassGEbinded <-NULL 
				for(i in 1:5) numGenesClassGEbinded<-rbind(numGenesClassGEbinded, numGenesClassGE[[i]]) # 5 = GE loops
				colsCV <- NULL
				for(colCVIter in colores) colsCV<-c(colsCV, rep(colCVIter, nIters))
				
				maxY <- ifelse(is.null(numGenesTPglobal), max(numGenesClassGEbinded), max(max(numGenesClassGEbinded), max(numGenesClass))) # Misma escala 
				barplot(numGenesClassGEbinded, beside=TRUE, main="Genes selected in the 5-fold CV for Generalization Error", density=50, col=colsCV, ylim=c(0, maxY), border=colsCV, xlab="Classes", ylab="Number of genes")

				# Selected genes
				numTrainGenesPlot <-numGenesClassGEbinded
				ind <- 0
				indexes <- NULL
				for(iter in 1:5) 
				{
					indexes <- rbind(indexes, cbind(rep(0, numClasses), 1:numClasses))
					for(cl in 1:numClasses)
					{			
						indexes[(iter-1)*numClasses+cl, 1] <- ind + which(numTrainGenesPlot[(ind+1):(ind+nIters), cl] == numTrainGenes[iter,cl])[1]
					}
					ind <- ind+nIters
				}
				
				# Hacer 0...
				temp <- matrix(nrow=nrow(numTrainGenesPlot), ncol=ncol(numTrainGenesPlot), data=0) 
				temp[indexes] <- numTrainGenesPlot[indexes]
				numTrainGenesPlot <- temp
				
				barplot(numTrainGenesPlot, beside=TRUE,add=TRUE, col=colsCV, border=colsCV)
			}
			
			dev.off()
		}