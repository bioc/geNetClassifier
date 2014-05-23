########################################
###
### CLASS:      GeNetClassifierReturn
###
########################################
### METHODS:
###
### (initialize)
### (show)
### plot
### overview
### gSlots. Alias: names
###
######################################## 

setClass(Class="GeNetClassifierReturn",
        representation=representation(
            call="call",
            classifier="list",
            classificationGenes="GenesRanking", 
            generalizationError="GeneralizationError",
            
            genesRanking="GenesRanking",
                genesRankingType="character",     # "significant" (only genes over lpThreshold) or "all" (full genes ranking)
            genesNetwork="list",                             # List of GenesNetwork for each class
                genesNetworkType="character"      # "classification" (Network within the classification genes) or "topGenes" (network of the top numGenesNetworkPlot genes)
            )
        )
        
    # To initialize with content, all slots should be provided.
    # test <- new("GeNetClassifierReturn",accuracy=ret$globalAccuracy , sensitivitySpecificity=ret$sensitivitySpecificity ,classificationGenes.stats=ret$genesStats , classificationGenes.num=ret$numTrainGenes.stats ,confMatrix= ret$confussionMatrix , querySummary=ret$predictionStats )
    setMethod("initialize", "GeNetClassifierReturn", function(.Object, call=NULL, classifier=NULL, classificationGenes=NULL, generalizationError=NULL, genesRanking=NULL, genesRankingType=NULL, genesNetwork=NULL, genesNetworkType=NULL)  
    {
        ##TODO ADD TYPE CHECKING
        .Object@call <- call

        if(!is.null(classifier)) .Object@classifier             <- classifier
        if(!is.null(classifier) && !any(class(classifier) %in% "svm")) warning("The classifier provided is not a SVM.")
        if(!is.null(classificationGenes)) .Object@classificationGenes <- classificationGenes
        
        if(!is.null(generalizationError)) .Object@generalizationError             <- generalizationError
        
        if(!is.null(genesRanking)) .Object@genesRanking         <- genesRanking
        if(!is.null(genesRankingType)) .Object@genesRankingType     <- genesRankingType
        
        if(!is.null(genesNetwork)) .Object@genesNetwork         <- genesNetwork
        if(!is.null(genesNetworkType)) .Object@genesNetworkType     <- genesNetworkType
            
        .Object
    })
    
    # When typing the GenesRanking object, doesn't show the whole object content, but the ranking of the top10 genes of each class
    setMethod("show", "GeNetClassifierReturn", function(object)
    {
        cat("R object summary:\n")
        
        if (!is.null(object@classificationGenes) && (is.list(object@classifier) && length(object@classifier)>0))
        {
            classGenes<-apply(getRanking(object@classificationGenes, showGeneID=TRUE, showGeneLabels=FALSE)$geneID, 2, function(x) {length(which(!is.na(x)))})
            cat("Classifier trained with ",length(object@classifier$SVMclassifier$fitted)," samples. ", sep="")
            if(!is.null(object@call$removeCorrelations) && !is.null(object@call$removeInteractions))
            {
                if(eval(object@call$removeCorrelations) && eval(object@call$removeInteractions)) cat("\nCorrelations and interactions (mutual information) between genes were removed.")
            }else{
                if(!is.null(object@call$removeCorrelations) && eval(object@call$removeCorrelations)) cat("\nCorrelations between genes were removed.")
                if(!is.null(object@call$removeInteractions) && eval(object@call$removeInteractions)) cat("\nInteractions (mutual information) between genes were removed.")
            }
            cat("\nTotal number of genes included in the classifier: ",sum(classGenes),". \nNumber of genes per class: \n",sep="")
            print(classGenes)    
            cat("For classificationGenes details: genesDetails(EXAMPLE@classificationGenes)")
        }
        
        if (!is.null(object@generalizationError)&& !is.empty(object@generalizationError))
        {
            cat("\n\nGeneralization error and gene stats calculated through 5-fold cross-validation:\n")
            print(slotNames(object@generalizationError))
        }

        if(!is.null(object@genesRanking) && (nrow(object@classificationGenes@ord)!=0))
        {
            cat("\nThe ranking of ",object@genesRankingType," genes contains (genes per class):\n",sep="")
            print(numGenes(object@genesRanking))
        }
        
        if(length(object@genesNetwork) > 0)
        {
            cat("\nThe networks calculated for the",object@genesNetworkType,"genes of each class contain:\n",sep=" ")
            print( rbind("Number of genes"= sapply(object@genesNetwork, getNumNodes), "Number of relations"=sapply(object@genesNetwork, getNumEdges)))
        }
        
        cat("\nAvailable slots in this R object:\n")
        #slots <- slotNames(object)
        #lapply (slotNames(object), function(x) if (!is.null(eval...x)) print(x))
        print(names(object))
        cat('To see an overview of all available slots type "overview(EXAMPLE)"\n')

    })
    
    plot.GeNetClassifierReturn <- function (x, y="missing", fileName=NULL, lpThreshold=0.95, numGenesLpPlot=1000, numGenesNetworkPlot=100, geneLabels=NULL, verbose=TRUE, ...)
    { 
        plotGeNetClassifierReturn( geNetClassifierReturn=x, fileName=fileName, lpThreshold=lpThreshold, numGenesLpPlot=numGenesLpPlot, numGenesNetworkPlot=numGenesNetworkPlot, geneLabels=geneLabels, verbose=verbose)
    }
    
    if(!exists("overview")) setGeneric(name="overview", def=function(object) standardGeneric("overview"))
    setMethod("overview", signature="GeNetClassifierReturn", definition=function(object)
    {
        #TODO ADD if( !=NULL...)
        
        cat("Slot @call:\n") 
        print(object@call)
        
        if(!is.null(object@classifier) && (is.list(object@classifier) && length(object@classifier)>0)) 
        {
            cat("\nSlot @classifier:")
            print(object@classifier)
        }
        if(!is.null(object@classificationGenes)&&(nrow(object@classificationGenes@ord)!=0)) 
        {
            cat("\n\nSlot @classificationGenes:\n")
            print(object@classificationGenes)
        }
        if(!is.null(object@generalizationError) && !is.empty(object@generalizationError)) 
        {
            cat("\nSlot @generalizationError:\n")
            print(object@generalizationError)
        }
        if(!is.null(object@genesRanking) && (nrow(object@classificationGenes@ord)!=0)) 
        {
            cat("\nSlot @genesRanking (of ",object@genesRankingType," genes):\n", sep="")
            print(object@genesRanking)
        }
        if(length(object@genesNetwork) > 0)
        {
            cat("\nSlot @genesNetwork: \nList of genes networks of the ",object@genesNetworkType," genes for each of the classes: \n", sep="")
            print( rbind("Number of genes"= sapply(object@genesNetwork, getNumNodes), "Number of relationships"=sapply(object@genesNetwork, getNumEdges)))
        }
    
    })
    
    
    if(!exists("gSlots")) setGeneric(name="gSlots", def=function(x) standardGeneric("gSlots"))
    setMethod("gSlots", signature="GeNetClassifierReturn", definition = function(x)
    {
        object <- x
        slots<-c("call")
        if(!is.null(object@classifier) && (is.list(object@classifier) && length(object@classifier)>0)) slots <- c(slots,"classifier")
        if(!is.null(object@classificationGenes) && (nrow(object@classificationGenes@ord)!=0)) slots <- c(slots,"classificationGenes")
        if(!is.null(object@generalizationError) && !is.empty(object@generalizationError)) slots <- c(slots,"generalizationError")
        if(!is.null(object@genesRanking) && (nrow(object@genesRanking@ord)!=0)) slots <- c(slots,"genesRanking")
        if(length(object@genesRankingType) > 0) slots <- c(slots,"genesRankingType")
        if(length(object@genesNetwork) > 0) slots <- c(slots,"genesNetwork")
        if(length(object@genesNetworkType) > 0) slots <- c(slots,"genesNetworkType")
        
        return (slots)
    })
    
    setMethod("names", signature="GeNetClassifierReturn", definition = function(x)
    {        
        return (gSlots(x))
    })

    # Method for GenesNetwork. Signature: GeNetClassifierReturn - Network will be extracted
    if(!exists("getSubNetwork")) setGeneric(name="getSubNetwork", def=function(network, genes, showWarnings=TRUE) standardGeneric("getSubNetwork")) #Copy from class.GenesNetwork
    setMethod("getSubNetwork", signature="GeNetClassifierReturn", definition=function(network, genes=NULL, showWarnings=TRUE)
    {
        return(getSubNetwork(network@genesNetwork, genes, showWarnings=showWarnings))
    })
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
