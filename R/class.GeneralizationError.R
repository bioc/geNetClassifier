
    # test <- new("GeneralizationError",accuracy=ret$globalAccuracy , sensitivitySpecificity=ret$sensitivitySpecificity ,classificationGenes.stats=ret$genesStats , classificationGenes.num=ret$numTrainGenes.stats ,confMatrix= ret$confussionMatrix ,probMatrix=ret$probabilityMatrix ,querySummary=ret$predictionStats )
     

########################################
###
### CLASS:      GeneralizationError
###
########################################
### METHODS:
###
### (initialize)
### (show)
### overview
### is.empty
###
###
########################################


setClass(Class="GeneralizationError",
        representation=representation(
            accuracy="matrix",
            sensitivitySpecificity="matrix",
            confMatrix="matrix",
            probMatrix    ="matrix",
            querySummary="list",
            classificationGenes.stats="list",
            classificationGenes.num="matrix"
            )
        )
        
    # To initialize with content, all slots should be provided.
    setMethod("initialize", "GeneralizationError", function(.Object, accuracy=NULL, sensitivitySpecificity=NULL,classificationGenes.stats=NULL, classificationGenes.num=NULL,confMatrix=NULL,probMatrix=NULL,querySummary=NULL) 
    {
        ##TODO ADD TYPE CHECKING
        .Object@accuracy <- accuracy
        .Object@sensitivitySpecificity <- sensitivitySpecificity
        .Object@classificationGenes.stats <- classificationGenes.stats
        .Object@classificationGenes.num <- classificationGenes.num
        .Object@confMatrix <- confMatrix
        .Object@probMatrix <- probMatrix
        .Object@querySummary <- querySummary
            
        .Object
    })
    
    # When typing the GenesRanking object, doesn't show the whole object content, but the ranking of the top10 genes of each class
    setMethod("show", "GeneralizationError", function(object)
    {
        cat("Estimated accuracy, sensitivity and specificity for the classifier:\n")
        print(object@accuracy)
        print(object@sensitivitySpecificity)
        cat('\nTo see all available statistics type "overview(EXAMPLE@generalizationError)"\n')
    })
    
    if(!exists("overview")) setGeneric(name="overview", def=function(object) standardGeneric("overview"))
    setMethod("overview", signature="GeneralizationError", definition=function(object)
    {
        cat("Slot @accuracy:\n")
        print(object@accuracy)
        cat("\n\nSlot @sensitivitySpecificity:\n")
        print(object@sensitivitySpecificity)
        cat("\n\nSlot @confMatrix:\n")
        print(object@confMatrix)
        cat("\n\nSlot @probMatrix:\n")
        print(object@probMatrix)
        cat("\n\nSlot @querySummary:\n")
        print(object@querySummary)
        cat("\n\nSlot @classificationGenes.num:\n")
        print(object@classificationGenes.num)
        cat("\n\nlapply(...@classificationGenes.stats,function(x)x[1:4,])):\n")
        print(lapply(object@classificationGenes.stats,function(x) x[1:ifelse(dim(x)[1]>=4, 4, dim(x)[1]),]))
    
    })
    
    # Check whether the object contains any info
    if(!exists("is.empty")) setGeneric(name="is.empty", def=function(object) standardGeneric("is.empty"))
    setMethod("is.empty",  signature="GeneralizationError",  function(object)
    {        
        if(nrow(object@accuracy)!=0) return(FALSE)
        if(nrow(object@sensitivitySpecificity)!=0) return(FALSE)
        if(nrow(object@confMatrix)!=0) return(FALSE)

        if(nrow(object@probMatrix)!=0) return(FALSE)
        if(length(object@querySummary)!=0) return(FALSE)
        if(nrow(object@classificationGenes.num)!=0) return(FALSE)
        if(length(object@classificationGenes.stats)!=0) return(FALSE)
        else return(TRUE)
    })
    
    
