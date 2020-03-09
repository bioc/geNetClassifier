
########################################
###
### CLASS:      GenesNetwork
###
########################################
### METHODS:
###
### (initialize)
### (show)
### plot.GenesNetwork
### overview
### nodes         // getNodes
### numNodes     // getNumNodes
### edges         // getEdges
### numEdges     // getNumEdges
### getSubNetwork
### network2txt
########################################

setClass(Class="GenesNetwork",
        representation=representation(
            nodes="character",
            edges = "matrix"
            )
        )
        
        
    # To initialize with content, all slots should be provided.
    # new("GenesNetwork", nodes = buildGenes[which(buildGenes[,cl]!="NA"),cl],edges = genesNetwork[[cl]])
    setMethod("initialize", "GenesNetwork", function(.Object, nodes=NULL, edges=NULL) 
    {
        ## ADD TYPE CHECKING
        .Object@nodes <- nodes
        .Object@edges <- edges
            
        .Object
    })
                                            # Alias: new("GenesNetwork", nodes = nodes, edges = edges)
                                                # if(!exists("GenesNetwork")) setGeneric(name="GenesNetwork", def=function(nodes=NULL, edges=NULL) standardGeneric("character"))
                                                # setMethod("GenesNetwork", "GenesNetwork", definition=function(nodes=NULL, edges=NULL) 
                                                # {
                                                    # new("GeneNetwork", nodes = nodes, edges = edges)
                                                # })
                                                
    # When typing the GenesRanking object, doesn't show the whole object content, but the ranking of the top10 genes of each class
    setMethod("show", "GenesNetwork", function(object)
    {
        cat("Attribute summary of the GenesNetwork:")
        cat("\nNumber of nodes (genes): ")
        print(getNumNodes(object))
        cat("Number of edges (relationships): ")
        print(getNumEdges(object))
    })
    
    #setGeneric(name="plot", def=function(x, y, ...)  standardGeneric("plot"))
    # setMethod(f="plot", signature=c(x="GenesNetwork", y="missing"), definition=function (x, y, ...)
    # {
        # message("Plotting default network view. For more options use plotNetwork()")
        ##plotNetwork(genesNetwork, classificationGenes=NULL, genesRanking=NULL, genesInfo=NULL,geneLabels=NULL, returniGraphs=FALSE, plotType="dynamic", plotAllNodesNetwork=TRUE, plotOnlyConnectedNodesNetwork=FALSE,  plotClassifcationGenesNetwork=FALSE, labelSize=0.5, width=800, height=500, fileName="genesNetwork.pdf", verbose=TRUE
        # plotNetwork (genesNetwork=x)
    # })
    
    # plot.GenesNetwork <- function (x, y="missing", classificationGenes=NULL, genesRanking=NULL, genesInfo=NULL,geneLabels=NULL, returniGraphs=FALSE, plotType="dynamic", plotAllNodesNetwork=TRUE, plotOnlyConnectedNodesNetwork=FALSE,  plotClassifcationGenesNetwork=FALSE, labelSize=0.5, width=800, height=500, fileName="genesNetwork.pdf", verbose=TRUE,...)
    # { 
    #     plotNetwork(genesNetwork=x, classificationGenes=classificationGenes, genesRanking=genesRanking, genesInfo=genesInfo,geneLabels=geneLabels, returniGraphs=returniGraphs, plotType=plotType, plotAllNodesNetwork=plotAllNodesNetwork, plotOnlyConnectedNodesNetwork=plotOnlyConnectedNodesNetwork,  plotClassifcationGenesNetwork=plotClassifcationGenesNetwork, labelSize=labelSize, width=width, height=height, fileName=fileName, verbose=verbose)
    # }

    # Version 1.2:
    plot.GenesNetwork <- function (x, y="missing",...)
    { 
        plotNetwork(genesNetwork=x,...)
    }
    
    if(!exists("overview")) setGeneric(name="overview", def=function(object) standardGeneric("overview"))
    setMethod("overview", signature="GenesNetwork", definition=function(object)
    {
        nNodesShow <- 10
        nEdgesShow <- 5
        
        if(length(object@nodes) < nNodesShow) nNodesShow <- length(object@nodes)
        if(dim(object@edges)[1] < nEdgesShow) nEdgesShow <- dim(object@edges)[1] 
        
        if(nNodesShow>0)
        {
            cat("getNodes(...)[1:",nNodesShow,"]:\n", sep="")
            print(getNodes(object)[1:nNodesShow])
            if(length(object@nodes) > nNodesShow) cat("... (",length(getNodes(object))," nodes)", sep="")
        } else cat("The network doesn't contain any nodes.")
        
        if(nEdgesShow>0)
        {
            cat("\n\ngetEdges(...)[1:",nEdgesShow,",]:\n", sep="")
            print(getEdges(object)[1:nEdgesShow,])
            if(dim(object@edges)[1] > nEdgesShow) cat("... (",nrow(getEdges(object))," edges)\n", sep="")
        } else cat("\n\nThe network doesn't contain any edges.\n")
    })
    
    # Get nodes
    if(!exists("getNodes")) setGeneric(name="getNodes", def=function(object) standardGeneric("getNodes"))
    setMethod("getNodes", signature="GenesNetwork", definition=function(object)
    {
        object@nodes
    })
    if(!exists("nodes") || class(nodes)[1]!="standardGeneric") setGeneric(name="nodes", def=function(object) standardGeneric("nodes"))
    setMethod("nodes", signature="GenesNetwork", definition=function(object)
    {
        getNodes(object)
    })
    
    
    # Get number of nodes
    if(!exists("getNumNodes")) setGeneric(name="getNumNodes", def=function(object) standardGeneric("getNumNodes"))
    setMethod("getNumNodes", signature="GenesNetwork", definition=function(object)
    {
        length(object@nodes)
    })
    if(!exists("numNodes")) setGeneric(name="numNodes", def=function(object) standardGeneric("numNodes"))
    setMethod("numNodes", signature="GenesNetwork", definition=function(object)
    {
        getNumNodes(object)
    })
    
    
    # Get edges
    if(!exists("getEdges"))setGeneric(name="getEdges", def=function(object) standardGeneric("getEdges"))
    setMethod("getEdges", signature="GenesNetwork", definition=function(object)
    {
        object@edges
    })
    if(!exists("edges") || class(edges)[1]!="standardGeneric") setGeneric(name="edges", def=function(object) standardGeneric("edges"))
    setMethod("edges", signature="GenesNetwork", definition=function(object)
    {
        getEdges(object)
    })
    
    
    ### Get number of edges
    if(!exists("getNumEdges"))setGeneric(name="getNumEdges", def=function(object) standardGeneric("getNumEdges"))
    setMethod("getNumEdges", signature="GenesNetwork", definition=function(object)
    {
        dim(object@edges)[1]
    })
    
    if(!exists("numEdges"))setGeneric(name="numEdges", def=function(object) standardGeneric("numEdges"))
    setMethod("numEdges", signature="GenesNetwork", definition=function(object)
    {
        getNumEdges(object)
    })

    # Extract a sub network
    if(!exists("getSubNetwork")) setGeneric(name="getSubNetwork", def=function(network, genes, showWarnings=TRUE) standardGeneric("getSubNetwork"))
    setMethod("getSubNetwork", signature="GenesNetwork", definition=function(network, genes=NULL, showWarnings=TRUE)
    {
        if(!is.logical(showWarnings)) showWarnings <- TRUE
        if(class(network)[1]!="GenesNetwork")  stop("The 'genes network' should be the network of one class returned by the classifier. (i.e. EXAMPLE@genesNetwork[[1]]")  
        if(is.matrix(genes)){
            if(dim(genes)[2] > 1) warning("The genes were provided in matrix format with more than one column, they might belong to different classes.", immediate=TRUE)
            genes <- as.vector(genes)
        }
        if(!is.character(genes))  stop("Genes should be given as a character vector. (i.e. genes=getRanking(EXAMPLE@classificationGenes)[,1]")

        genes <- genes[which(genes!="NA")]
    
        # Test whether the given nodes are in the network
        if(length(genes)>0)
        {                
            if(!all(genes %in% network@nodes) )
            {
                missingGenes <- genes[which(!genes %in% network@nodes)]
                genes <- genes[which(genes %in% network@nodes)]
                if(showWarnings) warning(paste("The following genes are not in the network: ", paste(missingGenes,  collapse=", "), ".",sep=""))
            }
        } else 
        {
                if(showWarnings)     warning("The given genes are not in the network.")
        }
        
        if(length(genes)>0) 
        {
            gen1 <- which(network@edges[,"gene1"] %in% genes) 
            gen2 <- which(network@edges[,"gene2"] %in% genes)
            subNet <- new("GenesNetwork", nodes= genes, edges=rbind(NULL,network@edges[gen1[which(gen1 %in% gen2)],]))            
        }else
        {
            subNet <- NULL
        }
        return(subNet)
    })
    
    setMethod("getSubNetwork", signature="list", definition=function(network, genes=NULL, showWarnings=TRUE)
    {
        if(!is.logical(showWarnings)) showWarnings <- TRUE
        if(!is.list(network) && !class(network[[1]]) !="GenesNetwork") stop("'network' should be a GenesNetwork object or the genes network returned by geNetClassifier. (i.e. EXAMPLE@genesNetwork)")
        if(class(genes)[1] =="GenesRanking") genes <- getRanking(genes, showGeneID=TRUE)$geneID
        if(!is.matrix (genes))  stop("When providing a list of networks (...@genesNetwork), the genes should be given as GenesRanking or matrix with named colums.")
        if(!any(colnames(genes) %in% names(network))) stop("The classes (names) of the network list and the colnames of the genes matrix don't match")
        
        classes <- names(network)[which(names(network) %in% colnames(genes))]
        subNetList <- NULL
        for( cl in classes )
        {        
            subNetList <- c(subNetList, netw=list(getSubNetwork(network[[cl]], genes[,cl], showWarnings=showWarnings)))
            names(subNetList)[length(subNetList)] <- cl
            if( length(classes) == 2) break
        }
        return(subNetList)
    })
    
    setMethod("getSubNetwork", signature="NULL", definition=function(network, genes=NULL, showWarnings=TRUE)
    {
        return(NULL)
    })
    
    # Moved to GeNetClassifierReturn 
    #setMethod("getSubNetwork", signature="GeNetClassifierReturn", definition=function(network, genes=NULL)
    #{
    #    return(getSubNetwork(network@genesNetwork, genes))
    #})
    

# Exports the network as .txt
if(!exists("network2txt")) setGeneric(name="network2txt", def=function(network, filePrefix=NULL, nwClass=NULL) standardGeneric("network2txt"))
setMethod("network2txt", signature="GenesNetwork", definition=function(network, filePrefix=NULL, nwClass=NULL)
{ 
    if (!is.list(network) && class(network) != "GenesNetwork") stop("network should be a 'GenesNetwork' or a list of networks.")
    if(!is.character(filePrefix) && !is.null(filePrefix)) stop("filePrefix should be of type character.")
    if(is.null(filePrefix)) filePrefix <- "genesNetwork"
    if(!is.null(nwClass) && !is.character(nwClass)) stop("nwClass should be of type character.")
    
    utils::write.table(getEdges(network), file=paste(filePrefix,"_edges_",nwClass,".txt", sep=""), sep="\t", row.names=FALSE)
    utils::write.table(getNodes(network), file=paste(filePrefix,"_nodes_",nwClass,".txt", sep=""), sep="\t", col.names=FALSE, row.names=FALSE)
})

setMethod("network2txt", signature="list", definition=function(network, filePrefix=NULL)
{
    for(nwClass in names(network))
    {
        network2txt(network[[nwClass]], filePrefix=filePrefix, nwClass=nwClass)
    }
})

















    
