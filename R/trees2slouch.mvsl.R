## This file is part of mvSLOUCH

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this package or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .


'.ouch2slouch.mvsl'<-function (tree) 
{
## Function from the slouch package. Author : Jason Pienaar
    if (!inherits(tree, "ouchtree")) 
	stop(sQuote("tree"), " must be of class ", sQuote("ouchtree"))
    N <- length(tree@nodes)
    tmp <- as(tree, "data.frame")
    tmp$ancestors <- as.character(tmp$ancestors)
    tmp$ancestors <- as.numeric(tmp$ancestors)
    tmp$times <- as.character(tmp$times)
    tmp$times <- as.numeric(tmp$times)
    tmp$nodes <- as.character(tmp$nodes)
    tmp$nodes <- as.numeric(tmp$nodes)
    tmp$ancestors[1] <- 0
    slouch_node <- 1:N
    ancestor <- rep(NA, times = N)
    for (i in 1:N) {
        if (tmp$labels[i] == "" || is.na(tmp$labels[i])) 
		tmp$labels[i] = NA
    }
    rownames(tmp) <- 1:nrow(tmp)
    names(tmp) = c("nodes", "ancestor", "time", "species")
    return(tmp)
}

'.ape2slouch.mvsl'<-function (tree) 
{
## Modified ouch2slouch (author Jason Pienaar) function by Krzysztof Bartoszek
    if (!inherits(tree, "phylo")){stop(sQuote("tree"), " must be of class ", sQuote("phylo"))}
    if (is.null(tree$time.of.nodes)){node_heights<-ape::node.depth.edgelength(tree)}
    else{node_heights<-tree$time.of.nodes}

    N <- nrow(tree$edge)+1
    tmp<-as.data.frame(matrix(0,ncol=4,nrow=N))
    
    names(tmp) <- c("nodes", "ancestors", "times", "labels")
    
    for (i in 1:N){
	iedge<-which(tree$edge[,2]==i)
	if (length(iedge)==1){tmp$ancestors[i]<-tree$edge[iedge,1]}
    }
    tmp$ancestors <- as.character(tmp$ancestors)
    tmp$ancestors <- as.numeric(tmp$ancestors)
    
    tmp$times<-node_heights
    tmp$times <- as.character(tmp$times)
    tmp$times <- as.numeric(tmp$times)
    
    tmp$nodes <- 1:N
    tmp$nodes <- as.character(tmp$nodes)
    tmp$nodes <- as.numeric(tmp$nodes)

    vnodes_edgecol1<-sort(unique(tree$edge[,1]))
    vnodes_edgecol2<-sort(unique(tree$edge[,2]))
    vnodes_tips<-setdiff(vnodes_edgecol2,vnodes_edgecol1) ## only those nodes that have an ending
    vnodes_internal<-setdiff(1:N,vnodes_tips)

    tmp$labels<-NA
    tmp$labels[vnodes_tips]<-tree$tip.label
    tmp$labels[vnodes_internal]<-vnodes_internal

    root_index<-setdiff(vnodes_edgecol1,vnodes_edgecol2)
    if (length(root_index)>1){stop("Error in the phylogeny! There seems to be more than one root!")}
    
    Ntips<-length(vnodes_tips)
    Ninternal<-length(vnodes_internal)
    
    tmp$nodes[vnodes_tips]<-(-1)*tmp$nodes[vnodes_tips]
    tmp$nodes[root_index]<-1
    tmp$ancestors[which(tmp$ancestors==root_index)]<-1
    vnodes_internal<-setdiff(vnodes_internal,root_index)
    
    tmp$nodes[vnodes_internal]<-2:(length(vnodes_internal)+1)
    for (i in 1:N){
	ianc<-tmp$ancestors[i]
	if ((ianc!=0) && (ianc!=1)){## not root node nor root ancestor
	    tmp$ancestors[i]<-which(vnodes_internal == ianc)+1
	}
    }
    tmp$nodes[which(tmp$nodes<0)]<- (Ninternal+1):(Ninternal+Ntips)
    tmp<-tmp[with(tmp,order(nodes)),]  

    tmp$ancestors[1] <- 0

    slouch_node <- 1:N
    ancestor <- rep(NA, times = N)
    for (i in 1:N) {
        if (tmp$labels[i] == "" || is.na(tmp$labels[i])) 
		tmp$labels[i] = NA
    }
    rownames(tmp) <- 1:nrow(tmp)
    names(tmp) <- c("nodes", "ancestor", "time", "species")
    return(tmp)
}

