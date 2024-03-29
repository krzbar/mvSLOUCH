simulate_clustered_phylogeny<-function(v_sizeclusts,joining_branchlengths=NULL,f_simclustphyl="sim.bd.taxa_Yule1",joiningphyl=NULL,b_change_joining_branches=FALSE,...){
## The function simulates a phylogeny that contains a specified number of clades each of a specified size.
## The resulting object is an enhanced phylo class object. The resulting tree will nearly always be NON-ultrametric.
## Input parameters
## v_sizeclusts: a vector with the sizes of the clades/clusters
## joining_branchlengths: default NULL, if joiningphyl is NULL, then has to be provided. A vector of two numbers.
##	the first element are the lenghts of the branches of the cluster joining phylogeny leading to the clusters.
##	The second element will be the lengths of the "internal" branches of the cluster joining phylogeny. 
##	If only a single number is provided, then all the branches of the joining phylogeny will have their
##	lengths equal to this value. 
## f_simclustphyl: what function to use to simulate the phylogeny inside each cluster. 
##	The default value of "sim.bd.taxa_Yule1" corresponds to a pure birth tree generated by
##	ape::rphylo(n=clade_size,birth=1,death=0), without a root branch
##	otherwise the user should pass an object of class "function" and its parameters in place of the ... .
##	The first parameter must be the number of contemporary leaves and be called n.
##	The function has to return a valid phylo object.
## joiningphyl: By what phylogeny are the clades to be joined by. 
##	Either NULL (default), a phylo object, the character string "sim.bd.taxa_Yule1" or an object of class "function". 
##	If NULL, then they are joined by a 
##	caterpillar/comb/pectinate phylogeny with the branch lengths as provided by the joining_branchlengths parameter.
##	If it is a phylo object, then they will be joined by it. Importantly the number of tips of this phylogeny has to 
##	equal the number of clusters. If "sim.bd.taxa_Yule1", then the joining phylogeny is simulated as a pure birth tree
##	with tips equalling the number of clusters by ape::rphylo(). If it is a function, then this is used
##	and its parameters are passed through ... . The first parameter must be the number of contemporary leaves 
##	and be called n. The function has to return a valid phylo object.
## b_change_joining_branches: logical. If joining phylogeny (parameter joiningphyl) was provided or simulates, should its
##	branches be changed according to what was provided in joining_branchlengths (if it was not NULL). By default FALSE
##	and the branch lengths are not changed. 
## ... : parameters to be passed to user provided f_simclustphyl and joiningphyl functions. Unless one knows exactly
##	what one is doing they should be passed by name. If there is a conflict of names, then one should pass wrapper
##	functions around these functions where the names conflict is resolved.
## Value
## The resulting object is a clustered_phylo object which inherits from the  phylo class and enhances it. 
## Apart from the standard phylo fields it has two additional ones:
## edges_clusters: a named list with length equalling the number of clades/clusters plus 1. The first element
##	of the list is called joining_tree and contains the indices (row numbers of the edge matrix, indices of the edge_length vector) 
##	of the edges inside the subtree joining the clusters. Afterwords element (i+1) is named cluster_i and contains a numeric
##	vector with the indices of the edges inside clade i.
## tips_clusters: a named list with length equalling the number of clades/clusters. Each field of the list is a numeric
##	vector containing the indices of the tips inside the clade. The names of element i of the list is cluster_i.
## Example call
## We use wrapper function for illustration
## my_sim.bd.taxa<-function(n,...){
##   ape::rphylo(n=n,...)
## }
##	f_simClustPhyl(v_sizeclusts=c(5,5,5),f_simclustphyl=my_sim.bd.taxa,b_change_joining_branches=TRUE,joining_branchlengths=c(20,NA),joining=my_sim.bd.taxa,lambda=1,mu=0)

    if (is.null(v_sizeclusts)){stop("Size clusters not provided!")}
    numclusts<-length(v_sizeclusts)
    if ((inherits(joiningphyl,"phylo"))&&(numclusts!=(joiningphyl$Nnode+1))){stop("Number of tips of connecting phylogeny does not equal the number of clusters.")}
    
    
    if (numclusts==1){
	if ((is.character(f_simclustphyl))&&(f_simclustphyl=="sim.bd.taxa_Yule1")){
		p1<-ape::rphylo(n=v_sizeclusts[1],birth=1,death=0)
    		p_joined<-p1;p_joined$root.edge<-joining_branchlengths[1];return(p_joined)
    	}
    	else if(is.function(f_simclustphyl)){
    	    p1<-f_simclustphyl(n=v_sizeclusts[1],...)
    	    if (!inherits(p1,"phylo")){
    		stop("The simulation function passed through f_simclustphyl does not return a phylo object!")
    	    }    	    
    	    p_joined<-p1;p_joined$root.edge<-joining_branchlengths[1];return(p_joined)
    	}else{stop("Cannot simulate the phylogeny with the provided f_simclustphyl.")}
    }
    if (!is.null(joining_branchlengths)){
        internal_joining_branchlength<-NA
        if (length(joining_branchlengths)==2){
    	    internal_joining_branchlength<-joining_branchlengths[[2]]
	    joining_branchlengths<-joining_branchlengths[[1]]
	}else{
	    if (length(joining_branchlengths)==1){internal_joining_branchlength<-joining_branchlengths}
	    else{stop("At the moment joining_branchlengths cannot have length more than 2!")}
	}
    }
    if (is.null(joiningphyl)){
        ## if no joining phylogeny provided, then caterpillar phylogeny as joining
        ## create fully unbalanced caterpillar tree
        if (is.null(joining_branchlengths)){stop("If joining phylogeny is not provided, joining_branchlengths is required!")}
        joiningphyl<-list(Nnode=numclusts-1,edge=matrix(NA,nrow=2*numclusts-2,ncol=2),edge.length=rep(NA,2*numclusts-2),tip.label=paste(paste0("t",1:numclusts)))

        if (numclusts>2){
    	    joiningphyl$edge[1,]<-c(numclusts+2,1)
    	    joiningphyl$edge.length[1]<-joining_branchlengths
    	    curredge<-2
    	    currintnode<-numclusts+2
    	    for (i in 2:(numclusts-1)){
    		joiningphyl$edge[curredge,]<-c(currintnode,i)
		joiningphyl$edge.length[curredge]<-joining_branchlengths
		joiningphyl$edge[(curredge+1),]<-c(currintnode+1,currintnode)
		if (is.na(internal_joining_branchlength)){
		    stop("The second element of joining_branchlengths cannot be NA if the joining phylogeny is NOT provided.")
		}
		joiningphyl$edge.length[curredge+1]<-internal_joining_branchlength
    		curredge<-curredge+2
		currintnode<-currintnode+1		
	    }
	    joiningphyl$edge[curredge,]<-c(numclusts+1,numclusts)
	    joiningphyl$edge[which(joiningphyl$edge==(2*numclusts))]<-numclusts+1
	    joiningphyl$edge.length[curredge]<-joining_branchlengths
	}else{
	    joiningphyl$edge<-rbind(c(3,1),c(3,2))
    	    joiningphyl$edge.length<-rep(joining_branchlengths,2)    	    
	}
	
	class(joiningphyl)<-"phylo"
    }else{
        if (!inherits(joiningphyl,"phylo")){
    	    if((inherits(joiningphyl,"character")) && (joiningphyl=="sim.bd.taxa_Yule1")){
		joiningphyl<-ape::rphylo(n=numclusts,birth=1,death=0)
	    }
	    else if(is.function(joiningphyl)){
    		joiningphyl<-joiningphyl(n=numclusts,...)
    		if (!inherits(joiningphyl,"phylo")){
    		    stop("The simulation function passed through joiningphyl does not return a phylo object!")
    		}    	        		
    	    }else{stop("Cannot simulate the phylogeny with the provided joiningphyl.")}
	}
	joiningphyl$root.edge<-NULL
        if (b_change_joining_branches){
            if ((!is.null(joining_branchlengths))&&(is.numeric(joining_branchlengths))){
		v_tip_branches<-which(joiningphyl$edge[,2]<(numclusts+1))
		joiningphyl$edge.length[v_tip_branches]<-joining_branchlengths[[1]]
		if (!is.na(internal_joining_branchlength)){
		    v_non_tip_branches<-setdiff(1:length(joiningphyl$edge.length),v_tip_branches)
		    if (length(v_non_tip_branches)>0){
			    joiningphyl$edge.length[v_non_tip_branches]<-internal_joining_branchlength
		    }	
		}	
	    }else{message("Cannot change branches leading to clusters as joining_branchlengths not provided")}
	}
    }

    numtipscurr<-0
    p_joined<-joiningphyl
    numalltips<-sum(v_sizeclusts)
    p_joined$edge[which(p_joined$edge>(numclusts+1))]<-p_joined$edge[which(p_joined$edge>(numclusts+1))]+numalltips
    p_joined$edge[which(p_joined$edge==(numclusts+1))]<-numalltips+1 ## final root node
    p_joined$edge[which(p_joined$edge<=numclusts)]<-p_joined$edge[which(p_joined$edge<=numclusts)]+numalltips+1
    labelfirstnodeinsubtree<-numalltips+joiningphyl$Nnode+length(joiningphyl$tip.label)+1
    p_joined$edges_clusters<-vector("list",numclusts+1)
    p_joined$tips_clusters<-vector("list",numclusts)
    names(p_joined$edges_clusters)<-c("joining_tree",paste0("cluster_",1:numclusts))
    names(p_joined$tips_clusters)<-paste0("cluster_",1:numclusts)
    p_joined$edges_clusters$joining_tree<-1:nrow(p_joined$edge)
    
    for (i in 1:numclusts){
        if ((is.character(f_simclustphyl))&&(f_simclustphyl=="sim.bd.taxa_Yule1")){
    	    pnext<-ape::rphylo(n=v_sizeclusts[i],birth=1,death=0)
    	}
    	else if(is.function(f_simclustphyl)){
    	    pnext<-f_simclustphyl(n=v_sizeclusts[i],...)
    	    if (!inherits(pnext,"phylo")){
    		stop("The simulation function passed through f_simclustphyl does not return a phylo object!")
    	    }    	    
    	}else{stop("Cannot simulate the phylogeny with the provided f_simclustphyl.")}
        pnext$root.edge<-NULL    
	    
        ## p_joined has correct numbering of all nodes only need to renumber nodes in pnext
        ## renumber in pnext internal
        intnodes_torelabel<-which(pnext$edge>(v_sizeclusts[i]+1))
        if (length(intnodes_torelabel)>0){
    	    pnext$edge[intnodes_torelabel]<-pnext$edge[intnodes_torelabel]+labelfirstnodeinsubtree-(v_sizeclusts[i]+2)
    	    labelfirstnodeinsubtree<-labelfirstnodeinsubtree+v_sizeclusts[i]-2
	}
        ## renumber in pnext root
        pnext$edge[which(pnext$edge==(v_sizeclusts[i]+1))]<-numalltips+i+1
	    
        ## renumber tips of pnext
        p_joined$tips_clusters[[i]]<-sort(unique(pnext$edge[which(pnext$edge<=v_sizeclusts[i])]))+numtipscurr
        pnext$edge[which(pnext$edge<=v_sizeclusts[i])]<-pnext$edge[which(pnext$edge<=v_sizeclusts[i])]+numtipscurr

        numtipscurr<-numtipscurr+v_sizeclusts[i]
        start_edge_cluster<-nrow(p_joined$edge)+1
        p_joined$edge<-rbind(p_joined$edge,pnext$edge)
        end_edge_cluster<-nrow(p_joined$edge)
        p_joined$edges_clusters[[i+1]]<-start_edge_cluster:end_edge_cluster
        p_joined$edge.length<-c(p_joined$edge.length,pnext$edge.length)
    }		
    ## change p_joined$Nnode
    p_joined$Nnode<-length(unique(p_joined$edge[,1]))
    ## label tips correctly
    p_joined$tip.label<-paste0("t",1:numalltips)
    class(p_joined)<-c("clustered_phylo","phylo")
    p_joined
}

plot.clustered_phylo<-function(x,clust_cols=NULL,clust_edge.width=NULL,clust_edge.lty=NULL,clust_tip.color="black",joiningphylo_col="black",joiningphylo_edge.width=1,joiningphylo_edge.lty=1,...){
    if (is.null(joiningphylo_col)){joiningphylo_col<-"black"}
    vedgecol<-rep(joiningphylo_col,length(x$edge.length))
    if (!is.null(clust_cols)){
	if (length(clust_cols)!=(length(x$edges_clusters)-1)){
	## recycle colours
	    opt_warn<-options("warn")
	    options(warn=-1)
	    new_clust_cols<-rep(NA,length(x$edges_clusters)-1)
	    new_clust_cols[]<-clust_cols
	    clust_cols<-new_clust_cols
	    options(warn=opt_warn$warn)
	}
	for (i in 2:length(x$edges_clusters)){
	    vedgecol[x$edges_clusters[[i]]]<-clust_cols[i-1]
	}
    }
    if (is.null(joiningphylo_edge.width)){joiningphylo_edge.width<-1}
    vedgewidth<-rep(joiningphylo_edge.width,length(x$edge.length))
    if (!is.null(clust_edge.width)){
	if (length(clust_edge.width)!=(length(x$edges_clusters)-1)){
	## recycle colours
	    opt_warn<-options("warn")
	    options(warn=-1)
	    new_clust_edge.width<-rep(NA,length(x$edges_clusters)-1)
	    new_clust_edge.width[]<-clust_edge.width
	    clust_edge.width<-new_clust_edge.width
	    options(warn=opt_warn$warn)
	}
	for (i in 2:length(x$edges_clusters)){
	    vedgewidth[x$edges_clusters[[i]]]<-clust_edge.width[i-1]
	}
    }
    if (is.null(joiningphylo_edge.lty)){joiningphylo_edge.lty<-1}
    vedgelty<-rep(joiningphylo_edge.lty,length(x$edge.length))
    if (!is.null(clust_edge.lty)){
	if (length(clust_edge.lty)!=(length(x$edges_clusters)-1)){
	## recycle colours
	    opt_warn<-options("warn")
	    options(warn=-1)
	    new_clust_edge.lty<-rep(NA,length(x$edges_clusters)-1)
	    new_clust_edge.lty[]<-clust_edge.lty
	    clust_edge.lty<-new_clust_edge.lty
	    options(warn=opt_warn$warn)
	}
	for (i in 2:length(x$edges_clusters)){
	    vedgelty[x$edges_clusters[[i]]]<-clust_edge.lty[i-1]
	}
    }
    if (is.null(clust_tip.color)){clust_tip.color<-"black"}
    if (length(clust_tip.color)!=length(x$tips_clusters)){
	## recycle colours
	    opt_warn<-options("warn")
	    options(warn=-1)
	    new_clust_tip.color<-rep(NA,length(x$tips_clusters))
	    new_clust_tip.color[]<-clust_tip.color
	    clust_tip.color<-new_clust_tip.color
	    options(warn=opt_warn$warn)
    }
    vtipcolor<-rep(NA,length(x$tip.label))
    for (i in 1:length(x$tips_clusters)){
	    vtipcolor[x$tips_clusters[[i]]]<-clust_tip.color[i]
    }
    
    class(x)<-"phylo"
    ape::plot.phylo(x,edge.color=vedgecol,edge.width=vedgewidth,edge.lty=vedgelty,tip.color=vtipcolor,...)
}
