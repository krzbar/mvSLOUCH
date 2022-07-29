## This file is part of mvSLOUCH

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this package or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .

.internal_phyltree_paths_BL<-function(phyltree){
## called in regimes.R
    if (!is.element("edge.length",names(phyltree))){.my_stop("ERROR: phylogeny does not have branch lengths. mvSLOUCH assumes continuous time, hence requires a dated phylogeny",TRUE)}
    phyltree_paths(phyltree)
}

phyltree_paths<-function(phyltree){
## called in fitch.mvsl.R, phyltree_paths.R

## root edges are ignored and removed

    vnodes_edgecol1<-sort(unique(phyltree$edge[,1]))
    vnodes_edgecol2<-sort(unique(phyltree$edge[,2]))
    vnodes_tips<-setdiff(vnodes_edgecol2,vnodes_edgecol1) ## only those nodes that are an ending

    if(!is.element("Ntips",names(phyltree))){
    ## tree can be non-binary, so we cannot just count number of edges
	phyltree$Ntips<-length(vnodes_tips)
	if (!is.element("Nnode",names(phyltree))){phyltree$Nnode<-length(vnodes_edgecol1)}
    }
    
    if (is.element("Nnode",names(phyltree))){
        if (phyltree$Nnode!=length(vnodes_edgecol1)){.my_stop("Inconsistency in the number of nodes of the tree!",TRUE)}
    }
    
    vdescnum_internal<-table(phyltree$edge[,1])
    vindex_root_path<-which(vdescnum_internal==1)

    if (length(vindex_root_path)>0){
    ## we have a root edge!
	    if ((is.element("edge.length",names(phyltree)))&&(length(phyltree$edge.length)==nrow(phyltree$edge))){phyltree$edge.length<-phyltree$edge.length[-vindex_root_path]}
	    phyltree$edge<-phyltree$edge[-vindex_root_path,]
    }
    

    if ((is.element("edge.length",names(phyltree))&&(length(phyltree$edge.length)<=1))||(length(phyltree$edge)<=2)){.my_stop("The tree has only 1 tip species. Impossible to work with",TRUE)}
    if (is.element("edge.length",names(phyltree))&&(nrow(phyltree$edge)!=length(phyltree$edge.length))){.my_stop("There is some inconsistency in the edges of the  provided phylogeny!",TRUE)}

    if (!is.element("tip.label",names(phyltree))){phyltree$tip.label<-paste("species_",1:phyltree$Ntips,sep="");.my_warning("The tip species of the phylogeny do not have names! Creating generic ones!",TRUE,FALSE)}	
    
    if ((!is.element("tip_species_index",names(phyltree)))||(!is.element("internal_nodes_index",names(phyltree)))||(!is.element("root_index",names(phyltree)))){
	vnodes_edgecol1<-sort(unique(phyltree$edge[,1]))
	vnodes_edgecol2<-sort(unique(phyltree$edge[,2]))
	vnodes_tips<-setdiff(vnodes_edgecol2,vnodes_edgecol1) ## only those nodes that are an ending
	if (!is.element("tip_species_index",names(phyltree))){phyltree$tip_species_index<-vnodes_tips}
	if (!is.element("internal_nodes_index",names(phyltree))){phyltree$internal_nodes_index<-vnodes_edgecol1}
	if (!is.element("root_index",names(phyltree))){
	    root_index<-setdiff(vnodes_edgecol1,vnodes_edgecol2)
	    if (length(root_index)==1){phyltree$root_index<-root_index}
	    else{.my_stop("Error in the phylogeny! There seems to be more than one root!\n",TRUE)}
	}
	if (phyltree$root_index!=(phyltree$Ntips+1)){
	    .my_stop(paste("Error in phylogeny! Root does not have index Ntips+1=",phyltree$Ntips+1,", but ",phyltree$root_index,". Number of tips is Ntips=",phyltree$Ntips,".\n",sep=""),TRUE)
	}
    }
    
    if(!is.element("path.from.root",names(phyltree))){    
	nnodes<-phyltree$Ntips+phyltree$Nnode[[1]]
	phyltree$path.from.root<-vector("list",phyltree$Ntips) 
	phyltree$path.from.root<-sapply(phyltree$path.from.root,function(x){list(nodes=c(),edges=c())},simplify=FALSE)
	names(phyltree$path.from.root)<-phyltree$tip.label
	phyltree$time.of.nodes<-rep(0,nnodes);
	vnodesdone<-rep(NA,nnodes)


	for (i in phyltree$tip_species_index){
	## we get around here assuming that the numbering of  tip nodes is 1:n as described in ape
	    currnode<-i
	    edge.path<-c()
	    node.path<-c()
	    while(currnode!=(phyltree$Ntips+1)){ ## phyltree$Ntips+1 is the root in ape
		if (is.na(vnodesdone[currnode])){
		    vnodesdone[currnode]<-i
		    node.path<-c(node.path,currnode)
		    edgenum<-which(phyltree$edge[,2]==currnode)
		    edge.path<-c(edge.path,edgenum)
		    currnode<-phyltree$edge[edgenum,1]
		}else{
		    i_index<-vnodesdone[currnode]
		    index_in_i_index<-which(phyltree$path.from.root[[i_index]]$nodes==currnode)
		    node.path<-c(node.path,rev(phyltree$path.from.root[[i_index]]$nodes[1:index_in_i_index]))
		    edgenum<-which(phyltree$edge[,2]==currnode)
	    	    index_in_i_index<-which(phyltree$path.from.root[[i_index]]$edges==edgenum)
		    edge.path<-c(edge.path,rev(phyltree$path.from.root[[i_index]]$edges[1:index_in_i_index]))
		    break()
		}    
	    }
	    if(currnode==phyltree$Ntips+1){node.path<-c(node.path,currnode)}
	    node.path<-rev(node.path)
	    edge.path<-rev(edge.path)
	    phyltree$path.from.root[[i]]$nodes<-node.path
	    phyltree$path.from.root[[i]]$edges<-edge.path

	    if (is.element("edge.length",names(phyltree))){
		## root will have zero path length
		if (i!=phyltree$root_index){phyltree$time.of.nodes[i]<-sum(phyltree$edge.length[edge.path])}
		if (length(node.path)>1){
		    for (j in (length(node.path)-1):1){
			jnode<-node.path[j]
			if (jnode!=phyltree$root_index){## root does not need to be done!
			    if (phyltree$time.of.nodes[jnode]==0){
				if (jnode>1){## if jnode==1 we are at path to root anyway
				    phyltree$time.of.nodes[jnode]<-sum(phyltree$edge.length[edge.path[1:(j-1)]])
				}
			    }else{break()}
			}
		    }
		}
	    }
	}
        phyltree$tree_height<-NA
        phyltree$tree_height<-max(phyltree$time.of.nodes) 
    }
    if (!is.null(phyltree$root.edge)){phyltree$root.edge<-NULL}

    ## $bDone_formvSLOUCH only checked in regimes.R
    phyltree$bDone_formvSLOUCH<-TRUE
    if (!is.element("edge.length",names(phyltree))){phyltree$bDone_formvSLOUCH<-FALSE}
    phyltree    
}

.phyltree_remove_path_fields<-function(phyltree){
## called in loglik.R
    if (is.element("path.from.root",names(phyltree))){phyltree$path.from.root<-NULL}
    if (is.element("time.of.nodes",names(phyltree))){phyltree$time.of.nodes<-NULL}
    if (is.element("Ntips",names(phyltree))){phyltree$Ntips<-NULL}
    if (is.element("tree_height",names(phyltree))){phyltree$tree_height<-NULL}
    if (is.element("bDone_formvSLOUCH",names(phyltree))){phyltree$bDone_formvSLOUCH<-NULL}
    if (is.element("tip_species_index",names(phyltree))){phyltree$tip_species_index<-NULL}
    if (is.element("internal_nodes_index",names(phyltree))){phyltree$internal_nodes_index<-NULL}
    if (is.element("root_index",names(phyltree))){phyltree$root_index<-NULL}
    if (is.element("mTreeDist",names(phyltree))){phyltree$mTreeDist<-NULL}
    if (is.element("mAncestorTimes",names(phyltree))){phyltree$mAncestorTimes<-NULL}
    if (is.element("vSpeciesPairs",names(phyltree))){phyltree$vSpeciesPairs<-NULL}
    if (is.element("likFun_OU",names(phyltree))){phyltree$likFun_OU<-NULL}
    if (is.element("likFun_BM",names(phyltree))){phyltree$likFun_BM<-NULL}
    if (is.element("likFun_BM_kX",names(phyltree))){phyltree$likFun_BM_kX<-NULL}
    if (is.element("likFun_BM_all",names(phyltree))){phyltree$likFun_BM_all<-NULL}
    if (is.element("kX",names(phyltree))){phyltree$kX<-NULL}
    if (is.element("kYX",names(phyltree))){phyltree$kYX<-NULL}
    ## field node.label as from version 1.2.9 PCMBase requires unique node.labels
    ## mvSLOUCH makes no use of this field, while PCMBase permits it to be NULL
    ## hence as it is user provided it is easier just to remove it here
    ## before calling PCMBase's likelihood functions
    if (is.element("node.label",names(phyltree))){phyltree$node.label<-NULL}
    if (is.element("b_usePCMBaseCpp",names(phyltree))){phyltree$b_usePCMBaseCpp<-NULL}
    phyltree
}

.phyltree_remove_tips<-function(phyltree,vtips){
# called in loglik.R, regimes.R 
# function assumes vtips is numeric
    if (length(vtips>0)){
	vedge_regime<-NA
	vedge_names<-NA
	vedge_jump<-NA
	if (is.element("edge.regime",names(phyltree))){vedge_regime<-phyltree$edge.regime}
	if (is.element("edge.jump",names(phyltree))){vedge_jump<-phyltree$edge.jump}	
	if (!is.null(names(phyltree$edge.length))){vedge_names<-names(phyltree$edge.length)}
	
	names(phyltree$edge.length)<-sapply(1:length(phyltree$edge.length),function(i){paste("edge_",i,sep="")},simplify=TRUE)
	if (!is.na(vedge_regime[1])){names(vedge_regime)<-names(phyltree$edge.length)}
	if (!is.na(vedge_jump[1])){names(vedge_jump)<-names(phyltree$edge.length)}
	if (!is.na(vedge_names[1])){names(vedge_names)<-names(phyltree$edge.length)}
	
	phyltree<-ape::drop.tip(phyltree,vtips,trim.internal=FALSE,collapse.singles=FALSE)
	v_new_empty_tips<-which(phyltree$tip.label=="NA")## a new tip could be created if a whole clade is in vtips, e.g. a cherry
	while(length(v_new_empty_tips)>0){## this has to stop as we know that there are tips with measurements in the tree	
	    phyltree<-ape::drop.tip(phyltree,v_new_empty_tips,trim.internal=FALSE,collapse.singles=FALSE)
	    v_new_empty_tips<-which(phyltree$tip.label=="NA")## a new tip could be created if a whole clade is in vtips, e.g. a cherry
	}

	if (!is.na(vedge_regime[1])){phyltree$edge.regime<-vedge_regime[names(phyltree$edge.length)]}
	if (!is.na(vedge_jump[1])){phyltree$edge.jump<-vedge_jump[names(phyltree$edge.length)]}
	if (!is.na(vedge_names[1])){names(phyltree$edge.length)<-vedge_names[names(phyltree$edge.length)]}
	else{names(phyltree$edge.length)<-NULL}
    }
    phyltree
}
