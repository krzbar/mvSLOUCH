## This file is part of mvSLOUCH

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this package or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .


.InitialRegimeSetup<-function(phyltree,regimes,regimes.times,mData,kX,kYX=NULL,root.regime=NULL,M.error=NULL,bSave_in_phyltree=FALSE){
## called in evolmodelest.R PhyloSDE.R, modelparamsummary.R, simulaVasicekprocphyl.R, OUphylregression.R
    bTake_from_phyltree<-FALSE
    vregimes<-NA
    pcmbase_model_box<-NA

    ## we assume the mvslouchRegimeStructure is an intrinsic to mvSLOUCH field and will not be used outside
    ## mvSLOUCH's output does not save the tree so this field should not be provided from the outside

    if (is.element("mvslouchRegimeStructure",names(phyltree))){bTake_from_phyltree<-TRUE}
    if ((!is.element("bDone_formvSLOUCH",names(phyltree)))||(!phyltree$bDone_formvSLOUCH)){phyltree<-.internal_phyltree_paths_BL(phyltree)}
    
    if (!(is.null(names(phyltree$edge.length)) && is.null(phyltree$edge.regime))){
	if (is.null(names(phyltree$edge.length))){names(phyltree$edge.length)<-phyltree$edge.regime}
	if (is.null(phyltree$edge.regime)){phyltree$edge.regime<-names(phyltree$edge.length)}
    }


    bOKregimes<-TRUE    
    if (bTake_from_phyltree){
	regimes.times<-phyltree$mvslouchRegimeStructure$regimes.times
	regimes<-phyltree$mvslouchRegimeStructure$regimes
	regimes.types<-phyltree$mvslouchRegimeStructure$regimes.types
	root.regime<-phyltree$mvslouchRegimeStructure$root.regime
	regimes.types.orig<-phyltree$mvslouchRegimeStructure$regimes.types.orig
	pcmbase_model_box<-phyltree$mvslouchRegimeStructure$pcmbase_model_box
	bOKregimes<-phyltree$mvslouchRegimeStructure$bOKregimes
    }else{	
	if (is.null(regimes.times)){regimes.times<-sapply(phyltree$path.from.root,function(lnodes_lineage,vtime_of_nodes){
    	    vtime_of_nodes[lnodes_lineage$nodes]
	},vtime_of_nodes=phyltree$time.of.nodes,simplify=FALSE)}
    
	nullRegimes<-FALSE
	if (is.null(regimes)){
	    if ((is.null(phyltree$edge.regime))&&(is.null(names(phyltree$edge.length)))){
    		nullRegimes<-TRUE
    		## -1 as the first entry of each element of regimes.times is the root
    		## number of edges is one less than the nodes on the path
    		regimes<-sapply(regimes.times,function(reg){rep("reg.1",length(reg)-1)},simplify=FALSE)        
		names(phyltree$edge.length)<-rep("reg.1",length(phyltree$edge.length))
		phyltree$edge.regime<-names(phyltree$edge.length)
	    }
	    else{
		if (!is.null(names(phyltree$edge.length))){
		    regimes<-names(phyltree$edge.length)
		}
		if (!is.null(phyltree$edge.regime)){
		    if (is.null(regimes)){regimes<-phyltree$edge.regime}
		    else{
			if (!all(regimes==phyltree$edge.regime)){
			    .my_warning("WARNING: In the provided phylogeny the regimes names in field edge.regime do not match those in names of edge.length vector. Using the field edge.regime",TRUE,TRUE)
			    regimes<-phyltree$edge.regime
			    names(phyltree$edge.length)<-regimes
			}
		    }
		}
		
	    }
	}
    }
    
    ## by now regimes should not be NULL
    if (!is.list(regimes)){## The regimes are given as a vector in the phylo tree format, i.e. the index in the vector corresponds to the edge number --- Change them to a list
	## this vector should be IDENTICAL to phyltree$edge.regime and names(phyltree$edge.length)
	## this part of the code should not be run if the regimes were passed through the phylogeny
	if (bTake_from_phyltree){.my_warning("WARNING: error in regimes provided in phylogeny! Trying to correct them!",TRUE,TRUE)}
	if (!((is.vector(regimes))&&(length(regimes)==length(phyltree$edge.length)))){bOKregimes<-FALSE;.my_stop("Wrong data structure for regimes, either vector or list of length equalling number of tree edges.",TRUE)}
	else{
    	    if (is.null(phyltree$edge.regime)){phyltree$edge.regime<-regimes}
    	    else{if (!all(regimes==phyltree$edge.regime)){phyltree$edge.regime<-regimes;.my_warning("WARNING: provided regimes differ from those in phylogenetic tree (field edge.regime). Using provided regimes.",TRUE,TRUE)}}
	    if (is.null(names(phyltree$edge.length))){names(phyltree$edge.length)<-regimes}
	    else{if (!all(regimes==names(phyltree$edge.length))){names(phyltree$edge.length)<-regimes;.my_warning("WARNING: provided regimes differ from those in phylogenetic tree (names of edge.length field). Using provided regimes.",TRUE,TRUE)}}
    	    
    	    if (is.null(root.regime)){
    		v_root_branches<-which(phyltree$edge[,1]==phyltree$root_index)
    		v_cand_regs<-unique(regimes[v_root_branches])
    		if (length(v_cand_regs)>1){
    		    root.regime<-names(which.max(table(regimes)[as.character(v_cand_regs)]))
    		}else{
    		    root.regime<-v_cand_regs
    		}    		
    	    }
    	    if (is.null(names(phyltree$edge.length))){names(phyltree$edge.length)<-regimes}
    	    vregimes<-regimes
    	    regimes<-sapply(phyltree$path.from.root,function(lnodes_lineage,vregimes){
    	    ## unlike taking @lineages from ouch object here we do not need to rev
    		as.character(sapply(lnodes_lineage$edges,function(reg,vregimes){vregimes[reg]},vregimes=vregimes,simplify=TRUE))
    	    },vregimes=vregimes,simplify=FALSE)    	    
	}
    }else{  	  	
        if (is.null(phyltree$edge.regime)&&is.null(names(phyltree$edge.length))){
    	    vregimes<-rep(NA,length(phyltree$edge.length))	    
	    inumedgesnotdone<-length(phyltree$edge.length)	    
	    for(i in 1:length(phyltree$path.from.root)){
		vedgelineage<-rev(phyltree$path.from.root[[i]]$edges)
		inumedgesinlineage<-length(vedgelineage)
		inumdonenow<-length(which(is.na(vregimes[vedgelineage])))
		vregimes[vedgelineage]<-regimes[[i]] 
		
		inumedgesnotdone<-inumedgesnotdone-inumdonenow		
		if (inumedgesnotdone==0){break()} 
	    }
	}	    
    }

## taken from Venelin's PCMBase, regimes as names of $edge.length is a convention from phytools

    if (is.null(phyltree$edge.regime)|| is.null(names(phyltree$edge.length))){
	if ((any(is.na(vregimes))||(length(vregimes!=length(phyltree$edge.length))))&&(!(is.null(phyltree$edge.regime)&& is.null(names(phyltree$edge.length))))){
	    if (!is.null(phyltree$edge.regime)){
		vregimes<-phyltree$edge.regime
		if(is.null(names(phyltree$edge.length))){
		    names(phyltree$edge.length)<-phyltree$edge.regime
		}else{if (!all(vregimes==names(phyltree$edge.length))){names(phyltree$edge.length)<-vregimes;.my_warning("WARNING: provided in tree field $edge.regime differ from those in tree's names of edge.length field. Using provided regimes.",TRUE,TRUE)}}
	    }
	    if ((is.null(phyltree$edge.regime))&&(!is.null(names(phyltree$edge.length)))){
		vregimes<-names(phyltree$edge.length)
		phyltree$edge.regime<-names(phyltree$edge.length)
	    }
	}
	if ((any(is.na(vregimes))||(length(vregimes!=length(phyltree$edge.length))))&&((is.null(phyltree$edge.regime)&& is.null(names(phyltree$edge.length)))))
	{vregimes<-rep("reg.1",length(phyltree$edge.length))}## This line should never actually be done as this correction is done right at the start of the function
	if (is.null(phyltree$edge.regime)){phyltree$edge.regime<-vregimes}
	else{if (!all(vregimes==phyltree$edge.regime)){phyltree$edge.regime<-vregimes;.my_warning("WARNING: provided regimes differ from those in phylogenetic tree (field edge.regime). Using provided regimes.",TRUE,TRUE)}}
	if (is.null(names(phyltree$edge.length))){names(phyltree$edge.length)<-vregimes}
	else{if (!all(vregimes==names(phyltree$edge.length))){names(phyltree$edge.length)<-vregimes;.my_warning("WARNING: provided regimes differ from those in phylogenetic tree (names of edge.length field). Using provided regimes.",TRUE,TRUE)}}
    }
    
    
    if (!bTake_from_phyltree){
	## if they were already taken, then no need to have code running again
	bOKregimes<-.CheckSanityRegimes(phyltree,regimes, regimes.times)
	if (bOKregimes) {
	    regimes.types.orig<-c()
    	    for (i in 1: length(regimes)){regimes.types.orig<-c(regimes.types.orig,unique(regimes[[i]]))}
    	    regimes.types.orig<-sort(unique(regimes.types.orig)) ## regime names are in alphabetical order
    	    regimes.types<-1:length(regimes.types.orig)
    	    regimes<-sapply(regimes,function(vregs,regimes.types.orig){
        	sapply(vregs,function(orgreg,regimes.types.orig){which(regimes.types.orig==orgreg)},regimes.types.orig=regimes.types.orig)
    	    },regimes.types.orig=regimes.types.orig,simplify=FALSE)
	}else{.my_stop("ERROR: There is an error with the regimes. Please see the documentation how to set them up. One possibility is to provide a vector (regimes parameter) of length equalling the number of branches in the tree. Each vector entry should be a character regime name.",TRUE)}
    }
    
    if ((!is.null(M.error))||(!bTake_from_phyltree && !is.null(kYX))){
	lmodelsetup<-.ModelSetupPCMBase(phyltree,M.error,kYX)
	pcmbase_model_box<-lmodelsetup$pcmbase_model_box
	phyltree<-lmodelsetup$phyltree
    }
    
    ## create likelihood functions for faster evaluation
    if (!is.null(mData)){
	##model_OU<-PCMBase::PCM("OU",k=kYX,regimes=regimes.types.orig)
	model_OU<-PCMBase::PCM("OU",k=kYX,regimes=names(pcmbase_model_box$Sigma_x[1,1,]))
        model_BM_kX<-PCMBase::PCM("BM",k=kX,regimes=names(pcmbase_model_box$Sigma_x[1,1,]))
	model_BM_all<-PCMBase::PCM("BM",k=kYX,regimes=names(pcmbase_model_box$Sigma_x[1,1,]))
	phyltree_tmp<-.phyltree_remove_path_fields(phyltree)
	phyltree$likFun_OU <- NA
	phyltree$likFun_BM_kX <- NA
	phyltree$likFun_BM_all <- NA
	
	
	if (requireNamespace("PCMBaseCpp",quietly=TRUE)){
#	    metaICpp_OU <- PCMBaseCpp::PCMInfoCpp(X=t(mData), tree=phyltree_tmp, model=model_OU)
#	    phyltree$likFun_OU <- PCMBase::PCMCreateLikelihood(X=t(mData), tree=phyltree_tmp, model=model_OU,metaI = metaICpp_OU)
#	    metaICpp_BM_kx <- PCMBaseCpp::PCMInfoCpp(X=t(mData[,(kYX-kX+1):kYX,drop=FALSE]), tree=phyltree_tmp, model=model_BM_kX)
#	    phyltree$likFun_BM_kX <- PCMBase::PCMCreateLikelihood(X=t(mData[,(kYX-kX+1):kYX,drop=FALSE]), tree=phyltree_tmp, model=model_BM_kX,metaI = metaICpp_BM_kx)
#	    metaICpp_BM_all <- PCMBaseCpp::PCMInfoCpp(X=t(mData), tree=phyltree_tmp, model=model_BM_all)
#	    phyltree$likFun_BM_all <- PCMBase::PCMCreateLikelihood(X=t(mData), tree=phyltree_tmp, model=model_BM_all,metaI = metaICpp_BM_all)    

    	    vNArows<-which(apply(mData,1,function(x){all(is.na(x))}))
            if (length(vNArows)>0){ 	    
                ## here it is assumed that the order of rows corresponds to the order of tips
                phyltree_tmp<-.phyltree_remove_tips(phyltree_tmp,vNArows)
                mData<-mData[-vNArows,,drop=FALSE]
                mData<-mData[phyltree_tmp$tip.label,,drop=FALSE] ## need to reorder in case species order changed                  
        	## Cpp does not support rows that are all NA but non Cpp version does
        	## not removing rows as then one needs to reorder data in case tip order changed
        	## after dropping tips
	    	## phyltree$likFun_OU <- PCMBase::PCMCreateLikelihood(X=t(mData), tree=phyltree_tmp, model=model_OU)
		## phyltree$likFun_BM_all <- PCMBase::PCMCreateLikelihood(X=t(mData), tree=phyltree_tmp, model=model_BM_all)    
	    }
	    ##else{
		##phyltree$likFun_OU <- PCMBase::PCMCreateLikelihood(X=t(mData), tree=phyltree_tmp, model=model_OU,metaI = PCMBaseCpp::PCMInfoCpp)
		##phyltree$likFun_BM_all <- PCMBase::PCMCreateLikelihood(X=t(mData), tree=phyltree_tmp, model=model_BM_all,metaI = PCMBaseCpp::PCMInfoCpp)    
	    ##}
	    phyltree$likFun_OU <- PCMBase::PCMCreateLikelihood(X=t(mData), tree=phyltree_tmp, model=model_OU,metaI = PCMBaseCpp::PCMInfoCpp)
	    phyltree$likFun_BM_all <- PCMBase::PCMCreateLikelihood(X=t(mData), tree=phyltree_tmp, model=model_BM_all,metaI = PCMBaseCpp::PCMInfoCpp)    
    
	    mData_kx<-mData[,(kYX-kX+1):kYX,drop=FALSE]
	    vNArows<-which(apply(mData_kx,1,function(x){all(is.na(x))}))
            if (length(vNArows)>0){ 	    
                phyltree_tmp<-.phyltree_remove_tips(phyltree_tmp,vNArows)
                mData_kx<-mData_kx[-vNArows,,drop=FALSE]
                mData_kx<-mData_kx[phyltree_tmp$tip.label,,drop=FALSE] ## need to reorder in case species order changed                  
	    	## Cpp does not support rows that are all NA but non Cpp version does
        	## not removing rows as then one needs to reorder data in case tip order changed
        	## after dropping tips
		##phyltree$likFun_BM_kX <- PCMBase::PCMCreateLikelihood(X=t(mData_kx), tree=phyltree_tmp, model=model_BM_kX)
	    }
	    ##else{phyltree$likFun_BM_kX <- PCMBase::PCMCreateLikelihood(X=t(mData_kx), tree=phyltree_tmp, model=model_BM_kX,metaI = PCMBaseCpp::PCMInfoCpp)}
	    phyltree$likFun_BM_kX <- PCMBase::PCMCreateLikelihood(X=t(mData_kx), tree=phyltree_tmp, model=model_BM_kX,metaI = PCMBaseCpp::PCMInfoCpp)
	    phyltree$b_usePCMBaseCpp<-TRUE
	}
	else{
	    phyltree$likFun_OU <- PCMBase::PCMCreateLikelihood(X=t(mData), tree=phyltree_tmp, model=model_OU)
	    phyltree$likFun_BM_kX <- PCMBase::PCMCreateLikelihood(X=t(mData[,(kYX-kX+1):kYX,drop=FALSE]), tree=phyltree_tmp, model=model_BM_kX)
	    phyltree$likFun_BM_all <- PCMBase::PCMCreateLikelihood(X=t(mData), tree=phyltree_tmp, model=model_BM_all)    
	    phyltree$b_usePCMBaseCpp<-FALSE
	}
	## phyltree_tmp, mData is not used after this place in this function
	phyltree$kX<-kX
	phyltree$kYX<-kYX
    }
    ## ========================================================

    
    lres<-NA
    if (!bSave_in_phyltree){
	## remove any regime stuff if they are, i.e. do opposite of what was done previously
	if (is.element("mvslouchRegimeStructure",names(phyltree))){phyltree$mvslouchRegimeStructure<-NULL}
	lres<-list(regimes=regimes,regimes.times=regimes.times,regimes.types=regimes.types,root.regime=root.regime,regimes.types.orig=regimes.types.orig,phyltree=phyltree,pcmbase_model_box=pcmbase_model_box,bOKregimes=bOKregimes)
    }else{## this option is if we are using evolmodelest, i.e. we will be reusing the same things for many model calculations
	phyltree$mvslouchRegimeStructure<-vector("list",7)
	names(phyltree$mvslouchRegimeStructure)<-c("regimes","regimes.times","regimes.types","root.regime","regimes.types.orig","pcmbase_model_box","bOKregimes")
	phyltree$mvslouchRegimeStructure$regimes<-regimes
	phyltree$mvslouchRegimeStructure$regimes.times<-regimes.times
	phyltree$mvslouchRegimeStructure$regimes.types<-regimes.types
	phyltree$mvslouchRegimeStructure$root.regime<-root.regime
	phyltree$mvslouchRegimeStructure$regimes.types.orig<-regimes.types.orig
	phyltree$mvslouchRegimeStructure$pcmbase_model_box<-pcmbase_model_box
	phyltree$mvslouchRegimeStructure$bOKregimes<-bOKregimes
	## phyltree does not need to be copied as it is done already
	lres<-phyltree
    }
    
    lres
}

.ModelSetupPCMBase<-function(phyltree,M_error,kYX){
## called in regimes.R
    
    ## number of unique regimes
    unique_regimes<-unique(phyltree$edge.regime)
    num_regs<-length(unique_regimes)

    arH <- abind::abind(matrix(NA,kYX,kYX), along=3, new.names=list( x=NULL, y=NULL,regime=c(unique_regimes[1])))
    arTheta <- abind::abind(rep(NA,kYX), along=2, new.names=list(xy=NULL,regime=c(unique_regimes[1]) ))
    arSigma_x <- abind::abind(matrix(NA,kYX,kYX), along=3, new.names=list(x=NULL, y=NULL,regime=c(unique_regimes[1]) ))
    arSigmae_x <- abind::abind(matrix(0,kYX,kYX), along=3, new.names=list( x=NULL, y=NULL,regime=c(unique_regimes[1])))

    if (num_regs>1){
	for (i in 2:num_regs){
	    arH <- abind::abind(arH,matrix(NA,kYX,kYX), along=3, new.names=list( x=NULL, y=NULL,regime=c(unique_regimes[1:i])))
	    arTheta <- abind::abind(arTheta,rep(NA,kYX), along=2, new.names=list(xy=NULL,regime=c(unique_regimes[1:i]) ))
	    arSigma_x <- abind::abind(arSigma_x, matrix(NA,kYX,kYX), along=3,new.names=list(x=NULL, y=NULL,regime=c(unique_regimes[1:i]) ))
	    arSigmae_x <- abind::abind(arSigmae_x,matrix(0,kYX,kYX), along=3, new.names=list( x=NULL, y=NULL,regime=c(unique_regimes[1:i])))
	}
    }
    
    if (!is.null(M_error)){
    ## there is measurement error present, otherwise the matrix was just 0s
	if (is.matrix(M_error)){## just the same value for all regimes
	    M_error<-sapply(1:length(phyltree$tip.label),function(i,x){x},x=M_error,simplify=FALSE)
	}	
	vregimes_to_remove<-setdiff(unique(phyltree$edge.regime[which(phyltree$edge[,2]<phyltree$Ntips+1)]),unique(phyltree$edge.regime[which(phyltree$edge[,2]>=phyltree$Ntips+1)]))
    	vcurrArnames<-unique_regimes	    

	for(i in 1:length(M_error)){
	    ## in ape format tips are 1:n	    
	    reg<-phyltree$edge.regime[which(phyltree$edge[,2]==i)]
	    new_regname<-paste(reg,"_merrorregime_node",i,sep="")
	    phyltree$edge.regime[which(phyltree$edge[,2]==i)]<-new_regname
	    vcurrArnames<-c(vcurrArnames,new_regname)
	    arH <- abind::abind(arH,matrix(NA,kYX,kYX), along=3, new.names=list( x=NULL, y=NULL,regime=vcurrArnames ))
	    arTheta <- abind::abind(arTheta,rep(NA,kYX), along=2, new.names=list(xy=NULL,regime=vcurrArnames))
	    arSigma_x <- abind::abind(arSigma_x, matrix(NA,kYX,kYX), along=3,new.names=list( x=NULL, y=NULL,regime=vcurrArnames))
	    arSigmae_x <- abind::abind(arSigmae_x,.changeSigmatoSyy(M_error[[i]],"UpperTri",NULL,NULL,TRUE), along=3, new.names=list( x=NULL, y=NULL,regime=vcurrArnames))
	}
	names(phyltree$edge.length)<-phyltree$edge.regime
	if (length(vregimes_to_remove)>0){
	    ## remove no measurement error regimes that are present only on pendant edges
	    for (reg in vregimes_to_remove){
	        arH <- arH[,,-which(names(arH[1,1,])==reg),drop=FALSE]
	        arTheta <- arTheta[,-which(names(arTheta[1,])==reg),drop=FALSE] 
	        arSigma_x <- arSigma_x[,,-which(names(arSigma_x[1,1,])==reg),drop=FALSE]
	        arSigmae_x <- arSigmae_x[,,-which(names(arSigmae_x[1,1,])==reg),drop=FALSE]
	    }
	}	
    }

    pcmbase_model_box <- list(H=arH[,,,drop=FALSE],
                     Theta=arTheta[,,drop=FALSE],
                     Sigma_x=arSigma_x[,,,drop=FALSE],
                     Sigmae_x=arSigmae_x[,,,drop=FALSE],
                     X0=rep(NA,kYX))
    class(pcmbase_model_box) <- 'OU'
    list(pcmbase_model_box=pcmbase_model_box,phyltree=phyltree)
}


.CheckSanityRegimes<-function(phyltree,regimes, regimes.times){
    bOKregimes<-TRUE
    if ((length(regimes)!=length(regimes.times))||(length(regimes)!=phyltree$Ntips)){
        bOKregimes<-FALSE;.my_stop("ERROR in regimes: wrong number of regimes, regimes.times vectors",TRUE)
    }
    else{        
        vOKtimes<-sapply(1:length(regimes),function(i,regimes,epochs,regimes.times){
            bret<-TRUE
            if (length(regimes[[i]])!=(length(regimes.times[[i]])-1)){bret<-FALSE;.my_stop(paste("ERROR in regimes: problem with ",i,"th regimes or regimes.times",sep=""),TRUE)}
            if (length(regimes.times[[i]])<length(epochs$nodes[[i]])){bret<-FALSE;.my_stop(paste("ERROR in regimes: ",i,"th regimes.times does not agree with phylogenetic tree",sep=""),TRUE)}
            if (length(regimes[[i]])<(length(epochs$nodes[[i]])-1)){bret<-FALSE;.my_stop(paste("ERROR in regimes: ",i,"th regimes does not agree with phylogenetic tree",sep=""),TRUE)}
            bret            
        },regimes=regimes,epochs=phyltree$paths.from.root,regimes.times=regimes.times)
        if (length(vOKtimes[vOKtimes])!=length(vOKtimes)){bOKregimes<-FALSE}
    }
    bOKregimes
}


.createRegimes<-function(PhylTree,vRegimes){
    if (!is.element("path.from.root",names(PhylTree))){PhylTree<-phyltree_paths(PhylTree)}

    lRegimes<-vector("list",5)
    names(lRegimes)<-c("regimes","regimeTypes","regimeTimes","regimeCoding","regimeVect")
    lRegimes$regimeVect<-vRegimes

    lRegimes$regimeTimes<-sapply(PhylTree$path.from.root,function(lnodes_lineage,vtime_of_nodes){
	vtime_of_nodes[lnodes_lineage$nodes]
    },vtime_of_nodes=PhylTree$time.of.nodes,simplify=FALSE)


    vRegimeTypes<-unique(vRegimes)
    lRegimes$regimeTypes<-as.character(1:length(vRegimeTypes))
    lRegimes$regimeCoding<-cbind(vRegimeTypes,as.character(1:length(vRegimeTypes)))

    ## vRegimes now correspond to the edge number on which it is present
    lRegs<-sapply(PhylTree$path.from.root,function(epochObj,Regs){vLin<-epochObj$edges;vRegs<-rep(NA,length(vLin)-1);if (length(vLin)>1){for(i in 1:length(vRegs)){vRegs[i]<-Regs[vLin[i]]}};vRegs },Regs=as.character(vRegimes),simplify=FALSE)
    ## path from root is only for tips unlike in the ouch tree case    
    ##lRegs<-lRegs[(PhylTree@nnodes-PhylTree@nterm+1):PhylTree@nnodes]
    
    ## this is only the translation from human regimes to numeric
    lRegs<-sapply(lRegs,function(vReg,mCode){LinRegs<-sapply(vReg,function(r,mCode){mCode[which(mCode[,1]==r),2]},mCode=mCode,simplify=TRUE);names(LinRegs)<-NULL;LinRegs},mCode=lRegimes$regimeCoding,simplify=FALSE)
    lRegimes$regimes<-lRegs
    lRegimes
}

 
 
