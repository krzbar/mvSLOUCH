## This file is part of mvSLOUCH

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this package or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .


.calcESSanalytical<-function(phyltree,proc.params,evolmodel,mData,Merror=NULL,vNAs=NULL,ESS.method="reg"){
## called from evolmodelest.R for model selection purposes

## mX is not needed directly, but it is needed to get detV

    if ((evolmodel=="mvslouch") && (.is_det0(proc.params$A))){
	evolmodel<-"ouch"
	proc.params<-.mvslouch_to_ouch_model(proc.params)
    }

    kYX<-ncol(mData);kX<-kYX
    if (evolmodel=="mvslouch"){kX<-ncol(proc.params$B)}
    regimesList<-.InitialRegimeSetup(phyltree=phyltree,regimes=NULL,regimes.times=NULL,mData,kX=kX,kYX=kYX,root.regime=NULL,M.error=Merror)
    bOKregimes<-regimesList$bOKregimes
    if (!bOKregimes){.my_stop("The regimes, regimes.times and phyltree variables are not consistant. Cannot perform estimation. Please see manual on their format.",TRUE)}
    regimes<-regimesList$regimes
    regimes.times<-regimesList$regimes.times
    regimes.types<-regimesList$regimes.types
    root.regime<-regimesList$root.regime
    regimes.types.orig<-regimesList$regimes.types.orig
    phyltree<-regimesList$phyltree
    proc.params$pcmbase_model_box<-regimesList$pcmbase_model_box ## this is the model object in which parameters for all regimes sit in for PCMBase
    proc.params$pcmbase_model_box<-.update_pcmbase_box_params(proc.params,evolmodel)

    
    n<-phyltree$Ntips
    bReturnPhylTree<-FALSE
    if (is.null(Merror) || all(is.na(Merror))) {Merror<-0}
    if (evolmodel=="slouch"){evolmodel<-"mvslouch"}
    if ((evolmodel=="bm") || (evolmodel=="ouch") || (evolmodel=="mvslouch")){
        res<-.getMVSLphylCovMatrix(phyltree,proc.params,evolmodel,regimes.times,regimes,regimes.types,Merror)
        V<-res$mCovPhyl.Merror
	procdim<-res$procdim
	if (is.element("phyltree",names(res))){phyltree<-res$phyltree;bReturnPhylTree<-TRUE}
	rm(res);gc(verbose=FALSE)
    }

    logdetV<-.detV(mX=mData, phyltree=phyltree, pcmbase_model_box=proc.params$pcmbase_model_box,log=TRUE)
    vspeciesvars<-diag(V)
    vobsspeciesvars<-vspeciesvars
    lspeciesvars<-sapply(1:n,function(i,procdim,V){V[(((i-1)*procdim+1):(i*procdim)),(((i-1)*procdim+1):(i*procdim))]},procdim,V=V,simplify=FALSE)
    lobsspeciesvars<-lspeciesvars
    if ((!is.null(vNAs)) && (length(vNAs)>0)){
	vobsspeciesvars<-vobsspeciesvars[-vNAs];
	for (i in 1:nrow(mData)){
	    if (length(is.na(mData[i,]))>0){
		iNAs<-which(is.na(mData[i,]))
		lobsspeciesvars[[i]]<-lobsspeciesvars[[i]][-iNAs,-iNAs,drop=FALSE]
## later in the code we create a variable  full.NA.species to store which species are completely NA
	    }
	}	
    }
    
    orgV<-V
    if ((!is.null(vNAs)) && (length(vNAs)>0)){V<-V[-vNAs,-vNAs,drop=FALSE]}
    samplesize<-n*procdim
    if ((!is.null(vNAs)) && (length(vNAs)>0)){samplesize<-samplesize-length(vNAs)}
    full.NA.species<-which(sapply(1:n,function(i,d,vNAs){res<-FALSE;if(length(intersect((((i-1)*d+1):(i*d)),vNAs))==d){res<-TRUE};res},d=procdim,vNAs=vNAs,simplify=TRUE))
    n.noNAs<-n - length(full.NA.species)

## samplesize should equal nV if no NAs!
    
    sumdiagV<-sum(vobsspeciesvars)
    sumlogdiagV<-sum(log(vobsspeciesvars))
    
    sumoffdiagV<-sum(V[upper.tri(V,diag=FALSE)])
    sumoffdiagabsV<-sum(abs(V[upper.tri(V,diag=FALSE)]))
    sumabsV<-sum(abs(V))
    sumV<-sum(V)
    nV<-nrow(V)
    corrV<-cov2cor(V)    
    sumoffdiagabscorrV<-sum(abs(corrV[upper.tri(corrV,diag=FALSE)]))
    rhoabsV<-mean(abs(corrV[upper.tri(corrV,diag=FALSE)]),na.rm=TRUE)

    ESS.factor<-0
    ESS<-1
    ESS.factor.model.selection<-0
    ESS.model.selection<-1



    ESS.factor.MI<-0
    ESS.MI<-1
    ESS.factor.MI<-1/log(exp(1)+(sumlogdiagV-logdetV)/2)
    ESS.MI<-(samplesize-1)*ESS.factor.MI+1
    if (ESS.factor.MI<0){
	.my_warning("Warning : negative MI ESS factor.\n",TRUE,FALSE)
	ESS.factor.MI<-0
	ESS.MI<-1
    }
    if (ESS.method=="MI"){ESS.factor<-ESS.factor.MI;ESS<-ESS.MI}

    ESS.factor.mean<-0
    ESS.mean<-1
    
    glsmodel<-NA ## for now we have this as NA, maybe in the future for some conditional model
    pcmbase_model_box_mean0<-.set_mean0_pcmbase_model_box(proc.params$pcmbase_model_box,glsmodel=glsmodel)

    tryCatch({
        ## we use the transformation between a correlation and covariance matrix R=(solve(diag(sdvalues))%*%S%*%solve(diag(sdvalues)))
        mD1sd<-matrix(1,ncol=procdim,nrow=n)
	mD1sd<-mD1sd*matrix(sqrt(vspeciesvars),ncol=procdim,nrow=n,byrow=TRUE)
	ESS.mean<-.pcmbaseDphylGaussian_RSS(mD1sd,phyltree,pcmbase_model_box_mean0)
	ESS.factor.mean<-(ESS.mean-1)/(samplesize-1)
	if (ESS.factor.mean<0){
	    .my_warning("Warning : negative mean ESS factor.\n",TRUE,FALSE)
	    ESS.factor.mean<-0
	    ESS.mean<-1
	}
    },error=function(e){.my_message("Problem with inverting between species correlation matrix",TRUE);.my_message(e,TRUE);.my_message("\n",TRUE)})
    if (ESS.method=="mean"){ESS.factor<-ESS.factor.mean;ESS<-ESS.mean;ESS.factor.model.selection<-ESS.factor.mean;ESS.model.selection<-ESS.mean}

    ESS.factor.reg<-0
    ESS.reg<-1
    
##    Here we have to calculate V explicitely
    tryCatch({
	vr.calc<-0
	for (species_i in 1:n){for(trait_k in 1:procdim){
	    i<-(species_i-1)*procdim+trait_k
	    if ((is.null(vNAs))||((!is.null(vNAs))&&(!is.element(i,vNAs)))){		
    		mDCi<-matrix(orgV[i,],ncol=procdim,nrow=n,byrow=TRUE)
		mDCi[species_i,trait_k]<-NA
		mDCi[which(is.na(mData))]<-NA
		regcalcedval<-1-.pcmbaseDphylGaussian_RSS(mDCi,phyltree,pcmbase_model_box_mean0)/vspeciesvars[i]		
		if (regcalcedval<0){.my_warning(paste("Error in calculating regression ESS for species ",species_i," and trait ",trait_k," --- negative conditional variance (",regcalcedval,"), taking absolute value.\n"),TRUE,FALSE);regcalcedval<-abs(regcalcedval)}
		vr.calc<-vr.calc+regcalcedval
	    }
	}}
	ESS.factor.reg<-vr.calc/samplesize
	ESS.reg<-(samplesize-1)*ESS.factor.reg+1
	if (ESS.factor.reg<0){
	    .my_warning("Warning : negative regression ESS factor.\n",TRUE,FALSE)
	    ESS.factor.reg<-0
	    ESS.reg<-1	    
	}
    },error=function(e){.my_message("Problem with calculating regression ESS.",TRUE);.my_message(e,TRUE);.my_message("\n",TRUE)})
    if (ESS.method=="reg"){ESS.factor<-ESS.factor.reg;ESS<-ESS.reg;ESS.factor.model.selection<-ESS.factor.reg;ESS.model.selection<-ESS.reg;}


    ESS.factor.mvMI<-ESS.factor.MI ##0
    ESS.mvMI<-ESS.MI ##1
    if (procdim>1){
	sumlogdetVj<-0
	for (i in setdiff(1:n,full.NA.species)){
	    Vj<-lobsspeciesvars[[i]]
	    logdetVj<-log(det(Vj))
	    sumlogdetVj<-sumlogdetVj+logdetVj
	}
	ESS.mvMI<-1
	ESS.factor.mvMI<-1/log(exp(1)+(sumlogdetVj-logdetV)/2)
	ESS.mvMI<-(n.noNAs-1)*ESS.factor.mvMI+1
        if (ESS.method=="mvMI"){ESS.factor<-ESS.factor.mvMI;ESS<-ESS.mvMI;ESS.factor.model.selection<-ESS.factor.MI;ESS.model.selection<-ESS.MI}
    }

    ESS.factor.mvreg<-ESS.factor.reg ##0
    ESS.mvreg<-ESS.reg ##1
    if (procdim>1){
	mvvr.calc<-0
	for (species_i in 1:n){
	    if ((length(full.NA.species)==0)||((length(full.NA.species)>0)&&(!is.element(species_i,full.NA.species)))){
		## this is a species with at least 1 observation
		mCi<-orgV[,((species_i-1)*procdim+1):(species_i*procdim)]
		mCi[((species_i-1)*procdim+1):(species_i*procdim),]<-NA
		mCVC<-matrix(NA,procdim,procdim)
		lmCikk<-vector("list",procdim)
		for (trait_k in 1:procdim){
		    i<-(species_i-1)*procdim+trait_k
		    if ((is.null(vNAs))||(!is.element(i,vNAs))){				
			mCikk<-matrix(mCi[,trait_k],ncol=procdim,nrow=n,byrow=TRUE)
			mCikk[which(is.na(mData))]<-NA
			lmCikk[[trait_k]]<-mCikk
			mCVC[trait_k,trait_k]<-.pcmbaseDphylGaussian_RSS(mCikk,phyltree,pcmbase_model_box_mean0)
		    }		    
		}
		if (procdim>1){
		    for (trait_k1 in 1:(procdim-1)){
			i1<-(species_i-1)*procdim+trait_k1
			for(trait_k2 in (trait_k1+1):procdim){
    			    i2<-(species_i-1)*procdim+trait_k2
			    if ((is.null(vNAs))||((!is.element(i1,vNAs))&&(!is.element(i2,vNAs)))){				
				mCik1mk2<-lmCikk[[trait_k1]]-lmCikk[[trait_k2]]
				k1mk2val<-.pcmbaseDphylGaussian_RSS(mCik1mk2,phyltree,pcmbase_model_box_mean0)
	    			mCVC[trait_k1,trait_k2]<-(mCVC[trait_k1,trait_k1]+mCVC[trait_k2,trait_k2]-k1mk2val)/2
	    			mCVC[trait_k2,trait_k1]<-mCVC[trait_k1,trait_k2]
			    }			
			}
		    }	    
		}
		mVi<-lobsspeciesvars[[species_i]]
		if (length(which(is.na(mData[species_i,])))>0){
		    iNAs<-which(is.na(mData[species_i,]))
		    mCVC<-mCVC[-iNAs,-iNAs,drop=FALSE]
		}
		mRegESSi<-diag(1,nrow=nrow(mVi),ncol=ncol(mVi))-solve(mVi)%*%mCVC
		mRegESSi<-(t(mRegESSi)+mRegESSi)/2
		## symmetric is immediate, but check for no errors
		## this is a non-degenarate variance matrix so cannot be semi-definite!
		if ((!.matrixcalc_is.symmetric.matrix(mRegESSi)) || (!.matrixcalc_is.positive.definite(mRegESSi))){
		    .my_warning(paste("Error in calculating multivariate regression ESS for species ",species_i, "--- variance matrix not positive-definite. Adjusting it using Matrix::nearPD(). \n"),TRUE,FALSE)
		    mRegESSi<-as.matrix(Matrix::nearPD(mRegESSi)$mat)
		}
		regcalcedval<-det(mRegESSi)
		if (regcalcedval<0){.my_warning(paste("Error in calculating multivariate regression ESS for species ",species_i, "--- negative conditional generalized variance (",regcalcedval,"), taking absolute value. \n"),TRUE,FALSE);regcalcedval<-abs(regcalcedval)}
		mvvr.calc<-mvvr.calc+regcalcedval
	    }
	}
	ESS.factor.mvreg<-mvvr.calc/n.noNAs
	ESS.mvreg<-(n.noNAs-1)*ESS.factor.mvreg+1

	if (ESS.factor.mvreg<0){
	    .my_warning("Warning : negative multivariate regression ESS factor.\n",TRUE,FALSE)
	    ESS.factor.mvreg<-0
	    ESS.mvreg<-1	    
	}
	ESS.mvreg<-(n.noNAs-1)*ESS.factor.mvreg+1
        if (ESS.method=="mvreg"){ESS.factor<-ESS.factor.mvreg;ESS<-ESS.mvreg;ESS.factor.model.selection<-ESS.factor.reg;ESS.model.selection<-ESS.reg}
    }
    rhon<-matrix(0,procdim,procdim)
    if (procdim==1){rhon[1,1]<-mean(c(corrV[upper.tri(corrV,diag=FALSE)]),na.rm=TRUE)}
    else{
	for(i in 1:procdim){
	    for(j in 1:procdim){
		itraits<-which(((1:nrow(corrV))-1)%%procdim+1 == i)
		jtraits<-which(((1:nrow(corrV))-1)%%procdim+1 == j)
		corrVij<-corrV[itraits,jtraits]
		rhon[i,j]<-mean(c(corrVij[upper.tri(corrVij,diag=FALSE)],na.rm=TRUE))
	    }
	}
    }
    if (bReturnPhylTree){lres<-list(ESS.factor=ESS.factor,ESS=ESS,ESS.factor.model.selection=ESS.factor,ESS.model.selection=ESS,rhon=rhon,ESS.coeffs=list(ESS.MI=list(ESS.MI.factor=ESS.factor.MI,ESS.MI=ESS.MI,ESS.factor.mvMI=ESS.factor.mvMI,ESS.mvMI=ESS.mvMI),ESS.mean=list(ESS.mean.factor=ESS.factor.mean,ESS.mean=ESS.mean),ESS.reg=list(ESS.reg.factor=ESS.factor.reg,ESS.reg=ESS.reg,ESS.mvreg.factor=ESS.factor.mvreg,ESS.mvreg=ESS.mvreg)),phyltree=phyltree)}
    else{lres<-list(ESS.factor=ESS.factor,ESS=ESS,ESS.factor.model.selection=ESS.factor,ESS.model.selection=ESS,rhon=rhon,ESS.coeffs=list(ESS.MI=list(ESS.MI.factor=ESS.factor.MI,ESS.MI=ESS.MI,ESS.factor.mvMI=ESS.factor.mvMI,ESS.mvMI=ESS.mvMI),ESS.mean=list(ESS.mean.factor=ESS.factor.mean,ESS.mean=ESS.mean),ESS.reg=list(ESS.reg.factor=ESS.factor.reg,ESS.reg=ESS.reg,ESS.mvreg.factor=ESS.factor.mvreg,ESS.mvreg=ESS.mvreg)))}
    lres
}

.getESScriteria<-function(LogLik,dof,ESS,ESS.factor,rhon,RSS){
    lParamsSummary<-list()
    lParamsSummary$LogLik<-LogLik
    lParamsSummary$dof<- dof
    lParamsSummary$m2loglik<- -2*LogLik
    lParamsSummary$aic<- -2*LogLik+2*lParamsSummary$dof
    lParamsSummary$aic.c<- lParamsSummary$aic +2*lParamsSummary$dof*(lParamsSummary$dof+1)/(ESS-lParamsSummary$dof-1)
    lParamsSummary$sic<- lParamsSummary$m2loglik+log(ESS)*lParamsSummary$dof
    lParamsSummary$bic<-lParamsSummary$m2loglik+lParamsSummary$dof*log(ESS)
    lParamsSummary$RSS<-RSS
    lParamsSummary$ESS<-ESS
    lParamsSummary$ESS.factor<-ESS.factor
    lParamsSummary$rhon<-rhon
    lParamsSummary
}


.getMVSLproc_size<-function(proc.params,evolmodel){
    if (evolmodel=="slouch"){evolmodel<-"mvslouch"}
    
    kY<-0;kX<-0
    if (evolmodel=="bm"){kY<-0;kX<-nrow(proc.params$vX0)}
    if (evolmodel=="slouch"){kY<-1;kX<-ncol(proc.params$B)}
    if (evolmodel=="ouch"){kY<-ncol(proc.params$A);kX<-0}
    if (evolmodel=="mvslouch"){kY<-ncol(proc.params$A);kX<-ncol(proc.params$B)}
    
    list(procdim=kY+kX)
}

.getMVSLphylCovMatrix<-function(phyltree,proc.params,evolmodel,regimes.times,regimes,regimes.types,Merror=NULL){
    
    ## This is essentially a dummy object that might become useful in the future when B for the mvSLOUCH model is estimated by GLS
    ## now it is used to pass the times of species
    lPrecalculates<-list(vSpecies_times=NA)
    lPrecalculates$vSpecies_times<-phyltree$time.of.nodes[phyltree$tip_species_index] 

    n<-phyltree$Ntips
    params<-list()
    if (evolmodel=="slouch"){evolmodel<-"mvslouch"}
    params$EvolModel<-evolmodel 
    
    
    if (params$EvolModel=="bm"){kY<-0;kX<-nrow(proc.params$vX0)}
    if (params$EvolModel=="slouch"){kY<-1;kX<-ncol(proc.params$B)}
    if (params$EvolModel=="ouch"){kY<-ncol(proc.params$A);kX<-0}
    if (params$EvolModel=="mvslouch"){kY<-ncol(proc.params$A);kX<-ncol(proc.params$B)}
    procdim<-kY+kX
    params$EstimParams<-list()

    modelParams<-vector("list",4); names(modelParams)<-c("regimeTimes","regimes","regimeTypes","Merror")
    modelParams$regimeTimes<-regimes.times; modelParams$regimes<-regimes ; modelParams$regimeTypes<-regimes.types
    if (is.null(Merror) || all(is.na(Merror))) {Merror<-0}
    
    
    modelParams<-c(modelParams,proc.params)
    bCalcA<-TRUE ; lexptcalc<-TRUE; if (evolmodel=="bm"){bCalcA<-FALSE;lexptcalc<-FALSE;}
    modelParams$precalcMatrices<-.decompEigenA.S(modelParams,lPrecalculates,NULL,list(bCalcA=bCalcA,bCovCalc=TRUE,dzetacalc=FALSE,lexptcalc=lexptcalc,kappacalc=FALSE,interceptcalc=FALSE),NULL)
    bReturnPhylTree<-FALSE
    if ((!is.element("mTreeDist",names(phyltree)))||(!is.element("vSpeciesPairs",names(phyltree)))){
	phyltree<-.calculate_Tree_dists_forCovariance(phyltree)
	bReturnPhylTree<-TRUE
    }
    mCovPhyl.Merror<-NA
    if (evolmodel=="bm"){
	## AncestorTimes are not saved as estimation under BM will be done at most once
	mCovPhyl.Merror<-.calc.phyl.cov(evolmodel,modelParams,mAncestorTimes=matrix(phyltree$time.of.nodes[1:phyltree$Ntips],ncol=phyltree$Ntips,nrow=phyltree$Ntips,byrow=FALSE)-phyltree$mTreeDist)
    }else{
	mCovPhyl.Merror<-.calc.phyl.cov(evolmodel,modelParams,mTreeDist=phyltree$mTreeDist,vSpeciesTime=phyltree$time.of.nodes[1:phyltree$Ntips],vSpeciesPairs=phyltree$vSpeciesPairs)
    }

    if (!is.null(Merror)){
	Merror<-.createMeasurementError(Merror,n,procdim)
	if (is.matrix(Merror)){
	    for (i in 1:n){
		mCovPhyl.Merror[(((i-1)*procdim+1):(i*procdim)),(((i-1)*procdim+1):(i*procdim))]<-mCovPhyl.Merror[(((i-1)*procdim+1):(i*procdim)),(((i-1)*procdim+1):(i*procdim))]+Merror
	    }
	}else{
	    for (i in 1:n){
		mCovPhyl.Merror[(((i-1)*procdim+1):(i*procdim)),(((i-1)*procdim+1):(i*procdim))]<-mCovPhyl.Merror[(((i-1)*procdim+1):(i*procdim)),(((i-1)*procdim+1):(i*procdim))]+Merror[[i]]
	    }	
	}
    }
    if (bReturnPhylTree){lres<-list(mCovPhyl.Merror=mCovPhyl.Merror,procdim=procdim,phyltree=phyltree)}
    else{lres<-list(mCovPhyl.Merror=mCovPhyl.Merror,procdim=procdim)}
    lres
}

.calculate_Tree_dists_forCovariance<-function(phyltree){
## called from getESS.R
## mTreeDist : matrix where M[i,j]= distance from node i to the last common ancestor of i and j
## mAncestorTimes: matrix where M[i,j]= time (counting from the root) of the most recent common ancestor of nodes i and j
## mAncestorTimes is not actually calculated here as it will only once (in ESS) be used for a BM so no need to pass it around in memory
## vSpeciesPairs: a technical vector used by sapply so that one does not have to create the vector 1:n^2 [upper.tri all the time]
    n<-phyltree$Ntips
    if (!is.element("vSpeciesPairs",names(phyltree))){
	phyltree$vSpeciesPairs<-matrix(1:n^2,n,n)[upper.tri(matrix(1,n,n),diag=TRUE)]
    }

    if (!is.element("mTreeDist",names(phyltree))){
	phyltree$mTreeDist<-matrix(1:n^2,nrow=n,ncol=n) ## need to calculate all the times of ancestors from mTreeDist,mSpecDist
	phyltree$mTreeDist<-apply(phyltree$mTreeDist,c(1,2),function(ij,path_from_root,time_of_nodes,n){
                i<-(ij-1)%%n+1
                j<-(ij-1)%/%n+1
		vAncijnode<-intersect(path_from_root[[i]]$nodes,path_from_root[[j]]$nodes)
		iAncijnode<-vAncijnode[length(vAncijnode)]
		time_of_nodes[i]-time_of_nodes[iAncijnode]        
	},path_from_root=phyltree$path.from.root,time_of_nodes=phyltree$time.of.nodes,n=n)
	colnames(phyltree$mTreeDist)<-phyltree$tip.label
	rownames(phyltree$mTreeDist)<-phyltree$tip.label	
    }
    phyltree
}