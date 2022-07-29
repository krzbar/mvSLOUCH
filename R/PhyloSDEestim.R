## This file is part of mvSLOUCH

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this package or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .


.check_input_trait_data<-function(mData,n=NULL,vSpeciesLabels=NULL){
## called in PhyloSDE.R, evolmodelest.R, modelparamssummary.R, OUphylregression.R
    if (is.null(mData)){.my_stop("The observed trait data is NULL! It has to be a matrix!",TRUE)}
    if (is.matrix(mData)){
	if (!is.numeric(mData)){.my_stop("The observed trait data is not a numeric matrix!",TRUE)}
    }else{
	if (is.data.frame(mData)){
	    mData<-as.matrix(mData)
	    if (!is.numeric(mData)){.my_stop("The observed trait data is not a matrix but a data.frame! It is also not numeric!",TRUE)}
	    else{.my_stop("The observed trait data is not a matrix but a data.frame!",TRUE)}
	}else{
	    if (is.numeric(mData)){
		if (!is.null(names(mData))){vmData_names<-names(mData)}else{vmData_names<-vSpeciesLabels}
		mData<-matrix(mData,ncol=1)
		colnames(mData)<-"Observed_trait"
		rownames(mData)<-vmData_names
		.my_warning("The observed trait data is not a matrix! Changing it to a single column matrix",TRUE,TRUE)
	    }else{.my_stop("The observed trait data is not a numeric matrix!",TRUE)}
	}
    }   
    ## checking for constant columns found on
    ## https://stackoverflow.com/questions/15068981/removal-of-constant-columns-in-r
    ## but changed == to all.equal since we also do not wnat nearly constant 
    ## due to possible singularities
    if (any(apply(mData, MARGIN = 2, function(x){isTRUE(all.equal(max(x, na.rm = TRUE), min(x, na.rm = TRUE)))}))){
	.my_stop("The observed trait data contains at least one constant or nearly constant between the species trait! Please remove the column(s) as this will make estimation impossible due to actual (or numerical) singularities.",TRUE)
    } 
    if ((!is.null(n))&&(nrow(mData)!=n)){.my_stop(paste("The number of (rows) of the observed trait data does not equal the number of tip species=",n),TRUE)}
    if (is.null(colnames(mData))){
	.my_warning("WARNING: The columns of that observed trait data matrix do not have names. Creating generic ones.",TRUE,FALSE)
	colnames(mData)<-paste("trait_",1:ncol(mData),sep="")
    }
    if (is.null(rownames(mData))){
	.my_warning("WARNING: The rows of that observed trait data matrix are not named according to the phylogeny. We assume that the order of the rows corresponds to the indexing of the tips in the phylogeny object!",TRUE,FALSE)
	if(!is.null(vSpeciesLabels)){
	    rownames(mData)<-vSpeciesLabels
	}
    }else{
	if (!setequal(rownames(mData),vSpeciesLabels)){
	    .my_stop("The row names of the observed trait data matrix (rownames(mData)) do not match the species names in the phylogeny (phyltree$tip.label)!",TRUE)
	}else{
	    mData<-mData[vSpeciesLabels,,drop=FALSE]
	}
    }
    if (any(apply(mData,2,function(vx){all(is.na(vx))}))){
	.my_stop("Cannot have a completely unobserved trait! Remove the completely NA column from the data and perhaps modelling setup!",TRUE)    
    }    
    if (any(apply(mData,1,function(vx){all(is.na(vx))}))){
	.my_stop("Cannot have a completely unobserved species! Remove the completely NA row from the data, remove the species from the tree, correct regimes vector if present and perhaps other modelling setup objects!",TRUE)    
    }    

    mData
}

.PhyloSDEestim<-function(phyltree,mData,kY,regimes=NULL,regimes.times=NULL,root.regime=NULL,predictors=NULL,params=NULL,M.error=NULL,estimate.root.state=FALSE,maxiter=c(10,50,100),cBestim_method="ML"){
## predictors have to be given as column numbers
    
##    cBestim_method<-"ML" ## can be GLS
##    cBestim_method<-"GLS" ## can be GLS
    
    lEstimResults<-NA
    b_printmessage<-FALSE

## =========================================================================================================================
## checking if mData is a matrix    
## PCMBase requires it to be a matrix
## checking if tree is in ape format
    if (!inherits(phyltree,"phylo")){.my_stop("The phylogenetic tree has to be of the phylo format (i.e. ape).",TRUE)}
    mData<-.check_input_trait_data(mData,NULL,phyltree$tip.label)
    n<-nrow(mData)
    kYX<-ncol(mData)
## =========================================================================================================================

## =========================================================================================================================
## checking if M.error is of the correct type and making it appropriate for further use
    M.error<-.createMeasurementError(M.error,n,kYX)
## =========================================================================================================================
    
    if (is.null(params)){params<-list()}
    if (!is.element("EvolModel",names(params))){params$EvolModel<-"mvslouch"}
    if (params$EvolModel=="slouch"){params$EvolModel<-"mvslouch";.my_warning("Call to slouch not yet implemented, using mvslouch model. You might need to run again with corrected input structure.",TRUE,FALSE)}	

## =====================================================================================
## Setup the regimes
    tmpkX<-kYX;if (!is.null(kY)){tmpkX<-kYX-kY}
    regimesList<-.InitialRegimeSetup(phyltree=phyltree,regimes=regimes,regimes.times=regimes.times,mData=mData,kX=tmpkX,kYX=kYX,root.regime=root.regime,M.error=M.error)
    bOKregimes<-regimesList$bOKregimes
    if (!bOKregimes){.my_stop("The regimes, regimes.times and phyltree variables are not consistant. Cannot perform estimation. Please see manual on their format.",TRUE)}
    regimes<-regimesList$regimes
    regimes.times<-regimesList$regimes.times
    regimes.types<-regimesList$regimes.types
    root.regime<-regimesList$root.regime
    regimes.types.orig<-regimesList$regimes.types.orig ## regimes names are in alphabetical order
    phyltree<-regimesList$phyltree
    pcmbase_model_box<-regimesList$pcmbase_model_box ## this is the model object in which parameters for all regimes sit in for PCMBase
    mData<-.check_input_trait_data(mData,n=phyltree$Ntips,vSpeciesLabels=phyltree$tip.label)
## =====================================================================================    
    if(estimate.root.state){if ((is.null(regimes))||(length(unique(regimes))==1)){estimate.root.state<-FALSE}}

    if (params$EvolModel=="bm"){kY<-0}
    if (params$EvolModel=="slouch"){kY<-1}
    if ((params$EvolModel=="ouch")&&(is.null(kY)||is.na(kY)||(kY!=kYX))){
        .my_message(paste("Wrong or missing value of kY, ",kY," supplied, changing it to : ",kYX,sep=""),b_printmessage);
        kY<-kYX
    }
    kX<-kYX-kY

	
    if (!is.element("method",names(params))){params$method<-"glsgc"} 
    if (params$EvolModel=="bm"){params$method<-"maxlik"}  
	
    if (!is.element("bShouldPrint",names(params))){params$bShouldPrint<-FALSE}
    if (!is.element("tol",names(params))){params$tol<-.set.tol(params$method)}
    if (!is.element("maxIter",names(params))){params$maxIter<-.set.maxiter(params$method,maxiter)}
    if (!is.element("maxTries",names(params))){params$maxTries<-10}
    if (!is.element("minLogLik",names(params))){params$minLogLik<- -Inf}
    if (!is.element("pcmbase_model_box",names(params))){params$pcmbase_model_box<-pcmbase_model_box}    

    params$EstimParams<-.set.estimparams(params,kY,kX,length(regimes.types),estimate.root.state,phyltree$tree_height,cBestim_method)
	
## ============= post-hoc corrections to EstimaParams =======================================    
    if (params$EstimParams$designToEstim$y0AncState){
        if (!is.null(root.regime)){params$EstimParams$designToEstim$y0Regime<-which(regimes.types.orig==root.regime)}
        else{
    	    if (is.element("y0Regime",names(params$EstimParams$designToEstim))){params$EstimParams$designToEstim$y0Regime<-which(regimes.types.orig==params$EstimParams$designToEstim$y0Regime)}	
	    else{params$EstimParams$designToEstim$y0Regime<-regimes[[1]][1]} ## need to choose something anyway ...	
	}
    }

    if (!is.null(predictors)){
        predictors<-intersect(1:kYX,predictors)## if user provided column out of range
        params$EstimParams$predictors<-predictors
        predictors<-colnames(mData)[predictors]
    }	
    if (!is.element("TerminalLabels",names(params$EstimParams))){params$EstimParams$TerminalLabels<-phyltree$tip.label}
    
    ## this is so for non mvslouch models B is not attempted to be optimized over
    if (params$EvolModel!="mvslouch"){params$EstimParams$Fixed$B<-NA}
    
    if (!is.null(M.error)){
        params$EstimParams$M_error<-M.error
    }

    params$EstimParams$RegimeTypes<-regimes.types ## required in modelparamtransform if mPsi would be optimized over in ML
## ==============================================================================================	
	
    if (params$bShouldPrint){.my_message("Run time, start: \n",TRUE);.my_message(paste(Sys.time(),collapse=" "),TRUE);}

    lEstimResults<-.maxlik.estim(mData=mData,EvolModel=params$EvolModel,method=params$method,regimes=regimes,regimeTypes=regimes.types,EstimationParams=params$EstimParams,PhylTree=phyltree,regimeTimes=regimes.times,bShouldPrint=params$bShouldPrint,tol=params$tol,maxIter=params$maxIter,maxTries=params$maxTries,minLogLik=params$minLogLik,regimes.types.orig=regimes.types.orig)
    if (params$bShouldPrint){.my_message("Running time, end \n",TRUE);.my_message(paste(Sys.time(),collapse=", "),TRUE)}
    lEstimResults<-.correct.names(lEstimResults,regimes.types.orig,if(params$EvolModel!="bm"){colnames(mData)[1:kY]}else{NULL},if ((params$EvolModel=="mvslouch")||(params$EvolModel=="bm")||(params$EvolModel=="slouch")){colnames(mData)[(kY+1):(kY+kX)]}else{NULL},predictors,params$EvolModel)
    lEstimResults
}

.correct.names<-function(listobj,vregime.names,vY.names,vX.names,predictors=NULL,EvolModel="mvslouch"){
    if (is.element("num_extracted",names(listobj))){listobj$num_extracted<-NULL}
    if (is.element("pcmbase_model_box",names(listobj))){listobj$pcmbase_model_box<-NULL}
    if (is.element("M_error",names(listobj))){listobj$M_error<-NULL}
    if ((EvolModel=="ouch")||(EvolModel=="bm")){
	if (is.element("BrownResult",names(listobj))){listobj$BrownResult<-NULL}
	if (is.element("B",names(listobj))){listobj$B<-NULL}
        if (is.element("B.confidence.interval",names(listobj))){listobj$B.confidence.interval<-NULL}
        if (is.element("B.regression.confidence.interval",names(listobj))){listobj$B.regression.confidence.interval<-NULL}
	if (is.element("Syx",names(listobj))){listobj$Syx<-NULL}
	if (is.element("Sxy",names(listobj))){listobj$Sxy<-NULL}
	if (is.element("Syx.confidence.interval",names(listobj))){listobj$Syx.confidence.interval<-NULL}
	if (is.element("Sxy.confidence.interval",names(listobj))){listobj$Sxy.confidence.interval<-NULL}
	if (is.element("optimal.regression",names(listobj))){listobj$optimal.regression<-NULL}
	if (is.element("A1B",names(listobj))){listobj$A1B<-NULL}
	if (is.element("conditional.cov.matrix",names(listobj))){listobj$conditional.cov.matrix<-NULL}
	if (is.element("conditional.corr.matrix",names(listobj))){listobj$conditional.corr.matrix<-NULL}
	if (is.element("evolutionary.regression",names(listobj))){listobj$evolutionary.regression<-NULL}
	if ((is.element("cov.with.optima",names(listobj)))&&(!is.character(listobj$cov.with.optima))){if (length(vY.names)==1){listobj$cov.with.optima<-matrix(listobj$cov.with.optima,1,1)};colnames(listobj$cov.with.optima)<-vY.names;rownames(listobj$cov.with.optima)<-vY.names}
	if ((is.element("corr.with.optima",names(listobj)))&&(!is.character(listobj$corr.with.optima))){if (length(vY.names)==1){listobj$corr.with.optima<-matrix(listobj$corr.with.optima,1,1)};colnames(listobj$corr.with.optima)<-vY.names;rownames(listobj$corr.with.optima)<-vY.names}
	if (EvolModel=="bm"){
    	    if (is.element("A",names(listobj))){listobj$A<-NULL}
    	    if (is.element("A.confidence.interval",names(listobj))){listobj$A.confidence.interval<-NULL}
	    if (is.element("mPsi",names(listobj))){listobj$mPsi<-NULL}	
	    if (is.element("mPsi.rotated",names(listobj))){listobj$mPsi.rotated<-NULL}		
    	    if (is.element("expmtA",names(listobj))){listobj$A<-NULL}
	    if (is.element("mPsi.confidence.interval",names(listobj))){listobj$mPsi.confidence.interval<-NULL}	
	    if (is.element("mPsi.regression.confidence.interval",names(listobj))){listobj$mPsi.regression.confidence.interval<-NULL}	

    	    if (is.element("expmtA.confidence.interval",names(listobj))){listobj$A.confidence.interval<-NULL}
	}
	if (EvolModel=="ouch"){
	    if (is.element("vX0",names(listobj))){listobj$vX0<-NULL}
	    if (is.element("vX0.confidence.interval",names(listobj))){listobj$vX0.confidence.interval<-NULL}
	    if (is.element("vX0.regression.confidence.interval",names(listobj))){listobj$vX0.regression.confidence.interval<-NULL}
	    if (is.element("Sxx",names(listobj))){listobj$Sxx<-NULL}
	    if (is.element("Sxx.confidence.interval",names(listobj))){listobj$Sxx.confidence.interval<-NULL}
	    if (is.element("optima.correlations",names(listobj))){if (length(vY.names)==1){listobj$optima.correlations<-matrix(listobj$optima.correlations,1,1)};colnames(listobj$optima.correlations)<-vY.names;rownames(listobj$optima.correlations)<-vY.names}    
	}
	if (EvolModel=="bm"){
	    if (is.element("vY0",names(listobj))){listobj$vY0<-NULL}
	    if (is.element("vY0.confidence.interval",names(listobj))){listobj$vY0.confidence.interval<-NULL}
	    if (is.element("vY0.regression.confidence.interval",names(listobj))){listobj$vY0.regression.confidence.interval<-NULL}
	    if (is.element("Syy",names(listobj))){listobj$Syy<-NULL}
	    if (is.element("Syy.confidence.interval",names(listobj))){listobj$Syy.confidence.interval<-NULL}
	}
    }
    ## here we remove mPsi0 from output
    if (is.element("mPsi0",names(listobj))){listobj$mPsi0<-NULL}	
    if (is.element("mPsi0.rotated",names(listobj))){listobj$mPsi0.rotated<-NULL}	
    if (is.element("mPsi0.confidence.interval",names(listobj))){listobj$mPsi0.confidence.interval<-NULL}	
    if (is.element("mPsi0.regression.confidence.interval",names(listobj))){listobj$mPsi0.regression.confidence.interval<-NULL}	
    if (is.element("regression_covariance_matrix",names(listobj))){
	v_mPsi0colsrows<-grep("mPsi0+",colnames(listobj$regression_covariance_matrix))
	if (length(v_mPsi0colsrows)>0){
	    listobj$regression_covariance_matrix<-listobj$regression_covariance_matrix[-v_mPsi0colsrows,-v_mPsi0colsrows,drop=FALSE]	    
	}
    }
    # ============================
    
    if (is.element("A",names(listobj))){if (length(vY.names)==1){listobj$A<-matrix(listobj$A,1,1)};colnames(listobj$A)<-vY.names;rownames(listobj$A)<-vY.names}    
    if (is.element("A.confidence.interval",names(listobj))){colnames(listobj$A.confidence.interval$Lower.end)<-vY.names;rownames(listobj$A.confidence.interval$Lower.end)<-vY.names;colnames(listobj$A.confidence.interval$Estimated.Point)<-vY.names;rownames(listobj$A.confidence.interval$Estimated.Point)<-vY.names;colnames(listobj$A.confidence.interval$Upper.end)<-vY.names;rownames(listobj$A.confidence.interval$Upper.end)<-vY.names}    
    if (is.element("B",names(listobj))){if ((length(vY.names)==1)||(length(vX.names)==1)){listobj$B<-matrix(listobj$B,nrow=length(vY.names),ncol=length(vX.names))};colnames(listobj$B)<-vX.names;rownames(listobj$B)<-vY.names}
    if (is.element("B.confidence.interval",names(listobj))){colnames(listobj$B.confidence.interval$Lower.end)<-vX.names;rownames(listobj$B.confidence.interval$Lower.end)<-vY.names;colnames(listobj$B.confidence.interval$Estimated.Point)<-vX.names;rownames(listobj$B.confidence.interval$Estimated.Point)<-vY.names;colnames(listobj$B.confidence.interval$Upper.end)<-vX.names;rownames(listobj$B.confidence.interval$Upper.end)<-vY.names}    
    if (is.element("B.regression.confidence.interval",names(listobj))){
	if (length(vX.names)==1){rownames(listobj$B.regression.confidence.interval)<-vY.names}    
	else{colnames(listobj$B.regression.confidence.interval$Lower.end)<-vX.names;rownames(listobj$B.regression.confidence.interval$Lower.end)<-vY.names;colnames(listobj$B.regression.confidence.interval$Estimated.Point)<-vX.names;rownames(listobj$B.regression.confidence.interval$Estimated.Point)<-vY.names;colnames(listobj$B.regression.confidence.interval$Upper.end)<-vX.names;rownames(listobj$B.regression.confidence.interval$Upper.end)<-vY.names}
    }	
    if (is.element("mPsi",names(listobj))){if ((length(vY.names)==1)||(length(vregime.names)==1)){listobj$mPsi<-matrix(listobj$mPsi,nrow=length(vY.names),ncol=length(vregime.names))};colnames(listobj$mPsi)<-vregime.names;rownames(listobj$mPsi)<-vY.names}
    if (is.element("mPsi.confidence.interval",names(listobj))){colnames(listobj$mPsi.confidence.interval$Lower.end)<-vregime.names;rownames(listobj$mPsi.confidence.interval$Lower.end)<-vY.names;colnames(listobj$mPsi.confidence.interval$Estimated.Point)<-vregime.names;rownames(listobj$mPsi.confidence.interval$Estimated.Point)<-vY.names;colnames(listobj$mPsi.confidence.interval$Upper.end)<-vregime.names;rownames(listobj$mPsi.confidence.interval$Upper.end)<-vY.names}    
    if (is.element("mPsi.regression.confidence.interval",names(listobj))){
	if (length(vregime.names)==1){rownames(listobj$mPsi.regression.confidence.interval)<-vY.names}    
	else{colnames(listobj$mPsi.regression.confidence.interval$Lower.end)<-vregime.names;rownames(listobj$mPsi.regression.confidence.interval$Lower.end)<-vY.names;colnames(listobj$mPsi.regression.confidence.interval$Estimated.Point)<-vregime.names;rownames(listobj$mPsi.regression.confidence.interval$Estimated.Point)<-vY.names;colnames(listobj$mPsi.regression.confidence.interval$Upper.end)<-vregime.names;rownames(listobj$mPsi.regression.confidence.interval$Upper.end)<-vY.names}
    }
    #if (is.element("mPsi0",names(listobj))){if (length(vY.names)==1){listobj$mPsi0<-matrix(listobj$mPsi0,nrow=1,ncol=1)};rownames(listobj$mPsi0)<-vY.names}
    #if (is.element("mPsi0.confidence.interval",names(listobj))){rownames(listobj$mPsi0.confidence.interval$Lower.end)<-vY.names;rownames(listobj$mPsi0.confidence.interval$Estimated.Point)<-vY.names;rownames(listobj$mPsi0.confidence.interval$Upper.end)<-vY.names}    
    #if (is.element("mPsi0.regression.confidence.interval",names(listobj))){rownames(listobj$mPsi0.regression.confidence.interval)<-vY.names}   
    if (is.element("vY0",names(listobj))){if (length(vY.names)==1){listobj$vY0<-matrix(listobj$vY0,1,1)};rownames(listobj$vY0)<-vY.names}
    if (is.element("vY0.confidence.interval",names(listobj))){rownames(listobj$vY0.confidence.interval$Lower.end)<-vY.names;rownames(listobj$vY0.confidence.interval$Estimated.Point)<-vY.names;rownames(listobj$vY0.confidence.interval$Upper.end)<-vY.names}
    if (is.element("vY0.regression.confidence.interval",names(listobj))){rownames(listobj$vY0.regression.confidence.interval)<-vY.names}
    if (is.element("vX0",names(listobj))){if (length(vX.names)==1){listobj$vX0<-matrix(listobj$vX0,1,1)};rownames(listobj$vX0)<-vX.names}
    if (is.element("vX0.confidence.interval",names(listobj))){rownames(listobj$vX0.confidence.interval$Lower.end)<-vX.names;rownames(listobj$vX0.confidence.interval$Estimated.Point)<-vX.names;rownames(listobj$vX0.confidence.interval$Upper.end)<-vX.names}
    if (is.element("vX0.regression.confidence.interval",names(listobj))){rownames(listobj$vX0.regression.confidence.interval)<-vX.names}
    if (is.element("BX0.regression.confidence.interval",names(listobj))){rownames(listobj$BX0.regression.confidence.interval)<-vX.names}
    if (is.element("Syy",names(listobj))){if (length(vY.names)==1){listobj$Syy<-matrix(listobj$Syy,1,1)};colnames(listobj$Syy)<-vY.names;rownames(listobj$Syy)<-vY.names}
    if (is.element("Syy.confidence.interval",names(listobj))){colnames(listobj$Syy.confidence.interval$Lower.end)<-vY.names;rownames(listobj$Syy.confidence.interval$Lower.end)<-vY.names;colnames(listobj$Syy.confidence.interval$Estimated.Point)<-vY.names;rownames(listobj$Syy.confidence.interval$Estimated.Point)<-vY.names;colnames(listobj$Syy.confidence.interval$Upper.end)<-vY.names;rownames(listobj$Syy.confidence.interval$Upper.end)<-vY.names}    
    if (is.element("Syx",names(listobj))){if ((length(vY.names)==1)||(length(vX.names)==1)){listobj$Syx<-matrix(listobj$Syx,nrow=length(vY.names),ncol=length(vX.names))};colnames(listobj$Syx)<-vX.names;rownames(listobj$Syx)<-vY.names}
    if (is.element("Syx.confidence.interval",names(listobj))){colnames(listobj$Syx.confidence.interval$Lower.end)<-vX.names;rownames(listobj$Syx.confidence.interval$Lower.end)<-vY.names;colnames(listobj$Syx.confidence.interval$Estimated.Point)<-vX.names;rownames(listobj$Syx.confidence.interval$Estimated.Point)<-vY.names;colnames(listobj$Syx.confidence.interval$Upper.end)<-vX.names;rownames(listobj$Syx.confidence.interval$Upper.end)<-vY.names}    
    if (is.element("Sxy",names(listobj))){if ((length(vY.names)==1)||(length(vX.names)==1)){listobj$Sxy<-matrix(listobj$Sxy,ncol=length(vY.names),nrow=length(vX.names))};colnames(listobj$Sxy)<-vY.names;rownames(listobj$Sxy)<-vX.names}
    if (is.element("Sxy.confidence.interval",names(listobj))){colnames(listobj$Sxy.confidence.interval$Lower.end)<-vY.names;rownames(listobj$Sxy.confidence.interval$Lower.end)<-vX.names;colnames(listobj$Sxy.confidence.interval$Estimated.Point)<-vY.names;rownames(listobj$Sxy.confidence.interval$Estimated.Point)<-vX.names;colnames(listobj$Sxy.confidence.interval$Upper.end)<-vY.names;rownames(listobj$Sxy.confidence.interval$Upper.end)<-vX.names}    
    if (is.element("Sxx",names(listobj))){if (length(vX.names)==1){listobj$Sxx<-matrix(listobj$Sxx,1,1)};colnames(listobj$Sxx)<-vX.names;rownames(listobj$Sxx)<-vX.names}
    if (is.element("Sxx.confidence.interval",names(listobj))){colnames(listobj$Sxx.confidence.interval$Lower.end)<-vX.names;rownames(listobj$Sxx.confidence.interval$Lower.end)<-vX.names;colnames(listobj$Sxx.confidence.interval$Estimated.Point)<-vX.names;rownames(listobj$Sxx.confidence.interval$Estimated.Point)<-vX.names;colnames(listobj$Sxx.confidence.interval$Upper.end)<-vX.names;rownames(listobj$Sxx.confidence.interval$Upper.end)<-vX.names}    
    if (is.element("expmtA",names(listobj))){if (length(vY.names)==1){listobj$expmtA<-matrix(listobj$expmtA,1,1)};colnames(listobj$expmtA)<-vY.names;rownames(listobj$expmtA)<-vY.names}
    if (is.element("mPsi.rotated",names(listobj))){if ((length(vY.names)==1)||(length(vregime.names)==1)){listobj$mPsi.rotated<-matrix(listobj$mPsi.rotated,nrow=length(vY.names),ncol=length(vregime.names))};colnames(listobj$mPsi.rotated)<-vregime.names;rownames(listobj$mPsi.rotated)<-vY.names}
    #if (is.element("mPsi0.rotated",names(listobj))){if (length(vY.names)==1)   {listobj$mPsi0.rotated<-matrix(listobj$mPsi0.rotated,nrow=1,ncol=1)};rownames(listobj$mPsi0.rotated)<-vY.names}
    if (is.element("optimal.regression",names(listobj))){if ((length(vY.names)==1)||(length(vX.names)==1)){listobj$optimal.regression<-matrix(listobj$optimal.regression,nrow=length(vY.names),ncol=length(vX.names))};rownames(listobj$optimal.regression)<-vY.names;colnames(listobj$optimal.regression)<-vX.names}
    if (is.element("A1B",names(listobj))){if ((length(vY.names)==1)||(length(vX.names)==1)){listobj$A1B<-matrix(listobj$A1B,nrow=length(vY.names),ncol=length(vX.names))};rownames(listobj$A1B)<-vY.names;colnames(listobj$A1B)<-vX.names}
    if (is.element("cov.matrix",names(listobj))){colnames(listobj$cov.matrix)<-c(vY.names,vX.names);rownames(listobj$cov.matrix)<-c(vY.names,vX.names)}
    if (is.element("corr.matrix",names(listobj))){colnames(listobj$corr.matrix)<-c(vY.names,vX.names);rownames(listobj$corr.matrix)<-c(vY.names,vX.names)}
    if (is.element("conditional.cov.matrix",names(listobj))){if (length(vY.names)==1){listobj$conditional.cov.matrix<-matrix(listobj$conditional.cov.matrix,1,1)};colnames(listobj$conditional.cov.matrix)<-vY.names;rownames(listobj$conditional.cov.matrix)<-vY.names}
    if (is.element("conditional.corr.matrix",names(listobj))){if (length(vY.names)==1){listobj$conditional.corr.matrix<-matrix(listobj$conditional.corr.matrix,1,1)};colnames(listobj$conditional.corr.matrix)<-vY.names;rownames(listobj$conditional.corr.matrix)<-vY.names}
    if (is.element("evolutionary.regression",names(listobj))){
	if ((length(vY.names)==1)||(length(vX.names)==1)){listobj$evolutionary.regression<-matrix(listobj$evolutionary.regression,ncol=length(vX.names),nrow=length(vY.names))}
	if (is.null(predictors)){colnames(listobj$evolutionary.regression)<-vX.names;rownames(listobj$evolutionary.regression)<-vY.names}
	else{
	    if (length(setdiff(predictors,vX.names))>0){
		listobj$regression.Y.on.X<-listobj$evolutionary.regression
		colnames(listobj$regression.Y.on.X)<-vX.names;rownames(listobj$regression.Y.on.X)<-vY.names
		listobj$evolutionary.regression<-NULL
		tryCatch({
		    listobj$evolutionary.regression<-listobj$cov.matrix[setdiff(c(vY.names,vX.names),predictors),predictors]%*%solve(listobj$cov.matrix[predictors,predictors])		    		
		    if ((length(predictors)==1)||(length(setdiff(c(vY.names,vX.names),predictors))==1)){listobj$evolutionary.regression<-matrix(listobj$evolutionary.regression,ncol=length(predictors),nrow=length(setdiff(c(vY.names,vX.names),predictors)))}
		    colnames(listobj$evolutionary.regression)<-predictors;rownames(listobj$evolutionary.regression)<-setdiff(c(vY.names,vX.names),predictors)
		},error=function(e){.my_message(paste("Error in evolutionary regression calculation",e),TRUE);.my_message("\n",TRUE)})		
		if (is.element("conditional.cov.matrix",names(listobj))){
		    listobj$conditional.cov.Y.on.X<-listobj$conditional.cov.matrix
		    listobj$conditional.corr.Y.on.X<-listobj$conditional.corr.matrix
		}
	       listobj$conditional.cov.matrix<-NULL
	       listobj$conditional.corr.matrix<-NULL
		tryCatch({
		    listobj$conditional.cov.matrix<-listobj$cov.matrix[setdiff(c(vY.names,vX.names),predictors),setdiff(c(vY.names,vX.names),predictors)]-listobj$cov.matrix[setdiff(c(vY.names,vX.names),predictors),predictors]%*%solve(listobj$cov.matrix[predictors,predictors])%*%listobj$cov.matrix[predictors,setdiff(c(vY.names,vX.names),predictors)]
		    if (length(setdiff(c(vY.names,vX.names),predictors))==1){listobj$conditional.cov.matrix<-matrix(listobj$conditional.cov.matrix,ncol=1,nrow=1)}
		    listobj$conditional.corr.matrix<-.my_cov2cor(listobj$conditional.cov.matrix)
		    colnames(listobj$conditional.cov.matrix)<-setdiff(c(vY.names,vX.names),predictors);rownames(listobj$conditional.cov.matrix)<-setdiff(c(vY.names,vX.names),predictors)
		    colnames(listobj$conditional.corr.matrix)<-setdiff(c(vY.names,vX.names),predictors);rownames(listobj$conditional.corr.matrix)<-setdiff(c(vY.names,vX.names),predictors)
		},error=function(e){.my_message(paste("Error in conditional covariance calculation",e),TRUE);.my_message("\n",TRUE)})		
	    }	    
	    else{colnames(listobj$evolutionary.regression)<-vX.names;rownames(listobj$evolutionary.regression)<-vY.names}	
	}
    }else{
	if ((is.element("cov.matrix",names(listobj)))&&(!is.null(predictors))){
	    listobj$evolutionary.regression<-NULL
	    tryCatch({	
    		listobj$evolutionary.regression<-listobj$cov.matrix[setdiff(c(vY.names,vX.names),predictors),predictors]%*%solve(listobj$cov.matrix[predictors,predictors])
		    if ((length(predictors)==1)||(length(setdiff(c(vY.names,vX.names),predictors))==1)){listobj$evolutionary.regression<-matrix(listobj$evolutionary.regression,ncol=length(predictors),nrow=length(setdiff(c(vY.names,vX.names),predictors)))}
		    colnames(listobj$evolutionary.regression)<-predictors;rownames(listobj$evolutionary.regression)<-setdiff(c(vY.names,vX.names),predictors)
	    },error=function(e){.my_message(paste("Error in evolutionary regression calculation",e),TRUE);.my_message("\n",TRUE)})		
	    if (is.element("conditional.cov.matrix",names(listobj))){
		listobj$conditional.cov.Y.on.X<-listobj$conditional.cov.matrix
		listobj$conditional.corr.Y.on.X<-listobj$conditional.corr.matrix
	    }
	    listobj$conditional.cov.matrix<-NULL
	    listobj$conditional.corr.matrix<-NULL
	    tryCatch({
		listobj$conditional.cov.matrix<-listobj$cov.matrix[setdiff(c(vY.names,vX.names),predictors),setdiff(c(vY.names,vX.names),predictors)]-listobj$cov.matrix[setdiff(c(vY.names,vX.names),predictors),predictors]%*%solve(listobj$cov.matrix[predictors,predictors])%*%listobj$cov.matrix[predictors,setdiff(c(vY.names,vX.names),predictors)]
		if (length(setdiff(c(vY.names,vX.names),predictors))==1){listobj$conditional.cov.matrix<-matrix(listobj$conditional.cov.matrix,ncol=1,nrow=1)}
		listobj$conditional.corr.matrix<-.my_cov2cor(listobj$conditional.cov.matrix)
		
		colnames(listobj$conditional.cov.matrix)<-setdiff(c(vY.names,vX.names),predictors);rownames(listobj$conditional.cov.matrix)<-setdiff(c(vY.names,vX.names),predictors)
		colnames(listobj$conditional.corr.matrix)<-setdiff(c(vY.names,vX.names),predictors);rownames(listobj$conditional.corr.matrix)<-setdiff(c(vY.names,vX.names),predictors)
	    },error=function(e){.my_message(paste("Error in conditional covariance calculation",e),TRUE);.my_message("\n",TRUE)})		
	}
    }
    if (is.element("trait.regression",names(listobj))){
        listobj$trait.regression<-sapply(1:length(vY.names),function(i,vY.names,vX.names,trait.reg){
	    		    trait.reg<-matrix(trait.reg[[i]],nrow=1);colnames(trait.reg)<-c(vY.names[-i],vX.names);rownames(trait.reg)<-vY.names[i]
			    trait.reg
	},vY.names=vY.names,vX.names=vX.names,trait.reg=listobj$trait.regression,simplify=FALSE)
        if (length(setdiff(predictors,vX.names))>0){
	    listobj$Y.regression<-listobj$trait.regression
	    if (is.element("cov.matrix",names(listobj))){
		vtraits<-setdiff(c(vY.names,vX.names),predictors)
		vtraitsIndex<-sort(sapply(vtraits,function(trait,vNames){which(vNames==trait)},vNames=c(vY.names,vX.names),simplify=TRUE))
		listobj$trait.regression<-NULL
		tryCatch({		
		    listobj$trait.regression<-sapply(1:length(vtraits),function(i,vtraitsIndex,mCov){
				mCov[vtraitsIndex[i],-vtraitsIndex[i],drop=FALSE]%*%solve(mCov[-vtraitsIndex[i],-vtraitsIndex[i],drop=FALSE])    
		    },vtraitsIndex=vtraitsIndex,mCov=listobj$cov.matrix,simplify=FALSE)
		    listobj$trait.regression<-sapply(1:length(vtraits),function(i,vY.names,vX.names,vtraitsIndex,trait.reg){
			    trait.reg<-matrix(trait.reg[[i]],nrow=1);colnames(trait.reg)<-c(vY.names,vX.names)[-vtraitsIndex[i]];rownames(trait.reg)<-c(vY.names,vX.names)[vtraitsIndex[i]]
			    trait.reg
		    },vY.names=vY.names,vX.names=vX.names,vtraitsIndex=vtraitsIndex,trait.reg=listobj$trait.regression,simplify=FALSE)
		},error=function(e){.my_message(paste("Error in trait regression calculation",e),TRUE);.my_message("\n",TRUE)})		
	    }else{listobj$trait.regression<-NULL}
	}
    }
    if ((is.element("cov.with.optima",names(listobj)))&&(!is.character(listobj$cov.with.optima))){if (length(vY.names)==1){listobj$cov.with.optima<-matrix(listobj$cov.with.optima,1,1)};colnames(listobj$cov.with.optima)<-vY.names;rownames(listobj$cov.with.optima)<-vY.names}
    if ((is.element("corr.with.optima",names(listobj)))&&(!is.character(listobj$corr.with.optima))){if (length(vY.names)==1){listobj$corr.with.optima<-matrix(listobj$corr.with.optima,1,1)};colnames(listobj$corr.with.optima)<-vY.names;rownames(listobj$corr.with.optima)<-vY.names}
    if ((is.element("optima.cov.matrix",names(listobj)))&&(!is.character(listobj$optima.cov.matrix))){if (length(vY.names)==1){listobj$optima.cov.matrix<-matrix(listobj$optima.cov.matrix,1,1)};colnames(listobj$optima.cov.matrix)<-vY.names;rownames(listobj$optima.cov.matrix)<-vY.names}
    if ((is.element("optima.corr.matrix",names(listobj)))&&(!is.character(listobj$optima.corr.matrix))){if (length(vY.names)==1){listobj$optima.corr.matrix<-matrix(listobj$optima.corr.matrix,1,1)};colnames(listobj$optima.corr.matrix)<-vY.names;rownames(listobj$optima.corr.matrix)<-vY.names}
    if ((is.element("stationary.cov.matrix",names(listobj)))&&(!is.character(listobj$stationary.cov.matrix))){
	if (length(vY.names)==1){listobj$stationary.cov.matrix<-matrix(listobj$stationary.cov.matrix,1,1)};colnames(listobj$stationary.cov.matrix)<-vY.names;rownames(listobj$stationary.cov.matrix)<-vY.names
	if (length(vY.names)==1){listobj$stationary.corr.matrix<-matrix(listobj$stationary.corr.matrix,1,1)};colnames(listobj$stationary.corr.matrix)<-vY.names;rownames(listobj$stationary.corr.matrix)<-vY.names
	if ((length(vX.names)==0)&&(!is.null(predictors))){
	    vResps<-setdiff(vY.names,predictors)
	    tryCatch({
		if (!is.element("limiting.regression",names(listobj))){
		    listobj$limiting.regression<-listobj$stationary.cov.matrix[vResps,predictors]%*%solve(listobj$stationary.cov.matrix[predictors,predictors])
		    if ((length(vResps)==1)||(length(predictors)==1)){listobj$limiting.regression<-matrix(listobj$limiting.regression,nrow=length(vResps),ncol=length(predictors))};rownames(listobj$limiting.regression)<-vResps;colnames(listobj$limiting.regression)<-predictors
		}
		if (!is.element("limiting.relationship.cov.matrix",names(listobj))){
		    listobj$limiting.relationship.cov.matrix<-listobj$limiting.regression%*%listobj$stationary.cov.matrix[predictors,predictors]%*%t(listobj$limiting.regression)
		    if (length(vResps)==1){listobj$limiting.relationship.cov.matrix<-matrix(listobj$limiting.relationship.cov.matrix,1,1)};colnames(listobj$limiting.relationship.cov.matrix)<-vResps;rownames(listobj$limiting.relationship.cov.matrix)<-vResps    
		    listobj$limiting.relationship.corr.matrix<-.my_cov2cor(listobj$limiting.relationship.cov.matrix)		    

		}
		if (!is.element("cov.with.limit",names(listobj))){
		    listobj$cov.with.limit<-listobj$stationary.cov.matrix[vResps,predictors]%*%t(listobj$limiting.regression)
		    mcov.curr<-listobj$stationary.cov.matrix[vResps,vResps,drop=FALSE]
		    mcov.stat<-listobj$limiting.relationship.cov.matrix[vResps,vResps,drop=FALSE]
		    listobj$corr.with.limit<-apply(matrix(0:((length(vResps))^2-1),length(vResps),length(vResps),byrow=TRUE),c(1,2),function(ij,kY,mcov.curr,mcov.stat,mcov.with){i<-ij%/%kY+1;j<-ij%%kY+1;mcov.with[i,j]/(sqrt(mcov.curr[i,i]*mcov.stat[j,j]))},kY=length(vResps),mcov.curr=mcov.curr,mcov.stat=mcov.stat,mcov.with=listobj$cov.with.limit)
		    if (length(vResps)==1){listobj$cov.with.limit<-matrix(listobj$cov.with.limit,1,1)};colnames(listobj$cov.with.limit)<-vResps;rownames(listobj$cov.with.limit)<-vResps    
		}
	    },error=function(e){.my_message(paste("Error in limiting regression and covariance calculation",e),TRUE);.my_message("\n",TRUE)})		
	}
	tryCatch({    
	    if ((!is.element("limiting.trait.regression",names(listobj)))&&(length(vY.names)>1)){ ## no point in doing the regressions if there is only a single trait
		listobj$limiting.trait.regression<-sapply(1:length(vY.names),function(i,mCov,vY.names){
			    limiting.trait.reg<-mCov[i,-i,drop=FALSE]%*%solve(mCov[-i,-i,drop=FALSE])    
			    limiting.trait.reg<-matrix(limiting.trait.reg,nrow=1);colnames(limiting.trait.reg)<-vY.names[-i];rownames(limiting.trait.reg)<-vY.names[i]
			    limiting.trait.reg	
		},mCov=listobj$stationary.cov.matrix,vY.names=vY.names,simplify=FALSE)
	    }
	},error=function(e){.my_message(paste("Error in limiting trait regression calculation",e),TRUE);.my_message("\n",TRUE)})
    }   
    if (is.element("StS",names(listobj))){colnames(listobj$StS)<-c(vY.names,vX.names);rownames(listobj$StS)<-c(vY.names,vX.names)}
    if (is.element("lower.summary",names(listobj))){
        if (is.element("LogLik",names(listobj$lower.summary))){listobj$lower.summary$LogLik<-NULL}
	if (is.element("dof",names(listobj$lower.summary))){listobj$lower.summary$dof<-NULL}
	if (is.element("m2loglik",names(listobj$lower.summary))){listobj$lower.summary$m2loglik<-NULL}
	if (is.element("aic",names(listobj$lower.summary))){listobj$lower.summary$aic<-NULL}
	if (is.element("aic.c",names(listobj$lower.summary))){listobj$lower.summary$aic.c<-NULL}
	if (is.element("sic",names(listobj$lower.summary))){listobj$lower.summary$sic<-NULL}
	if (is.element("bic",names(listobj$lower.summary))){listobj$lower.summary$bic<-NULL}
	if (is.element("RSS",names(listobj$lower.summary))){listobj$lower.summary$RSS<-NULL}
	if (is.element("R2",names(listobj$lower.summary))){listobj$lower.summary$R2<-NULL}
	if (is.element("R2_phylaverage",names(listobj$lower.summary))){listobj$lower.summary$R2_phylaverage<-NULL}
	if (is.element("RSS_comment",names(listobj$lower.summary))){listobj$lower.summary$RSS_comment<-NULL}
	if (is.element("R2_comment",names(listobj$lower.summary))){listobj$lower.summary$R2_comment<-NULL}
	if (is.element("R2_phylaverage_comment",names(listobj$lower.summary))){listobj$lower.summary$R2_phylaverage_comment<-NULL}
	if (is.element("RSS_non_phylogenetic",names(listobj$lower.summary))){listobj$lower.summary$RSS_non_phylogenetic<-NULL}
	if (is.element("R2_non_phylogenetic",names(listobj$lower.summary))){listobj$lower.summary$R2_non_phylogenetic<-NULL}
	if (is.element("RSS_non_phylogenetic_comment",names(listobj$lower.summary))){listobj$lower.summary$RSS_non_phylogenetic_comment<-NULL}
	if (is.element("R2_non_phylogenetic_comment",names(listobj$lower.summary))){listobj$lower.summary$R2_non_phylogenetic_comment<-NULL}
	if (is.element("RSS_non_phylogenetic_conditional_on_predictors",names(listobj$lower.summary))){listobj$lower.summary$RSS_non_phylogenetic_conditional_on_predictors<-NULL}
	if (is.element("R2_non_phylogenetic_conditional_on_predictors",names(listobj$lower.summary))){listobj$lower.summary$R2_non_phylogenetic_conditional_on_predictors<-NULL}
	if (is.element("RSS_non_phylogenetic_conditional_on_predictors_comment",names(listobj$lower.summary))){listobj$lower.summary$RSS_non_phylogenetic_conditional_on_predictors_comment<-NULL}
	if (is.element("R2_non_phylogenetic_conditional_on_predictors_comment",names(listobj$lower.summary))){listobj$lower.summary$R2_non_phylogenetic_conditional_on_predictors_comment<-NULL}
	if (is.element("RSS_conditional_on_predictors",names(listobj$lower.summary))){listobj$lower.summary$RSS_conditional_on_predictors<-NULL}
	if (is.element("R2_conditional_on_predictors",names(listobj$lower.summary))){listobj$lower.summary$R2_conditional_on_predictors<-NULL}
	if (is.element("RSS_conditional_on_predictors_comment",names(listobj$lower.summary))){listobj$lower.summary$RSS_conditional_on_predictors_comment<-NULL}
	if (is.element("R2_conditional_on_predictors_comment",names(listobj$lower.summary))){listobj$lower.summary$R2_conditional_on_predictors_comment<-NULL}
    }
    if (is.element("upper.summary",names(listobj))){
        if (is.element("LogLik",names(listobj$upper.summary))){listobj$upper.summary$LogLik<-NULL}
	if (is.element("dof",names(listobj$upper.summary))){listobj$upper.summary$dof<-NULL}
	if (is.element("m2loglik",names(listobj$upper.summary))){listobj$upper.summary$m2loglik<-NULL}
	if (is.element("aic",names(listobj$upper.summary))){listobj$upper.summary$aic<-NULL}
	if (is.element("aic.c",names(listobj$upper.summary))){listobj$upper.summary$aic.c<-NULL}
	if (is.element("sic",names(listobj$upper.summary))){listobj$upper.summary$sic<-NULL}
	if (is.element("bic",names(listobj$upper.summary))){listobj$upper.summary$bic<-NULL}
	if (is.element("RSS",names(listobj$upper.summary))){listobj$upper.summary$RSS<-NULL}
	if (is.element("R2",names(listobj$upper.summary))){listobj$upper.summary$R2<-NULL}
	if (is.element("R2_phylaverage",names(listobj$upper.summary))){listobj$upper.summary$R2_phylaverage<-NULL}
	if (is.element("RSS_comment",names(listobj$upper.summary))){listobj$upper.summary$RSS_comment<-NULL}
	if (is.element("R2_comment",names(listobj$upper.summary))){listobj$upper.summary$R2_comment<-NULL}
	if (is.element("R2_phylaverage_comment",names(listobj$upper.summary))){listobj$upper.summary$R2_phylaverage_comment<-NULL}
	if (is.element("RSS_non_phylogenetic",names(listobj$upper.summary))){listobj$upper.summary$RSS_non_phylogenetic<-NULL}
	if (is.element("R2_non_phylogenetic",names(listobj$upper.summary))){listobj$upper.summary$R2_non_phylogenetic<-NULL}
	if (is.element("RSS_non_phylogenetic_comment",names(listobj$upper.summary))){listobj$upper.summary$RSS_non_phylogenetic_comment<-NULL}
	if (is.element("R2_non_phylogenetic_comment",names(listobj$upper.summary))){listobj$upper.summary$R2_non_phylogenetic_comment<-NULL}
	if (is.element("RSS_non_phylogenetic_conditional_on_predictors",names(listobj$upper.summary))){listobj$upper.summary$RSS_non_phylogenetic_conditional_on_predictors<-NULL}
	if (is.element("R2_non_phylogenetic_conditional_on_predictors",names(listobj$upper.summary))){listobj$upper.summary$R2_non_phylogenetic_conditional_on_predictors<-NULL}
	if (is.element("RSS_non_phylogenetic_conditional_on_predictors_comment",names(listobj$upper.summary))){listobj$upper.summary$RSS_non_phylogenetic_conditional_on_predictors_comment<-NULL}
	if (is.element("R2_non_phylogenetic_conditional_on_predictors_comment",names(listobj$upper.summary))){listobj$upper.summary$R2_non_phylogenetic_conditional_on_predictors_comment<-NULL}
	if (is.element("RSS_conditional_on_predictors",names(listobj$upper.summary))){listobj$upper.summary$RSS_conditional_on_predictors<-NULL}
	if (is.element("R2_conditional_on_predictors",names(listobj$upper.summary))){listobj$upper.summary$R2_conditional_on_predictors<-NULL}
	if (is.element("RSS_conditional_on_predictors_comment",names(listobj$upper.summary))){listobj$upper.summary$RSS_conditional_on_predictors_comment<-NULL}
	if (is.element("R2_conditional_on_predictors_comment",names(listobj$upper.summary))){listobj$upper.summary$R2_conditional_on_predictors_comment<-NULL}
    }
    if (is.element("regressCovar",names(listobj))){listobj$regressCovar<-NULL}    
    if (is.element("parameter_signs",names(listobj))){## remove trailing parameter_signs field if it was NULL, i.e. not used
	if (is.null(listobj$parameter_signs)){listobj$parameter_signs<-NULL}
    }
    sapply(listobj,function(obj,vregime.names,vY.names,vX.names,predictors,EvolModel){if (is.list(obj)){.correct.names(obj,vregime.names,vY.names,vX.names,predictors,EvolModel)}else{obj}},vregime.names=vregime.names,vY.names=vY.names,vX.names=vX.names,predictors=predictors,EvolModel=EvolModel,simplify=FALSE)
}

.set.tol<-function(method){
    tol=switch(method,
	glsgc=c(0.01,0.01),
	maxlik=NA
    )
    tol
}

.set.maxiter<-function(method,maxiter=c(10,50,100)){
## make this a function of numtips and remove not needed methods
    if ((!is.vector(maxiter)) || (!is.numeric(maxiter)) || (length(maxiter)!=3)){
	maxiter<-c(10,50,100);
	.my_warning("WARNING: maxiter passed in a wrong way setting it to default of c(10,50,100)",TRUE,FALSE)
    }
    maxiter=switch(method,
	glsgc=maxiter,
	maxlik=NA
    )
    maxiter
}

.set.estimparams<-function(params,kY,kX,numregs,estimate.root.state=FALSE,tree_height=NULL,cBestim_method="ML"){
    if (!is.element("EstimParams",names(params))){EstimParams<-list()}
    else{EstimParams<-params$EstimParams}
    EstimParams$vVars<-NULL
    EstimParams$conditional<-FALSE

## =====================================================================
## at the moment we do not allow for control of these    
    EstimParams$signsvX0<-NULL
    EstimParams$signsSxx<-NULL
    EstimParams$signsSxy<-NULL
    EstimParams$signsSyx<-NULL
    EstimParams$signsmPsi0<-NULL
## =====================================================================

    if((is.element("pcmbase_model_box",names(params)))&&(!is.element("pcmbase_model_box",names(EstimParams)))){EstimParams$pcmbase_model_box<-params$pcmbase_model_box}

    if (!is.element("cCalc.type",names(EstimParams))){EstimParams$cCalc.type<-"PCMBase"}
    if (!is.element("estimate.root.state",names(EstimParams))){EstimParams$estimate.root.state<-estimate.root.state}
    
    if (!is.element("kY",names(EstimParams))){EstimParams$kY<-kY}
    if (!is.element("kX",names(EstimParams))){EstimParams$kX<-kX}
    if (!is.element("Atype",names(EstimParams))){EstimParams$Atype<-"Invertible"}
    ## when B is done by GLS, line below is commented out 
    #if (!is.element("Btype",names(EstimParams))){EstimParams$Btype<-"Any"}
    ## when B is done by ML line below is needed
    if (cBestim_method=="ML"){if (!is.element("Btype",names(EstimParams))){EstimParams$Btype<-"Any"}}
    ## =====================================
    if (!is.element("diagA",names(EstimParams))){EstimParams$diagA<-NULL}
    if (!is.element("Syytype",names(EstimParams))){EstimParams$Syytype<-"UpperTri"}
    if (!is.element("diagSyy",names(EstimParams))){EstimParams$diagSyy<-"Positive"}

    if (!is.element("maxAabsval",names(EstimParams))){if (is.null(tree_height)){EstimParams$maxAabsval<-100}else{EstimParams$maxAabsval<-log(2)/(0.005*tree_height)}}
    if (!is.element("maxSyyabsval",names(EstimParams))){EstimParams$maxSyyabsval<-5*EstimParams$maxAabsval}


    if (!is.element("signsA",names(EstimParams))){EstimParams$signsA<-NULL}
    ## B is done by ML ONLY (but same line for GLS)
    if ((cBestim_method=="GLS") || (cBestim_method=="ML")){
        if (!is.element("signsB",names(EstimParams))){EstimParams$signsB<-NULL}
    }
    ## =====================================
    if (!is.element("signsmPsi",names(EstimParams))){EstimParams$signsmPsi<-NULL}
    if (!is.element("signsmPsi0",names(EstimParams))){EstimParams$signsmPsi0<-NULL}
    if (!is.element("signsvY0",names(EstimParams))){EstimParams$vY0<-NULL}
    if (!is.element("signsvX0",names(EstimParams))){EstimParams$vX0<-NULL}
    if (!is.element("signsSyy",names(EstimParams))){EstimParams$Syy<-NULL}    
    if (!is.element("signsSyx",names(EstimParams))){EstimParams$Syx<-NULL}
    if (!is.element("signsSxy",names(EstimParams))){EstimParams$Sxy<-NULL}
    if (!is.element("signsSxx",names(EstimParams))){EstimParams$Sxx<-NULL}        
    if (!is.element("maximMethod",names(EstimParams))){EstimParams$maximMethod<-"optim"}
    if (!is.element("conf.level",names(EstimParams))){EstimParams$conf.level<-0.95}
    if (!is.element("calcCI",names(EstimParams))){EstimParams$calcCI<-FALSE}
    if (!is.element("designToEstim",names(EstimParams))){EstimParams$designToEstim<-list()}
    if (!is.element("y0",names(EstimParams$designToEstim))){EstimParams$designToEstim$y0<-TRUE}
    if (!is.element("psi",names(EstimParams$designToEstim))){EstimParams$designToEstim$psi<-TRUE}
    if (!is.element("psi0",names(EstimParams$designToEstim))){EstimParams$designToEstim$psi0<-FALSE}
    if ((!is.element("X0",names(EstimParams$designToEstim)))&&(params$EvolModel=="mvslouch")){EstimParams$designToEstim$X0<-FALSE}
    if (cBestim_method=="GLS"){
	## B is done by GLS 
	if ((!is.element("B",names(EstimParams$designToEstim)))&&(params$EvolModel=="mvslouch")){EstimParams$designToEstim$B<-TRUE}
    }
    if (cBestim_method=="ML"){
	## B is done by ML 
        if ((!is.element("B",names(EstimParams$designToEstim)))&&(params$EvolModel=="mvslouch")){EstimParams$designToEstim$B<-FALSE}    
    }
    ## ==============================================
    
    if ((!is.element("BX0",names(EstimParams$designToEstim)))&&(params$EvolModel=="mvslouch")){EstimParams$designToEstim$BX0<-FALSE}
    if ((!is.element("UseX0",names(EstimParams$designToEstim)))&&(params$EvolModel=="mvslouch")){EstimParams$designToEstim$UseX0<-TRUE}	
    if ((!is.element("y0AncState",names(EstimParams$designToEstim))&&(!EstimParams$estimate.root.state))||(numregs==1)){EstimParams$designToEstim$y0AncState<-TRUE}
    if (!is.element("y0AncState",names(EstimParams$designToEstim))&&(numregs>1)&&(EstimParams$estimate.root.state)){EstimParams$designToEstim$y0AncState<-FALSE}
    if ((!is.element("y0OnlyFixed",names(EstimParams$designToEstim)))&&(params$EvolModel=="mvslouch")){EstimParams$designToEstim$y0OnlyFixed<-FALSE}
    if ((cBestim_method=="GLS") || (cBestim_method=="ML")){
	## should phylogenetic correlations be taken into account when setting up the regression design matrix (conditional expectation E[Y|X]) (FALSE)
	## or not, i.e. do a simple regression (TRUE)
	## FALSE means that one caculates from E[Y|X] Cov(Y,X)%*%Var(X)%*%(X-EX), where Y, X are the response and predictor traits of 
	## all the sepcies stacked on each other and the phylogeny sits instide Cov(Y,X) and Var(X)
	## However, this means that the estimation procedure looses the O(n) fast likelihood evaluation property, as it is impossible to 
	## get the Cov(Y,X)%*%Var(X) part of the design matrix without first calculating Var((Y,X)) in O(n^2) time
	## TRUE means that we assume that Var((Y,X)) is diagonal, hence we can do all operations in O(n) as we just calculate
	## the variance blocks on the diagonal. Tip heights are taken into account, so some phylogenetic information on the species is retained.
	if (!is.element("SimpReg",names(EstimParams$designToEstim))){EstimParams$designToEstim$SimpReg<-TRUE}    	
    }
    ## ================================================================
    if ((!is.element("iRegLin",names(EstimParams$designToEstim)))&&(params$EvolModel=="mvslouch")){EstimParams$designToEstim$iRegLin<-TRUE}
    if ((!is.element("optimMethod",names(EstimParams$designToEstim)))&&(params$EvolModel=="mvslouch")){EstimParams$designToEstim$optimMethod<-"optim"}
    if ((!is.element("BFullXNA",names(EstimParams$designToEstim)))&&(params$EvolModel=="mvslouch")){EstimParams$designToEstim$BFullXNA<-TRUE}	
    if ((!is.element("FullNAYX",names(EstimParams$designToEstim)))&&(params$EvolModel=="mvslouch")){EstimParams$designToEstim$FullNAYX<-TRUE}	
    if (!is.element("FullNAY",names(EstimParams$designToEstim))){EstimParams$designToEstim$FullNAY<-TRUE}	
    if (!is.element("sigmaRule",names(EstimParams$designToEstim))){EstimParams$designToEstim$sigmaRule<-3}	
    if (cBestim_method=="GLS"){
	## B is done by GLS 
	if (!is.element("YnonCondX",names(EstimParams$designToEstim))){EstimParams$designToEstim$YnonCondX<-FALSE}	
    }
    if (cBestim_method=="ML"){
	## B is done by ML
	if (!is.element("YnonCondX",names(EstimParams$designToEstim))){EstimParams$designToEstim$YnonCondX<-TRUE}	
	##if (!is.element("YnonCondX",names(EstimParams$designToEstim))){EstimParams$designToEstim$YnonCondX<-FALSE}	
    }
    ## ================================================================
    EstimParams$designToEstim$Atype<-EstimParams$Atype
    EstimParams$designToEstim$Btype<-EstimParams$Btype
    EstimParams$designToEstim$Syytype<-EstimParams$Syytype
    
## =====================================================================================
## for parameters estimated by regression we do not allow fixing of sign
    if (is.element("signsmPsi",names(EstimParams))){
	if (EstimParams$designToEstim$psi){
	    EstimParams$signsmPsi[which(EstimParams$signsmPsi=="-")]<-NA
    	    EstimParams$signsmPsi[which(EstimParams$signsmPsi=="+")]<-NA
    	    EstimParams$signsmPsi<-as.numeric(EstimParams$signsmPsi)
    	    EstimParams$designToEstim$signsmPsi<-EstimParams$signsmPsi ## This is needed as a copy for phylgls
        }else{
    	    EstimParams$designToEstim$signsmPsinonGLS<-EstimParams$signsmPsi
        }
    }
    if (is.element("signsmPsi0",names(EstimParams))){
	if (EstimParams$designToEstim$psi0){
	    EstimParams$signsmPsi0[which(EstimParams$signsmPsi0=="-")]<-NA
    	    EstimParams$signsmPsi0[which(EstimParams$signsmPsi0=="+")]<-NA
    	    EstimParams$signsmPsi0<-as.numeric(EstimParams$signsmPsi0)
    	    EstimParams$designToEstim$signsmPsi0<-EstimParams$signsmPsi0 ## This is needed as a copy for phylgls
	}else{
	    EstimParams$designToEstim$signsmPsi0nonGLS<-EstimParams$signsmPsi0
	}
    }
    if (is.element("signsvY0",names(EstimParams))){
	if (EstimParams$designToEstim$y0){
	    EstimParams$signsvY0[which(EstimParams$signsvY0=="-")]<-NA
    	    EstimParams$signsvY0[which(EstimParams$signsvY0=="+")]<-NA
    	    EstimParams$signsvY0<-as.numeric(EstimParams$signsvY0)
    	    EstimParams$designToEstim$signsvY0<-EstimParams$signsvY0 ## This is needed as a copy for phylgls
    	}else{
    	    EstimParams$designToEstim$signsvY0nonGLS<-EstimParams$signsvY0
    	}
    }
    if (is.element("signsB",names(EstimParams))){
	if (EstimParams$designToEstim$B){
	    EstimParams$signsB[which(EstimParams$signsB=="-")]<-NA
    	    EstimParams$signsB[which(EstimParams$signsB=="+")]<-NA
    	    EstimParams$signsB<-as.numeric(EstimParams$signsB)
    	    EstimParams$designToEstim$signsB<-EstimParams$signsB ## This is needed as a copy for phylgls
    	}else{
    	    EstimParams$designToEstim$signsBnonGLS<-EstimParams$signsB 
    	}
    }
    if (is.element("signsA",names(EstimParams))){
    	    EstimParams$designToEstim$signsAnonGLS<-EstimParams$signsA    	
    }
    if (is.element("signsSyy",names(EstimParams))){
    	    EstimParams$designToEstim$signsSyynonGLS<-EstimParams$signsSyy
    }
    if (is.element("signsSxy",names(EstimParams))){
    	    EstimParams$designToEstim$signsSXynonGLS<-EstimParams$signsSxy
    }
    if (is.element("signsSyx",names(EstimParams))){
    	    EstimParams$designToEstim$signsSyxnonGLS<-EstimParams$signsSyx
    }
##    if (is.element("signsSxx",names(EstimParams))){
##    ## This part of code should never be done at the moment
##    	    EstimParams$designToEstim$signsSxxnonGLS<-EstimParams$signsSxx
##    }
##    if (is.element("signsvX0",names(EstimParams))){
##    ## This part of code should never be done at the moment
##    	    EstimParams$designToEstim$signsvX0nonGLS<-EstimParams$signsvX0
##    }

## =======================================================================================
    
    if (!is.element("Fixed",names(EstimParams))){EstimParams$Fixed<-list()}
    if (!is.element("mPsi",names(EstimParams$Fixed))){EstimParams$Fixed$mPsi<-NA}
    if (!is.element("mPsi0",names(EstimParams$Fixed))){	if ((kY>0)&&((numregs==1)||(!EstimParams$designToEstim$psi0)||(!EstimParams$designToEstim$y0AncState))){EstimParams$Fixed$mPsi0<-matrix(0,ncol=1,nrow=kY)}else{EstimParams$Fixed$mPsi0<-NA}}
    ## No need to correct for estimating both mPsi0 and Y0 if there is only one regime as this mPsi0 is immediately zeroed in one regime case
    if (!is.element("vY0",names(EstimParams$Fixed))){EstimParams$Fixed$vY0<-NA}
    if (!is.element("vX0",names(EstimParams$Fixed))){EstimParams$Fixed$vX0<-NA}
    if (!is.element("Sxx",names(EstimParams$Fixed))){EstimParams$Fixed$Sxx<-NA}
    ## B GLS needs this so that there is no attempts to optize over B
    if (cBestim_method=="GLS"){
	if (!is.element("B",names(EstimParams$Fixed))){EstimParams$Fixed$B<-NA}
    }
    ## ================================================
    if (!is.element("Sxy",names(EstimParams$Fixed))){if((kX>0)&&(kY>0)){EstimParams$Fixed$Sxy<-matrix(0,nrow=kX,ncol=kY)}else{EstimParams$Fixed$Sxy<-NA}}    
    if (!is.element("Syx",names(EstimParams$Fixed))){if((kX>0)&&(kY>0)){EstimParams$Fixed$Syx<-matrix(0,nrow=kY,ncol=kX)}else{EstimParams$Fixed$Syx<-NA}}        
    if (!is.element("KnownParams",names(EstimParams))){
	vKnownNames<-c()
	EstimParams$KnownParams<-list()
	j<-1	
	if (length(EstimParams$Fixed)>0){#&&(length(EstimParams$KnownParams)>0)){
	    for (i in 1:length(EstimParams$Fixed)){
		if (!is.na(EstimParams$Fixed[[i]][1])){
		    EstimParams$KnownParams[[j]]<-EstimParams$Fixed[[i]]
		    vKnownNames<-c(vKnownNames,names(EstimParams$Fixed)[i])
		    j<-j+1
		}	    
	    }
	}
	if (j>1){names(EstimParams$KnownParams)<-vKnownNames}
    }else{
	vListNull<-c()
	if ((length(EstimParams$Fixed)>0)&&(length(EstimParams$KnownParams)>0)){	
	    for (j in 1:length(EstimParams$KnownParams)){
		i<-which(names(EstimParams$Fixed)==names(EstimParams$KnownParams)[j])
		if (is.na(EstimParams$Fixed[[i]][1])){vListNull<-c(vListNull,names(EstimParams$KnowParams)[j])}
	    }	
	}
	if (length(vListNull)>0){for(par in vListNull){EstimParams$KnownParams[[which(names(EstimParams$KnownParams)==par)]]<-NULL}}
    }
    
    
    if ((!is.element("optim_parscale",names(EstimParams)))||(is.null(EstimParams$optim_parscale))){
	EstimParams$optim_parscale<-c("parscale_A"=1,"logparscale_A"=4,"logparscale_other"=2) ## (3,5,1) is also a good option from experiments
    }else{
	if ((!is.element("parscale_A",names(EstimParams$optim_parscale)))||(is.null(EstimParams$optim_parscale["parscale_A"]))||(is.na(EstimParams$optim_parscale["parscale_A"]))){
	    EstimParams$optim_parscale["parscale_A"]<-1
	}
	if ((!is.element("logparscale_A",names(EstimParams$optim_parscale)))||(is.null(EstimParams$optim_parscale["logparscale_other"]))||(is.na(EstimParams$optim_parscale["logparscale_other"]))){
	    EstimParams$optim_parscale["logparscale_A"]<-4
	}
	if ((!is.element("logparscale_other",names(EstimParams$optim_parscale)))||(is.null(EstimParams$optim_parscale["logparscale_other"]))||(is.na(EstimParams$optim_parscale["logparscale_other"]))){
	    EstimParams$optim_parscale["logparscale_other"]<-2
	}
    }
    

    setparres<-.set.paramatrizationnames(EstimParams,params$EvolModel,kY,kX,numregs)    
    parnames<-setparres$parnames
    EstimParams$StartPoint<-rnorm(length(parnames),sd=3) ## we are drawing parameters here that might be overwritten in the next step, but this is negligable w.r.t time optimality
    names(EstimParams$StartPoint)<-parnames

    if (!is.null(EstimParams$lStartPoint)){
	vprovided_params<-.par.inv.transform(EstimParams$lStartPoint,EstimParams)
	EstimParams$StartPoint[names(vprovided_params)]<-vprovided_params
	EstimParams$lStartPoint<-NULL
    }
    EstimParams$parscale<-setparres$parscale ##.set_parscale(parnames)
    
    EstimParams
}


.set_parscale<-function(parnames){
    v_parscale<-rep(1,length(parnames))
    v_parscale[which(sapply(parnames,function(x){substring(x, 1, 1)=="A"},simplify=TRUE))]<- 3
    v_parscale
}

.correct_for_diagonalSigns_parscale<-function(parscale,EstimParams,paramname,logelement_scale,loglogelement_scale,scalarelement_scale,kY,kX=NULL){
	if (is.null(kX)){kX<-kY}
	if (is.element(paste("diag",paramname,sep=""),names(EstimParams))){
	    diagtype<-EstimParams[[which(names(EstimParams)==paste("diag",paramname,sep=""))]]
	    bislogdiag<-FALSE
	    if (!is.null(diagtype)){
		bislogdiag<-is.element(diagtype,c("Positive","Negative"))
	    }
	    
	    matrixtype<-EstimParams[[which(names(EstimParams)==paste(paramname,"type",sep=""))]]
	    vdiag_els<-switch(matrixtype,
		     SingleValueDiagonal={1},
                     Diagonal={1:kY},
                     UpperTri={diag(.par.transform.uppertri.matrix(1:((kY+1)*kY/2),kY))},
                     LowerTri={diag(.par.transform.lowertri.matrix(1:((kY+1)*kY/2),kY))},
                     TwoByTwo={c(1,3)},
                     DecomposablePositive={1:kY},
                     DecomposableNegative={1:kY},
                     Invertible={(kY^2-kY+1):(kY^2)},
                     Any={diag(matrix(1:(kY*kX),kY,kX,byrow=TRUE))},
                     0
	    )
	    
	    if(is.element(paste("signs",paramname,sep=""),names(EstimParams))){
		matrixsigns<-EstimParams[[which(names(EstimParams)==paste("signs",paramname,sep=""))]]
                vnotNA<-which(!is.na(matrixsigns))
		tologscale<-c(which(matrixsigns=="+"),which(matrixsigns=="-"))
		matrixsigns[tologscale]<-NA
		ondiag<-intersect(tologscale,vdiag_els)
                if (matrixtype=="Any"){
            	    if (length(tologscale)>0){
            		if (length(vnotNA)==0)
            		    parscale[tologscale]<-logelement_scale                	    
                	    if ((length(ondiag)>0)&&(bislogdiag)){parscale[ondiag]<-loglogelement_scale}
                	}else{
                	    tmpparscale<-matrix(scalarelement_scale,nrow=nrow(matrixsigns),ncol=ncol(matrixsigns))
                	    tmpparscale[tologscale]<-logelement_scale
                	    tmpparscale[ondiag]<-loglogelement_scale
                	    tmpparscale[which(is.na(matrixsigns))]<-NA
                	    parscale<-c(t(tmpparscale))
                	    parscale<-setdiff(parscale,NA)                	    
                	}
            	    }
            	if ((matrixtype=="SingleValueDiagonal")||(matrixtype=="Diagonal")){
            	    if ((length(ondiag)>0)&&(bislogdiag)){parscale[ondiag]<-loglogelement_scale}
            	    if ((length(ondiag)>0)&&(!bislogdiag)){parscale[ondiag]<-logelement_scale}
            	}
            }else{if ((vdiag_els[1]!=0)&&(bislogdiag)){parscale[vdiag_els]<-logelement_scale}}
	}
	parscale
}

.set.paramatrizationnames<-function(EstimParams,EvolModel,kY,kX,numregs){
    parnames<-c()
    parscale<-c()
    
    scalarelementA_scale<-EstimParams$optim_parscale[["parscale_A"]]
    logelementA_scale<-EstimParams$optim_parscale[["logparscale_A"]]
    loglogelementA_scale<-2*EstimParams$optim_parscale[["logparscale_A"]]

    logelement_scale<-1
    loglogelement_scale<-EstimParams$optim_parscale[["logparscale_other"]]
    scalarelement_scale<-2*EstimParams$optim_parscale[["logparscale_other"]]
    
    if (EvolModel=="bm"){
    	if (!is.element("vY0",names(EstimParams$Fixed))){parnames<-c(parnames,.generatenames("vY0",kY))}
    	if (!is.element("Sxx",names(EstimParams$Fixed))){
	    if (is.element("Sxxtype",names(EstimParams))){
		Sxxlength=switch(EstimParams$Sxxtype,
		     SingleValueDiagonal={1},
                     Diagonal={kY},
                     Symmetric={(kY+1)*kY/2},
                     UpperTri={(kY+1)*kY/2},
                     LowerTri={(kY+1)*kY/2},
                     Any={kY^2}
		)		
		parnames<-c(parnames,.generatenames("Sxx",Sxxlength))
		parscale<-c(parscale,.correct_for_diagonalSigns_parscale(rep(scalarelement_scale,Sxxlength),EstimParams,"Sxx",logelement_scale,loglogelement_scale,scalarelement_scale,kY))
	    }
	}	

    }
    if (EvolModel=="ouch"){
    	if (!is.element("A",names(EstimParams$Fixed))){
	    if (is.element("Atype",names(EstimParams))){
		Alength=switch(EstimParams$Atype,
		     SingleValueDiagonal={1},
                     Diagonal={kY},
                     SymmetricPositiveDefinite={(kY+1)*kY/2},
                     Symmetric={(kY+1)*kY/2},
                     TwoByTwo={4},
                     UpperTri={(kY+1)*kY/2},
                     LowerTri={(kY+1)*kY/2},
                     DecomposablePositive={kY^2},
                     DecomposableNegative={kY^2},
                     DecomposableReal={kY^2},
                     Invertible={kY^2},
                     Any={
                        num_vals<-kY^2
                        if (is.element("signsA",names(EstimParams))){
                    	    vToNA<-c(which(EstimParams$signsA=="+"),which(EstimParams$signsA=="-"))
                    	    if (length(vToNA)>0){EstimParams$signsA[vToNA]<-NA}
                    	    num_vals<-num_vals-length(which(!is.na(EstimParams$signsA)))
                        }
                        num_vals
                     }
		)
		parnames<-c(parnames,.generatenames("A",Alength))
		parscale<-c(parscale,.correct_for_diagonalSigns_parscale(rep(scalarelementA_scale,Alength),EstimParams,"A",logelementA_scale,loglogelementA_scale,scalarelementA_scale,kY))
	    }
	}
	if (!is.element("vY0",names(EstimParams$Fixed))){parnames<-c(parnames,.generatenames("vY0",kY))}
	if (!is.element("mPsi",names(EstimParams$Fixed))){
	    if (is.element("mPsitype",names(EstimParams))){
		psilength=switch(EstimParams$mPsitype,
	    	     Global={kY},
	    	     Regimes={kY*numregs}
		)
		if (is.element("signsmPsi",names(EstimParams))){
                    	vToNA<-c(which(EstimParams$signsmPsi=="+"),which(EstimParams$signsmPsi=="-"))
                    	if (length(vToNA)>0){EstimParams$signsmPsi[vToNA]<-NA}
                    	psilength<-psilength-length(which(!is.na(EstimParams$signsmPsi)))
                }                    
		parnames<-c(parnames,.generatenames("Psi",psilength))
	    }
	}
	if (!is.element("mPsi0",names(EstimParams$Fixed))){parnames<-c(parnames,.generatenames("Psi0",kY))}
	if (!is.element("Syy",names(EstimParams$Fixed))){
	    if (is.element("Syytype",names(EstimParams))){
		Syylength=switch(EstimParams$Syytype,
		     SingleValueDiagonal={1},
                     Diagonal={kY},
                     Symmetric={(kY+1)*kY/2},
                     UpperTri={
                        (kY+1)*kY/2                        
                     },
                     LowerTri={
                        (kY+1)*kY/2
                     },
                     Any={
                        num_vals<-kY^2
                        if (is.element("signsSyy",names(EstimParams))){
                    	    vToNA<-c(which(EstimParams$signsSyy=="+"),which(EstimParams$signsSyy=="-"))
                    	    if (length(vToNA)>0){EstimParams$signsSyy[vToNA]<-NA}
                    	    num_vals<-num_vals-length(which(!is.na(EstimParams$signsSyy)))
                        }
                        num_vals                        
                     }
		)
		parnames<-c(parnames,.generatenames("Syy",Syylength))
		parscale<-c(parscale,.correct_for_diagonalSigns_parscale(rep(scalarelement_scale,Syylength),EstimParams,"Syy",logelement_scale,loglogelement_scale,scalarelement_scale,kY))
	    }
	}	   
    }
    if (EvolModel=="slouch"){}
    if (EvolModel=="mvslouch"){
	if (!is.element("A",names(EstimParams$Fixed))){
	    if (is.element("Atype",names(EstimParams))){
		Alength=switch(EstimParams$Atype,
		     SingleValueDiagonal={1},
                     Diagonal={kY},
                     SymmetricPositiveDefinite={(kY+1)*kY/2},
                     Symmetric={(kY+1)*kY/2},
                     TwoByTwo={4},
                     UpperTri={(kY+1)*kY/2},
                     LowerTri={(kY+1)*kY/2},
                     DecomposablePositive={kY^2},
                     DecomposableNegative={kY^2},
                     DecomposableReal={kY^2},
                     Invertible={kY^2},
                     Any={
                        num_vals<-kY^2
                        if (is.element("signsA",names(EstimParams))){
                    	    vToNA<-c(which(EstimParams$signsA=="+"),which(EstimParams$signsA=="-"))
                    	    if (length(vToNA)>0){EstimParams$signsA[vToNA]<-NA}
                    	    num_vals<-num_vals-length(which(!is.na(EstimParams$signsA)))
                        }
                        num_vals
                     }
		)
		parnames<-c(parnames,.generatenames("A",Alength))
		parscale<-c(parscale,.correct_for_diagonalSigns_parscale(rep(scalarelementA_scale,Alength),EstimParams,"A",logelementA_scale,loglogelementA_scale,scalarelementA_scale,kY))
	    }
	}
	if (!is.element("B",names(EstimParams$Fixed))){
	    if (is.element("Btype",names(EstimParams))){
		Blength=switch(EstimParams$Btype,
	    	     MinusA={0},
	    	     SingleValue={1},
	    	     SingleValueDiagonal={1},
                     Diagonal={min(kY,kX)},
                     Symmetric={min(kY,kX)*(min(kY,kX)+1)/2+kY*kX-min(kY,kX)^2},
                     Any={
                        num_vals<-kY*kX
                        if (is.element("signsB",names(EstimParams))){
                    	    vToNA<-c(which(EstimParams$signsB=="+"),which(EstimParams$signsB=="-"))
                    	    if (length(vToNA)>0){EstimParams$signsB[vToNA]<-NA}
                    	    num_vals<-num_vals-length(which(!is.na(EstimParams$signsB)))
                        }
                        num_vals
                     }
		)
		parnames<-c(parnames,.generatenames("B",Blength))
		parscale<-c(parscale,.correct_for_diagonalSigns_parscale(rep(scalarelement_scale,Blength),EstimParams,"B",logelement_scale,loglogelement_scale,scalarelement_scale,kY,kX))

	    }
	}
	if (!is.element("vX0",names(EstimParams$Fixed))){parnames<-c(parnames,.generatenames("vX0",kX))}	
	if (!is.element("vY0",names(EstimParams$Fixed))){parnames<-c(parnames,.generatenames("vY0",kY))}
	if (!is.element("mPsi",names(EstimParams$Fixed))){
	    if (is.element("mPsitype",names(EstimParams))){
		psilength=switch(EstimParams$mPsitype,
	    	     Global={kY},
	    	     Regimes={kY*numregs}
		)
		if (is.element("signsmPsi",names(EstimParams))){
                    	vToNA<-c(which(EstimParams$signsmPsi=="+"),which(EstimParams$signsmPsi=="-"))
                    	if (length(vToNA)>0){EstimParams$signsmPsi[vToNA]<-NA}
                    	psilength<-psilength-length(which(!is.na(EstimParams$signsmPsi)))
                }
		parnames<-c(parnames,.generatenames("Psi",psilength))
	    }
	}
	if (!is.element("mPsi0",names(EstimParams$Fixed))){parnames<-c(parnames,.generatenames("mPsi0",kY))}
	if (!is.element("Syy",names(EstimParams$Fixed))){
	    if (is.element("Syytype",names(EstimParams))){
		Syylength=switch(EstimParams$Syytype,
		     SingleValueDiagonal={1},
                     Diagonal={kY},
                     Symmetric={(kY+1)*kY/2},
                     UpperTri={
                        (kY+1)*kY/2                        
                     },
                     LowerTri={
                      (kY+1)*kY/2
                    },
                     Any={
                        num_vals<-kY^2
                        if (is.element("signsSyy",names(EstimParams))){
                    	    vToNA<-c(which(EstimParams$signsSyy=="+"),which(EstimParams$signsSyy=="-"))
                    	    if (length(vToNA)>0){EstimParams$signsSyy[vToNA]<-NA}
                    	    num_vals<-num_vals-length(which(!is.na(EstimParams$signsSyy)))
                        }
			num_vals

                    }
		)
		parnames<-c(parnames,.generatenames("Syy",Syylength))
		parscale<-c(parscale,.correct_for_diagonalSigns_parscale(rep(scalarelement_scale,Syylength),EstimParams,"Syy",logelement_scale,loglogelement_scale,scalarelement_scale,kY))
	    }
	}	
	if (!is.element("Syx",names(EstimParams$Fixed))){
	    num_vals<-kY*kX
	    if (is.element("signsSyx",names(EstimParams))){
                    	    vToNA<-c(which(EstimParams$signsSyx=="+"),which(EstimParams$signsSyx=="-"))
                    	    if (length(vToNA)>0){EstimParams$signsSyx[vToNA]<-NA}
                    	    num_vals<-num_vals-length(which(!is.na(EstimParams$signsSyx)))
            }
	    parnames<-c(parnames,.generatenames("Syx",num_vals))
	}	
	if (!is.element("Sxy",names(EstimParams$Fixed))){
	    num_vals<-kY*kX
	    if (is.element("signsSxy",names(EstimParams))){
                    	    vToNA<-c(which(EstimParams$signsSxy=="+"),which(EstimParams$signsSxy=="-"))
                    	    if (length(vToNA)>0){EstimParams$signsSxy[vToNA]<-NA}
                    	    num_vals<-num_vals-length(which(!is.na(EstimParams$signsSxy)))
            }
	    parnames<-c(parnames,.generatenames("Sxy",num_vals))
	}	
	if (!is.element("Sxx",names(EstimParams$Fixed))){
	    if (is.element("Sxxtype",names(EstimParams))){
		Sxxlength=switch(EstimParams$Sxxtype,
		     OneSigmaDiagonal={1},
                     Diagonal={kX},
                     Symmetric={(kX+1)*kX/2},
                     Any={kX^2}
		)
		parnames<-c(parnames,.generatenames("Sxx",Sxxlength))
	    	parscale<-c(parscale,.correct_for_diagonalSigns_parscale(rep(scalarelement_scale,Sxxlength),EstimParams,"Sxx",logelement_scale,loglogelement_scale,scalarelement_scale,kY))
	    }
	}	
    }    
    list(parnames=parnames,parscale=parscale)
}

.generatenames<-function(prefix,len){
    gennednames<-c()	
    if (len==1){gennednames<-c(prefix)}
    if (len==2){gennednames<-c(paste(prefix,"start",sep=""),paste(prefix,"end",sep=""))}
    if (len>2){
	gennednames<-sapply(1:len,function(x,prefix){paste(prefix,"_",x,sep="")},prefix=prefix)
	gennednames[1]<-paste(prefix,"start",sep="")
	gennednames[len]<-paste(prefix,"end",sep="")	                             
    }	    	    
    gennednames
}

.createMeasurementError<-function(M_error,n,kYX){
    bM_errorOK<-FALSE
    if (is.null(M_error)){M_error<-NULL;bM_errorOK<-TRUE}
    if (is.vector(M_error,mode="numeric")){
	if (((length(M_error)==1)|| (length(M_error)==kYX)) && all(M_error>=0)){
	    M_error<-diag(M_error,nrow=kYX,ncol=kYX)
	    M_error<-sapply(1:n,function(i,x){x},x=M_error,simplify=FALSE)
	    bM_errorOK<-TRUE
	}
    }
    if (is.matrix(M_error)){
##	if ((all(dim(M_error)==kYX))&&(matrixcalc::is.symmetric.matrix(M_error))&&(matrixcalc::is.positive.semi.definite(M_error))){
	if ((all(dim(M_error)==kYX))&&(.matrixcalc_is.symmetric.matrix(M_error))&&(.matrixcalc_is.positive.semi.definite(M_error))){
	    M_error<-sapply(1:n,function(i,x){x},x=M_error,simplify=FALSE)    
	    bM_errorOK<-TRUE
	}    
    }
    if ((is.list(M_error))&&(length(M_error)==n)){
	bM_errorOK<-TRUE
	i<-1
	while(bM_errorOK && (i<=n)){
	    if (is.vector(M_error[[i]],mode="numeric")){
		if (((length(M_error[[i]])==1)|| (length(M_error[[i]])==kYX)) && all(M_error[[i]]>=0)){
		    M_error[[i]]<-diag(M_error[[i]],nrow=kYX,ncol=kYX)
		}
		else{bM_errorOK<-FALSE}
	    }
    	    if (is.matrix(M_error[[i]])){
		if (!((all(dim(M_error[[i]])==kYX))&&(.matrixcalc_is.symmetric.matrix(M_error[[i]]))&&(.matrixcalc_is.positive.semi.definite(M_error[[i]])))){bM_errorOK<-FALSE}    
	    }    
	    i<-i+1
	}
    }
    if (!bM_errorOK){
	.my_stop(paste("Measurement error is provided incorrectly. The admissable types are: a single non-negative number (all traits and species have same measurement error), a non-negative vector of length ",kYX," (all species have the same diagonal measurement error matrix), a ",kYX," by ",kYX," symmetric semi-positive definite matrix (all species have the same measurement error matrix) or a list of length ",n, " (species specific measurement error, each list entry a single non-negative number, non-negative vector of length ", kYX," or ",kYX," by ",kYX," symmetric semi-positive definite matrix). The different entries of the list can be of different types (i.e. single number, vector or matrix). Notice that the measurement errors have to be independent between species!",sep=""),TRUE)
    }
    M_error
}
