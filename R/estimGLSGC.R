## This file is part of mvSLOUCH

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this package or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .

.glsgc.estim<-function(mData,EvolModel,PhylTree,EstimationParams,modelParams,lPrecalculates,tol=c(0.01,0.01),maxIter=c(50,50,100),bShouldPrint=FALSE,maxTries=15,minLogLik=-Inf){
    ## lPrecalculates should NEVER be NULL, we need the times of species to be passed to .decompEigenA.S()
    if ((is.null(lPrecalculates))||(!is.list(lPrecalculates))){
        lPrecalculates<-list(vSpecies_times=NA)
    }
    if (!is.element("vSpecies_times",names(lPrecalculates))){   
        lPrecalculates$vSpecies_times<-PhylTree$time.of.nodes[PhylTree$tip_species_index] 
    }

    mY<-mData[,1:EstimationParams$kY,drop=FALSE]
    .my_message("beginEstim",bShouldPrint)

    attempt2<-1
    doneOK<-FALSE

    while(!doneOK && (attempt2<=maxTries)){
	tryCatch({
	    ## Some error was generated, change starting parameters
	    if (attempt2>1){EstimationParams$StartPoint[1:length(EstimationParams$StartPoint)]<-rnorm(length(EstimationParams$StartPoint))}

	    EstimationParams<-.beginEstimationParams(EvolModel,EstimationParams,mData,PhylTree)

	    vEstim<-EstimationParams$StartPoint
	    vNames<-names(EstimationParams$StartPoint)
        
	    ## Calculate model for starting point
	    attempt<-1
	    bShouldTry<-TRUE
	    vEstim1<-vEstim
	    .my_message(paste("Starting point of heuristic search procedure : "),bShouldPrint)
	    .my_message(paste(vEstim,collapse=", "),bShouldPrint)
	    
	    while(bShouldTry&&(attempt<=maxTries)){
		.my_message("evalPoint \n",bShouldPrint)
		modelParams<-.par_transform_withPCMBase(vEstim,EstimationParams,EvolModel,modelParams)		
		tmpEvaluatedPoint<-.EvaluatePoint(EvolModel,PhylTree,mData,mY,modelParams,lPrecalculates,EstimationParams,tol[2],maxIter[2],bShouldPrint,TRUE,NULL,TRUE,FALSE,minLogLik=minLogLik)
		gc(verbose=bShouldPrint)
		modelParams<-tmpEvaluatedPoint$modelParams
		MaxLogLik<-tmpEvaluatedPoint$LogLik
		if (is.infinite(MaxLogLik)||is.nan(MaxLogLik)||(is.na(MaxLogLik))||(isTRUE(all.equal(MaxLogLik,minLogLik)))){
		    MaxLogLik<- minLogLik
		    vEstim[1:length(vEstim)]<-rnorm(length(vEstim))/10
		    .my_message("Point generated error, randomly changing it. New start point : ",bShouldPrint)
		    .my_message(paste(vEstim,collapse=", "),bShouldPrint)		    
		}
		else{
		    bShouldTry<-FALSE
		    MaxSearchPoint<-modelParams
		    MaxSearchPointParams<-vEstim
		    vEstim1<-vEstim*100
		    LogLik<-MaxLogLik
		}
		attempt<-attempt+1
	    }
	    doneOK<- TRUE
	},error=function(e){.my_message(paste("Restarting from new point, caught error:",e),FALSE)})
#	},error=function(e){.my_message(paste("Restarting from new point, caught error:",e),TRUE)})
	attempt2<-attempt2+1
    }
    if (attempt2 > maxTries){MaxLogLik<- NA}
    iter<-1
    if (is.infinite(MaxLogLik)||is.nan(MaxLogLik)||(is.na(MaxLogLik))||(isTRUE(all.equal(MaxLogLik,minLogLik)))){
	.my_message("Cannot start search from this point errors generated",TRUE)
	EstimationParams<-.beginEstimationParams(EvolModel,EstimationParams,mData,PhylTree)
	vEstim<-EstimationParams$StartPoint
	.my_message(paste(vEstim,collapse=", "),TRUE)
	iter<-maxIter[1]+1
	LogLik<- minLogLik
    }

    if (length(vEstim)==0){iter<-maxIter[1]+1} ## e.g. the user provided all parameters to be optimized over so there is no point in running the optimization
## run HeuristicSearch algorithm
    while((.calc.vec.dist(vEstim,vEstim1)>=tol[1])&&(iter<=maxIter[1])){
	## maximize loglikelihood over params
	vEstim1<-vEstim
	LogLik<- minLogLik
	if (EstimationParams$maximMethod=="optim"){## maximize likelihood by optim
	    tryCatch({
		.my_message("optim",bShouldPrint)
		c_optim_method<-"Nelder-Mead"
    		if (length(vEstim)==1){c_optim_method<-"BFGS"}
		max_optim_iter<-maxIter[3]
		OptimRes<-optim(
		    par=vEstim,
	    	    fn=.MinusPhylLogLikFunc,method=c_optim_method,phyltree=PhylTree,EvolModel=EvolModel,modelParams=modelParams,EstimationParams=EstimationParams,lPrecalculates=lPrecalculates,mData=mData,vNames=vNames,minLogLik=minLogLik,
	    	    control=list(maxit=max_optim_iter,parscale=EstimationParams$parscale)
		)
		gc(verbose=bShouldPrint)
		LogLik<-(-1)*OptimRes$value
		if (is.na(LogLik)||is.nan(LogLik)){LogLik<- minLogLik}
		vEstim<-OptimRes$par	   
		names(vEstim)<-vNames
		if (LogLik>MaxLogLik){
		    MaxLogLik<-LogLik
		    MaxSearchPoint<-modelParams
		    MaxSearchPointParams<-vEstim
		    MaxSearchPoint<-.EvaluatePoint(EvolModel,PhylTree,mData,mY,MaxSearchPoint,lPrecalculates,EstimationParams,NA,NA,NA,FALSE,lPoint=.par_transform_withPCMBase(vEstim,EstimationParams,EvolModel),TRUE,FALSE,minLogLik=minLogLik)$modelParams
		}
	    },error=function(e){.my_message(paste("Error in optim",e),FALSE);.my_message("\n",FALSE)})
	}
	.my_message(paste(EstimationParams$maximMethod,"found"),bShouldPrint)
	.my_message(paste(c(vEstim,"LogLik"=LogLik),collapse=", "),bShouldPrint);
	.my_message(paste("in iteration",iter,"of search.\n"),bShouldPrint)
	
	attempt<-1
	bShouldTry<-TRUE
	vEstimOrg<-vEstim
	while(bShouldTry&&(attempt<=maxTries)){
	    LogLik<- NA
	    tryCatch({
		modelParams<-.par_transform_withPCMBase(vEstim,EstimationParams,EvolModel,modelParams)
    		.my_message("evalPoint\n",bShouldPrint)
 		tmpEvaluatedPoint<-.EvaluatePoint(EvolModel,PhylTree,mData,mY,modelParams,lPrecalculates,EstimationParams,tol[2],maxIter[2],bShouldPrint,TRUE,NULL,TRUE,FALSE,minLogLik=minLogLik)
		gc(verbose=bShouldPrint)
		modelParams<-tmpEvaluatedPoint$modelParams    	    
		LogLik<-tmpEvaluatedPoint$LogLik
	    },error=function(e){.my_message(e,FALSE);.my_message("\n",FALSE)})	    
	    if (is.na(LogLik)||is.nan(LogLik)){LogLik<- minLogLik}
	    if (is.infinite(LogLik)||is.nan(LogLik)||(is.na(LogLik))||(isTRUE(all.equal(LogLik,minLogLik)))){
	        LogLik<- minLogLik
		vEstim<-jitter(vEstim)
		.my_message("Point generated error, randomly changing it. Point : ",bShouldPrint);
		.my_message(paste(vEstim,collapse=", "),bShouldPrint);
	    }
	    else{
		bShouldTry<-FALSE
		if (LogLik>MaxLogLik){
	    	    MaxLogLik<-LogLik
	    	    MaxSearchPoint<-modelParams
		    MaxSearchPointParams<-vEstim
		}
	    }
	    attempt<-attempt+1
	}
	if (is.infinite(LogLik)||is.nan(LogLik)||(is.na(LogLik))||(isTRUE(all.equal(LogLik,minLogLik)))){
	    attempt<-1
	    bShouldTry<-TRUE	
	    while(bShouldTry&&(attempt<=maxTries)){	    
	    	if (attempt==maxTries-1){vEstim[1:length(vEstim)]<-rnorm(length(vEstim))/10}else{vEstim<-jitter(MaxSearchPointParams)}		
		LogLik<-NA
		tryCatch({
		    modelParams<-.par_transform_withPCMBase(vEstim,EstimationParams,EvolModel,modelParams)
		    tmpEvaluatedPoint<-.EvaluatePoint(EvolModel,PhylTree,mData,mY,modelParams,lPrecalculates,EstimationParams,tol[2],maxIter[2],bShouldPrint,TRUE,NULL,TRUE,FALSE,minLogLik=minLogLik)
    		    gc(verbose=bShouldPrint)
    		    modelParams<-tmpEvaluatedPoint$modelParams
		    LogLik<-tmpEvaluatedPoint$LogLik
		},error=function(e){.my_message(e,FALSE);.my_message("\n",FALSE)})	    
		if (is.na(LogLik)||is.nan(LogLik)){LogLik<- minLogLik}
		if (is.infinite(LogLik)||is.nan(LogLik)||(is.na(LogLik))||(isTRUE(all.equal(LogLik,minLogLik)))){
		    .my_message("Point generated error, randomly changing it. Point : ",bShouldPrint)
		    .my_message(paste(vEstim,collapse=", "),bShouldPrint);
    		    LogLik<- minLogLik
		}
		else{
		    bShouldTry<-FALSE
		    .my_message("Restarting search from jittered best found estimate.",bShouldPrint)
		    if (LogLik>MaxLogLik){
	    		MaxLogLik<-LogLik
	    		MaxSearchPoint<-modelParams
			MaxSearchPointParams<-vEstim
		    }
	        }	    
	    attempt<-attempt+1
	    }	    
	}
	if (is.na(LogLik)||is.nan(LogLik)){LogLik<- minLogLik}
	if (is.infinite(LogLik)||is.nan(LogLik)||(is.na(LogLik))||(isTRUE(all.equal(LogLik,minLogLik)))){
	    .my_message("Cannot search from this point errors generated.",TRUE)
	    vEstim<-vEstimOrg
	    .my_message(vEstim,TRUE);
	    iter<-maxIter[1]+1
    	    LogLik<- minLogLik
	}		
	iter<-iter+1	
    }    
    if (bShouldPrint){
        if (.calc.vec.dist(vEstim,vEstim1)>=tol[1]){.my_message("Heuristic search procedure did not reach convergence at point.",bShouldPrint)}
        else{.my_message("Search algorithm converged.",bShouldPrint)}
    }
    
## Evaluate result of heuristic search algorithm    
    MaxLikHeuristicSearch<-vector("list",2)
    names(MaxLikHeuristicSearch)<-c("FinalFound","MaxLikFound")
    MaxLikHeuristicSearch$MaxLikFound<-"Same as final found"
    HeuristicSearchFinalFind<-vector("list",4)
    names(HeuristicSearchFinalFind)<-c("HeuristicSearchPointFinalFind","ParamsInModel","ParamSummary","LogLik")
    HeuristicSearchFinalFind$HeuristicSearchPointFinalFind<-c(vEstim,LogLik)
    names(HeuristicSearchFinalFind$HeuristicSearchPointFinalFind)<-c(names(EstimationParams$StartPoint),"LogLik")
    .my_message("Found estimate : ",bShouldPrint);.my_message(HeuristicSearchFinalFind$HeuristicSearchPointFinalFind,bShouldPrint);.my_message("\n",bShouldPrint)
    HeuristicSearchFinalFind$ParamsInModel<-.par_transform_withPCMBase(HeuristicSearchFinalFind$HeuristicSearchPointFinalFind[-length(HeuristicSearchFinalFind$HeuristicSearchPointFinalFind)],EstimationParams,EvolModel,modelParams)
    tree_height<-1
    if ((is.element("tree_height",names(PhylTree)))&&(!is.na(PhylTree$tree_height))){tree.height<-PhylTree$tree_height}
    tryCatch({	
	tmpEvaluatedPoint<-.EvaluatePoint(EvolModel,PhylTree,mData,mY,HeuristicSearchFinalFind$ParamsInModel,lPrecalculates,EstimationParams,tol[2],maxIter[2],bShouldPrint,TRUE,NULL,TRUE,FALSE,minLogLik=minLogLik)
	
	HeuristicSearchFinalFind$ParamsInModel<-tmpEvaluatedPoint$modelParams	
	HeuristicSearchFinalFind$LogLik<-tmpEvaluatedPoint$LogLik    	
	vVars2<-NULL
        n<-nrow(mData)                
        if ((EvolModel=="mvslouch")||(EvolModel=="ouch")||(EvolModel=="slouch")){vVars2<-c()}
	if ((is.null(EstimationParams$predictors)&&((EvolModel=="mvslouch")||(EvolModel=="ouch")||(EvolModel=="slouch")))){
	    if ((EvolModel=="mvslouch")||((EvolModel=="ouch"))){vVars2<-1:EstimationParams$kY}
	}else{
	    if (!is.null(EstimationParams$predictors)){
		NumVar<-length(mData)/n
	        vVars2<-setdiff(1:NumVar,EstimationParams$predictors)
	    }
        }
        orgEvolModel<-EvolModel
        if ((EvolModel=="mvslouch") && .is_det0(tmpEvaluatedPoint$modelParams$A,NULL)){
    	    tmpEvaluatedPoint$modelParams<-.mvslouch_to_ouch_model(tmpEvaluatedPoint$modelParams)
    	    tmpEvaluatedPoint$modelParams$pcmbase_model_box<-.update_pcmbase_box_params(tmpEvaluatedPoint$modelParams,"ouch",vDo=c("H"=TRUE,"Theta"=TRUE,"Sigma_x"=FALSE,"X0"=TRUE))
    	    tmpEvaluatedPoint$modelParams$precalcMatrices<-.decompEigenA.S(tmpEvaluatedPoint$modelParams,lPrecalculates,EstimationParams$designToEstim,list(bCalcA=TRUE,bCovCalc=TRUE,dzetacalc=FALSE,lexptcalc=TRUE,kappacalc=FALSE,interceptcalc=FALSE),NULL)
    	    EvolModel<-"ouch"
        }
        
    	RSS<-.calc.phyl.LogLik.traits(PhylTree,mData,lPrecalculates=lPrecalculates,EvolModel,modelParams=tmpEvaluatedPoint$modelParams,vVars=EstimationParams$vVars,conditional=EstimationParams$conditional,TRUE,minLogLik=minLogLik,vVars2=vVars2)
	EvolModel<-orgEvolModel
	model_params_to_use<-HeuristicSearchFinalFind$ParamsInModel
	if ((EvolModel=="mvslouch") && .is_det0(model_params_to_use$A,NULL)){
	    EvolModel<-"mvslouchtoouchdetA0"
	}	
	tryCatch({
	    HeuristicSearchFinalFind$ParamSummary<-.Params.summary(PhylTree,model_params_to_use,EvolModel,EstimationParams$designToEstim,mData,tree_height,HeuristicSearchFinalFind$LogLik,length(vEstim),RSS,KnownParams=EstimationParams$KnownParams,conf.level=EstimationParams$conf.level,vVars=EstimationParams$vVars,conditional=EstimationParams$conditional,minLogLik=minLogLik,bfullCI=EstimationParams$calcCI)
	},error=function(e){.my_message(paste("Cannot summarize final found point",e),TRUE);.my_message("\n",TRUE)})
	
	EvolModel<-orgEvolModel
	HeuristicSearchFinalFind$ParamsInModel<-.cleanUpModelParams(HeuristicSearchFinalFind$ParamsInModel)
    },error=function(e){.my_message(paste("Cannot evaluate final found point",e),TRUE);.my_message("\n",TRUE)})
    MaxLikHeuristicSearch$FinalFound<-HeuristicSearchFinalFind
    if (LogLik<MaxLogLik){
	HeuristicSearchMaxFind<-vector("list",4)
	names(HeuristicSearchMaxFind)<-c("HeuristicSearchPointMaxLik","ParamsInModel","ParamSummary","LogLik")
	HeuristicSearchMaxFind$HeuristicSearchPointMaxLik<-c(MaxSearchPointParams,MaxLogLik)
	names(HeuristicSearchMaxFind$HeuristicSearchPointMaxLik)<-c(names(EstimationParams$StartPoint),"LogLik")
	.my_message("Maxmimum likelihood found estimate : ",bShouldPrint);.my_message(HeuristicSearchMaxFind$HeuristicSearchPointMaxLik,bShouldPrint);.my_message("\n",bShouldPrint)
        HeuristicSearchMaxFind$ParamsInModel<-MaxSearchPoint
	HeuristicSearchMaxFind$LogLik<-MaxLogLik
	
	orgEvolModel<-EvolModel
        model_params_to_use<-HeuristicSearchMaxFind$ParamsInModel
        if ((EvolModel=="mvslouch") && .is_det0(tmpEvaluatedPoint$modelParams$A,NULL)){
    	    model_params_to_use<-.mvslouch_to_ouch_model(model_params_to_use)
    	    model_params_to_use$pcmbase_model_box<-.update_pcmbase_box_params(model_params_to_use,"ouch",vDo=c("H"=TRUE,"Theta"=TRUE,"Sigma_x"=FALSE,"X0"=TRUE))
    	    model_params_to_use$precalcMatrices<-.decompEigenA.S(model_params_to_use,lPrecalculates,EstimationParams$designToEstim,list(bCalcA=TRUE,bCovCalc=TRUE,dzetacalc=FALSE,lexptcalc=TRUE,kappacalc=FALSE,interceptcalc=FALSE),NULL)
    	    EvolModel<-"ouch"
        }
	RSS<-.calc.phyl.LogLik.traits(PhylTree,mData,lPrecalculates=lPrecalculates,EvolModel,modelParams=model_params_to_use,vVars=EstimationParams$vVars,conditional=EstimationParams$conditional,TRUE,minLogLik=minLogLik,vVars2=vVars2)
	EvolModel<-orgEvolModel
	model_params_to_use<-HeuristicSearchMaxFind$ParamsInModel
	if ((EvolModel=="mvslouch") && .is_det0(model_params_to_use$A,NULL)){
	    EvolModel<-"mvslouchtoouchdetA0"
	}
	tryCatch({	
	    HeuristicSearchMaxFind$ParamSummary<-.Params.summary(PhylTree,model_params_to_use,EvolModel,EstimationParams$designToEstim,mData,tree_height,HeuristicSearchMaxFind$LogLik,length(vEstim),RSS,KnownParams=EstimationParams$KnownParams,conf.level=EstimationParams$conf.level,vVars=EstimationParams$vVars,conditional=EstimationParams$conditional,minLogLik=minLogLik,bfullCI=EstimationParams$calcCI)
	},error=function(e){.my_message(paste("Cannot summarize maximum likelihood found point",e),TRUE);.my_message("\n",TRUE)})

	EvolModel<-orgEvolModel
	HeuristicSearchMaxFind$ParamsInModel<-.cleanUpModelParams(HeuristicSearchMaxFind$ParamsInModel)
	MaxLikHeuristicSearch$MaxLikFound<-HeuristicSearchMaxFind
    }                        
    MaxLikHeuristicSearch
}
