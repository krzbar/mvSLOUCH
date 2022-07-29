## This file is part of mvSLOUCH

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this package or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .


.beginEstimationParams<-function(EvolModel,EstimationParams,mData,PhylTree){
## called in estimGLSGC.R
    if ((EvolModel=="mvslouch")||(EvolModel=="ouch")){
	if (EvolModel=="mvslouch"){
	    vBMvars<-(EstimationParams$kY+1):(EstimationParams$kY+EstimationParams$kX)
	    mX<-mData[,vBMvars,drop=FALSE]

	    bmEstim<-.bm.estim(mX,PhylTree,EstimationParams$pcmbase_model_box,NULL,FALSE,vBMvars)
    	    mVxx<-bmEstim$StS
    	    if (is.element("Sxy",names(EstimationParams$Fixed$Sxy))){
    		## at least at the start Sxy should not sit anywhere but in the fixed
    		mVxx<-mVxx+EstimationParams$Fixed$Sxy%*%t(EstimationParams$Fixed$Sxy)
    	    }

	    for (cpcmbase_reg in names(EstimationParams$pcmbase_model_box$Sigma_x[1,1,])){
    		## Brownian motion with drift can also be considered here
        	EstimationParams$pcmbase_model_box$Sigma_x[vBMvars,vBMvars,which(names(EstimationParams$pcmbase_model_box$Sigma_x[1,1,])==cpcmbase_reg)]<-.changeSigmatoSyy(mVxx,"UpperTri",NULL,NULL,FALSE)
	    }
	    EstimationParams$pcmbase_model_box<-.set_pcmbase_model_box_X0(EstimationParams$pcmbase_model_box,bmEstim$vX0[,1])

	    EstimationParams$modelParams$pcmbase_model_box<-EstimationParams$pcmbase_model_box
	    EstimationParams$Fixed$Sxx<-bmEstim$Sxx
	    EstimationParams$Fixed$vX0<-bmEstim$vX0


	    mXmX0<-t(mX) ## in previous version of mvSLOUCH it was transposed in this way!
	    if (EstimationParams$designToEstim$UseX0){
		mXmX0<-mXmX0-matrix(EstimationParams$Fixed$vX0,ncol=nrow(mX),nrow=length(EstimationParams$Fixed$vX0),byrow=FALSE) 
	    }  

	    EstimationParams$Data<-vector("list",2)
	    names(EstimationParams$Data)<-c("mX","mXmX0")
	    EstimationParams$Data$mX<-mX
	    EstimationParams$Data$mXmX0<-mXmX0	    
	}
	if (EvolModel=="ouch"){}
    }
    EstimationParams    
}

.EvaluatePoint<-function(EvolModel,phyltree,mData,mY,modelParams,lPrecalculates,EstimationParams,tol,maxIter,bShouldPrint,bFull,lPoint=NULL,calcLogLik=TRUE,bSummarizePoint=FALSE,minLogLik=-Inf,t=1){
## bShouldPrint parameter is not used at the moment, but will be for mvSLOUCH GLS
## mData can be passed as NULL this just means that one should not calcualte the loglikelihood and return NA for it

    LogLik<- minLogLik
    if (((EvolModel=="mvslouch")||(EvolModel=="ouch"))&&(!bFull)&&(!is.null(lPoint))){
        modelParams$A<-lPoint$A
        modelParams$Syy<-lPoint$Syy
    }
    if ((EvolModel=="mvslouch")||(EvolModel=="ouch")){
	## lPrecalculates should never be NULL    
	## it is used to pass the times of species
	if (is.null(lPrecalculates)||(!is.list(lPrecalculates))){
    	    lPrecalculates<-list(vSpecies_times=NA)
	}
	if (!is.element("vSpecies_times",names(lPrecalculates))){   
    	    lPrecalculates$vSpecies_times<-phyltree$time.of.nodes[phyltree$tip_species_index] 
	}
    }
    if (EvolModel=="mvslouch"){
    	    tmpvY0<-NA
    	    if (!EstimationParams$designToEstim$y0){tmpvY0<-EstimationParams$Fixed$vY0}
    	    vAncPsi<-NA
    	    vX0<-NA
	    mPsi<-NA
	    mPsi0<-NA
    	    if (!EstimationParams$designToEstim$psi){mPsi<-EstimationParams$Fixed$mPsi}
	    if (is.element("mPsi0",names(EstimationParams$Fixed))){mPsi0<-EstimationParams$Fixed$mPsi0}
	    if (EstimationParams$designToEstim$y0AncState){
		if (!EstimationParams$designToEstim$psi){
		    vAncPsi<-EstimationParams$Fixed$mPsi[,EstimationParams$designToEstim$y0Regime]
		    if (!EstimationParams$designToEstim$psi0 && is.element("mPsi0",names(EstimationParams$Fixed))){vAncPsi<-vAncPsi+EstimationParams$Fixed$mPsi0}
		}
		if (!EstimationParams$designToEstim$y0OnlyFixed && !EstimationParams$designToEstim$B && EstimationParams$designToEstim$UseX0){vX0<-EstimationParams$Fixed$vX0;y0OnlyFixed<-FALSE}
	    }
    	    bDzeta<-FALSE;bKappa<-FALSE ## these are not present in non-conditional calculations
	    if (EstimationParams$designToEstim$B && !EstimationParams$designToEstim$YnonCondX){bDzeta<-TRUE;bKappa<-TRUE}
    	    modelParams$precalcMatrices<-.decompEigenA.S(modelParams,lPrecalculates,EstimationParams$designToEstim,list(bCalcA=TRUE,bCovCalc=TRUE,dzetacalc=bDzeta,lexptcalc=TRUE,kappacalc=bKappa,interceptcalc=TRUE),EstimationParams$Data$mXmX0)
    }
    if ((EvolModel=="ouch")||(EvolModel=="bm")){bDzeta<-FALSE;bKappa<-FALSE}
    if (EvolModel=="ouch"){
	    tmpvY0<-NA
	    if (!EstimationParams$designToEstim$y0){tmpvY0<-EstimationParams$Fixed$vY0}
    	    if (!EstimationParams$designToEstim$psi){mPsi<-EstimationParams$Fixed$mPsi}
    	    if (is.element("mPsi0",names(EstimationParams$Fixed))){mPsi0<-EstimationParams$Fixed$mPsi0}
    	    modelParams$precalcMatrices<-.decompEigenA.S(modelParams,lPrecalculates,EstimationParams$designToEstim,list(bCalcA=TRUE,bCovCalc=TRUE,dzetacalc=FALSE,lexptcalc=TRUE,kappacalc=FALSE,interceptcalc=TRUE),NULL)
    }

## bFull is internal sometimes we already have parameters calculated and we are only trying to obtain the summary
    if (bFull && ((EvolModel=="ouch")||(EvolModel=="mvslouch"))){
	## updated pcmbase model box is in modelparams
	modelParams<-.do_phylGLSestimation(EvolModel,phyltree,mY,EstimationParams$designToEstim,EstimationParams$conditional,bShouldPrint,modelParams,EstimationParams$Data$mX,maxIter=maxIter,tol=tol)
    }

    if ((calcLogLik)||(!is.null(mData))){
	orgEvolModel<-EvolModel
    	model_params_to_use<-modelParams
        if ((EvolModel=="mvslouch") && .is_det0(model_params_to_use$A,NULL)){
            model_params_to_use<-.mvslouch_to_ouch_model(model_params_to_use)
            model_params_to_use$pcmbase_model_box<-.update_pcmbase_box_params(model_params_to_use,"ouch",vDo=c("H"=TRUE,"Theta"=TRUE,"Sigma_x"=FALSE,"X0"=TRUE))
            model_params_to_use$precalcMatrices<-.decompEigenA.S(model_params_to_use,lPrecalculates,EstimationParams$designToEstim,list(bCalcA=TRUE,bCovCalc=TRUE,dzetacalc=FALSE,lexptcalc=TRUE,kappacalc=FALSE,interceptcalc=FALSE),NULL)
    	    EvolModel<-"ouch"
        }
	LogLik<-.calc.phyl.LogLik.traits(phyltree,mData,lPrecalculates=lPrecalculates,EvolModel=EvolModel,modelParams=model_params_to_use,vVars=EstimationParams$vVars,conditional=EstimationParams$conditional,RSS=FALSE,minLogLik=minLogLik)
	EvolModel<-orgEvolModel
    }
    else{if(calcLogLik){LogLik<-NA}}
    
    
    lPointSummary<-NULL
    if (bSummarizePoint){
        n<-phyltree$Ntips
	if ((EvolModel=="mvslouch")||(EvolModel=="ouch")){kY<-nrow(modelParams$A)}
        if ((EvolModel=="mvslouch")||(EvolModel=="bm")||(EvolModel=="slouch")){kX<-ncol(modelParams$Sxx)}
        if (EvolModel=="slouch"){kY<-1;}
	vVars2<-NULL
	if ((EvolModel=="mvslouch")||(EvolModel=="ouch")||(EvolModel=="slouch")){vVars2<-c()}
	if ((is.null(EstimationParams$predictors)&&((EvolModel=="mvslouch")||(EvolModel=="ouch")||(EvolModel=="slouch")))){
	    if ((EvolModel=="mvslouch")||(EvolModel=="ouch")){vVars2<-1:kY}
	}else{
	    if (!is.null(EstimationParams$predictors)){
		NumVar<-length(mData)/n
		vVars2<-setdiff(1:NumVar,EstimationParams$predictors)
	    }	
	}

	if(!is.null(mData)){
	    orgEvolModel<-EvolModel
    	    model_params_to_use<-modelParams
    	    if ((EvolModel=="mvslouch") && .is_det0(model_params_to_use$A,NULL)){
        	model_params_to_use<-.mvslouch_to_ouch_model(model_params_to_use)
        	model_params_to_use$pcmbase_model_box<-.update_pcmbase_box_params(model_params_to_use,"ouch",vDo=c("H"=TRUE,"Theta"=TRUE,"Sigma_x"=FALSE,"X0"=TRUE))
        	model_params_to_use$precalcMatrices<-.decompEigenA.S(model_params_to_use,lPrecalculates,EstimationParams$designToEstim,list(bCalcA=TRUE,bCovCalc=TRUE,dzetacalc=FALSE,lexptcalc=TRUE,kappacalc=FALSE,interceptcalc=FALSE),NULL)
    		EvolModel<-"ouch"
    	    }    	    
	    RSS<-"Error in calculating RSS in parameter summary."
	    tryCatch({
		RSS<-.calc.phyl.LogLik.traits(phyltree,mData,lPrecalculates=lPrecalculates,EvolModel=EvolModel,modelParams=model_params_to_use,vVars=EstimationParams$vVars,conditional=EstimationParams$conditional,RSS=TRUE,minLogLik=minLogLik,vVars2=vVars2)
	    },error=function(e){.my_message("Error in RSS calculation in summarizing evaluated point.",FALSE);.my_message(e,FALSE);.my_message("\n",FALSE)})
	    EvolModel<-orgEvolModel
	}
	else{RSS<-NA}

        modelParamsTmp<-modelParams
        modelParamsTmp$B<-NA
         if (EvolModel=="mvslouch"){
	    if (EstimationParams$designToEstim$B && !EstimationParams$designToEstim$YnonCondX){
		modelParams$precalcMatrices[[4]]$lDzetaKappa<-.decompEigenA.S(modelParamsTmp,lPrecalculates,EstimationParams$designToEstim,list(bCalcA=TRUE,bCovCalc=TRUE,dzetacalc=bDzeta,lexptcalc=TRUE,kappacalc=bKappa,interceptcalc=TRUE),EstimationParams$Data$mXmX0)[[4]]$lDzetaKappa
    	    }
        }
	if ((EvolModel=="ouch")||(EvolModel=="mvslouch")){ 
	    modelParams$regressCovar<-.do_phylGLSestimation(EvolModel,phyltree,mY,EstimationParams$designToEstim,EstimationParams$conditional,bShouldPrint,modelParams,EstimationParams$Data$mX,maxIter=maxIter,tol=tol)$minvDV1D
	}
	lPointSummary<-.Params.summary(phyltree,modelParams,EvolModel,EstimationParams$designToEstim,mData=mData,t=t,LogLik=LogLik,npar0=length(modelParams$vPoint),RSS=RSS,KnownParams=EstimationParams$KnownParams,conf.level=EstimationParams$conf.level,vVars=EstimationParams$vVars,conditional=EstimationParams$conditional,minLogLik=minLogLik,bfullCI=EstimationParams$calcCI)
    }
    list(modelParams=modelParams,LogLik=LogLik,PointSummary=lPointSummary)
}

