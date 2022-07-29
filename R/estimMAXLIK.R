## This file is part of mvSLOUCH

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this package or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .


.maxlik.estim<-function(mData,PhylTree,EvolModel,EstimationParams,regimeTimes=NULL,regimes=NULL,regimeTypes=NULL,method="igls",tol=0.0001,maxIter=50,bShouldPrint=FALSE,maxTries=10,minLogLik=-Inf,regimes.types.orig=regimes.types.orig){
    MaxLikEstim<-NA
    
    ## This is essentially a dummy object that might become useful in the future when B for the mvSLOUCH model is estimated by GLS
    ## now it is used to pass the times of species
    lPrecalculates<-list(vSpecies_times=NA)
    lPrecalculates$vSpecies_times<-PhylTree$time.of.nodes[PhylTree$tip_species_index] 
    
    if (is.null(EstimationParams$vVars)){EstimationParams$vVars<-NULL}
    if (is.null(EstimationParams$conditional)){EstimationParams$conditional<-FALSE}
    if (is.null(EstimationParams$Atype)){EstimationParams$Atype<-"DecomposableReal"}
    if (is.null(EstimationParams$Btype)){EstimationParams$Btype<-"Any"}
    if (is.null(EstimationParams$mPsitype)){EstimationParams$mPsitype<-"Global"}
    if (is.null(EstimationParams$Syytype)){EstimationParams$Syytype<-"Symmetric"}
    if (is.null(EstimationParams$Sxxtype)){EstimationParams$Sxxtype<-"Symmetric"}
    

    ## .calculate.Tree.dists is deprecated with ape and pcmbase
    ## it might only be used or called when B is estimated by GLS (with full phylogenetic covariance and residuals) and not max lik
    
    modelParams<-vector("list",5)
    names(modelParams)<-c("regimeTimes","regimes","regimeTypes","pcmbase_model_regimes","regimes.types.orig")
    modelParams$regimeTimes<-regimeTimes
    modelParams$regimes<-regimes
    modelParams$regimeTypes<-regimeTypes		
    modelParams$pcmbase_model_box<-EstimationParams$pcmbase_model_box ## in EstimationParams we have a place holder without estimates!
    modelParams$regimes.types.orig<-regimes.types.orig ## this has to be passed as PCMBase will use original regime names while in GLS we just use the regime indices
    if (!is.null(EstimationParams$M_error)){
	modelParams$M_error<-EstimationParams$M_error
    }
    if ((method=="maxlik")&&(EvolModel=="bm")){
	MaxLikEstim<-vector("list",3)
	names(MaxLikEstim)<-c("BrownResult","ParamsInModel","ParamSummary")
	MaxLikEstim$BrownResult<-.bm.estim(mData,PhylTree,modelParams$pcmbase_model_box,modelParams$regimes.types.orig,minLogLik=minLogLik)
	MaxLikEstim$ParamsInModel<-list("Sxx"=MaxLikEstim$BrownResult$Sxx,"vX0"=MaxLikEstim$BrownResult$vX0,regressCovar=MaxLikEstim$BrownResult$regressCovar)
	MaxLikEstim$ParamSummary<-.Params.summary(PhylTree,MaxLikEstim$ParamsInModel,"bm",list(X0=TRUE),mData,NULL,MaxLikEstim$BrownResult$LogLik,NULL,MaxLikEstim$BrownResult$RSS,KnownParams=EstimationParams$KnownParams,conf.level=EstimationParams$conf.level,vVars=EstimationParams$vVars,conditional=EstimationParams$conditional,minLogLik=minLogLik,bfullCI=EstimationParams$calcCI)
##	MaxLikEstim$ParamSummary<-.Params.summary(PhylTree,MaxLikEstim$ParamsInModel,"bm",NULL,mData,NULL,MaxLikEstim$BrownResult$LogLik,NULL,MaxLikEstim$BrownResult$RSS,KnownParams=EstimationParams$KnownParams,conf.level=EstimationParams$conf.level,vVars=EstimationParams$vVars,conditional=EstimationParams$conditional,minLogLik=minLogLik,bfullCI=EstimationParams$calcCI)
    }    
    if (method=="glsgc"){
        MaxLikEstim<-.glsgc.estim(mData=mData,EvolModel=EvolModel,PhylTree=PhylTree,EstimationParams=EstimationParams,modelParams=modelParams,lPrecalculates=lPrecalculates,tol=tol,maxIter=maxIter,bShouldPrint=bShouldPrint,maxTries=maxTries,minLogLik=minLogLik)
    }
    MaxLikEstim
}

