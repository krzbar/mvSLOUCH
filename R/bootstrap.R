## This file is part of mvSLOUCH

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this package or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .

parametric.bootstrap<-function(estimated.model,phyltree,values.to.bootstrap=NULL,regimes=NULL,root.regime=NULL,M.error=NULL,predictors=NULL,kY=NULL,numboot=100,Atype=NULL,Syytype=NULL,diagA=NULL,parameter_signs=NULL,start_point_for_optim=NULL,parscale=NULL,min_bl=0.0003,maxiter=c(10,50,100),estimateBmethod="ML"){
    model.components<-.boot.extract.model.components(estimated.model,Atype,Syytype,diagA,parameter_signs=parameter_signs,start_point_for_optim=start_point_for_optim)
    evolmodel<-model.components$evolmodel
    Atype<-model.components$Atype
    Syytype<-model.components$Syytype
    modelParams<-model.components$modelParams
    diagA<-model.components$diagA
    kY<-model.components$kY
    parameter_signs<-model.components$parameter_signs
    start_point_for_optim<-model.components$start_point_for_optim
    parscale<-model.components$parscale
    
    
    res<-sapply(1:numboot,function(i,evolmodel,phyltree,modelParams,M.error,predictors,regimes,regimes.times,root.regime,kY,Atype,Syytype,diagA,values.to.bootstrap,parameter_signs,start_point_for_optim,parscale){
	.my_message(paste("Doing bootstrap iteration: ",i,". Start time: ",Sys.time(),sep=""),TRUE)
	if (!is.null(start_point_for_optim)){
	    orgnames<-names(start_point_for_optim)
	    start_point_for_optim<-sapply(start_point_for_optim,function(x){
		y<-matrix(jitter(c(x)),nrow=nrow(x),ncol=ncol(x)) ## R does everything by column, so will be consitent here
		colnames(y)<-colnames(x)
		rownames(y)<-rownames(x)
		y
	    },simplify=FALSE)
	    names(start_point_for_optim)<-orgnames
	}
	estres<-switch(evolmodel,
	    bm=.bm.sim.est(phyltree,modelParams,M.error,predictors,min_bl),
	    ouch=.ouou.sim.est(phyltree,modelParams,regimes,regimes.times,root.regime,Atype,Syytype,diagA,M.error,predictors,parameter_signs,start_point_for_optim,parscale,min_bl,maxiter[c(1,3)]),
	    slouch=.oubm.sim.est(phyltree,modelParams,regimes,regimes.times,root.regime,Atype,Syytype,diagA,M.error,predictors,kY,parameter_signs,start_point_for_optim,min_bl,maxiter,estimateBmethod=estimateBmethod),
	    mvslouch=.oubm.sim.est(phyltree,modelParams,regimes,regimes.times,root.regime,Atype,Syytype,diagA,M.error,predictors,kY,parameter_signs,start_point_for_optim,parscale,min_bl,maxiter,estimateBmethod=estimateBmethod)
	)
	.boot.extract(estres,values.to.bootstrap)
    },evolmodel=evolmodel,phyltree=phyltree,modelParams=modelParams,M.error=M.error,predictors=predictors,regimes=regimes,regimes.times=NULL,root.regime=root.regime,kY=kY,Atype=Atype,Syytype=Syytype,diagA=diagA,values.to.bootstrap=values.to.bootstrap,parameter_signs=parameter_signs,start_point_for_optim=start_point_for_optim,parscale=parscale,simplify=FALSE)
    
    if (!is.null(values.to.bootstrap)){
	res<-list(paramatric.bootstrap.estimation.replicates=res,bootstrapped.parameters=NULL)
	res$bootstrapped.parameters<-sapply(values.to.bootstrap,function(x,res){sapply(res,function(bootex,x){bootex$bootstrapped.values[[which(names(bootex$bootstrapped.values)==x)]]},x=x,simplify=FALSE)},res=res$paramatric.bootstrap.estimation.replicates,simplify=FALSE)
	names(res$bootstrapped.parameters)<-values.to.bootstrap
    }
    res
}

.bm.sim.est<-function(phyltree,modelParams,M.error,predictors,min_bl=0.0003){
    simres<-simulBMProcPhylTree(phyltree,X0=modelParams$vX0,Sigma=modelParams$Sxx,dropInternal=TRUE,M.error=M.error)
    estres<-BrownianMotionModel(phyltree,mData=simres,predictors=predictors,M.error=M.error,min_bl=min_bl)
    if (!is.element("data",names(estres))){estres$data<-simres}
    estres
}

.ouou.sim.est<-function(phyltree,modelParams,regimes,regimes.times,root.regime,Atype,Syytype,diagA,M.error,predictors,parameter_signs,start_point_for_optim,parscale,min_bl=0.0003,maxiter=c(10,100)){    
    simres<-simulOUCHProcPhylTree(phyltree,modelParams=modelParams,regimes=regimes,regimes.times=regimes.times,dropInternal=TRUE,M.error=M.error)
    estres<-ouchModel(phyltree,mData=simres,regimes=regimes,regimes.times=regimes.times,root.regime=root.regime,predictors=predictors,M.error=M.error,Atype=Atype,Syytype=Syytype,diagA=diagA,parameter_signs=parameter_signs,start_point_for_optim=start_point_for_optim,parscale=parscale,min_bl=min_bl,maxiter=maxiter)
    if (!is.element("data",names(estres))){estres$data<-simres}
    estres

}

.oubm.sim.est<-function(phyltree,modelParams,regimes,regimes.times,root.regime,Atype,Syytype,diagA,M.error,predictors,kY,parameter_signs,start_point_for_optim,parscale,min_bl=0.0003,maxiter=c(10,50,100),estimateBmethod="ML"){    
    simres<-simulMVSLOUCHProcPhylTree(phyltree,modelParams=modelParams,regimes=regimes,regimes.times=regimes.times,dropInternal=TRUE, M.error=M.error)
    estres<-mvslouchModel(phyltree,mData=simres,kY=kY,regimes=regimes,regimes.times=regimes.times,root.regime=root.regime,predictors=predictors,M.error=M.error,Atype=Atype,Syytype=Syytype,diagA=diagA,parameter_signs=parameter_signs,start_point_for_optim=start_point_for_optim,parscale=parscale,min_bl=min_bl,maxiter=maxiter,estimateBmethod=estimateBmethod)
    if (!is.element("data",names(estres))){estres$data<-simres}
    estres

}

.boot.extract.model.components<-function(estimated.model,Atype=NULL,Syytype=NULL,diagA=NULL,modelParams=NULL,parameter_signs=NULL,start_point_for_optim=NULL,parscale=NULL){
    res<-list(evolmodel=NA,Atype=NA,Syytype=NA,modelParams=NA,diagA=NA,kY=NA)
    if ((is.element("BestModel",names(estimated.model)))&&(is.list(estimated.model$BestModel))){estimated.model<-estimated.model$BestModel}    
    else{
	if ((is.element("MaxLikFound",names(estimated.model)))&&(is.list(estimated.model$MaxLikFound))){estimated.model<-estimated.model$MaxLikFound}
	else{
	    if ((is.element("FinalFound",names(estimated.model)))&&(is.list(estimated.model$FinalFound))){estimated.model<-estimated.model$FinalFound}
	}
    }    
    if (is.null(modelParams)){modelParams<-.boot.getval("ParamsInModel",estimated.model)}
    if (is.na(modelParams[1])){.my_warning("Cannot extract estimated parameters and do bootstrapping!",TRUE,TRUE)}
    else{
	evolmodel<-.boot.getval("evolmodel",estimated.model)
	if (is.null(Atype)){Atype<-.boot.getval("Atype",estimated.model)}
	if (is.null(Syytype)){Syytype<-.boot.getval("Syytype",estimated.model)}
	if (is.null(diagA)){diagA<-.boot.getval("diagA",estimated.model)}
	kY<-.boot.getval("kY",estimated.model)
    
	if (is.na(evolmodel)){
    	    if (is.element("B",names(modelParams))){evolmodel<-"mvslouch"}
    	    else{
    		if (is.element("A",names(modelParams))){evolmodel<-"ouch"}
    		else{if (is.element("Sxx",names(modelParams))){evolmodel<-"bm"}}    		
    	    }	
	}
    
	if ((evolmodel=="ouch")||(evolmodel=="mvslouch")){
	## try to guess the structure if it is NOT provided!
	    if (is.na(Atype)){
		Atype<-"Invertible"
		if (all(modelParams$A[!diag(nrow(modelParams$A))] == 0)){
		    Atype<-"Diagonal"
		    if (length(unique(diag(modelParams$A)))==1){Atype<-"SingleValueDiagonal"}
		}else{
		    if (all(modelParams$A[lower.tri(modelParams$A)] == 0)){Atype<-"UpperTri"}
		    else{
			if (all(modelParams$A[upper.tri(modelParams$A)] == 0)){Atype<-"LowerTri"}
			else{
			    eigval<-eigen(modelParams$A)$values
			    if (all(Im(eigval)==0)){
				if (all(eigval>0)){Atype<-"DecomposablePositive"}
				else{
				    if (all(eigval<0)){Atype<-"DecomposableNegative"}
				    else{Atype<-"DecomposableReal"}
				}
			    }
			}
		    }
		}		
	    }
	    if (is.na(diagA)){
		diagA<-NULL
		if (all(diag(modelParams$A)>0)){diagA<-"Positive"}
		if (all(diag(modelParams$A)<0)){diagA<-"Negative"}
	    }
	    if (is.na(Syytype)){
		Syytype<-"UpperTri"
		if (all(modelParams$Syytype[!diag(nrow(modelParams$Syytype))] == 0)){
		    Syytypetype<-"Diagonal"
		    if (length(unique(diag(modelParams$Syy)))==1){Syytype<-"SingleValueDiagonal"}
		}
		else{if (all(modelParams$Syytype[upper.tri(modelParams$Syytype)] == 0)){Syytypetype<-"LowerTri"}}    	    
    	    }
	    if (is.na(kY)){
		if (evolmodel=="mvslouch"){kY<-nrow(modelParams$A)}
	    }
	    
	    if (is.null(parameter_signs)){
		parameter_signs<-.boot.getval("parameter_signs",estimated.model)
	        if ((length(parameter_signs)==0)||((!is.list(parameter_signs))&&(length(parameter_signs)==1)&&(is.na(parameter_signs[1])))){
		    parameter_signs<-NULL
		}
	    }
	    
	    if (is.null(parscale)){
		parscale<-.boot.getval("parscale",estimated.model)
	        if ((length(parscale)==0)||(is.na(parscale[1]))){
		    parscale<-NULL
		}
	    }
	    
	    if (is.null(start_point_for_optim)){
		start_point_for_optim<-.boot.getval("start_point_for_optim",estimated.model)
	    	if ((length(start_point_for_optim)==0)||(any(is.na(start_point_for_optim)))){
		    start_point_for_optim<-NULL
		}	    
	    }
	}
	res<-list(evolmodel=evolmodel,Atype=Atype,Syytype=Syytype,modelParams=modelParams,diagA=diagA,kY=kY,parameter_signs=parameter_signs,start_point_for_optim=start_point_for_optim,parscale=parscale)
    }
    res
}

.boot.extract<-function(estres,values.to.bootstrap){
    estres<-.boot.cleanup(estres)
    bootval<-NA
    if (!is.null(values.to.bootstrap)){
	bootval<-sapply(values.to.bootstrap,.boot.getval,estres=estres,simplify=FALSE)
	names(bootval)<-values.to.bootstrap
    }
    list(bootstrapped.values=bootval,estimation.results=estres)
}

.boot.cleanup<-function(estres){
    if (is.element("MaxLikFound",names(estres))){
	if (is.list(estres$MaxLikFound)){estres<-estres$MaxLikFound}
	else{estres<-estres$FinalFound}
    }
    estres
}

.boot.getval<-function(stattoget,estres){
    res<-NA
    if (is.element(stattoget,names(estres))){res<-estres[[which(names(estres)==stattoget)]]}
    else{
	if (is.list(estres)){
	    lextracted<-sapply(estres,function(elestres,stattoget){res<-NA;if(is.list(elestres)){res<-.boot.getval(stattoget,elestres)};res},stattoget=stattoget,simplify=FALSE)
	    are.extracted.na<-sapply(lextracted,function(x){(is.null(x)||is.na(x[1]))},simplify=TRUE)
	    if (!all(are.extracted.na)){
		## get the first element - there should be only 1 
		resindex<-which(!are.extracted.na)[1] 
		res<-lextracted[[resindex]]
	    }
	}
    }

    res
}
