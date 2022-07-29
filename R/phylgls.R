## This file is part of mvSLOUCH

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this package or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .

.do_phylGLSestimation<-function(evolmodel,phyltree,mY,designToEstim,bConditional,bShouldPrint, model_params, mX,maxIter=50,tol=0.01){
## function called in modelparams.R
    iter<-1
    if ((!is.element("iRegLin",names(designToEstim)))||(!designToEstim$iRegLin)||(evolmodel!="mvslouch")||(!is.element("B",names(designToEstim)))){iter<-maxIter}
    diff_ests<-100    
    mYorg<-mY
    prevGLSest<-NA
    ## an NA in invA indicates problems, i.e. a bad point for the estimator, try to get out quickly
    ## however it could also be handled as a singular ouch model, PCMBase can handle these
    ## if(is.element("invA",model_params$precalcMatrices[[1]]) && (any(is.na(model_params$precalcMatrices[[1]]$invA)))){iter<-maxIter}
    while((iter<=maxIter)&&(diff_ests>tol)){
	mD<-.design_matrix_construction(evolmodel,n=nrow(mY),model_params=model_params,designToEstim=designToEstim,mX=mX)
	if (ncol(mD)>1){
	    ## if the user set all parameters to be estimated by regression as known, then there is nothing to do here
	    vIntercept<-NA
	    if ((length(model_params$precalcMatrices)>=5)&&(is.element("intercept",names(model_params$precalcMatrices[[5]])))){
		vIntercept<-model_params$precalcMatrices[[5]]$intercept
	    }else{vIntercept<-rep(0,ncol(mY)*nrow(mY))} ## dummy operation for the moment if intercept unknown
	    ##if (!designToEstim$YnonCondX){vDintercept<-mD[,1]}else{vDintercept<-mD[,1]}
	    vfullintercept<-.get_fullGLSintercept(evolmodel,vIntercept,mD[,1],designToEstim$YnonCondX,model_params)
	    mY<-.correct_phylGLS_response_by_intercept(mY,vfullintercept)
	    mD<-mD[,-1,drop=FALSE]
	    glsmodel<-NA
	    if ((bConditional)){
		if (bShouldPrint){.my_message("Conditional GLS estimation not implemented yet. Doing unconditional GLS. \n",TRUE)}
		glsmodel<-NA
	    }
	
	    if (evolmodel=="mvslouch"){
		mY<-cbind(mY,matrix(0,nrow=nrow(mY),ncol=ncol(model_params$Sxx)))
	    }
	    lGLSres<-.pcmbaseDphylOU_GLS(mY,mD,phyltree,model_params$pcmbase_model_box,glsmodel=glsmodel)
	    if (!is.null(lGLSres$vGLSest)){
		vNAsGLS<-which(is.na(lGLSres$vGLSest))
		if (length(vNAsGLS)>0){		    
		    if(length(prevGLSest)==length(lGLSres$vGLSest)){
			lGLSres$vGLSest[vNAsGLS]<-prevGLSest[vNAsGLS]
		    }
		    vNAsGLS<-which(is.na(lGLSres$vGLSest))
		    if (length(vNAsGLS)>0){lGLSres$vGLSest[vNAsGLS]<-0}
		    .my_warning("NAs produced in phylogenetic GLS estimation part of MLE. Setting parameters to 0.",FALSE)
    		}
    	    }
	    prevGLSest<-lGLSres$vGLSest
	    if ((is.element("B",names(designToEstim)))&&(is.element("B",names(model_params)))&&(designToEstim$B)){
		prevB<-model_params$B
	    }
	    model_params<-.extract_GLS_results(evolmodel,lGLSres,model_params,designToEstim)
	    if ((is.element("B",names(designToEstim)))&&(is.element("B",names(model_params)))&&(designToEstim$B)){
		if (!is.na(prevB[1])){diff_ests<-.calc.vec.dist(c(prevB),c(model_params$B))}
	    }else{diff_ests<-0}
	}
	else{
	    model_params<-.extract_from_signs(evolmodel,model_params,designToEstim)
	}
	mY<-mYorg
	iter<-iter+1
    }

    model_params
}


.correct_phylGLS_response_by_intercept<-function(mY,intercept){
## called in OUphylregression.R, phylgls.R
## intercept is already calculated in precalcs.R
## the intercept is stored as a vector for consitency with the old slow version of mvSLOUCH at the moment
    mY<-mY-matrix(intercept,nrow=nrow(mY),ncol=ncol(mY),byrow=TRUE)
    mY
}

.extract_GLS_results<-function(evolmodel,lGLSres,model_params,designToEstim){
    lGLSres$vGLSest[which(abs(lGLSres$vGLSest)<1e-15)]<-0
    lGLSres$minvDV1D[which(abs(lGLSres$minvDV1D)<1e-15)]<-0

    model_params=switch(evolmodel,
                bm=.extract_GLS_results_bm(lGLSres$vGLSest,model_params,designToEstim),
                ouch=.extract_GLS_results_ouch(lGLSres$vGLSest,model_params,designToEstim),
#               slouch=.extract_GLS_results_slouch(lGLSres$vGLSest,model_params,designToEstim),
                mvslouch=.extract_GLS_results_mvslouch(lGLSres$vGLSest,model_params,designToEstim)
            )
    model_params$regressCovar<-lGLSres$minvDV1D    
    model_params$regressCovar<-model_params$regressCovar[(1:model_params$num_extracted),(1:model_params$num_extracted),drop=FALSE]
    model_params$num_extracted<-NULL
    model_params
}


.get_fullGLSintercept<-function(evolmodel,vIntercept,mDintercept,YnonCondX,model_params){
    switch(evolmodel,
	bm=vIntercept+mDintercept,
	ouch=vIntercept+mDintercept,
#	slouch=.get_fullGLSintercept_slouch(vIntercept,mDintercept,YnonCondX,model_params)
	mvslouch=.get_fullGLSintercept_mvslouch(vIntercept,mDintercept,YnonCondX,model_params)
    )
}

.get_fullGLSintercept_mvslouch<-function(vIntercept,mDintercept,YnonCondX,model_params){
## called in OUphylregression.R, phylgls.R
    vInterfrommD<-rep(0,length(vIntercept))
# design matrix has the dummy estimate of X0-X0 regardless of YnonCondX or not
#    if (YnonCondX){vInterfrommD<-mDintercept}
#    else{
#	kY<-ncol(model_params$A)
#	kX<-ncol(model_params$Sxx)
#	n<-length(model_params$regimes)
#	for (i in 1:n){ ## we also estimate X0-X0 - corrected in intercept so the estimate here should be of the vector 0
#	    vInterfrommD[((i-1)*(kY)+1):((i-1)*(kY)+kY)]<- mDintercept[((i-1)*(kY+kX)+1):((i-1)*(kY+kX)+kY)]
#	}
#    }
    kY<-ncol(model_params$A)
    kX<-ncol(model_params$Sxx)
    n<-length(model_params$regimes)
    for (i in 1:n){ ## we also estimate X0-X0 - corrected in intercept so the estimate here should be of the vector 0
        vInterfrommD[((i-1)*(kY)+1):((i-1)*(kY)+kY)]<- mDintercept[((i-1)*(kY+kX)+1):((i-1)*(kY+kX)+kY)]
    }
    if (length(vInterfrommD)!=length(vIntercept)){.my_stop("Error in calculating intercept from known GLS parameters!")}
    vIntercept+vInterfrommD
}

.extract_GLS_results_bm<-function(vestim_params,model_params,designToEstim){
    vDoParamsUpdate<-c("H"=FALSE,"Theta"=FALSE,"Sigma_x"=FALSE,"X0"=FALSE)
    kX<-nrow(model_params$Sxx)

    CurrPos<-1
    ## get X0    
    if (designToEstim$X0){
	if ((is.element("signsvX0",names(designToEstim)))&&(length(which(!is.na(designToEstim$signsvX0)))>0)){
	    vToExtractvX0<-which(is.na(designToEstim$signsvX0))
	    mEstGLSvX0<-matrix(0,ncol=1,nrow=kX)
	    if (length(vToExtractvX0)<length(mEstGLSvX0)){
		if (length(vToExtractvX0)>0){
		    mEstGLSvX0[-vToExtractvX0]<-designToEstim$signsvX0[-vToExtractvX0]
		}else{
		    mEstGLSvX0<-designToEstim$signsvX0
		}		
	    }
	    if (length(vToExtractvX0)>0){
		mEstGLSvX0[vToExtractvX0]<-vestim_params[CurrPos:(CurrPos+length(vToExtractvX0)-1)]
	    }
	    updateCurrPos<-length(vToExtractvX0)
	}else{
	    mEstGLSvX0<-matrix(vestim_params[CurrPos:(CurrPos+kX-1)],nrow=kX,ncol=1)
	    updateCurrPos<-kX
	}
	model_params$vX0<-mEstGLSvX0;CurrPos<-CurrPos+updateCurrPos;vDoParamsUpdate["X0"]<-TRUE
    }
    model_params$pcmbase_model_box<-.update_pcmbase_box_params_bm(model_params,vDo=vDoParamsUpdate)
    model_params$num_extracted<-CurrPos-1
    model_params
}

.extract_GLS_results_ouch<-function(vestim_params,model_params,designToEstim){
## called in OUphylregression.R, phylgls.R
    vDoParamsUpdate<-c("H"=FALSE,"Theta"=FALSE,"Sigma_x"=FALSE,"X0"=FALSE)
    kY<-nrow(model_params$A)

    CurrPos<-1    

    if (designToEstim$y0 && !designToEstim$y0AncState){
	if ((is.element("signsvY0",names(designToEstim)))&&(length(which(!is.na(designToEstim$signsvY0)))>0)){
	    vToExtractvY0<-which(is.na(designToEstim$signsvY0))
	    mEstGLSvY0<-matrix(0,ncol=1,nrow=kY)
	    if (length(vToExtractvY0)<length(mEstGLSvY0)){
		if (length(vToExtractvY0)>0){
		    mEstGLSvY0[-vToExtractvY0]<-designToEstim$signsvY0[-vToExtractvY0]
		}else{
		    mEstGLSvY0<-designToEstim$signsvY0
		}
	    }
	    if (length(vToExtractvY0)>0){
		mEstGLSvY0[vToExtractvY0]<-vestim_params[1:length(vToExtractvY0)]
	    }
	    updateCurrPos<-length(vToExtractvY0)+1
	}else{
	    mEstGLSvY0<-matrix(vestim_params[1:kY],ncol=1)
	    updateCurrPos<-kY+1	    
	}
	model_params$vY0<-mEstGLSvY0;CurrPos<-updateCurrPos;if(updateCurrPos>1){vDoParamsUpdate["X0"]<-TRUE}
    }    

    if (designToEstim$psi){
	if ((is.element("signsmPsi",names(designToEstim)))&&(length(which(!is.na(designToEstim$signsmPsi)))>0)){
	    vToExtractmPsi<-which(is.na(designToEstim$signsmPsi))
	    mEstGLSmPsi<-matrix(0,nrow=kY,ncol=length(model_params$regimeTypes))
	    if (length(vToExtractmPsi)<length(mEstGLSmPsi)){
		if (length(vToExtractmPsi)>0){
		    mEstGLSmPsi[-vToExtractmPsi]<-designToEstim$signsmPsi[-vToExtractmPsi]
		}else{
		    mEstGLSmPsi<-designToEstim$signsmPsi
		}
	    }
	    if (length(vToExtractmPsi)>0){
		mEstGLSmPsi[vToExtractmPsi]<-vestim_params[CurrPos:(CurrPos+length(vToExtractmPsi)-1)]
	    }
	    updateCurrPos<-length(vToExtractmPsi)
	}else{
	    mEstGLSmPsi<-matrix(vestim_params[CurrPos:(CurrPos+length(model_params$regimeTypes)*kY-1)],nrow=kY,ncol=length(model_params$regimeTypes),byrow=FALSE)
	    updateCurrPos<-length(model_params$regimeTypes)*kY
	}
	model_params$mPsi<-mEstGLSmPsi;CurrPos<-CurrPos+updateCurrPos;if(updateCurrPos>0){vDoParamsUpdate["Theta"]<-TRUE}
    }
    if (designToEstim$psi0){
        if ((is.element("signsvmPsi0",names(designToEstim)))&&(length(which(!is.na(designToEstim$signsmPsi0)))>0)){
	    vToExtractmPsi0<-which(is.na(designToEstim$signsmPsi0))
	    mEstGLSmPsi0<-matrix(0,ncol=1,nrow=kY)
	    if (length(vToExtractmPsi0)<length(mEstGLSmPsi0)){
		if (length(vToExtractmPsi0)>0){
		    mEstGLSmPsi0[-vToExtractmPsi0]<-designToEstim$signsmPsi0[-vToExtractmPsi0]
		}else{
		    mEstGLSmPsi0<-designToEstim$signsmPsi0
		}
	    }
	    if (length(vToExtractmPsi0)>0){
		mEstGLSmPsi0[vToExtractmPsi0]<-vestim_params[CurrPos:(CurrPos+length(vToExtractmPsi0)-1)]
	    }
	    updateCurrPos<-length(vToExtractmPsi0)
	}else{
	    mEstGLSmPsi0<-matrix(vestim_params[CurrPos:(CurrPos+kY-1)],nrow=kY,ncol=1)
	    updateCurrPos<-kY
	}
	model_params$mPsi0<-mEstGLSmPsi0;CurrPos<-CurrPos+updateCurrPos
    }
    if (designToEstim$y0 && designToEstim$y0AncState){
	model_params$vY0<-matrix(model_params$mPsi[,designToEstim$y0Regime,drop=FALSE],ncol=1,nrow=kY)
	if (is.element("mPsi0",names(model_params))&&(!is.null(model_params$mPsi0))&&(!is.na(model_params$mPsi0[1]))){
	    model_params$vY0<-matrix(model_params$vY0+model_params$mPsi0,ncol=1,nrow=kY)
	}
	vDoParamsUpdate["X0"]<-TRUE
    }
	
    model_params$pcmbase_model_box<-.update_pcmbase_box_params_ouch(model_params,vDo=vDoParamsUpdate)
    model_params$num_extracted<-CurrPos-1
    model_params
}

.extract_GLS_results_mvslouch<-function(vestim_params,model_params,designToEstim){
## there is no need for information on conditionalYcT or not
## as information if B was estimated via GLS sits in designToEstim$B
## while at this point we do not care how the values were calculated we only want to extract them

    vDoParamsUpdate<-c("H"=FALSE,"Theta"=FALSE,"Sigma_x"=FALSE,"X0"=FALSE)
    kY<-nrow(model_params$A)
    kX<-nrow(model_params$Sxx)

    CurrPos<-1

    if (designToEstim$y0 && !designToEstim$y0AncState){
	if ((is.element("signsvY0",names(designToEstim)))&&(length(which(!is.na(designToEstim$signsvY0)))>0)){
	    vToExtractvY0<-which(is.na(designToEstim$signsvY0))
	    mEstGLSvY0<-matrix(0,ncol=1,nrow=kY)
	    if (length(vToExtractvY0)<length(mEstGLSvY0)){
		if (length(vToExtractvY0)>0){
		    mEstGLSvY0[-vToExtractvY0]<-designToEstim$signsvY0[-vToExtractvY0]
		}else{
		    mEstGLSvY0<-designToEstim$signsvY0
		}
	    }
	    if (length(vToExtractvY0)>0){
		mEstGLSvY0[vToExtractvY0]<-vestim_params[1:length(vToExtractvY0)]
	    }
	    updateCurrPos<-length(vToExtractvY0)+1
	}else{
	    mEstGLSvY0<-matrix(vestim_params[1:kY],ncol=1)
	    updateCurrPos<-kY+1	    
	}
	model_params$vY0<-mEstGLSvY0;CurrPos<-updateCurrPos;if(updateCurrPos>1){vDoParamsUpdate["X0"]<-TRUE}
    }
    
    ## we get y0 last as we might be mapping it back to optimal ancestral state
    if (designToEstim$psi){
	if ((is.element("signsmPsi",names(designToEstim)))&&(length(which(!is.na(designToEstim$signsmPsi)))>0)){
	    vToExtractmPsi<-which(is.na(designToEstim$signsmPsi))
	    mEstGLSmPsi<-matrix(0,nrow=kY,ncol=length(model_params$regimeTypes))
	    if (length(vToExtractmPsi)<length(mEstGLSmPsi)){
		if (length(vToExtractmPsi)>0){
		    mEstGLSmPsi[-vToExtractmPsi]<-designToEstim$signsmPsi[-vToExtractmPsi]
		}else{
		    mEstGLSmPsi<-designToEstim$signsmPsi
		}
	    }
	    if (length(vToExtractmPsi)>0){
		mEstGLSmPsi[vToExtractmPsi]<-vestim_params[CurrPos:(CurrPos+length(vToExtractmPsi)-1)]
	    }
	    updateCurrPos<-length(vToExtractmPsi)
	}else{
	    mEstGLSmPsi<-matrix(vestim_params[CurrPos:(CurrPos+length(model_params$regimeTypes)*kY-1)],nrow=kY,ncol=length(model_params$regimeTypes),byrow=FALSE)
	    updateCurrPos<-length(model_params$regimeTypes)*kY
	}
	model_params$mPsi<-mEstGLSmPsi;CurrPos<-CurrPos+updateCurrPos;if(updateCurrPos>0){vDoParamsUpdate["Theta"]<-TRUE}

    }
    ## get psi0
    if (designToEstim$psi0){
	if ((is.element("signsmPsi0",names(designToEstim)))&&(length(which(!is.na(designToEstim$signsmPsi0)))>0)){
	    vToExtractmPsi0<-which(is.na(designToEstim$signsmPsi0))
	    mEstGLSmPsi0<-matrix(0,ncol=1,nrow=kY)
	    if (length(vToExtractmPsi0)<length(mEstGLSmPsi0)){
		if (length(vToExtractmPsi0)>0){
		    mEstGLSmPsi0[-vToExtractmPsi0]<-designToEstim$signsmPsi0[-vToExtractmPsi0]
		}else{
		    mEstGLSmPsi0<-designToEstim$signsmPsi0
		}
	    }
	    if (length(vToExtractmPsi0)>0){
		mEstGLSmPsi0[vToExtractmPsi0]<-vestim_params[CurrPos:(CurrPos+length(vToExtractmPsi0)-1)]
	    }
	    updateCurrPos<-length(vToExtractmPsi0)
	}else{
	    mEstGLSmPsi0<-matrix(vestim_params[CurrPos:(CurrPos+kY-1)],nrow=kY,ncol=1)
	    updateCurrPos<-kY
	}
	model_params$mPsi0<-mEstGLSmPsi0;CurrPos<-CurrPos+updateCurrPos

    }
    ## get B
    if (designToEstim$B){
	if ((is.element("signsB",names(designToEstim)))&&(length(which(!is.na(designToEstim$signsB)))>0)){
	    vToExtractB<-which(is.na(designToEstim$signsB))
	    mEstGLSB<-matrix(0,nrow=kY,ncol=kX)
	    mTmpB<-t(mEstGLSB)
	    ## R fills in by default by column
	    if (length(vToExtractB)>0){
		mTmpB[vToExtractB]<-vestim_params[CurrPos:(CurrPos+length(vToExtractB)-1)]
	    }
	    mEstGLSB<-t(mTmpB)
	    if (length(vToExtractB)<length(mEstGLSB)){
		if (length(vToExtractB)>0){
		    mEstGLSB[-vToExtractB]<-designToEstim$signsB[-vToExtractB]
		}else{mEstGLSB<-designToEstim$signsB}
	    }
	    updateCurrPos<-length(vToExtractB)
	}else{
	    mEstGLSB<-matrix(vestim_params[CurrPos:(CurrPos+kX*kY-1)],nrow=kY,ncol=kX,byrow=TRUE)
	    updateCurrPos<-kX*kY
	}
	model_params$B<-mEstGLSB;CurrPos<-CurrPos+updateCurrPos;if(updateCurrPos>0){vDoParamsUpdate["H"]<-TRUE}
	model_params$precalcMatrices[[1]]$A1B<-model_params$precalcMatrices[[1]]$invA%*%model_params$B
    }
    
    ## get X0    
    if (designToEstim$X0){
	if ((is.element("signsvX0",names(designToEstim)))&&(length(which(!is.na(designToEstim$signsvX0)))>0)){
	    vToExtractvX0<-which(is.na(designToEstim$signsvX0))
	    mEstGLSvX0<-matrix(0,ncol=1,nrow=kX)
	    if (length(vToExtractvX0)<length(mEstGLSvX0)){
	        if (length(vToExtractvX0)>0){
	    	    mEstGLSvX0[-vToExtractvX0]<-designToEstim$signsvX0[-vToExtractvX0]
	        }else{
	    	    mEstGLSvX0<-designToEstim$signsvX0
	        }
	    }
	    if (length(vToExtractvX0)>0){
		mEstGLSvX0[vToExtractvX0]<-vestim_params[CurrPos:(CurrPos+length(vToExtractvX0)-1)]
	    }
	    updateCurrPos<-length(vToExtractvX0)
	}else{
	    mEstGLSvX0<-matrix(vestim_params[CurrPos:(CurrPos+kX-1)],nrow=kX,ncol=1)
	    updateCurrPos<-kX
	}
	model_params$vX0<-mEstGLSvX0;CurrPos<-CurrPos+updateCurrPos;vDoParamsUpdate["X0"]<-TRUE
    }
    ## get y0
    if (designToEstim$y0 && designToEstim$y0AncState){
	model_params$vY0<-matrix(model_params$mPsi[,designToEstim$y0Regime,drop=FALSE],ncol=1,nrow=kY)
	if (is.element("mPsi0",names(model_params))&&(!is.null(model_params$mPsi0))&&(!is.na(model_params$mPsi0[1]))){
	    model_params$vY0<-matrix(model_params$vY0+model_params$mPsi0,ncol=1,nrow=kY)
	}
        ## how to deal with this for det(A)=0?
        if ((!designToEstim$y0OnlyFixed)&&(!is.na(model_params$precalcMatrices[[1]]$A1B[1]))){model_params$vY0<-model_params$vY0-model_params$precalcMatrices[[1]]$A1B%*%model_params$vX0}
	model_params$vY0<-matrix(model_params$vY0,ncol=1)
	vDoParamsUpdate["X0"]<-TRUE
    }    
    model_params$pcmbase_model_box<-.update_pcmbase_box_params_mvslouch(model_params,vDo=vDoParamsUpdate)
    model_params$num_extracted<-CurrPos-1
    model_params
}



.extract_from_signs<-function(evolmodel,model_params,designToEstim){
## there is no need for information on conditionalYcT or not
## as information if B was estimated via GLS sits in designToEstim$B
## while at this point we do not care how the values were calculated we only want to extract them

    vDoParamsUpdate<-c("H"=FALSE,"Theta"=FALSE,"Sigma_x"=FALSE,"X0"=FALSE)
    if (evolmodel=="mvslouch" || evolmodel=="ouch"){kY<-nrow(model_params$A)}else{kY<-0}
    if (evolmodel=="mvslouch" || evolmodel=="bm"){kX<-nrow(model_params$Sxx)}else{kX<-0}

    CurrPos<-1
    if (evolmodel=="mvslouch" || evolmodel=="ouch"){
	if (designToEstim$y0 && !designToEstim$y0AncState){
	    if ((is.element("signsvY0",names(designToEstim)))&&(length(which(!is.na(designToEstim$signsvY0)))>0)){
		vToExtractvY0<-which(!is.na(designToEstim$signsvY0))
		mEstGLSvY0<-matrix(0,ncol=1,nrow=kY)
		if (length(vToExtractvY0)>0){
		    mEstGLSvY0[vToExtractvY0]<-designToEstim$signsvY0[vToExtractvY0]
		    updateCurrPos<-length(vToExtractvY0)+1
		    model_params$vY0<-mEstGLSvY0;if(updateCurrPos>1){vDoParamsUpdate["X0"]<-TRUE}
		}	    
	    }	
	}
    
    
	## we get y0 last as we might be mapping it back to optimal ancestral state
	if (designToEstim$psi){
	    if ((is.element("signsmPsi",names(designToEstim)))&&(length(which(!is.na(designToEstim$signsmPsi)))>0)){
		vToExtractmPsi<-which(!is.na(designToEstim$signsmPsi))
		mEstGLSmPsi<-matrix(0,nrow=kY,ncol=length(model_params$regimeTypes))
		if (length(vToExtractmPsi)>0){
		    mEstGLSmPsi[vToExtractmPsi]<-designToEstim$signsmPsi[vToExtractmPsi]
		    updateCurrPos<-length(vToExtractmPsi)
	    	    model_params$mPsi<-mEstGLSmPsi;if(updateCurrPos>0){vDoParamsUpdate["Theta"]<-TRUE}
		}	    
	    }
	}
	## get psi0
	if (designToEstim$psi0){
	    if ((is.element("signsmPsi0",names(designToEstim)))&&(length(which(!is.na(designToEstim$signsmPsi0)))>0)){
		vToExtractmPsi0<-which(!is.na(designToEstim$signsmPsi0))
	        mEstGLSmPsi0<-matrix(0,ncol=1,nrow=kY)
		if (length(vToExtractmPsi0)>0){
		    mEstGLSmPsi0[vToExtractmPsi0]<-designToEstim$signsmPsi0[vToExtractmPsi0]
	    	    updateCurrPos<-length(vToExtractmPsi0)
		    model_params$mPsi0<-mEstGLSmPsi0
		}
	    }
	}
    }
    if (evolmodel=="mvslouch"){
	## get B
	if (is.element("B",names(designToEstim)) && designToEstim$B){
	    if ((is.element("signsB",names(designToEstim)))&&(length(which(!is.na(designToEstim$signsB)))>0)){
		vToExtractB<-which(!is.na(designToEstim$signsB))
		mEstGLSB<-matrix(0,nrow=kY,ncol=kX)
		if (length(vToExtractB)>0){
		    mEstGLSB[vToExtractB]<-designToEstim$signsB[vToExtractB]
		    updateCurrPos<-length(vToExtractB)
		    model_params$B<-mEstGLSB;if(updateCurrPos>0){vDoParamsUpdate["H"]<-TRUE}
		    model_params$precalcMatrices[[1]]$A1B<-model_params$precalcMatrices[[1]]$invA%*%model_params$B
		}
	    }	
	}
    }
    
    if (evolmodel=="mvslouch" || evolmodel=="bm"){
	## get X0        
	if (is.element("X0",names(designToEstim)) &&designToEstim$X0){
	    if ((is.element("signsvX0",names(designToEstim)))&&(length(which(!is.na(designToEstim$signsvX0)))>0)){
		vToExtractvX0<-which(!is.na(designToEstim$signsvX0))
		mEstGLSvX0<-matrix(0,ncol=1,nrow=kX)
		if (length(vToExtractvX0)>0){
	    	    mEstGLSvX0[-vToExtractvX0]<-designToEstim$signsvX0[-vToExtractvX0]
	    	    updateCurrPos<-length(vToExtractvX0)
	    	    model_params$vX0<-mEstGLSvX0;vDoParamsUpdate["X0"]<-TRUE
		}	
	    }
	}
    }

    if (evolmodel=="mvslouch" || evolmodel=="ouch"){
	## get y0
	if (designToEstim$y0 && designToEstim$y0AncState){
	    model_params$vY0<-matrix(model_params$mPsi[,designToEstim$y0Regime,drop=FALSE],ncol=1,nrow=kY)
	    if (is.element("mPsi0",names(model_params))&&(!is.null(model_params$mPsi0))&&(!is.na(model_params$mPsi0[1]))){
		model_params$vY0<-matrix(model_params$vY0+model_params$mPsi0,ncol=1,nrow=kY)
	    }
    	    ## how to deal with this for det(A)=0?
    	    if ((is.element("y0OnlyFixed",designToEstim))&&(!designToEstim$y0OnlyFixed)&&(!is.na(model_params$precalcMatrices[[1]]$A1B[1]))){model_params$vY0<-model_params$vY0-model_params$precalcMatrices[[1]]$A1B%*%model_params$vX0}
	    model_params$vY0<-matrix(model_params$vY0,ncol=1)
    	    vDoParamsUpdate["X0"]<-TRUE
	}    
    }

    model_params$pcmbase_model_box=switch(evolmodel,
	bm=.update_pcmbase_box_params_bm(model_params,vDo=vDoParamsUpdate),
	ouch=.update_pcmbase_box_params_ouch(model_params,vDo=vDoParamsUpdate),
	mvslouch=.update_pcmbase_box_params_mvslouch(model_params,vDo=vDoParamsUpdate)
    )
    model_params
}


.design_matrix_construction<-function(evolmodel,n,model_params=NULL,designToEstim=NULL,kYX=NULL,mX=NULL){
## called in OUphylregression.R phylgls.R
    mD=switch(evolmodel,
                bm=.design_matrix_bm(kYX=kYX,n=n),
                ouch=.design_matrix_ouch(model_params,designToEstim,model_params$precalcMatrices[[3]]$lexpmtA,model_params$precalcMatrices[[3]]$lexptjA),
#               slouch=.design_matrix_slouch(),
                mvslouch=.design_matrix_mvslouch(model_params,designToEstim,model_params$precalcMatrices[[3]]$lexpmtA,model_params$precalcMatrices[[3]]$lexptjA,mX)
            )
    mD
}

.design_matrix_bm<-function(kYX,n){
## we have the  response as species on top of species
## hence we need copies of the Id matrix on top of each other
    mD<-rep(1,n)%x%diag(1,kYX,kYX)
    vinterceptcorrect<-rep(0,kYX*n)
    mD<-cbind(vinterceptcorrect,mD)
    mD
}

.set_mean0_pcmbase_model_box<-function(pcmbase_model_box,glsmodel=NA){
## called in OUphylregression.R, phylgls.R
## So far this is for OU (BM, OUOU, OUBM) type models

    pcmbase_model_box_mean0<-pcmbase_model_box
    bisTheta<-is.element("Theta",names(pcmbase_model_box_mean0))
    bisX0<-is.element("X0",names(pcmbase_model_box_mean0))
    for (cpcmbase_reg in names(pcmbase_model_box_mean0$Sigma_x[1,1,])){
        if (bisTheta){
    	    k<-length(pcmbase_model_box_mean0$Theta[,which(names(pcmbase_model_box_mean0$Theta[1,])==cpcmbase_reg)])
    	    pcmbase_model_box_mean0$Theta[,which(names(pcmbase_model_box_mean0$Theta[1,])==cpcmbase_reg)] <-rep(0,k)
    	}
    	if (bisX0){
    	    k<-length(pcmbase_model_box_mean0$X0)
    	    pcmbase_model_box_mean0<-.set_pcmbase_model_box_X0(pcmbase_model_box_mean0,rep(0,k))
	}
    }
## this is for mvslouch if we will need to do the conditional model
    if (!is.na(glsmodel)){class(pcmbase_model_box_mean0)<-glsmodel}
    pcmbase_model_box_mean0
}


.pcmbaseDphylOU_GLS<-function(mY,mD,phyltree,pcmbase_model_box,glsmodel=NA){
## called in OUphylregression.R phylgls.R
## intercept needs to be corrected with PRIOR to GLS
## glsmodel is for mvslouch if we will need to do the conditional model
## we cannot have a completely unobserved trait

## y = mD%*% beta + eps; 
## y = mY with rows stacked
## y = c(t(mY))
## beta=(solve(t(mD)%*%solve(V)%*%mD))%*%t(mD)%*%solve(V)%*%y
    kX<-ncol(mD) ## number of predictor variables
    kY<-ncol(mY) ## number of response variables
    vY<-c(t(mY)) ## the response vectorized
    ## mY is our data matrix, the rows are species, columns are traits
    ## c() stacks columns
    ## t(mY) makes each species be a column, then c(t(mY)) makes the rows stacked onto each other

    mD_org<-mD
    vNAys<-which(is.na(vY))
    if (length(vNAys)>0){
	## we have NA responses
	## we need to make NA those rows of mD that correspond to these NAs
	## so that it will be taken into account in RSS calculations, i.e. V rows/columns also removed
	mD[vNAys,]<-NA
    }

    pcmbase_model_box_mean0<-.set_mean0_pcmbase_model_box(pcmbase_model_box,glsmodel=glsmodel)
    mDV1D<-matrix(NA,kX,kX)

## first find t(mD)%*%solve(V)%*%mD
#    for (i in 1:kX){
#	vDi<-mD[,i]
#	mvD<-matrix(vDi,ncol=kY,byrow=TRUE)
#	mDV1D[i,i]<-.pcmbaseDphylGaussian_RSS(mvD,phyltree,pcmbase_model_box_mean0,glsmodel)		
#    }
    diag(mDV1D)<-sapply(1:kX,function(i,mD,phyltree,pcmbase_model_box_mean0,glsmodel){
	vDi<-mD[,i]
	mvD<-matrix(vDi,ncol=kY,byrow=TRUE)
	.pcmbaseDphylGaussian_RSS(mvD,phyltree,pcmbase_model_box_mean0,glsmodel)		    
    },mD=mD,phyltree=phyltree,pcmbase_model_box_mean0=pcmbase_model_box_mean0,glsmodel=glsmodel,simplify=TRUE)
    if (kX>1){## double for loop cannot be done if kX is 1
##	for (i in 1:(kX-1)){
##    	    for(j in (i+1):kX){ ## upper triangle
##		vDi<-mD[,i]
##		vDj<-mD[,j]
##		mvDij<-matrix(vDi-vDj,ncol=kY,byrow=TRUE)
##		rss_calc_ij<-.pcmbaseDphylGaussian_RSS(mvDij,phyltree,pcmbase_model_box_mean0,glsmodel)		
##		mDV1D[i,j]<-(mDV1D[i,i]+mDV1D[j,j]-rss_calc_ij)/2
##		mDV1D[j,i]<-mDV1D[i,j]
##	    }
##	}
	mIJ<-matrix(1:(kX*kX),kX,kX,byrow=FALSE)
	vUppTri<-mIJ[upper.tri(mIJ,diag=FALSE)]
	vvals<-sapply(vUppTri,function(ij,kX,mD,kY,phyltree,pcmbase_model_box_mean0,glsmodel){
	    i<-((ij-1)%%kX+1)
	    j<-((ij-1)%/%kX+1)
	    vDi<-mD[,i]
	    vDj<-mD[,j]
	    mvDij<-matrix(vDi-vDj,ncol=kY,byrow=TRUE)
	    rss_calc_ij<-.pcmbaseDphylGaussian_RSS(mvDij,phyltree,pcmbase_model_box_mean0,glsmodel)		
	    (mDV1D[i,i]+mDV1D[j,j]-rss_calc_ij)/2
	},kX=kX,mD=mD,kY=kY,phyltree=phyltree,pcmbase_model_box_mean0=pcmbase_model_box_mean0,glsmodel=glsmodel,simplify=TRUE)
	mDV1D[vUppTri]<-vvals
	mDV1D[lower.tri(mDV1D,diag=FALSE)]<-t(mDV1D)[lower.tri(mDV1D,diag=FALSE)]
    }
## we now have t(mD)%*%solve(V)%*%mD

## now we cacluate t(mD)%*%solve(V)%*%y
## it will be similar except that instead of mD we have y
    vDV1y<-rep(NA,kX)
    yV1y<-.pcmbaseDphylGaussian_RSS(mY,phyltree,pcmbase_model_box_mean0,glsmodel)    

    
##    for (i in 1:kX){
##        vDi<-mD[,i]	
##	mvDiY<-matrix(vDi-vY,ncol=kY,byrow=TRUE)
##	rss_calc_iY<-.pcmbaseDphylGaussian_RSS(mvDiY,phyltree,pcmbase_model_box_mean0,glsmodel)		
##	vDV1y[i]<-(mDV1D[i,i]+yV1y-rss_calc_iY)/2
##    }

    vDV1y<-sapply(1:kX,function(i,mD,kY,yV1y,vY,phyltree,pcmbase_model_box_mean0,glsmodel){
        vDi<-mD[,i]	
	mvDiY<-matrix(vDi-vY,ncol=kY,byrow=TRUE)
	rss_calc_iY<-.pcmbaseDphylGaussian_RSS(mvDiY,phyltree,pcmbase_model_box_mean0,glsmodel)		
	(mDV1D[i,i]+yV1y-rss_calc_iY)/2
    },mD=mD,kY=kY,yV1y=yV1y,vY=vY,phyltree=phyltree,pcmbase_model_box_mean0=pcmbase_model_box_mean0,glsmodel=glsmodel,simplify=TRUE)

## and finally calculate the resulting GLS estimate
    minvDV1D<-solve(mDV1D)
    vGLSest<-minvDV1D%*%vDV1y
## the resulting objects DO NOT have any column, row names. The order of the variables has to be used
    if (is.null(vGLSest)){
	numelGLS<-0
	if (is.matrix(minvDV1D)){
	    numelGLS<-nrow(minvDV1D)
	}
	if (numelGLS>0){
	    vGLSest<-rep(0,numelGLS)
	}else{
	    .my_stop("Cannot do phylogenetic GLS estimation part of MLE.",TRUE)
	}
    }
    list(vGLSest=vGLSest, minvDV1D=minvDV1D)
}


.pcmbaseDphylGaussian_RSS<-function(mX,phyltree,pcmbase_model_box,glsmodel=NA){
    ## precreated likelihood function cannot be used as we are using different data and different model setup, dimension
    loglik<-.callPCMBase_mvlik(mX,phyltree, pcmbase_model_box,b_islog=TRUE,minLogLik=-1e16,b_useLikFunc=FALSE)
    logdetV<-.detV(mX=mX,phyltree = phyltree, pcmbase_model_box = pcmbase_model_box,log=TRUE,glsmodel=glsmodel)
    numobsvalues<-length(which(!is.na(mX)))
    RSS<-(-2)*(loglik+0.5*logdetV+ numobsvalues*log(2*pi)/2)    
    RSS
}


.design_matrix_ouch<-function(modelParams,designToEstim,lexpmtA,lexptjA){
## order of columns in mX, has to be the same as order of species in phylogeny
    A<-modelParams$A
    kY<-nrow(A) ## number of traits
    n<-length(modelParams$regimes) ## number of species
    vinterceptcorrect<-rep(0,kY*n)

    iNumToEstim<-0
    if (designToEstim$y0 && !designToEstim$y0AncState){iNumToEstim<-iNumToEstim+kY} ## original trait
    if (designToEstim$psi){ iNumToEstim<-iNumToEstim+length(modelParams$regimeTypes)*kY }
    if (designToEstim$psi0){ iNumToEstim<-iNumToEstim+kY }
    mD<-matrix(0,nrow=n*kY,ncol=iNumToEstim)

    currXcol<-1
    
    ## setup Y0 column ---------------------------------------------------------------------------
    if (designToEstim$y0 && !designToEstim$y0AncState){ 
	for (i in 1:n){## for each species
	    mD[((i-1)*kY+1):(i*kY),1:kY]<-lexpmtA[[i]]##sapply(lexpmtA,function(x){x},simplify=TRUE) 
	}
	updateXcol<-kY
	if (is.element("signsvY0",names(designToEstim))){
	    vY0notGLS<-which(!is.na(designToEstim$signsvY0))
	    vtmpy0<-designToEstim$signsvY0
	    vtmpy0[which(is.na(vtmpy0))]<-0
	    vinterceptcorrect<-vinterceptcorrect+mD%*%vtmpy0
	    if (length(vY0notGLS)>0){
		mD<-mD[,-c(vY0notGLS+currXcol-1),drop=FALSE]
		updateXcol<-updateXcol-length(vY0notGLS);iNumToEstim<-iNumToEstim-length(vY0notGLS)
		if (ncol(mD)==0){mD<-matrix(0,nrow=n*kY,ncol=iNumToEstim)}
	    }
	}
	currXcol<-currXcol+updateXcol
    }
    #---------------------------------------------------------------------------------------------
    ## setup psi ---------------------------------------------------------------------------------
    ## regime names assumed to be numbers 1 .. max and according to this constructed
    if (designToEstim$psi){ 
	for (i in 1:n){## for each species
    	    mRegTmp<-matrix(0,ncol=length(modelParams$regimeTypes)*kY,nrow=kY)
    	    for (j in 1:length(modelParams$regimes[[i]])){
		jType<-modelParams$regimes[[i]][j]
		mRegTmp[,((jType-1)*kY+1):(jType*kY)]<-mRegTmp[,((jType-1)*kY+1):(jType*kY)]+(lexptjA[[i]][[j+1]]-lexptjA[[i]][[j]])
    	    }
    	    if (designToEstim$y0 && designToEstim$y0AncState){ 
		jType=designToEstim$y0Regime
		mRegTmp[,((jType-1)*kY+1):(jType*kY)]<-mRegTmp[,((jType-1)*kY+1):(jType*kY)]+lexpmtA[[i]]
	    }
    	    mD[((i-1)*kY+1):(i*kY),(currXcol):(currXcol+length(modelParams$regimeTypes)*kY-1)]<-mRegTmp
	}
	updateXcol<-length(modelParams$regimeTypes)*kY
	if (is.element("signsmPsi",names(designToEstim))){
	    vmPsinotGLS<-which(!is.na(designToEstim$signsmPsi))
	    vtmppsi<-designToEstim$signsmPsi
	    vtmppsi[which(is.na(vtmppsi))]<-0
	    vinterceptcorrect<-vinterceptcorrect+mD[,(currXcol):(currXcol+length(modelParams$regimeTypes)*kY-1),drop=FALSE]%*%vtmppsi
	    if (length(vmPsinotGLS)>0){
	    ## should work as all seems to be columnwise
		mD<-mD[,-c(vmPsinotGLS+currXcol-1),drop=FALSE]
		updateXcol<-updateXcol-length(vmPsinotGLS);iNumToEstim<-iNumToEstim-length(vmPsinotGLS)
		if (ncol(mD)==0){mD<-matrix(0,nrow=n*kY,ncol=iNumToEstim)}
	    }
	}
	currXcol<-currXcol+updateXcol
    }
    #---------------------------------------------------------------------------------------------    
    ## setup Psi0 column ---------------------------------------------------------------------------
    if (designToEstim$psi0){ 
	## mD[,currXcol:(currXcol+kY-1)]<-sapply(lexpmtA,function(x){diag(1,ncol(x),nrow(x))-x},simplify=TRUE)
	for (i in 1:n){## for each species
	    x<-lexpmtA[[i]]
	    mD[((i-1)*kY+1):(i*kY),currXcol:(currXcol+kY-1)]<-diag(1,ncol(x),nrow(x))-x
	}
	updateXcol<-kY
	if (is.element("signsmPsi0",names(designToEstim))){
	    mPsi0notGLS<-which(!is.na(designToEstim$signsmPsi0))
	    vtmppsi0<-designToEstim$signsmPsi0
	    vtmppsi0[which(is.na(vtmppsi0))]<-0
	    vinterceptcorrect<-vinterceptcorrect+mD[,currXcol:(currXcol+kY-1),drop=FALSE]%*%vtmppsi0
	    if (length(mPsi0notGLS)>0){
		mD<-mD[,-c(mPsi0notGLS+currXcol-1),drop=FALSE]
		updateXcol<-updateXcol-length(mPsi0notGLS);iNumToEstim<-iNumToEstim-length(mPsi0notGLS)
		if (ncol(mD)==0){mD<-matrix(0,nrow=n*kY,ncol=iNumToEstim)}
	    }
	}
	currXcol<-currXcol+updateXcol
    }
    #---------------------------------------------------------------------------------------------
    if (ncol(mD)>0){
	mD<-cbind(vinterceptcorrect,mD)	
    }
    else{mD<-matrix(vinterceptcorrect,ncol=1)}
    mD[which(abs(mD)<1e-15)]<-0
    mD
}

.calc.vec.dist<-function(v1,v2,method="euclid",params=list(p=2)){
    dist<-0
    if (method=="euclid"){dist<-(sum((v1-v2)^params$p))^(1/params$p)}
    dist
}


.design_matrix_mvslouch<-function(modelParams,designToEstim,lexpmtA,lexptjA,mX){
    mD<-NA
    mDtmp<-.design_matrix_ouch(modelParams,designToEstim,lexpmtA,lexptjA)

    kY<-nrow(modelParams$A)
    kX<-ncol(modelParams$Sxx)
    n<-length(modelParams$regimes) ## number of species
    vinterceptcorrect<-mDtmp[,1];mDtmp<-mDtmp[,-1,drop=FALSE]
    if (!designToEstim$YnonCondX){	
	## seutp B -----------------------------------------------------------------------------------
	if (designToEstim$B){
    	    invA<-modelParams$precalcMatrices[[1]]$invA
	    A<-modelParams$A
	    X0<-modelParams$vX0[,1]
	    mCovXX<-modelParams$precalcMatrices[[2]]$S22 
	    mCovTX<-modelParams$precalcMatrices[[2]]$S12 
	    currXcol<-ncol(mDtmp)+1
	    mDtmp<-cbind(mDtmp,matrix(0,ncol=kY*kX,nrow=n*kY))
	    lDzeta<-modelParams$precalcMatrices[[4]]$lDzeta
	    lDzetaKappa<-modelParams$precalcMatrices[[4]]$lDzetaKappa
	    iNumToEstim<-ncol(mDtmp)+kY*kX
	    updateXcol<-kY*kX
	    
    	    for (i in 1:n){## for each species
        	if (designToEstim$UseX0){mDtmp[((i-1)*kY+1):(i*kY),(currXcol):(currXcol+kY*kX-1)]<-lDzeta[[i]]%x%matrix(X0,nrow=1)}
        	if (!designToEstim$SimpReg){.my_message("Full phylogenetic regression for estimating B not implemented yet!",FALSE);designToEstim$SimpReg<-TRUE}
        	if (designToEstim$SimpReg){
            	    vToAdd<-rep(0,kX*kY)
            	    viXmX0<-mX[i,]-X0
            	    vNAX<-which(is.na(viXmX0))
            	    if (length(vNAX)>0){
                	if ((length(vNAX)<kX)&&designToEstim$BFullXNA){
                    	    viXmX0<-viXmX0[-vNAX]
                    	    tmpimDzetaKappa<-(lDzetaKappa[[i]]%*%solve(mCovXX))
                    	    if(nrow(tmpimDzetaKappa)>1){imDzetaKappa<-tmpimDzetaKappa[,-vNAX]}else{imDzetaKappa<-matrix(tmpimDzetaKappa[,-vNAX],nrow=1)}
                    	    vToAdd<-imDzetaKappa%x%matrix(viXmX0,nrow=1) 
                	}
            	    }else{vToAdd<-lDzetaKappa[[i]]%x%matrix(viXmX0,nrow=1)}                     	    
            	    mDtmp[((i-1)*kY+1):(i*kY),(currXcol):(currXcol+kY*kX-1)]<-mDtmp[((i-1)*kY+1):(i*kY),(currXcol):(currXcol+kY*kX-1)]+vToAdd
        	}
        	if (designToEstim$y0 && designToEstim$y0AncState && !designToEstim$y0OnlyFixed && designToEstim$UseX0){mDtmp[((i-1)*kY+1):(i*kY),(currXcol):(currXcol+kY*kX-1)]<-mDtmp[((i-1)*kY+1):(i*kY),(currXcol):(currXcol+kY*kX-1)]-(lexpmtA[[i]]%*%invA)%x%matrix(X0,nrow=1)}    
    	    }
    	    if (is.element("signsB",names(designToEstim))){
		vBnotGLS<-which(!is.na(designToEstim$signsB))
		vtmpB<-designToEstim$signsB
		vtmpB[which(is.na(vtmpB))]<-0
		vinterceptcorrect<-vinterceptcorrect+mDtmp[,(currXcol):(currXcol+kY*kX-1),drop=FALSE]%*%vtmpB
		if (length(vBnotGLS)>0){
		    mDtmp<-mDtmp[,-c(vBnotGLS+currXcol-1),drop=FALSE]
		    updateXcol<-updateXcol-length(vBnotGLS);iNumToEstim<-iNumToEstim-length(vBnotGLS)
		    if (ncol(mDtmp)==0){mDtmp<-matrix(0,nrow=n*kY,ncol=iNumToEstim)}
		}
	    }
	}
	#---------------------------------------------------------------------------------------------              
    }
    mD<-matrix(0,nrow=(kY+kX)*n,ncol=ncol(mDtmp)+kX)
    for (i in 1:n){ ## we also estimate X0-X0 - corrected in intercept so the estimate here should be of the vector 0
	if (ncol(mDtmp)>0){mD[((i-1)*(kY+kX)+1):((i-1)*(kY+kX)+kY),1:ncol(mDtmp)]<-mDtmp[((i-1)*kY+1):(i*kY),]}
	mD[((i-1)*(kY+kX)+kY+1):(i*(kY+kX)),(ncol(mDtmp)+1):(ncol(mDtmp)+kX)]<-diag(1,kX,kX)
    }
    if (ncol(mD)>0){mD<-cbind(c(vinterceptcorrect,rep(0,n*kX)),mD)}
    else{mD<-matrix(c(vinterceptcorrect,rep(0,n*kX)),ncol=1)}
    mD[which(abs(mD)<1e-15)]<-0
    mD
}


.detV<-function(mX, phyltree, pcmbase_model_box,log=TRUE,glsmodel=NA){
## called in getESS.R, phylgls.R
    resdet<-0    
    tryCatch({
	pcmbase_model_box_mean0<-.set_mean0_pcmbase_model_box(pcmbase_model_box,glsmodel=glsmodel)
	vNANAN<-which(is.na(c(mX)))
	if (length(vNANAN)>0){mX[setdiff(1:(ncol(mX)*nrow(mX)),vNANAN)]<-0}
	else{mX<-matrix(0,nrow=nrow(mX),ncol=ncol(mX))}

	## precreated likelihood function cannot be used as we are using data with all entries equalling 0, also could be for different dimensions of data
	loglik<-.callPCMBase_mvlik(mX,phyltree, pcmbase_model_box_mean0,b_islog=TRUE,minLogLik=-1e16,b_useLikFunc=FALSE)
	resdet<-(-2)*(((nrow(mX)*ncol(mX)-length(vNANAN))/2)*log(2*pi)+loglik)},error=function(e){.my_message(paste(".detV: Error in calculating between-species-between-traits variance-covariance matrix: ",e),FALSE)})
    if (!log){resdet<-exp(resdet)}

    resdet
}
