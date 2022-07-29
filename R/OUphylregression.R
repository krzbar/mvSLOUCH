## This file is part of mvSLOUCH

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this package or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .


OU_phylreg<-function(mY,mD,phyltree,modelParams,regimes=NULL,kY=NULL,M.error=NULL,signif_level=0.05,regimes.times=NULL,root.regime=NULL,b_GLSB=FALSE,b_GLSX0=FALSE,signsB=NULL,signsvX0=NULL,estimate.root.state=FALSE){
    evolmodel<-.id_evolmodel(modelParams,b_throwerror=TRUE)
    .OU_phylreg_internal(mY=mY,mD=mD,phyltree=phyltree,modelParams=modelParams,evolmodel=evolmodel,regimes=regimes,kY=kY,M.error=M.error,signif_level=signif_level,regimes.times=regimes.times,root.regime=root.regime,b_GLSB=b_GLSB,b_GLSX0=b_GLSX0,signsB=signsB,signsvX0=signsvX0,estimate.root.state=estimate.root.state,b_calcmeanRSS=TRUE)
}

.OU_phylreg_internal<-function(mY,mD,phyltree,modelParams,evolmodel,regimes=NULL,kY=NULL,M.error=NULL,signif_level=0.05,regimes.times=NULL,root.regime=NULL,b_GLSB=FALSE,b_GLSX0=FALSE,signsB=NULL,signsvX0=NULL,estimate.root.state=FALSE,b_calcmeanRSS=TRUE){
    mY<-.check_input_trait_data(mY,n=phyltree$Ntips,vSpeciesLabels=phyltree$tip.label)
    mY_orig<-mY

    if (!all(is.na(mD))){
        if (mD[1]!="phylaverage"){
            if ((!is.matrix(mD))||(nrow(mD)!=(nrow(mY)*ncol(mY)))){
                .my_stop(paste0('Design matrix has to be either NA, "phylaverage" or be a matrix wth number of rows equalling the number of elements of mY: ',nrow(mY)*ncol(mY),'=',nrow(mY),'*',ncol(mY),'.'))
            }
        }
    }
    M.error<-.createMeasurementError(M.error,nrow(mY),ncol(mY))
    
## =====================================================================================
## Setup the regimes, pcmabase_model_box code from PhyloSDEestim.R
    kYX<-0
    kX<-NULL
    if (evolmodel=="bm"){kYX<-nrow(modelParams$Sxx)}
    if (evolmodel=="ouch"){kYX<-ncol(mY)}
    if (evolmodel=="mvslouch"){kX<-ncol(modelParams$Sxx);kYX<-ncol(modelParams$A)+kX}    
    tmpkX<-kYX;if (!is.null(kY)){tmpkX<-kYX-kY} ## tmpkX does not correspond to kX (!)
#    if (evolmodel=="bm"){kY<-nrow(modelParams$Sxx);kX<-kY}
    if (evolmodel=="ouch"){if(is.null(kY)){kY<-ncol(modelParams$A)}}
    if (evolmodel=="mvslouch"){if (is.null(kY)){kY<-ncol(modelParams$A)};}  

    regimesList<-.InitialRegimeSetup(phyltree=phyltree,regimes=regimes,regimes.times=regimes.times,mData=mY,kX=tmpkX,kYX=kYX,root.regime=root.regime,M.error=M.error)
    bOKregimes<-regimesList$bOKregimes
    if (!bOKregimes){.my_stop("The regimes, regimes.times and phyltree variables are not consistant. Cannot perform the phylogenetic regression. Please see manual on their format.",TRUE)}
    mY<-.check_input_trait_data(mY,n=phyltree$Ntips,vSpeciesLabels=phyltree$tip.label)
    modelParams$pcmbase_model_box<-regimesList$pcmbase_model_box ## this is the model object in which parameters for all regimes sit in for PCMBase
    modelParams$pcmbase_model_box<-.update_pcmbase_box_params(modelParams,evolmodel,vDo=c("H"=TRUE,"Theta"=TRUE,"Sigma_x"=TRUE,"X0"=TRUE))
    pcmbase_model_box<-modelParams$pcmbase_model_box
    
    regimes<-regimesList$regimes
    regimes.times<-regimesList$regimes.times
    regimes.types<-regimesList$regimes.types
    root.regime<-regimesList$root.regime
    regimes.types.orig<-regimesList$regimes.types.orig ## regimes names are in alphabetical order
    phyltree<-regimesList$phyltree
    if(estimate.root.state){if ((is.null(regimes))||(length(regimes.types.orig)==1)){estimate.root.state<-FALSE}}
## ===================================================================================== 
    
    b_interceptcalc<-FALSE
    b_estimateOUpars<-FALSE
    b_signspassed<-FALSE
    b_dzetacalc<-FALSE
    b_kappacalc<-FALSE
    if (all(is.na(mD))){
        b_estimateOUpars<-TRUE
	if ((is.null(kX))&&(!is.null(kY))){kX<-kYX-kY}
	## code from estimMAXLIK.R
	modelParams$regimeTimes<-regimes.times
	modelParams$regimes<-regimes
	modelParams$regimeTypes<-regimes.types                
	modelParams$regimes.types.orig<-regimes.types.orig ## this has to be passed as PCMBase will use original regime names while in GLS we just use the regime indices

        ## code taken from simulVasicekprocphyl.R
	params<-list()
	params$EvolModel<-evolmodel
	lPrecalculates<-list(vSpecies_times=NA)
        lPrecalculates$vSpecies_times<-phyltree$time.of.nodes[phyltree$tip_species_index] 
	EstimationParams<-.set.estimparams(params,kY,kX,length(regimes.types),estimate.root.state,cBestim_method="GLS") 
        designToEstim<-EstimationParams$designToEstim
	b_createtmpNAB<-FALSE
	if (evolmodel=="mvslouch"){
	    if (b_GLSB){
		designToEstim$B<-TRUE

		if (!is.null(signsB)){
		    if ((!is.matrix(signsB))||(ncol(signsB)!=kX)||(nrow(signsB)!=kY)){
			.my_stop(paste0("Wrong form of signsB, it should be a ",kY," by ",kX," matrix!"))
		    }
		    designToEstim$signsB<-signsB
		    b_signspassed<-TRUE
		}
		if ((is.matrix(modelParams$B))&&(ncol(modelParams$B)==kX)&&(nrow(modelParams$B)==kY)){
		    vNAs<-which(is.na(modelParams$B))
		    if(length(vNAs)>0){modelParams$B[vNAs]<-0;b_createtmpNAB<-TRUE;modelParams$pcmbase_model_box<-.update_pcmbase_box_params(modelParams,evolmodel,vDo=c("H"=TRUE,"Theta"=FALSE,"Sigma_x"=FALSE,"X0"=FALSE));pcmbase_model_box<-modelParams$pcmbase_model_box} ##starting conditions!
		}else{
		    modelParams$B<-matrix(0,nrow=kY,ncol=kX)
		    b_createtmpNAB<-TRUE
		    modelParams$pcmbase_model_box<-.update_pcmbase_box_params(modelParams,evolmodel,vDo=c("H"=TRUE,"Theta"=FALSE,"Sigma_x"=FALSE,"X0"=FALSE))
		    pcmbase_model_box<-modelParams$pcmbase_model_box
		}
	    }else{
		designToEstim$B<-FALSE   
		if ((!is.matrix(modelParams$B))||(ncol(modelParams$B)!=kX)||(nrow(modelParams$B)!=kY)||(any(is.na(modelParams$B)))||(any(!is.numeric(modelParams$B)))){
		    .my_stop(paste0("Since b_GLSB is FALSE, B should be a ",kY," by ",kX," matrix without any NAs!"))
		}
	    }	
	}
	if ((evolmodel=="mvslouch")||(evolmodel=="bm")){
	    if (b_GLSX0){
		designToEstim$X0<-TRUE
		if (!is.null(signsvX0)){
		    if ((!is.matrix(signsvX0))||(ncol(signsvX0)!=1)||(nrow(signsvX0)!=kX)){
			.my_stop(paste0("Wrong form of signsvX0, it should be a ",kX," by 1 matrix!"))
		    }
		    designToEstim$signsvX0<-signsvX0
		    b_signspassed<-TRUE
		}
		if ((is.matrix(modelParams$vX0))&&(ncol(modelParams$vX0)==1)&&(nrow(modelParams$vX0)==kX)){
		    vNAs<-which(is.na(modelParams$vX0))
		    if(length(vNAs)>0){modelParams$vX0[vNAs]<-0;modelParams$pcmbase_model_box<-.update_pcmbase_box_params(modelParams,evolmodel,vDo=c("H"=FALSE,"Theta"=FALSE,"Sigma_x"=FALSE,"X0"=TRUE));pcmbase_model_box<-modelParams$pcmbase_model_box} ##starting conditions!
		}else{
		    if (evolmodel=="mvslouch"){mX<-mY[,(kY+1):(kX+kY),drop=FALSE]}
		    else{mX<-mY}
		    modelParams$vX0<-matrix(c(apply(mX,2,mean)),nrow=kX,ncol=1)
		    modelParams$pcmbase_model_box<-.update_pcmbase_box_params(modelParams,evolmodel,vDo=c("H"=FALSE,"Theta"=FALSE,"Sigma_x"=FALSE,"X0"=TRUE))
		    pcmbase_model_box<-modelParams$pcmbase_model_box
		}
	    }else{
		designToEstim$X0<-FALSE   
		if ((!is.matrix(modelParams$vX0))||(ncol(modelParams$vX0)!=1)||(nrow(modelParams$vX0)!=kX)||(any(is.na(modelParams$vX0)))||(any(!is.numeric(modelParams$vX0)))){
		    .my_stop(paste0("Since b_GLSX0 is FALSE, vX0 should be a ",kX," by 1 matrix without any NAs!"))
		}
	    }		    
	}
    	if ((evolmodel=="mvslouch")||(evolmodel=="ouch")){
	    if (is.element("mPsi",names(modelParams))){
		if (any(is.na(modelParams$mPsi))||(any(!is.numeric(modelParams$mPsi)))){
		    designToEstim$psi<-TRUE
		    if(!all(is.na(modelParams$mPsi))){
			if((is.matrix(modelParams$mPsi))&&(nrow(modelParams$mPsi)==kY)&&(ncol(modelParams$mPsi)==length(regimes.types.orig))){
			    designToEstim$signsmPsi<-modelParams$mPsi
			    b_signspassed<-TRUE
			}else{
			    .my_stop(paste0("mPsi should be a ",kY," by ",length(regimes.types.orig)," matrix!"))
			}
		    }
		}else{
		    designToEstim$psi<-FALSE
		}    	
	    }else{
		designToEstim$psi<-TRUE
		modelParams$mPsi<-matrix(0,nrow=kY,ncol=length(regimes.types.orig))
		colnames(modelParams$mPsi)<-regimes.types.orig
	    }
	    if (is.element("vY0",names(modelParams))){
		    if (any(is.na(modelParams$vY0))||(any(!is.numeric(modelParams$vY0)))){
			designToEstim$y0<-TRUE
			if(!all(is.na(modelParams$vY0))){
			    if((is.matrix(modelParams$vY0)&&(nrow(modelParams$vY0)==kY)&&(ncol(modelParams$vY0)==1))){
				designToEstim$signsvY0<-modelParams$vY0
				b_signspassed<-TRUE
			    }else{
				.my_stop(paste0("vY0 should be a ",kY," by 1 matrix!"))
			    }
			}
		    }else{
			designToEstim$y0<-FALSE
		    }    	
	    }else{
	        designToEstim$y0<-TRUE
	        modelParams$vY0<-matrix(0,nrow=kY,ncol=1)
	    }	    	
	}
        if (designToEstim$y0AncState){
    	    if (!is.null(root.regime)){designToEstim$y0Regime<-which(regimes.types.orig==root.regime)}
    	    else{
    		if (is.element("y0Regime",names(designToEstim))){designToEstim$y0Regime<-which(regimes.types.orig==designToEstim$y0Regime)}        
        	else{designToEstim$y0Regime<-regimes[[1]][1]} ## need to choose something anyway ...     
    	    }
	}

	b_interceptcalc<-TRUE
	mXmX0<-NULL
	if (evolmodel=="mvslouch"){
    	    mX<-mY[,(kY+1):(kX+kY),drop=FALSE]    
    	    mY<-mY[,1:kY,drop=FALSE]
    	    b_dzetacalc<-(!designToEstim$YnonCondX)&&(designToEstim$B)
    	    b_kappacalc<-(!designToEstim$YnonCondX)&&(designToEstim$B)
    	    mXmX0<-t(mX) ## in previous version of mvSLOUCH it was transposed in this way!
	    # if X0 is unknown it is taken as the arithmetic average of the observations
    	    if (designToEstim$UseX0){
                mXmX0<-mXmX0-matrix(modelParams$vX0,ncol=nrow(mX),nrow=kX,byrow=FALSE) 
    	    }  
	}    

	modelParams$precalcMatrices<-.decompEigenA.S(modelParams,lPrecalculates,designToEstim,list(bCalcA=TRUE,bCovCalc=TRUE,dzetacalc=b_dzetacalc,lexptcalc=TRUE,kappacalc=b_kappacalc,interceptcalc=b_interceptcalc),mXmX0) 
	if (b_createtmpNAB){
	    modelParams_tmp<-modelParams
	    modelParams_tmp$B<-NA
	    modelParams$precalcMatrices[[4]]$lDzetaKappa<-.decompEigenA.S(modelParams_tmp,lPrecalculates,designToEstim,list(bCalcA=TRUE,bCovCalc=TRUE,dzetacalc=b_dzetacalc,lexptcalc=TRUE,kappacalc=b_kappacalc,interceptcalc=b_interceptcalc),mXmX0)[[4]]$lDzetaKappa
	}
	mD<-.design_matrix_construction(evolmodel,n=nrow(mY),model_params=modelParams,designToEstim=designToEstim,mX=mX)
## from phylgls.R
	 if (ncol(mD)>1){
            if (b_interceptcalc){
        	## if the user set all parameters to be estimated by regression as known, then there is nothing to do here
        	vIntercept<-NA
        	if ((length(modelParams$precalcMatrices)>=5)&&(is.element("intercept",names(modelParams$precalcMatrices[[5]])))){
            	    vIntercept<-modelParams$precalcMatrices[[5]]$intercept
        	}else{vIntercept<-rep(0,ncol(mY)*nrow(mY))} ## dummy operation for the moment if intercept unknown
        	vfullintercept<-.get_fullGLSintercept(evolmodel,vIntercept,mD[,1],designToEstim$YnonCondX,modelParams)
        	mY<-.correct_phylGLS_response_by_intercept(mY,vfullintercept)
        	mD<-mD[,-1,drop=FALSE]
        	glsmodel<-NA        
	    }    
            if (evolmodel=="mvslouch"){
                mY<-cbind(mY,matrix(0,nrow=nrow(mY),ncol=ncol(modelParams$Sxx)))
            }
        }else{
    	    .my_stop("Error in producing the design matrix, it has 0 columns!")
        }
    }
    if (mD[1]=="phylaverage"){mD<-rep(1,nrow(mY))%x%diag(ncol(mY))}
    lresult_regression<-.pcmbaseDphylOU_GLS(mY,mD,phyltree,pcmbase_model_box,glsmodel=NA)
    vGLSest<-lresult_regression$vGLSest
    minvDV1D<-lresult_regression$minvDV1D

    if (b_estimateOUpars) {
	modelParams<-.extract_GLS_results(evolmodel,lresult_regression,modelParams,designToEstim)
	pcmbase_model_box<-modelParams$pcmbase_model_box
	if (is.element("mPsi",names(modelParams))){
	    if (is.null(colnames(modelParams$mPsi))){
		colnames(modelParams$mPsi)<-modelParams$regimes.types.orig
	    }
	}
    }

    if ((b_estimateOUpars)&&(!b_signspassed)&&(!is.null(signif_level))) {
	regCIs<-.calcCI(modelParams,designToEstim,signif_level,minvDV1D)

    }else{
	regCIs<-matrix(NA,ncol=3,nrow=length(vGLSest))
	regCIs[,2]<-vGLSest
	vRegCIs<-qnorm(1-(1-signif_level)/2)*sqrt(diag(minvDV1D))
	regCIs[,1]<-vGLSest-vRegCIs
	regCIs[,3]<-vGLSest+vRegCIs
    }

    mYmmDb<-mY
    if (b_estimateOUpars){
	pcmbase_model_box_mean0<-modelParams$pcmbase_model_box
    }else{
        pcmbase_model_box_mean0<-.set_mean0_pcmbase_model_box(pcmbase_model_box,glsmodel=NA)
	mYmmDb<-mY-matrix(mD%*%vGLSest,ncol=ncol(mY),nrow=nrow(mY),byrow=TRUE)
    }
    RSS<- .pcmbaseDphylGaussian_RSS(mYmmDb,phyltree,pcmbase_model_box_mean0,glsmodel=NA) 
    RSS_phylaverage<-NA
    RSS_average<-NA
    if (b_calcmeanRSS){
	mDmean<-rep(1,nrow(mY))%x%diag(ncol(mY))        
	vGLSmeanest<-.pcmbaseDphylOU_GLS(mY,mD=mDmean,phyltree,pcmbase_model_box,glsmodel=NA)$vGLSest
	mYmGLSmean<-mY-matrix(vGLSmeanest,nrow=nrow(mY),ncol=ncol(mY),byrow=TRUE)
	RSS_phylaverage<- .pcmbaseDphylGaussian_RSS(mYmGLSmean,phyltree,pcmbase_model_box_mean0,glsmodel=NA) 
	
	mYmmean<-mY
	for (i in 1:ncol(mYmmean)){
            mYmmean[,i]<-mY[,i]-mean(mY[,i])
	}
	RSS_average<- .pcmbaseDphylGaussian_RSS(mYmmean,phyltree,pcmbase_model_box_mean0,glsmodel=NA) 
    }    
    R2_phylaverage<-1-RSS/RSS_phylaverage
    R2_average<-1-RSS/RSS_average
    modelParams<-.cleanUpModelParams(modelParams)    
    modelParams$pcmbase_model_box<-NULL
    list(vGLSest=vGLSest, regression.covariance.matrix=minvDV1D,regression.confidence.intervals=regCIs,modelParams=modelParams,mD=mD,RSS=RSS,R2_average=R2_average,R2_phylaverage=R2_phylaverage,RSS_average=RSS_average,RSS_phylaverage=RSS_phylaverage,phyltree=phyltree)    
}

OU_RSS<-function(mY,phyltree,modelParams,M.error=NULL,do_centre=NA,regimes=NULL,regimes.times=NULL,root.regime=NULL){
    evolmodel<-.id_evolmodel(modelParams,b_throwerror=TRUE)
    mY<-.check_input_trait_data(mY,n=phyltree$Ntips,vSpeciesLabels=phyltree$tip.label)
    mY_orig<-mY
    M.error<-.createMeasurementError(M.error,nrow(mY),ncol(mY))
    
## =====================================================================================
## Setup the regimes, pcmabase_model_box code from PhyloSDEestim.R
    kYX<-ncol(mY)
    kY<-NULL
    tmpkX<-kYX;if (!is.null(kY)){tmpkX<-kYX-kY} ## tmpkX does not correspond to kX (!)
    kY<-ncol(mY)
    if (is.element("A",names(modelParams))){kY<-ncol(modelParams$A)}
    regimesList<-.InitialRegimeSetup(phyltree=phyltree,regimes=regimes,regimes.times=regimes.times,mData=mY,kX=tmpkX,kYX=kYX,root.regime=root.regime,M.error=M.error)
    bOKregimes<-regimesList$bOKregimes
    if (!bOKregimes){.my_stop("The regimes, regimes.times and phyltree variables are not consistant. Cannot perform the phylogenetic regression. Please see manual on their format.",TRUE)}
    mY<-.check_input_trait_data(mY,n=phyltree$Ntips,vSpeciesLabels=phyltree$tip.label)
    modelParams$pcmbase_model_box<-regimesList$pcmbase_model_box ## this is the model object in which parameters for all regimes sit in for PCMBase
    modelParams$pcmbase_model_box<-.update_pcmbase_box_params(modelParams,evolmodel,vDo=c("H"=TRUE,"Theta"=TRUE,"Sigma_x"=TRUE,"X0"=TRUE))
    pcmbase_model_box<-modelParams$pcmbase_model_box
    phyltree<-regimesList$phyltree    

    if (!is.na(do_centre)){
	if (do_centre=="average"){
	    mY<-mY-matrix(apply(mY,2,mean),nrow=nrow(mY),ncol=ncol(mY),byrow=TRUE)
	    pcmbase_model_box_mean0<-.set_mean0_pcmbase_model_box(pcmbase_model_box,glsmodel=NA)
	}else if (do_centre=="evolutionary_model"){
#	    if (do_centre!=evolmodel){
#		.my_stop('Wrong value of do_centre, should be either NA, "average", "phylaverage", "bm", "ouch" or "mvslouch" and the same as evolmodel.')
#    	    }
        ## centring by mean is done inside phylgls calcs using PCMBase's likelihood calculations!
        }else if (do_centre=="phylaverage"){
            mDmean<-rep(1,nrow(mY))%x%diag(ncol(mY))        
            vGLSmeanest<-.pcmbaseDphylOU_GLS(mY,mD=mDmean,phyltree,pcmbase_model_box,glsmodel=NA)$vGLSest
            pcmbase_model_box<-.set_mean0_pcmbase_model_box(pcmbase_model_box,glsmodel=NA)
            mY<-mY-matrix(vGLSmeanest,nrow=nrow(mY),ncol=ncol(mY),byrow=TRUE)   
	}else{
	    .my_stop('Wrong value of do_centre, should be either NA, "average", "phylaverage" or "evolutionary_model".')
	}
    }else{
	pcmbase_model_box<-.set_mean0_pcmbase_model_box(pcmbase_model_box,glsmodel=NA)
    }
    RSS<-.pcmbaseDphylGaussian_RSS(mY,phyltree,pcmbase_model_box,glsmodel=NA)
    list(RSS=RSS,phyltree=phyltree)
}



OU_xVz<-function(mX,mZ,phyltree,modelParams,M.error=NULL,do_centre=NA,regimes=NULL,regimes.times=NULL,root.regime=NULL){
    evolmodel<-.id_evolmodel(modelParams,b_throwerror=TRUE)
    mX-.check_input_trait_data(mX,n=phyltree$Ntips,vSpeciesLabels=phyltree$tip.label)
    mZ<-.check_input_trait_data(mZ,n=phyltree$Ntips,vSpeciesLabels=phyltree$tip.label)
    M.error<-.createMeasurementError(M.error,nrow(mX),ncol(mX))
    
## =====================================================================================
## Setup the regimes, pcmabase_model_box code from PhyloSDEestim.R
    kYX<-ncol(mX)
    kY<-NULL
    tmpkX<-kYX;if (!is.null(kY)){tmpkX<-kYX-kY} ## tmpkX does not correspond to kX (!)
    kY<-ncol(mX)
    if (is.element("A",names(modelParams))){kY<-ncol(modelParams$A)}

    regimesList<-.InitialRegimeSetup(phyltree=phyltree,regimes=regimes,regimes.times=regimes.times,mData=mX,kX=tmpkX,kYX=kYX,root.regime=root.regime,M.error=M.error)
    bOKregimes<-regimesList$bOKregimes
    if (!bOKregimes){.my_stop("The regimes, regimes.times and phyltree variables are not consistant. Cannot perform the phylogenetic regression. Please see manual on their format.",TRUE)}
    mX<-.check_input_trait_data(mX,n=phyltree$Ntips,vSpeciesLabels=phyltree$tip.label)
    mZ<-.check_input_trait_data(mZ,n=phyltree$Ntips,vSpeciesLabels=phyltree$tip.label)    
    modelParams$pcmbase_model_box<-regimesList$pcmbase_model_box ## this is the model object in which parameters for all regimes sit in for PCMBase
    modelParams$pcmbase_model_box<-.update_pcmbase_box_params(modelParams,evolmodel,vDo=c("H"=TRUE,"Theta"=TRUE,"Sigma_x"=TRUE,"X0"=TRUE))
    pcmbase_model_box<-modelParams$pcmbase_model_box
    phyltree<-regimesList$phyltree
    if (!is.na(do_centre)){
	if (do_centre=="average"){
    	    mX<-mX-matrix(apply(mX,2,mean),nrow=nrow(mX),ncol=ncol(mX),byrow=TRUE)
    	    mZ<-mZ-matrix(apply(mZ,2,mean),nrow=nrow(mZ),ncol=ncol(mZ),byrow=TRUE)
    	    pcmbase_model_box<-.set_mean0_pcmbase_model_box(pcmbase_model_box,glsmodel=NA)
	}else if (do_centre=="evolutionary_model"){
#	else if ((do_centre=="bm")||(do_centre=="ouch")||(do_centre=="mvslouch")){
#    	  if (do_centre!=evolmodel){
#                .my_stop('Wrong value of do_centre, should be either NA, "average", "phylaverage", "bm", "ouch" or "mvslouch" and the same as evolmodel.')
#            }
        ## centring by mean is done inside phylgls calcs using PCMBase's likelihood calculations!
        }else if (do_centre=="phylaverage"){
            mDmean<-rep(1,nrow(mX))%x%diag(ncol(mX))        
            pcmbase_model_box<-.set_mean0_pcmbase_model_box(pcmbase_model_box,glsmodel=NA)
            vGLSmeanest_X<-.pcmbaseDphylOU_GLS(mX,mD=mDmean,phyltree,pcmbase_model_box,glsmodel=NA)$vGLSest
            mX<-mX-matrix(vGLSmeanest_X,nrow=nrow(mX),ncol=ncol(mX),byrow=TRUE)         
            vGLSmeanest_Z<-.pcmbaseDphylOU_GLS(mZ,mD=mDmean,phyltree,pcmbase_model_box,glsmodel=NA)$vGLSest
            mZ<-mZ-matrix(vGLSmeanest_Z,nrow=nrow(mZ),ncol=ncol(mZ),byrow=TRUE)         
        }else{
            .my_stop('Wrong value of do_centre, should be either NA, "average", "phylaverage" or "evolutionary_model".')
        }
    }else{
	pcmbase_model_box<-.set_mean0_pcmbase_model_box(pcmbase_model_box,glsmodel=NA)
    }    
    XZVXZ<-.pcmbaseDphylGaussian_RSS(mX-mZ,phyltree,pcmbase_model_box,glsmodel=NA)
    XVX<-.pcmbaseDphylGaussian_RSS(mX,phyltree,pcmbase_model_box,glsmodel=NA)
    ZVZ<-.pcmbaseDphylGaussian_RSS(mZ,phyltree,pcmbase_model_box,glsmodel=NA)    
    list(xVz=(XVX+ZVZ-XZVXZ)/2,phyltree=phyltree)
}

.get_phyl_mean<-function(phyltree,modelParams,evolmodel,M.error=NULL,regimes=NULL,regimes.times=NULL,root.regime=NULL){
## this function is untested and unused, but kept for the sake of the code
## =====================================================================================
## Setup the regimes, pcmabase_model_box code from PhyloSDEestim.R
    kY<-NULL
    kYX<-0
    if (evolmodel=="bm"){kYX<-nrow(modelParams$vX0)}
    if (evolmodel=="ouch"){kY<-ncol(modelParams$A);kYX<-kY}
    if (evolmodel=="mvslouch"){kY<-ncol(modelParams$A);kYX<-kY+ncol(modelParams$B);}    
    tmpkX<-kYX;if (!is.null(kY)){tmpkX<-kYX-kY} ## tmpkX does not correspond to kX (!)
    M.error<-.createMeasurementError(M.error,length(phyltree$Ntips),kYX)
    mY<-matrix(NA,ncol=kYX,nrow=length(phyltree$Ntips))
    
    regimesList<-.InitialRegimeSetup(phyltree=phyltree,regimes=regimes,regimes.times=regimes.times,mData=mY,kX=tmpkX,kYX=kYX,root.regime=root.regime,M.error=M.error)
    bOKregimes<-regimesList$bOKregimes
    if (!bOKregimes){.my_stop("The regimes, regimes.times and phyltree variables are not consistant. Cannot perform the phylogenetic regression. Please see manual on their format.",TRUE)}
    modelParams$pcmbase_model_box<-regimesList$pcmbase_model_box ## this is the model object in which parameters for all regimes sit in for PCMBase
    modelParams$pcmbase_model_box<-.update_pcmbase_box_params(modelParams,evolmodel,vDo=c("H"=TRUE,"Theta"=TRUE,"Sigma_x"=TRUE,"X0"=TRUE))
    pcmbase_model_box<-modelParams$pcmbase_model_box
    phyltree<-regimesList$phyltree    

    regimes<-regimesList$regimes
    regimes.times<-regimesList$regimes.times
    regimes.types<-regimesList$regimes.types
    root.regime<-regimesList$root.regime
    regimes.types.orig<-regimesList$regimes.types.orig ## regimes names are in alphabetical order

## ===================================================================================== 
    
    ## code from estimMAXLIK.R
    modelParams$regimeTimes<-regimes.times
    modelParams$regimes<-regimes
    modelParams$regimeTypes<-regimes.types                
    modelParams$regimes.types.orig<-regimes.types.orig ## this has to be passed as PCMBase will use original regime names while in GLS we just use the regime indices

    ## code taken from simulVasicekprocphyl.R
    lPrecalculates<-list(vSpecies_times=NA)
    lPrecalculates$vSpecies_times<-phyltree$time.of.nodes[phyltree$tip_species_index] 
    
    designToEstim<-list(B=FALSE,YnonCondX=TRUE) ## just to have a FALSE condition inside .decompEigenA.S but should not be necessary at all, i.e. designToEstim could be NULL
    modelParams$precalcMatrices<-.decompEigenA.S(modelParams,lPrecalculates,designToEstim,list(bCalcA=TRUE,bCovCalc=FALSE,dzetacalc=FALSE,lexptcalc=TRUE,kappacalc=FALSE,interceptcalc=FALSE),NULL) 
    modelParams$mPsi0<-NULL
    vPhylMean<-.calc.phyl.mean(lPrecalculates$vSpecies_times,evolmodel,modelParams)
    ## vPhylMean is filled in species by species
    matrix(vPhylMean,nrow=phyltree$Ntips,ncol=length(vPhylMean)/phyltree$Ntips,byrow=TRUE)	
}

.id_evolmodel<-function(modelParams,b_throwerror=TRUE){
    evolmodel<-"error"
    if (setequal(c("vX0","Sxx"),names(modelParams))){evolmodel<-"bm"}
    else{
	if ((setequal(c("vY0","Syy","mPsi","A"),names(modelParams)))||(setequal(c("vY0","Syy","mPsi","A","mPsi0"),names(modelParams)))){evolmodel<-"ouch"}
	    else{
		if ((setequal(c("vY0","Syy","mPsi","A","vX0","Sxx","B"),names(modelParams)))||(setequal(c("vY0","Syy","mPsi","A","mPsi0","vX0","Sxx","B"),names(modelParams)))||(setequal(c("vY0","Syy","mPsi","A","vX0","Sxx","B","Syx","Sxy"),names(modelParams)))||(setequal(c("vY0","Syy","mPsi","A","mPsi0","vX0","Sxx","B","Syx","Sxy"),names(modelParams)))){evolmodel<-"mvslouch"}
	    }
    }
    if ((evolmodel=="error")&&(b_throwerror)){
	.my_stop("Supplied model parameters are inconsistent with a BM, OUBM or OUOU model of evolution.")
    }
    evolmodel
}
