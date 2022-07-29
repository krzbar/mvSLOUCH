## This file is part of mvSLOUCH

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this package or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .


.MinusPhylLogLikFunc<-function(par,phyltree,EvolModel,modelParams,EstimationParams,lPrecalculates,mData,vNames,minLogLik=-Inf){
    LogLik<- (-1)*minLogLik ## optim minimizes!
    tryCatch({
	    if (EvolModel=="ouch"){
		LogLik<-.MinusPhylLogLikFuncouch(par,phyltree,modelParams,EstimationParams,lPrecalculates,mData,vNames,minLogLik)
	    }
	    if (EvolModel=="mvslouch"){
		LogLik<-.MinusPhylLogLikFuncMVslouch(par,phyltree,modelParams,EstimationParams,lPrecalculates,mData,vNames,minLogLik)
	    }
	},error=function(e){.my_message(paste("Caught:",e),FALSE);.my_message("\n",FALSE)})
    LogLik
}

.MinusPhylLogLikFuncMVslouch<-function(par,phyltree,modelParamsHeuristicSearch,EstimationParams,lPrecalculates,mData,vNames,minLogLik=-Inf){
## lPrecalculates should never be NULL    
## it is used to pass the times of species
    if (is.null(lPrecalculates)||(!is.list(lPrecalculates))){
	lPrecalculates<-list(vSpecies_times=NA)
    }
    if (!is.element("vSpecies_times",names(lPrecalculates))){	
	lPrecalculates$vSpecies_times<-phyltree$time.of.nodes[phyltree$tip_species_index] 
    }
    
    names(par)<-vNames
    modelParams<-.par_transform_withPCMBase(par,EstimationParams,"mvslouch",modelParamsHeuristicSearch)

    if (EstimationParams$designToEstim$B){modelParams$B<-modelParamsHeuristicSearch$B}
    modelParams$mPsi<-modelParamsHeuristicSearch$mPsi
    modelParams$mPsi0<-modelParamsHeuristicSearch$mPsi0
    modelParams$vY0<-modelParamsHeuristicSearch$vY0
    
    if (.is_det0(modelParams$A,NULL)){
        modelParams<-.mvslouch_to_ouch_model(modelParams)
        modelParams$pcmbase_model_box<-.update_pcmbase_box_params(modelParams,"ouch",vDo=c("H"=TRUE,"Theta"=TRUE,"Sigma_x"=FALSE,"X0"=TRUE))
	modelParams$precalcMatrices<-.decompEigenA.S(modelParams,lPrecalculates,EstimationParams$designToEstim,list(bCalcA=TRUE,bCovCalc=TRUE,dzetacalc=FALSE,lexptcalc=TRUE,kappacalc=FALSE,interceptcalc=FALSE),NULL)
	res<-(-1)*.calc.phyl.LogLik.traits(phyltree,mData,lPrecalculates=lPrecalculates,"ouch",modelParams=modelParams,vVars=EstimationParams$vVars,conditional=EstimationParams$conditional,FALSE,minLogLik=minLogLik)
    }else{
        modelParams$pcmbase_model_box<-.update_pcmbase_box_params(modelParams,"mvslouch",vDo=c("H"=TRUE,"Theta"=TRUE,"Sigma_x"=FALSE,"X0"=TRUE))
	modelParams$precalcMatrices<-.decompEigenA.S(modelParams,lPrecalculates,EstimationParams$designToEstim,list(bCalcA=TRUE,bCovCalc=TRUE,dzetacalc=FALSE,lexptcalc=TRUE,kappacalc=FALSE,interceptcalc=FALSE),NULL)
	res<-(-1)*.calc.phyl.LogLik.traits(phyltree,mData,lPrecalculates=lPrecalculates,"mvslouch",modelParams=modelParams,vVars=EstimationParams$vVars,conditional=EstimationParams$conditional,FALSE,minLogLik=minLogLik)
    }
    res
}

.MinusPhylLogLikFuncouch<-function(par,phyltree,modelParamsHeuristicSearch,EstimationParams,lPrecalculates,mData,vNames,minLogLik=-Inf){
## lPrecalculates should never be NULL    
## it is used to pass the times of species
    if (is.null(lPrecalculates)||(!is.list(lPrecalculates))){
	lPrecalculates<-list(vSpecies_times=NA)
    }
    if (!is.element("vSpecies_times",names(lPrecalculates))){	
	lPrecalculates$vSpecies_times<-phyltree$time.of.nodes[phyltree$tip_species_index] 
    }


    names(par)<-vNames
    modelParams<-.par_transform_withPCMBase(par,EstimationParams,"ouch",modelParamsHeuristicSearch)

    modelParams$mPsi<-modelParamsHeuristicSearch$mPsi
    modelParams$mPsi0<-modelParamsHeuristicSearch$mPsi0
    modelParams$vY0<-modelParamsHeuristicSearch$vY0
    modelParams$pcmbase_model_box<-.update_pcmbase_box_params(modelParams,"ouch",vDo=c("H"=FALSE,"Theta"=TRUE,"Sigma_x"=FALSE,"X0"=TRUE))

    modelParams$precalcMatrices<-.decompEigenA.S(modelParams,lPrecalculates,EstimationParams$designToEstim,list(bCalcA=TRUE,bCovCalc=TRUE,dzetacalc=FALSE,lexptcalc=TRUE,kappacalc=FALSE,interceptcalc=FALSE),NULL)

    (-1)*.calc.phyl.LogLik.traits(phyltree,mData,lPrecalculates=lPrecalculates,"ouch",modelParams=modelParams,vVars=EstimationParams$vVars,conditional=EstimationParams$conditional,FALSE,minLogLik=minLogLik)                
}

.calc.phyl.LogLik.traits<-function(phyltree,mData,lPrecalculates,EvolModel,modelParams,vVars=NULL,conditional=FALSE,RSS=FALSE,minLogLik=-Inf,vVars2=NULL){
## vVars are used if one wants to work with everything marginally or conditionally
## vVars2 are used if one just wants to calculate a conditional RSS
## be careful! numbers in vVars2 correspond to AFTER vVars selection was done!

    phylLogLik<-minLogLik  
    { 
	tryCatch({
		if (!(is.null(vVars))){
		    kYkX<-ncol(mData)
		    if ((max(vVars)>kYkX)||(min(vVars)<1)){.my_warning(paste("Provided variables out of the bound 1:",kYkX," correcting to this",sep=""),TRUE,FALSE)}
		    vVars<-intersect(unique(vVars),1:kYkX) ## we are not really interested in degenerate disitributions as same likelihood and singularity in dmvnorm
		    if (length(vVars)< kYkX){## no point if we are taking the same variables since we don't consider degenerate ones
			if (conditional){## conditional on the other variables distribution
			    .my_warning("WARNING: Maximizing over conditional likelihood is not implemented yet. Considering the full likelihood. \n",TRUE,FALSE)
			    vVars<-1:kYkX
			}else{## just the marginal distribution of the variables
			    vVarsRemove<-setdiff(1:kYkX,vVars)
			    mData[,vVarsRemove]<-NA ## in PCMBase if things are NA, then this is equivalent to marginalizing
			}
		    }
		}
		if (!RSS){
		    phylLogLik<-.callPCMBase_mvlik(mData,phyltree, modelParams$pcmbase_model_box,b_islog=TRUE,minLogLik=minLogLik)
		    if (phylLogLik<minLogLik){phylLogLik<-minLogLik}
		}else{
		    if (!is.null(vVars2)){
			n<-nrow(mData)
			kYkX<-ncol(mData)
			vVars_resp<-intersect(1:kYkX,unique(vVars2)) ## we are not really interested in degenerate disitributions as same likelihood and singularity in dmvnorm
			vVars_pred<-setdiff(1:kYkX,vVars_resp)
			kY<-length(vVars2)
			
			phylLogLik<-vector("list",3)
			names(phylLogLik)<-c("RSS","R2","R2_phylaverage")
			phylLogLik$RSS<-"Not calculated due to error"
			phylLogLik$R2<-"Not calculated due to error"			
			phylLogLik$R2_phylaverage<-"Not calculated due to error"			
			
			if (kY<kYkX){
			    ## currently we calculate a non-phylogenetic RSS, i.e. we do not take into account the correlation between the species, only between the traits but inside a species
			    phylLogLik<-c(phylLogLik,vector("list",4))
			    names(phylLogLik)<-c(names(phylLogLik)[1:3],"RSS_conditional_on_predictors","R2_conditional_on_predictors","RSS_non_phylogenetic_conditional_on_predictors","R2_non_phylogenetic_conditional_on_predictors")
			    phylLogLik$RSS_conditional_on_predictors<-"Not implemented yet"
			    phylLogLik$R2_conditional_on_predictors<-"Not implemented yet"
			    phylLogLik$RSS_non_phylogenetic_conditional_on_predictors<-"Not calculated due to error"
			    phylLogLik$R2_non_phylogenetic_conditional_on_predictors<-"Not calculated due to error"
			    
			    vIntercept_non_phyl<-apply(mData,2,mean,na.rm=TRUE)

			    vX0<-NULL;vY0<-NULL;mPsi<-NULL;A1B<-NULL;lexpmtA<-NULL;lexpmtjA<-NULL;regimes<-modelParams$regimes;mPsi0<-NULL
			    lAcalcs<-NULL;lScalcs<-NULL
			    vspecies_times<-phyltree$time.of.nodes[1:n]
			    
			    if (is.element("vX0",names(modelParams))){vX0<-modelParams$vX0}
			    if (is.element("vY0",names(modelParams))){vY0<-modelParams$vY0}
			    if (is.element("mPsi",names(modelParams))){mPsi<-modelParams$mPsi}
			    if (is.element("regimes",names(modelParams))){regimes<-modelParams$regimes}
			    if (is.element("precalcMatrices",names(modelParams))&&is.list(modelParams$precalcMatrices)&&(is.element("A1B",names(modelParams$precalcMatrices[[1]])))){A1B<-modelParams$precalcMatrices[[1]]$A1B}
			    if (is.element("precalcMatrices",names(modelParams))&&is.list(modelParams$precalcMatrices)&&(length(modelParams$precalcMatrices)>2)&&is.element("lexpmtA",names(modelParams$precalcMatrices[[3]]))){lexpmtA<-modelParams$precalcMatrices[[3]]$lexpmtA}
			    if (is.element("precalcMatrices",names(modelParams))&&is.list(modelParams$precalcMatrices)&&(length(modelParams$precalcMatrices)>2)&&is.element("lexpmtjA",names(modelParams$precalcMatrices[[3]]))){lexpmtjA<-modelParams$precalcMatrices[[3]]$lexpmtjA}
			    if (is.element("precalcMatrices",names(modelParams))&&is.list(modelParams$precalcMatrices)&&(length(modelParams$precalcMatrices)>1)){lAcalcs<-modelParams$precalcMatrices[[1]];lScalcs<-modelParams$precalcMatrices[[2]]}

			    if (is.element("mPsi0",names(modelParams))){mPsi0<-modelParams$mPsi0}
			    			    
			    tryCatch({
				mRSS_nonphyl_values<-t(sapply(1:n,function(i,vVars_resp,vVars_pred,EvolModel,mData,vIntercept_non_phyl,vX0,vY0,mPsi,A1B,lexpmtA,lexpmtjA,regimes,mPsi0){
				    RSS_non_phyl<-0;RSS_non_phyl_null<-0
				    s_warn<-0
				    vNA_resp<-which(is.na(mData[i,vVars_resp]));  if (length(vNA_resp)>0){vVars_resp<-vVars_resp[-vNA_resp]}
				    if (length(vVars_resp)>0){## we could have in principle a completely NA observation
					vMeani<-.calc_species_mean(EvolModel,lexpmtA[[i]],lexpmtjA[[i]],regimes[[i]],vX0=vX0,vY0=vY0,mPsi=mPsi,A1B=A1B,mPsi0=mPsi0)
					bis_non_degen<-FALSE
					bis_predNA<-FALSE
					mCovi<-.calc_species_covariance(EvolModel,vspecies_times[i],lAcalcs,lScalcs,method="plus.v",b_correct_forPD=FALSE)
					if (!.is0_Merror(modelParams$pcmbase_model_box)){
					    mCovi<-mCovi+modelParams$M_error[[i]]
					}
					if ((length(mCovi)>0) &&(!is.na(mCovi[1]))){
					    if ((!.matrixcalc_is.symmetric.matrix(mCovi)) || (!.matrixcalc_is.positive.definite(mCovi))){
        					s_warn<- -1 
        					bis_non_degen<-FALSE
    						tryCatch({mCovi<-as.matrix(Matrix::nearPD(mCovi)$mat);bis_non_degen<-TRUE},error=function(e){.my_message(e,FALSE)},warning=function(e){.my_message(e,FALSE)})
					    }else{bis_non_degen<-TRUE}				   
					    if (bis_non_degen){
						vNA_pred<-which(is.na(mData[i,vVars_pred]));  if (length(vNA_pred)>0){vVars_pred<-vVars_pred[-vNA_pred]}
						if (length(vVars_pred)>0){
						    vMeaniCond<-vMeani[vVars_resp]
						    mCoviCond<-mCovi[vVars_resp,vVars_resp,drop=FALSE]
						    if ((!.matrixcalc_is.symmetric.matrix(mCoviCond)) || (!.matrixcalc_is.positive.definite(mCoviCond))){
        				    		s_warn<- -1 
        				    		bis_non_degen<-FALSE
    							tryCatch({mCoviCond<-as.matrix(Matrix::nearPD(mCoviCond)$mat);bis_non_degen<-TRUE},error=function(e){.my_message(e,FALSE)})
						    }
						    if (bis_non_degen){
							tryCatch({bis_non_degen<-(!.is_det0(mCoviCond))&&(!.is_det0(mCovi[vVars_pred,vVars_pred,drop=FALSE]))},error=function(e){.my_message(paste("Problem with covariance of species ",i,". It has infinite or missing values hence we cannot calculate conditional on predictors RSS.",sep=""),TRUE)})
						    }
						}else{
						    bis_predNA<-TRUE
						}
					    }
					}
					if (bis_non_degen){
					    if (length(vVars_pred)>0){	
						b_no_error<-FALSE				    
						tryCatch({
						    mRegRespPred<-mCovi[vVars_resp,vVars_pred,drop=FALSE]%*%solve(mCovi[vVars_pred,vVars_pred,drop=FALSE])
						    vMeaniCond<-vMeaniCond+mRegRespPred%*%(mData[i,vVars_pred]-vMeani[vVars_pred])
						    mCoviCond<-mCoviCond-mRegRespPred%*%mCovi[vVars_pred,vVars_resp,drop=FALSE]
						    b_no_error<-TRUE
						},error=function(e){.my_message(e,FALSE)})
					    
						b_no_error2<-FALSE				    
						if ((!.matrixcalc_is.symmetric.matrix(mCoviCond)) || (!.matrixcalc_is.positive.definite(mCoviCond))){
        					    s_warn<- -1 
    						    tryCatch({mCoviCond<-as.matrix(Matrix::nearPD(mCoviCond)$mat);b_no_error2<-TRUE},error=function(e){.my_message(e,FALSE)})						
						}else{b_no_error2<-TRUE}
						if(b_no_error2){
						    b_no_error2<-FALSE
						    tryCatch({			
							invmCoviCond<-solve(mCoviCond)
							RSS_non_phyl<-t((mData[i,vVars_resp]-vMeaniCond))%*%invmCoviCond%*%(mData[i,vVars_resp]-vMeaniCond)
							## vIntercep_non_phyl is the expectation calculated for ALL the columns of mData
							RSS_non_phyl_null<-t(mData[i,vVars_resp]-vIntercept_non_phyl[vVars_resp])%*%invmCoviCond%*%(mData[i,vVars_resp]-vIntercept_non_phyl[vVars_resp])
							b_no_error2<-TRUE
						    },error=function(e){.my_message(e,FALSE)})
						}
						if (!b_no_error || !b_no_error2){s_warn<- -1}
					    }
					}else{
					    RSS_non_phyl<-0;RSS_non_phyl_null<-0;
					    if (!bis_predNA){
					    	s_warn<- -1
						.my_message(paste("Cannot calculate conditional on predictors RSS for species ",i,". ",sep=""),FALSE)
					    }
					}				    	    
				    }
				    c(RSS_non_phyl,RSS_non_phyl_null,s_warn)
				},vVars_resp=vVars_resp,vVars_pred=vVars_pred,EvolModel=EvolModel,mData=mData,vIntercept_non_phyl=vIntercept_non_phyl,vX0=vX0,vY0=vY0,mPsi=mPsi,A1B=A1B,lexpmtA=lexpmtA,lexpmtjA=lexpmtjA,regimes=regimes,mPsi0=mPsi0,simplify=TRUE))

				vRSS_nonphyl_values<-apply(mRSS_nonphyl_values[,1:2],2,sum)
				v_negRSS<-which(mRSS_nonphyl_values[,3]== -1)
				if (length(v_negRSS)>0){
				    phylLogLik$RSS_non_phylogenetic_conditional_on_predictors_comment<-paste("WARNING: the covariance matrix for species ",paste(v_negRSS,collapse=",")," is not symmetric positive definite! Trying to correct with Matrix::nearPD. Non-phylogenetic RSS calculations are highly doubtful!  This is probably the result of the optimizer getting stuck at some extreme value, you are advised to rerun the optimization. A possible explanation is that the search procedure finished in a very bad local maximum of the likelihood surface. One suggested solution is to re-run on a simpler, i.e. fewer free parameters and then use the output as a starting point for the considered here model. A simpler model can be set using the parameters Atype (from experience this is could be the key parameter to focus on) and Syytype. To use the output as a starting point for the optimization, use the start_point_for_optim parameter provided to the estimation function(s). See mvSLOUCH's manual.",sep="")
				}
				
				phylLogLik$RSS_non_phylogenetic_conditional_on_predictors<-vRSS_nonphyl_values[1]
				RSS_non_phylogenetic_null_model<-vRSS_nonphyl_values[2]
				phylLogLik$R2_non_phylogenetic_conditional_on_predictors<-1-phylLogLik$RSS_non_phylogenetic_conditional_on_predictors/RSS_non_phylogenetic_null_model
				if (!is.element("R2_non_phylogenetic_conditional_on_predictors",names(phylLogLik))||is.na(phylLogLik$R2_non_phylogenetic_conditional_on_predictors)||(phylLogLik$R2_non_phylogenetic_conditional_on_predictors<0)){
				    phylLogLik$R2_non_phylogenetic_conditional_on_predictors_comment<-"R2 is negative, consider rerunning estimation or a different model of evolution. A possible explanation is that the search procedure finished in a very bad local maximum of the likelihood surface. One suggested solution is to re-run on a simpler, i.e. fewer free parameters and then use the output as a starting point for the considered here model. A simpler model can be set using the parameters Atype (from experience this is could be the key parameter to focus on) and Syytype. To use the output as a starting point for the optimization, use the start_point_for_optim parameter provided to the estimation function(s). See mvSLOUCH's manual."
				}
			    },error=function(e){.my_message("Error in calculating conditional RSS: ",TRUE);.my_message(e,TRUE);.my_message("\n",TRUE)})
			}#else{
		    	#    phylLogLik<-vector("list",2)
			#    names(phylLogLik)<-c("RSS","R2")
			#    phylLogLik$RSS<-"Not calculated due to error"
			#    phylLogLik$R2<-"Not calculated due to error"
			##Always calculate the RSS and R2 for the complete model even in the presence of predictors
			tryCatch({
			    vIntercp<-apply(mData,2,mean,na.rm=TRUE)
			    mInterceptCentredData<-mData-matrix(c(vIntercp),nrow=n,ncol=length(vIntercp),byrow=TRUE)
			    #RSS_null_model<-sum((mInterceptCentredData)^2,na.rm=TRUE)
			    pcmbase_model_box_mean0<-.set_mean0_pcmbase_model_box(modelParams$pcmbase_model_box,glsmodel=NA)
			    RSS_null_model<-.pcmbaseDphylGaussian_RSS(mInterceptCentredData,phyltree,pcmbase_model_box_mean0)
			    ##RSS_null_model<-.pcmbaseDphylGaussian_RSS(mInterceptCentredData,phyltree,modelParams$pcmbase_model_box)
			    phylLogLik$RSS<-.pcmbaseDphylGaussian_RSS(mData,phyltree,modelParams$pcmbase_model_box)
			    phylLogLik$R2<-1-phylLogLik$RSS/RSS_null_model
			    
			    mDglsmean<-rep(1,nrow(mData))%x%diag(ncol(mData))        
    			    vGLSmeanest<-.pcmbaseDphylOU_GLS(mData,mD=mDglsmean,phyltree,modelParams$pcmbase_model_box,glsmodel=NA)$vGLSest
    			    mDatamGLSmean<-mData-matrix(vGLSmeanest,nrow=nrow(mData),ncol=ncol(mData),byrow=TRUE)
    			    RSS_phylaverage<- .pcmbaseDphylGaussian_RSS(mDatamGLSmean,phyltree,pcmbase_model_box_mean0,glsmodel=NA) 
			    phylLogLik$R2_phylaverage<-1-phylLogLik$RSS/RSS_phylaverage
				
			    if (!is.element("RSS",names(phylLogLik))||is.na(phylLogLik$RSS)||(phylLogLik$RSS<0)){
			        phylLogLik$RSS_comment<-"RSS is negative, consider rerunning estimation or a different (also non-phylogenetic) model of evolution. A possible explanation is that the search procedure finished in a very bad local maximum of the likelihood surface. One suggested solution is to re-run on a simpler, i.e. fewer free parameters and then use the output as a starting point for the considered here model. A simpler model can be set using the parameters Atype (from experience this is could be the key parameter to focus on) and Syytype. To use the output as a starting point for the optimization, use the start_point_for_optim parameter provided to the estimation function(s). See mvSLOUCH's manual."
			    }
			    if (!is.element("R2",names(phylLogLik))||is.na(phylLogLik$R2)||(phylLogLik$R2<0)){
			        phylLogLik$R2_comment<-"R2 is negative, consider rerunning estimation or a different (also non-phylogenetic) model of evolution. A possible explanation is that the search procedure finished in a very bad local maximum of the likelihood surface. One suggested solution is to re-run on a simpler, i.e. fewer free parameters and then use the output as a starting point for the considered here model. A simpler model can be set using the parameters Atype (from experience this is could be the key parameter to focus on) and Syytype. To use the output as a starting point for the optimization, use the start_point_for_optim parameter provided to the estimation function(s). See mvSLOUCH's manual."
			    }
			    if (!is.element("R2_phylaverage",names(phylLogLik))||is.na(phylLogLik$R2_phylaverage)||(phylLogLik$R2_phylaverage<0)){
			        phylLogLik$R2_phylaverage_comment<-"R2 based on the phylogenetic RSS average is negative, consider rerunning estimation or a different (also non-phylogenetic) model of evolution. A possible explanation is that the search procedure finished in a very bad local maximum of the likelihood surface. One suggested solution is to re-run on a simpler, i.e. fewer free parameters and then use the output as a starting point for the considered here model. A simpler model can be set using the parameters Atype (from experience this is could be the key parameter to focus on) and Syytype. To use the output as a starting point for the optimization, use the start_point_for_optim parameter provided to the estimation function(s). See mvSLOUCH's manual."
			    }
			},error=function(e){.my_message("Error in calculating RSS: ",TRUE);.my_message(e,TRUE);.my_message("\n",TRUE)})
			#}
		    }
		    else{phylLogLik<-.pcmbaseDphylGaussian_RSS(mData,phyltree,modelParams$pcmbase_model_box)}
		}
	    },warning=function(w){.my_warning(paste("Warning in likelihood/RSS calculations",w,"\n"),FALSE,FALSE)},error=function(e){.my_message(e,FALSE);.my_message("\n",FALSE)}
	)
    }
    phylLogLik
}


.callPCMBase_mvlik<-function(mData,phyltree, pcmbase_model_box,b_islog=TRUE,minLogLik=-Inf,b_useLikFunc=TRUE){
## called in estimBM.R, loglik.R, phylgls.R
    phyltree_org<-phyltree
    
    phyltree<-.phyltree_remove_path_fields(phyltree)
    phylLogLik<- minLogLik
    ## mData has to be a matrix we do not correct it to as.matrix() at this stage
    ## such a change should be done earlier

    modelclass<-intersect(c("OU","BM"),class(pcmbase_model_box))
    class(pcmbase_model_box)<-"list"
    pcmbase_model_box_params<-pcmbase_model_box
    
    
    vNArows<-which(apply(mData,1,function(x){all(is.na(x))}))
    if (length(vNArows)>0){	
	bremtipnames<-FALSE
	if (!is.element("tip.label",names(phyltree))){numtips<-length(unique(c(phyltree$edge)))-phyltree$Nnode;bremtipnames<-TRUE}
	brem_rownames<-FALSE
	if (is.null(rownames(mData))){rownames(mData)<-phyltree$tip.label;brem_rownames<-TRUE}
	
	## The error/warning below is not reported as it is caused by A computationally singular
	## causing all entries in evaluated data to be NA
	## so the optimizer should just get out of this place 
	## through an error later on due to empty tree and NA data
	## in likelihood evaluation or other calculations
	tryCatch({
		phyltree<-.phyltree_remove_tips(phyltree,vNArows)
	    },warning=function(m){
                .my_message("Warning raised by ape::drop_tip : ",FALSE);.my_message(m,FALSE);.my_message("\n",FALSE);
            },error=function(m){
                .my_message("Error raised by ape::drop_tip : ",FALSE);.my_message(m,FALSE);.my_message("\n",FALSE);
        })
	mData<-mData[-vNArows,,drop=FALSE]
	mData<-mData[phyltree$tip.label,,drop=FALSE] ## need to reorder in case species order changed
	if (bremtipnames){phyltree$tip.label<-NULL}
	if (brem_rownames){rownames(mData)<-NULL}
	b_useLikFunc<-FALSE
    }
    pcmbase_model_box<-PCMBase::PCM(model=modelclass,k=ncol(mData),regimes=names(pcmbase_model_box_params$Sigma_x[1,1,]),params=pcmbase_model_box_params)
    
    bno_error<-FALSE
    tryCatch({
        if (!b_useLikFunc){	
    	    if (isTRUE(phyltree_org$b_usePCMBaseCpp)){
    		phylLogLik<-PCMBase::PCMLik(X= t(mData),tree = phyltree, model = pcmbase_model_box,log=b_islog, metaI = PCMBaseCpp::PCMInfoCpp);
    	    }else{
    		phylLogLik<-PCMBase::PCMLik(X= t(mData),tree = phyltree, model = pcmbase_model_box,log=b_islog);
    	    }
    	}else{
	    v_pcmbasemodelparams <- double(PCMBase::PCMParamCount(o=pcmbase_model_box))
	    PCMBase::PCMParamLoadOrStore(o=pcmbase_model_box, vecParams=v_pcmbasemodelparams, offset=0, load=FALSE)
	    if (modelclass=="OU"){
		phylLogLik<-phyltree_org$likFun_OU(v_pcmbasemodelparams)
	    }
	    if (modelclass=="BM"){
		if (ncol(mData)==phyltree_org$kX){
		    phylLogLik<-phyltree_org$likFun_BM_kX(v_pcmbasemodelparams)
		}else{
		    phylLogLik<-phyltree_org$likFun_BM_all(v_pcmbasemodelparams)
		}
	    }
	}
	bno_error<-TRUE
    }
    ,"logic_error"=function(e){.my_message(paste("Error in call to PCMBaseCpp : \n",e,"\n",sep=""),FALSE)}
    ,"std::logic_error"=function(e){.my_message(paste("Error in call to PCMBaseCpp : \n",e,"\n",sep=""),FALSE)}
    ,"C++Error"=function(e){.my_message(paste("Error in call to PCMBaseCpp : \n",e,"\n",sep=""),FALSE)}
    ,error=function(e){.my_message(paste("Error in call to PCMBase::PCMLik : \n",e,"\n",sep=""),FALSE)},warning=function(e){.my_message(paste("Error in call to PCMBase::PCMLik : \n",e,"\n",sep=""),FALSE)})
    if (bno_error){
	attributes(phylLogLik)<-NULL
    }else{phylLogLik<- minLogLik}
    
    phylLogLik
}

.my_stop<-function(m,bdoPrint=TRUE){
## called in OUphylregression.R
    class(m)<-"character" 
    if (bdoPrint){stop(m,call. = FALSE)}
    else{stop(call. = FALSE)}
    NA
}

.my_message<-function(m,bdoPrint=TRUE){
    class(m)<-"character" 
    if (bdoPrint){message(m)}
    NA
}

.my_warning<-function(m,bdoPrint=TRUE,baswarning=TRUE){
    class(m)<-"character" 
    if (baswarning){
	if (bdoPrint){warning(m, call. = FALSE)}
	else{warning(call. = FALSE)}
    }else{
	if (bdoPrint){message(m)}
    }
    NA
}

