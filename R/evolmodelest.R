## This file is part of mvSLOUCH

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this package or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .


estimate.evolutionary.model<-function(phyltree,mData,regimes=NULL,root.regime=NULL,M.error=NULL,repeats=5,model.setups=NULL,predictors=NULL,kY=NULL,doPrint=FALSE,pESS=NULL,estimate.root.state=FALSE,min_bl=0.0003,maxiter=c(10,50,100)){
## checking if data is correctly passed is done heavily in .InitialRegimeSetup
    
    testedModels<-list()
    j<-1
        
## List below defines all possible interesting models
    if (is.null(model.setups)){
	model.setups<-.generate.basic.model.setups()    
    }else{
	if (!is.list(model.setups)){
	    model.setups<-switch(model.setups,
		univariate=.generate.univ.model.setups(),
		basic=.generate.basic.model.setups(),
		fundamental=.generate.fund.model.setups(),
		extended=.generate.ext.model.setups(),
		all=.generate.all.model.setups(),
		.my_stop("Incorrect model list name provided.",TRUE)
	   )
	}
    }

## ====================================================================================
## Initialization of regimes and tree pre-calculations, no point of doing this
## expensive calculation for every tested model   
    mData<-.check_input_trait_data(mData,NULL,phyltree$tip.label)
    M.error<-.createMeasurementError(M.error,nrow(mData),ncol(mData))
    kYX<-ncol(mData)
    kX<-kYX
    if (!is.null(kY)){
	tmpkX<-kYX-kY
    }
    else{
	##if (!is.null(predictors)){tmpkX<-kYX-length(predictors)}
	if (!is.null(predictors)){tmpkX<-length(predictors)}
	else{tmpkX<-kYX-1}
    }
    if (tmpkX<=0){
	.my_stop('Cannot have all variables as responses (i.e. kY>=ncol(mData)). Either set kY to NULL or choose a smaller (than number of columns in mData) number of responses ("OU traits" in OUBM model).',TRUE)
    }
    phyltree<-.InitialRegimeSetup(phyltree,regimes,regimes.times=NULL,mData=mData,kX=tmpkX,kYX=kYX,root.regime=root.regime,M.error=M.error,bSave_in_phyltree=TRUE)
    mData<-.check_input_trait_data(mData,n=phyltree$Ntips,vSpeciesLabels=phyltree$tip.label)
## =====================================================================================

    if ((!is.vector(maxiter))|| (!is.numeric(maxiter)) || (length(maxiter)!=3)){
        maxiter<-c(10,50,100);
        .my_warning("WARNING: maxiter passed in a wrong way, setting it to default of c(10,50,100)",TRUE,FALSE)
    }

    
    BestModel<-list(BestModel=NA,aic.c=10000,bic=10000,i=NA,model=NA,evolmodel=NA)
    if (!is.null(pESS)){
	phyltreeESS<-phyltree
	if (pESS!="only_calculate"){BestModelESS<-list(BestModel=NA,aic.c=10000,bic=10000,i=NA,model=NA,ESScrit=NA)}
	vNAs<-which(is.na(c(t(mData))))
	if (length(vNAs)==0){vNAs<-NULL}
	vNAs_forESS<-vNAs
	M_error_for_ESS<-M.error
    }else{phyltreeESS<-NULL;BestModelESS<-NULL;vNAs_forESS<-NULL;M_error_for_ESS<-NULL}
    
    used_model_setups<-list()
    for (i in 1:repeats){
	for (k in 1:length(model.setups)){
	    if ((model.setups[[k]]$evolmodel=="bm") && (i==1)){
	    ## no point in doing BM more than once as it will always give the same estimation result
		used_model_setups<-c(used_model_setups,model.setups[k])
		if (doPrint){.my_message("Doing estimation for BM model.\n",TRUE)}
		BMres<-NULL
		tryCatch({
		    BMres<-.internal_BrownianMotionModel(phyltree=phyltree,mData=mData,predictors=predictors,M.error=M.error,min_bl=min_bl)
		},error=function(e){.my_message(e,TRUE);.my_message("\n",TRUE)})
		testedModels[[j]]<-list()
		testedModels[[j]]$result<-BMres;testedModels[[j]]$aic.c<-NA;testedModels[[j]]$bic<-NA;testedModels[[j]]$model<-model.setups[[k]];j<-j+1	    
		
		if (!is.null(pESS)){
    		    mData_forESS<-mData
    		    vNAs<-which(is.na(c(t(mData))))
		    if (length(vNAs)==0){vNAs<-NULL}
		    vNAs_forESS<-vNAs
    		}

		l_estres_postproc<-.postproc_estres(BMres, BestModel, "bm", model.setups[[k]], j-1, pESS, BestModelESS=BestModelESS, phyltreeESS=phyltreeESS, mData=mData_forESS, M.error=M_error_for_ESS,vNAs=vNAs_forESS)
		BestModel<-l_estres_postproc$BestModel
		testedModels[[j-1]]$aic.c<-l_estres_postproc$aic.c
		testedModels[[j-1]]$bic<-l_estres_postproc$bic
		if (!is.null(pESS)){
		    testedModels[[j-1]]$ESScalcs<-l_estres_postproc$calcESS
		}
		phyltreeESS<-l_estres_postproc$phyltreeESS
		BestModelESS<-l_estres_postproc$BestModelESS    		
	    }
	    if (model.setups[[k]]$evolmodel=="ouch"){
		used_model_setups<-c(used_model_setups,model.setups[k])
		if(doPrint){.my_message(paste("Doing estimation for ouch model with A: ",model.setups[[k]]$Atype," with diagonal: ",model.setups[[k]]$diagA," Syy: ",model.setups[[k]]$Syytype,"\n",sep=""),TRUE)}
		OUres<-NULL
		lStartPoint<-NULL
		
		bdoanalytical_start<-TRUE
		if (is.element("start_point_for_optim",names(model.setups[[k]]))){
		    if ((i==2)||(repeats==1)){lStartPoint<-model.setups[[k]]$start_point_for_optim;bdoanalytical_start<-FALSE}
		}
		if ((bdoanalytical_start)&&(i==1)){
		    lStartPoint<-.createStartPointsASyyB(mData,phyltree$tree_height,model.setups[[k]],ncol(mData),FALSE)		    
		}
		tryCatch({
		    OUres<-.internal_ouchModel(phyltree=phyltree,mData=mData,regimes=regimes,regimes.times=NULL,root.regime=root.regime,predictors=predictors,M.error=M.error,Atype=model.setups[[k]]$Atype,Syytype=model.setups[[k]]$Syytype,diagA=model.setups[[k]]$diagA,estimate.root.state=estimate.root.state,parameter_signs=model.setups[[k]]$parameter_signs,lStartPoint=lStartPoint,parscale=model.setups[[k]]$parscale,min_bl=min_bl,maxiter=c(maxiter[1],maxiter[3]))
		},error=function(e){.my_message(e,TRUE);.my_message("\n",TRUE)})

		testedModels[[j]]<-list()
		testedModels[[j]]$result<-OUres;testedModels[[j]]$aic.c<-NA;testedModels[[j]]$bic<-NA;testedModels[[j]]$model<-model.setups[[k]];j<-j+1
    		
    		if (!is.null(pESS)){
    		    mData_forESS<-mData
    		    vNAs<-which(is.na(c(t(mData))))
		    if (length(vNAs)==0){vNAs<-NULL}
		    vNAs_forESS<-vNAs
    		}
    		l_estres_postproc<-.postproc_estres(OUres, BestModel, "ouch", model.setups[[k]], j-1, pESS, BestModelESS=BestModelESS, phyltreeESS=phyltreeESS, mData=mData_forESS, M.error=M_error_for_ESS,vNAs=vNAs_forESS)
		BestModel<-l_estres_postproc$BestModel
		testedModels[[j-1]]$aic.c<-l_estres_postproc$aic.c
		testedModels[[j-1]]$bic<-l_estres_postproc$bic
		if (!is.null(pESS)){
		    testedModels[[j-1]]$ESScalcs<-l_estres_postproc$calcESS
		}
		phyltreeESS<-l_estres_postproc$phyltreeESS
		BestModelESS<-l_estres_postproc$BestModelESS    		
	    }
	    if (model.setups[[k]]$evolmodel=="mvslouch"){
		used_model_setups<-c(used_model_setups,model.setups[k])
		if (is.null(kY)){
		    if (!is.null(predictors)){
			mData.mvsl<-mData[,c(setdiff(1:ncol(mData),predictors),predictors)]
			kY<-ncol(mData)-predictors
		    }else{kY<-1;mData.mvsl<-mData}
		}else{mData.mvsl<-mData}
		if(doPrint){.my_message(paste("Doing estimation for mvslouch model with A: ",model.setups[[k]]$Atype," with diagonal: ",model.setups[[k]]$diagA," Syy: ",model.setups[[k]]$Syytype,"\n",sep=""),TRUE)}
	    	mvslres<-NULL
		
		lStartPoint<-NULL
		
		bdoanalytical_start<-TRUE
		if (is.element("start_point_for_optim",names(model.setups[[k]]))){
		    if ((i==2)||(repeats==1)){lStartPoint<-model.setups[[k]]$start_point_for_optim;bdoanalytical_start<-FALSE}
		}
		if (!is.element("estimateBmethod",names(model.setups[[k]]))){
		    estimateBmethod<-"ML"		    
		}else{estimateBmethod<-model.setups[[k]]$estimateBmethod}

		if ((bdoanalytical_start)&&(i==1)){
		    bdoB<-TRUE
		    if (is.element("estimateBmethod",names(model.setups[[k]])) && (model.setups[[k]]$estimateBmethod=="GLS")){bdoB<-FALSE}
		    lStartPoint<-.createStartPointsASyyB(mData.mvsl,phyltree$tree_height,model.setups[[k]],kY,bdoB)		    
		}
		
	    	tryCatch({
	    	    mvslres<-.internal_mvslouchModel(phyltree=phyltree,mData=mData.mvsl,kY,regimes=regimes,regimes.times=NULL,root.regime=root.regime,predictors=predictors,M.error=M.error,Atype=model.setups[[k]]$Atype,Syytype=model.setups[[k]]$Syytype,diagA=model.setups[[k]]$diagA,estimate.root.state=estimate.root.state,parameter_signs=model.setups[[k]]$parameter_signs,lStartPoint=lStartPoint,parscale=model.setups[[k]]$parscale,min_bl=min_bl,maxiter=maxiter,estimateBmethod=estimateBmethod)
		},error=function(e){.my_message(e,TRUE);.my_message("\n",TRUE)})
				
		testedModels[[j]]<-list()
		testedModels[[j]]$result<-mvslres;testedModels[[j]]$aic.c<-NA;testedModels[[j]]$bic<-NA;testedModels[[j]]$model<-model.setups[[k]];j<-j+1

    		if (!is.null(pESS)){
    		    mData_forESS<-mData.mvsl
    		    vNAs<-which(is.na(c(t(mData.mvsl))))
		    if (length(vNAs)==0){vNAs<-NULL}
		    vNAs_forESS<-vNAs
    		}
    		l_estres_postproc<-.postproc_estres(mvslres, BestModel, "mvslouch", model.setups[[k]], j-1, pESS, BestModelESS=BestModelESS, phyltreeESS=phyltreeESS, mData=mData_forESS, M.error=M_error_for_ESS,vNAs=vNAs_forESS)
		BestModel<-l_estres_postproc$BestModel
		testedModels[[j-1]]$aic.c<-l_estres_postproc$aic.c
		testedModels[[j-1]]$bic<-l_estres_postproc$bic
		if (!is.null(pESS)){
		    testedModels[[j-1]]$ESScalcs<-l_estres_postproc$calcESS
		}
		phyltreeESS<-l_estres_postproc$phyltreeESS
		BestModelESS<-l_estres_postproc$BestModelESS	    
	    }
	}
    }   
    tryCatch({BestModel<-.describe.best.model(BestModel)},error=function(e){.my_message("Cannot describe best model, please go into returned object: ",TRUE);.my_message(e,TRUE);.my_message("\n",TRUE)})
    model.setups<-used_model_setups
    ## remove empty trailing parameter_signs
    model.setups<-sapply(model.setups, function(x){if((is.element("parameter_signs",names(x)))&&(length(x$parameter_signs)==0)){x$parameter_signs<-NULL};x},simplify=FALSE)
    ##model.setups<-rep(model.setups,repeats)
    res<-NA
    if ((!is.null(pESS))&&(pESS!="only_calculate")){
	tryCatch({BestModelESS<-.describe.best.model(BestModelESS)},error=function(e){.my_message("Cannot describe best ESS model, please go into returned object: ",TRUE);.my_message(e,TRUE);.my_message("\n",TRUE)})
	res<-list(BestModelESS=BestModelESS,BestModel=BestModel,testedModels=testedModels,model.setups=model.setups,repeats=repeats)
    }else{
    	res<-list(BestModel=BestModel,testedModels=testedModels,model.setups=model.setups,repeats=repeats)
    }
    res
}


.postproc_estres<-function(estres, BestModel, evolmodel, model_setup, i, pESS, BestModelESS=NULL, phyltreeESS=NULL, mData=NULL, M.error=NULL,vNAs=NULL) {
## phyltree and phyltreeESS are the same thing!
    loutput<-list(BestModel=BestModel,aic.c=NA,bic=NA,calcESS=NULL,phyltreeESS=phyltreeESS,BestModelESS=BestModelESS)
    if (!is.null(estres)){
	modname<-"BM"
	if (evolmodel=="mvslouch"){modname<-"OUBM"}
	if (evolmodel=="ouch"){modname<-"OUOU"}
	if (evolmodel=="bm"){model_call<-"bm"}else{
	        model_call<-paste(modname,": ",evolmodel," model with A: ",model_setup$Atype," with diagonal: ",model_setup$diagA," Syy: ",model_setup$Syytype,sep="")
	}
	if( (evolmodel=="mvslouch")||(evolmodel=="ouch")){
	    if (is.list(estres$MaxLikFound)){estres<-estres$MaxLikFound}
	    else{estres<-estres$FinalFound}
	}
	tryCatch({BestModel<-.return_best_model(BestModel,estres,evolmodel,model_call,model_setup,i,is_ESS=FALSE,calcESS=NULL)},error=function(e){.my_message("Cannot compare models",FALSE);.my_message(e,FALSE);.my_message("\n",FALSE)})
	aic.c_value<- Inf
	bic_value<- Inf
	if (is.element("ParamSummary",names(estres))){
	    if ((is.element("aic.c",names(estres$ParamSummary)))&&(is.element("bic",names(estres$ParamSummary)))){
		aic.c_value<-estres$ParamSummary$aic.c ## if some error and not present, then these are NULL
		bic_value<-estres$ParamSummary$bic
	    }
	}
	calcESS<-NULL
	if ((!is.null(pESS))&&(!is.null(phyltreeESS))){
	    tryCatch({
		calcESS<-.calcESSanalytical(phyltreeESS,proc.params=estres$ParamsInModel,evolmodel=evolmodel,mData=mData,Merror=M.error,vNAs=vNAs,ESS.method=pESS)
		if (is.element("phyltreeESS",names(calcESS))){
		    phyltreeESS<-calcESS$phyltree
		    calcESS$phyltree<-NULL
		}
		if (pESS!="only_calculate"){
		    tryCatch({BestModelESS<-.return_best_model(BestModelESS,estres,evolmodel,model_call,model_setup,i,is_ESS=TRUE,calcESS=calcESS)},error=function(e){.my_message("Cannot compare ESS models",FALSE);.my_message(e,FALSE);.my_message("\n",FALSE)})
		}
	    },error=function(e){.my_message("Error: cannot calculate ESS!",TRUE);.my_message(e,TRUE);.my_message("\n",TRUE)})
	}
	loutput<-list(BestModel=BestModel,aic.c=aic.c_value,bic=bic_value,calcESS=calcESS,phyltreeESS=phyltreeESS,BestModelESS=BestModelESS)

    }    
    loutput
}

.return_best_model<-function(currbest,candidate_model,evolmodel,model_call,model_setup,i,is_ESS=FALSE,calcESS=NULL){
    bisbetter<-c(FALSE,FALSE,FALSE,TRUE)
    if (is_ESS){
	if (is.null(calcESS)){
	    .my_message("No ESS object cannot compare models",FALSE)
	}else{
	    tryCatch({
		crit_to_cf<-.getESScriteria(candidate_model$ParamSummary$LogLik,candidate_model$ParamSummary$dof,calcESS$ESS.model.selection,calcESS$ESS.factor.model.selection,calcESS$rhon,candidate_model$ParamSummary$RSS)
	    },error=function(e){.my_message("No ESS object cannot compare models",FALSE);.my_message(e,FALSE);.my_message("\n",FALSE)})
	}
    }else{crit_to_cf<-candidate_model$ParamSummary}
    tryCatch({bisbetter<-.check_is_better(currbest,crit_to_cf)},error=function(e){.my_message("Cannot compare models",FALSE);.my_message(e,FALSE);.my_message("\n",FALSE)})
    
    if (bisbetter[1]){
	if(is_ESS){currbest$ESScrit<-crit_to_cf}
    	currbest$BestModel<-candidate_model
	currbest$aic.c<-crit_to_cf$aic.c
	currbest$bic<-crit_to_cf$bic				
	currbest$i<-i
	currbest$model.call<-model_call
	currbest$model<-model_setup
	currbest$evolmodel<-evolmodel
	if (bisbetter[2]){.my_warning(paste("WARNING: Infinite AICc (used BIC for comparison) - model is saturated, this indicates too few observations! Model: ",currbest$model.call),TRUE,FALSE)}	
    	if (bisbetter[3]){.my_warning(paste("WARNING: Missing AICc value, using BIC to compare Model: ",currbest$model.call),TRUE,FALSE)}
    }else{
	if (bisbetter[4]){.my_warning(paste("WARNING: Cannot calculate information criteria for model: ",candidate_model$model.call),TRUE,FALSE)}
    }
    currbest
}

.check_is_better<-function(currbest,candidate_model){
    bisbetter<-TRUE
    bwarn_inf_aic.c<-FALSE
    bwarn_na_aic.c<-FALSE
    bwarn_no_infcrit<-FALSE
    curr_aic.c<-currbest$aic.c
    curr_bic<-currbest$bic
    cand_aic.c<- Inf
    cand_bic<- Inf
    
    if (!(is.null(currbest)&&is.null(candidate_model))){
	if (!is.null(currbest)){
	    if (is.null(candidate_model)){bisbetter<-FALSE}
	    if (bisbetter && (!is.element("aic.c",names(candidate_model)))){cand_aic.c<-NA}
	    if (bisbetter && (!is.na(cand_aic.c)) && (is.na(candidate_model$aic.c) || is.infinite(candidate_model$aic.c) )){
		if (is.infinite(candidate_model$aic.c)){bwarn_inf_aic.c<-TRUE}
		cand_aic.c<-NA
	    }
	    if (bisbetter && (!is.na(cand_aic.c))) {cand_aic.c<-candidate_model$aic.c}
	    if (bisbetter && is.na(cand_aic.c) && (!is.element("bic",names(candidate_model)))){	cand_bic<-NA}
	    if (bisbetter && is.na(cand_aic.c) && (!is.na(cand_bic)) && (is.na(candidate_model$bic) || is.infinite(candidate_model$bic) )){cand_bic<-NA}
	    if (bisbetter && is.na(cand_aic.c) && (!is.na(cand_bic))) {cand_bic<-candidate_model$bic}
	    if (bisbetter){
		bisbetter<-FALSE
		if (!is.na(cand_aic.c)){
		    if (cand_aic.c< curr_aic.c){bisbetter<-TRUE}
		}else{
		    if (!is.na(cand_bic)){
			if (!bwarn_inf_aic.c){bwarn_na_aic.c<-TRUE}
			if (cand_bic< curr_bic){bisbetter<-TRUE}
		    }else{bwarn_no_infcrit<-TRUE}
		}
	    }
	}
    }else{bwarn_no_infcrit<-TRUE}
    c(bisbetter,bwarn_inf_aic.c,bwarn_na_aic.c,bwarn_no_infcrit)
}

.describe.best.model<-function(BestModel){
    ballPosEig<-FALSE 
    numtraits<-1
    
    if (is.element("evolmodel",names(BestModel)) && is.element(BestModel$evolmodel,c("bm","ouch","mvslouch"))){
	if ((BestModel$evolmodel=="ouch")||(BestModel$evolmodel=="mvslouch")){
	    if (length(which(Re(BestModel$BestModel$ParamSummary$phyl.halflife$halflives["halflife",])>0)) == length(BestModel$BestModel$ParamSummary$phyl.halflife$halflives["halflife",])){ballPosEig<-TRUE}
	    numtraits<-ncol(BestModel$BestModel$ParamsInModel$A)
	}
	if (BestModel$evolmodel=="bm"){numtraits<-ncol(BestModel$BestModel$ParamsInModel$Sxx)}    
	BestModel$model.description<-.generate.model.description(BestModel$evolmodel,BestModel$model,numtraits,ballPosEig)
	BestModel$key.properties<-.extract.model.key.properties(BestModel$evolmodel,BestModel$BestModel,k=numtraits)
        BestModel$parameter.SE<-.generate.model.se(BestModel$evolmodel)
	if (is.element("ESScrit",names(BestModel))){
	    BestModel<-BestModel[c("model.description","key.properties","ESScrit","parameter.SE","aic.c","evolmodel","model","model.call","BestModel","i")]
	}else{
	    BestModel<-BestModel[c("model.description","key.properties","parameter.SE","aic.c","evolmodel","model","model.call","BestModel","i")]
	}
	if ((BestModel$evolmodel=="ouch")||(BestModel$evolmodel=="mvslouch")){
	    BestModel$BestModel<-BestModel$BestModel[c("ParamsInModel","ParamSummary","LogLik","HeuristicSearchPointFinalFind")]
	}
	if (BestModel$evolmodel=="bm"){
	    BestModel$BestModel<-BestModel$BestModel[c("ParamsInModel","ParamSummary")]
	    BestModel$BestModel$LogLik<-BestModel$BestModel$ParamSummary$LogLik
	    BestModel$BestModel$HeuristicSearchPointFinalFind<-NA
	}
    }else{
	.my_warning("Found best model seems empty, please check final outputted object!",TRUE,TRUE)
    }
    BestModel
}

.generate.list.of.model.setups<-function(vevolmodels,vAtypes,vSyytypes,vdiagA){
    msetups<-expand.grid(vevolmodels,vAtypes,vSyytypes,vdiagA, stringsAsFactors = FALSE)   
   msetups<-as.matrix(msetups) 
    lres<-sapply(1:nrow(msetups),function(i,msetups){x<-msetups[i,];x<-unname(x);list(evolmodel=x[1],Atype=x[2],Syytype=x[3],diagA=x[4],parameter_signs=list())},msetups=msetups,simplify=FALSE)
    if (!is.element("bm",vevolmodels)){
	lres[[length(lres)+1]]<-list(evolmodel="bm")
    }
    lres
}
                
.generate.basic.model.setups<-function(){
    vevolmodels<-c("ouch","mvslouch")
    vAtypes<-c("Diagonal","UpperTri","LowerTri","DecomposablePositive","DecomposableReal")
    vSyytypes<-c("Diagonal","UpperTri")
    vdiagA<-c("Positive")
    .generate.list.of.model.setups(vevolmodels,vAtypes,vSyytypes,vdiagA)
}

.generate.univ.model.setups<-function(){
    vevolmodels<-c("ouch")
    vAtypes<-c("Diagonal")
    vSyytypes<-c("Diagonal")
    vdiagA<-c("Positive","Negative")
    .generate.list.of.model.setups(vevolmodels,vAtypes,vSyytypes,vdiagA)
}
                                    
.generate.fund.model.setups<-function(){
    vevolmodels<-c("ouch","mvslouch")
    vAtypes<-c("Diagonal","UpperTri","LowerTri","SymmetricPositiveDefinite","DecomposablePositive","DecomposableReal","Invertible")
    vSyytypes<-c("Diagonal","UpperTri")
    vdiagA<-c("Positive",NULL)
    .generate.list.of.model.setups(vevolmodels,vAtypes,vSyytypes,vdiagA)
}
                                                        
.generate.ext.model.setups<-function(){
    generate.model.setups()
}
                                                            
.generate.all.model.setups<-function(){
    .my_message("WARNING: all allowed model setups will be analyzed. This will take a very long time!\n",TRUE)
    vevolmodels<-c("ouch","mvslouch")
    vAtypes<-c("SingleValueDiagonal","Diagonal","UpperTri","LowerTri","SymmetricPositiveDefinite","Symmetric","DecomposablePositive","DecomposableNegative","DecomposableReal","Invertible")
    vSyytypes<-c("SingleValueDiagonal","Diagonal","UpperTri","Symmetric")
    vdiagA<-c("Positive",NULL,"Negative")
    .generate.list.of.model.setups(vevolmodels,vAtypes,vSyytypes,vdiagA)
}


generate.model.setups<-function(){
## This function should be really hidden but is made available so that the user can create a personalized version of it to speed up estimation.
	list(
	    list(evolmodel="bm"),
	    list(evolmodel="ouch",Atype="Invertible",Syytype="UpperTri",diagA=NULL),
	    list(evolmodel="ouch",Atype="DecomposableReal",Syytype="UpperTri",diagA=NULL),
	    list(evolmodel="ouch",Atype="DecomposablePositive",Syytype="UpperTri",diagA=NULL),
	    list(evolmodel="ouch",Atype="SymmetricPositiveDefinite",Syytype="UpperTri",diagA=NULL),
	    list(evolmodel="ouch",Atype="UpperTri",Syytype="UpperTri",diagA=NULL),
	    list(evolmodel="ouch",Atype="LowerTri",Syytype="UpperTri",diagA=NULL),
	    list(evolmodel="ouch",Atype="Diagonal",Syytype="UpperTri",diagA=NULL),
	    list(evolmodel="ouch",Atype="SingleValueDiagonal",Syytype="UpperTri",diagA=NULL),
	    list(evolmodel="ouch",Atype="Invertible",Syytype="UpperTri",diagA="Positive"),
	    list(evolmodel="ouch",Atype="DecomposableReal",Syytype="UpperTri",diagA="Positive"),
	    list(evolmodel="ouch",Atype="DecomposablePositive",Syytype="UpperTri",diagA="Positive"),
	    list(evolmodel="ouch",Atype="UpperTri",Syytype="UpperTri",diagA="Positive"),
	    list(evolmodel="ouch",Atype="LowerTri",Syytype="UpperTri",diagA="Positive"),
	    list(evolmodel="ouch",Atype="Diagonal",Syytype="UpperTri",diagA="Positive"),
	    list(evolmodel="ouch",Atype="SingleValueDiagonal",Syytype="UpperTri",diagA="Positive"),
	    list(evolmodel="ouch",Atype="Invertible",Syytype="Diagonal",diagA=NULL),
	    list(evolmodel="ouch",Atype="DecomposableReal",Syytype="Diagonal",diagA=NULL),
	    list(evolmodel="ouch",Atype="DecomposablePositive",Syytype="Diagonal",diagA=NULL),
	    list(evolmodel="ouch",Atype="SymmetricPositiveDefinite",Syytype="Diagonal",diagA=NULL),
	    list(evolmodel="ouch",Atype="UpperTri",Syytype="Diagonal",diagA=NULL),
	    list(evolmodel="ouch",Atype="LowerTri",Syytype="Diagonal",diagA=NULL),
	    list(evolmodel="ouch",Atype="Diagonal",Syytype="Diagonal",diagA=NULL),
	    list(evolmodel="ouch",Atype="SingleValueDiagonal",Syytype="Diagonal",diagA=NULL),
	    list(evolmodel="ouch",Atype="Invertible",Syytype="Diagonal",diagA="Positive"),
	    list(evolmodel="ouch",Atype="DecomposableReal",Syytype="Diagonal",diagA="Positive"),
	    list(evolmodel="ouch",Atype="DecomposablePositive",Syytype="Diagonal",diagA="Positive"),
	    list(evolmodel="ouch",Atype="UpperTri",Syytype="Diagonal",diagA="Positive"),
	    list(evolmodel="ouch",Atype="LowerTri",Syytype="Diagonal",diagA="Positive"),
	    list(evolmodel="ouch",Atype="Diagonal",Syytype="Diagonal",diagA="Positive"),
	    list(evolmodel="ouch",Atype="SingleValueDiagonal",Syytype="Diagonal",diagA="Positive"),
	    list(evolmodel="ouch",Atype="Invertible",Syytype="SingleValueDiagonal",diagA=NULL),
	    list(evolmodel="ouch",Atype="DecomposableReal",Syytype="SingleValueDiagonal",diagA=NULL),
	    list(evolmodel="ouch",Atype="DecomposablePositive",Syytype="SingleValueDiagonal",diagA=NULL),
	    list(evolmodel="ouch",Atype="SymmetricPositiveDefinite",Syytype="SingleValueDiagonal",diagA=NULL),
	    list(evolmodel="ouch",Atype="UpperTri",Syytype="SingleValueDiagonal",diagA=NULL),
	    list(evolmodel="ouch",Atype="LowerTri",Syytype="SingleValueDiagonal",diagA=NULL),
	    list(evolmodel="ouch",Atype="Diagonal",Syytype="SingleValueDiagonal",diagA=NULL),
	    list(evolmodel="ouch",Atype="SingleValueDiagonal",Syytype="SingleValueDiagonal",diagA=NULL),
	    list(evolmodel="ouch",Atype="Invertible",Syytype="SingleValueDiagonal",diagA="Positive"),
	    list(evolmodel="ouch",Atype="DecomposableReal",Syytype="SingleValueDiagonal",diagA="Positive"),
	    list(evolmodel="ouch",Atype="DecomposablePositive",Syytype="SingleValueDiagonal",diagA="Positive"),
	    list(evolmodel="ouch",Atype="UpperTri",Syytype="SingleValueDiagonal",diagA="Positive"),
	    list(evolmodel="ouch",Atype="LowerTri",Syytype="SingleValueDiagonal",diagA="Positive"),
	    list(evolmodel="ouch",Atype="Diagonal",Syytype="SingleValueDiagonal",diagA="Positive"),
	    list(evolmodel="ouch",Atype="SingleValueDiagonal",Syytype="SingleValueDiagonal",diagA="Positive"),
	    list(evolmodel="mvslouch",Atype="Invertible",Syytype="UpperTri",diagA=NULL),
	    list(evolmodel="mvslouch",Atype="DecomposableReal",Syytype="UpperTri",diagA=NULL),
	    list(evolmodel="mvslouch",Atype="DecomposablePositive",Syytype="UpperTri",diagA=NULL),
	    list(evolmodel="mvslouch",Atype="SymmetricPositiveDefinite",Syytype="UpperTri",diagA=NULL),
	    list(evolmodel="mvslouch",Atype="UpperTri",Syytype="UpperTri",diagA=NULL),
	    list(evolmodel="mvslouch",Atype="LowerTri",Syytype="UpperTri",diagA=NULL),
	    list(evolmodel="mvslouch",Atype="Diagonal",Syytype="UpperTri",diagA=NULL),
	    list(evolmodel="mvslouch",Atype="SingleValueDiagonal",Syytype="UpperTri",diagA=NULL),
	    list(evolmodel="mvslouch",Atype="Invertible",Syytype="UpperTri",diagA="Positive"),
	    list(evolmodel="mvslouch",Atype="DecomposableReal",Syytype="UpperTri",diagA="Positive"),
	    list(evolmodel="mvslouch",Atype="DecomposablePositive",Syytype="UpperTri",diagA="Positive"),
	    list(evolmodel="mvslouch",Atype="UpperTri",Syytype="UpperTri",diagA="Positive"),
	    list(evolmodel="mvslouch",Atype="LowerTri",Syytype="UpperTri",diagA="Positive"),
	    list(evolmodel="mvslouch",Atype="Diagonal",Syytype="UpperTri",diagA="Positive"),
	    list(evolmodel="mvslouch",Atype="SingleValueDiagonal",Syytype="UpperTri",diagA="Positive"),
	    list(evolmodel="mvslouch",Atype="Invertible",Syytype="Diagonal",diagA=NULL),
	    list(evolmodel="mvslouch",Atype="DecomposableReal",Syytype="Diagonal",diagA=NULL),
	    list(evolmodel="mvslouch",Atype="DecomposablePositive",Syytype="Diagonal",diagA=NULL),
	    list(evolmodel="mvslouch",Atype="SymmetricPositiveDefinite",Syytype="Diagonal",diagA=NULL),
	    list(evolmodel="mvslouch",Atype="UpperTri",Syytype="Diagonal",diagA=NULL),
	    list(evolmodel="mvslouch",Atype="LowerTri",Syytype="Diagonal",diagA=NULL),
	    list(evolmodel="mvslouch",Atype="Diagonal",Syytype="Diagonal",diagA=NULL),
	    list(evolmodel="mvslouch",Atype="SingleValueDiagonal",Syytype="Diagonal",diagA=NULL),
	    list(evolmodel="mvslouch",Atype="Invertible",Syytype="Diagonal",diagA="Positive"),
	    list(evolmodel="mvslouch",Atype="DecomposableReal",Syytype="Diagonal",diagA="Positive"),
	    list(evolmodel="mvslouch",Atype="DecomposablePositive",Syytype="Diagonal",diagA="Positive"),
	    list(evolmodel="mvslouch",Atype="UpperTri",Syytype="Diagonal",diagA="Positive"),
	    list(evolmodel="mvslouch",Atype="LowerTri",Syytype="Diagonal",diagA="Positive"),
	    list(evolmodel="mvslouch",Atype="Diagonal",Syytype="Diagonal",diagA="Positive"),
	    list(evolmodel="mvslouch",Atype="SingleValueDiagonal",Syytype="Diagonal",diagA="Positive"),
	    list(evolmodel="mvslouch",Atype="Invertible",Syytype="SingleValueDiagonal",diagA=NULL),
	    list(evolmodel="mvslouch",Atype="DecomposableReal",Syytype="SingleValueDiagonal",diagA=NULL),
	    list(evolmodel="mvslouch",Atype="DecomposablePositive",Syytype="SingleValueDiagonal",diagA=NULL),
	    list(evolmodel="mvslouch",Atype="SymmetricPositiveDefinite",Syytype="SingleValueDiagonal",diagA=NULL),
	    list(evolmodel="mvslouch",Atype="UpperTri",Syytype="SingleValueDiagonal",diagA=NULL),
	    list(evolmodel="mvslouch",Atype="LowerTri",Syytype="SingleValueDiagonal",diagA=NULL),
	    list(evolmodel="mvslouch",Atype="Diagonal",Syytype="SingleValueDiagonal",diagA=NULL),
	    list(evolmodel="mvslouch",Atype="SingleValueDiagonal",Syytype="SingleValueDiagonal",diagA=NULL),
	    list(evolmodel="mvslouch",Atype="Invertible",Syytype="SingleValueDiagonal",diagA="Positive"),
	    list(evolmodel="mvslouch",Atype="DecomposableReal",Syytype="SingleValueDiagonal",diagA="Positive"),
	    list(evolmodel="mvslouch",Atype="DecomposablePositive",Syytype="SingleValueDiagonal",diagA="Positive"),
	    list(evolmodel="mvslouch",Atype="UpperTri",Syytype="SingleValueDiagonal",diagA="Positive"),
	    list(evolmodel="mvslouch",Atype="LowerTri",Syytype="SingleValueDiagonal",diagA="Positive"),
	    list(evolmodel="mvslouch",Atype="Diagonal",Syytype="SingleValueDiagonal",diagA="Positive"),
	    list(evolmodel="mvslouch",Atype="SingleValueDiagonal",Syytype="SingleValueDiagonal",diagA="Positive")
	)
}

.generate.model.description<-function(evolmodel,model.setup,numtraits,ballPosEig){
    chdesc<-"This is a qualitative description of the best found model by the estimation package. See also the field $key.properties for quantative descriptions. "
    if (evolmodel=="bm"){	
	chdesc<-paste(chdesc,"BM model. Brownian motion model, no stabilizing selection observed. ",sep="")
	if (numtraits==1){chdesc<-paste(chdesc,"Phenotype randomly oscillates ",sep="")}
	else{chdesc<-paste(chdesc,"Phenotypes randomly oscillate ",sep="")}
	chdesc<-paste(chdesc,"around the initial root state accumulating more and more variance with time.",sep="")
    }
    if ((evolmodel=="ouch") || (evolmodel=="mvslouch")){    
	if (numtraits==1){
	    if ((evolmodel=="ouch") && ballPosEig){chdesc<-paste(chdesc,"OU model. Single trait under stabilizing selection, adapting to a deterministic optimum. ",sep="")}
	    if ((evolmodel=="ouch") && !ballPosEig){chdesc<-paste(chdesc,"OU model. Single trait diverging. Results needs careful interpretation ",sep="")}
	    if ((evolmodel=="mvslouch") && ballPosEig){chdesc<-paste(chdesc,"OUBM model. Single trait under stabilizing selection, adapting to a randomly evolving optimum. ",sep="")}
	    if ((evolmodel=="mvslouch") && !ballPosEig){chdesc<-paste(chdesc,"OUBM model. Single trait diverging from a randomly evolving value. Results needs careful interpretation ",sep="")}
	}else{
	    if (!ballPosEig){
		if (evolmodel=="ouch") {chdesc<-paste(chdesc,"OUOU model. Multiple traits diverging. Results needs careful interpretation ",sep="")}
		if (evolmodel=="mvslouch") {chdesc<-paste(chdesc,"OUBM model. Multiple traits diverging from a randomly evolving value. Results needs careful interpretation ",sep="")}
	    }else{
		if (evolmodel=="ouch") {chdesc<-paste(chdesc,"OUOU model. Multiple traits under stabilizing selection, adapting to a deterministic optimum. ",sep="")}
		if (evolmodel=="mvslouch") {chdesc<-paste(chdesc,"OUBM model. Multiple traits under stabilizing selection, adapting to s randomly evolving optimum. ",sep="")}
		if ((model.setup$Atype== "Diagonal")|| (model.setup$Atype=="SingleValueDiagonal")){
		    if ((model.setup$Syytype== "Diagonal")|| (model.setup$Syytype=="SingleValueDiagonal")){		
			chdesc<-paste(chdesc,"All traits are evolving independently. They do not effect each others' adaptation to the primary optimum nor are their random perturbations correlated. ",sep="")
			if (model.setup$Atype=="SingleValueDiagonal"){chdesc<-paste(chdesc,"All traits have the same rate of adaptation to the primary optimum. ",sep="")}
			if (model.setup$Syytype=="SingleValueDiagonal"){chdesc<-paste(chdesc,"All traits experiance the same magnitude of random perturbations. ",sep="")}		    
		    }else{
			chdesc<-paste(chdesc,"All traits are adapting independently to the primary optimum. All dependencies between traits are due to the random perturbations. ",sep="")
			if (model.setup$Atype=="SingleValueDiagonal"){chdesc<-paste(chdesc,"All traits have the same rate of adaptation to the primary optimum. ",sep="")}		    
		    }		
		}
		else{
		    if ((model.setup$Atype== "UpperTri")|| (model.setup$Atype=="LowerTri")){
			chdesc<-paste(chdesc,"We have a clear causality pattern in traits' adaptation and primary optimum. ",sep="")
			if (model.setup$Atype== "UpperTri"){
			    chdesc<-paste(chdesc,"The ''bottom'' trait in the A matrix is effecting the primary optimum of all other traits. The second last one all the traits' except the bottom one's and so on. The ''top'' trait's primary optimum is effected by all other traits. ",sep="")
			}else{
			    chdesc<-paste(chdesc,"The ''top'' trait in the A matrix is effecting the primary optimum of all other traits. The second one all the traits' except the top one's and so on. The ''bottom'' trait's primary optimum is effected by all other traits. ",sep="")
			}
		    }else{
			chdesc<-paste(chdesc,"There is no clear causality pattern in traits' adaptation. All traits seem to effect each other's primary optimum. However look in the A matrix for 0 or very small values. These will indicate bivariate causality effects. ",sep="")
		    }
		    if ((model.setup$Syytype== "Diagonal")|| (model.setup$Syytype=="SingleValueDiagonal")){		
			chdesc<-paste(chdesc,"The random perturbations are not correlated between different traits. All dependencies between the traits come from their co-adaptation. ",sep="")
			if (model.setup$Syytype=="SingleValueDiagonal"){chdesc<-paste(chdesc,"All traits experiance the same magnitude of random perturbations. ",sep="")}		    
		    }else{
			chdesc<-paste(chdesc,"The random perturbations are correlated between different traits. The traits are dependent both through their co-adaptation and correlated random perturbations. ",sep="")		    
		    }
		}
	    }
	}
    }    
    chdesc
}

.generate.model.se<-function(evolmodel){
    chsedesc<-"The function did not calculate any confidence intervals as this takes a long time. You have to do it manually. Here is code to call: "
    if (evolmodel=="bm"){
        chsedesc<-paste(chsedesc,'parametric.bootstrap(estimated.model=your_estimated_model,phyltree=your_phylogenetic_tree,values.to.bootstrap=c(vector of parameter names, e.g.: "Sxx", "StS"),numboot=desired number of bootstrap replicates, and see other technical parameters in manual) ',sep='')
    }
    if (evolmodel=="ouch"){
	chsedesc<-paste(chsedesc,',parametric.bootstrap(estimated.model=your_estimated_model,phyltree=your_phylogenetic_tree,values.to.bootstrap=c(vector of parameter names, e.g.: "evolutionary.regression", "optimal.regression"),regimes=your_regimes,root.regime=your_root_regime (does not need to be provided),M.error=your_measurement_error,predictors=your_predictors (if you want conditional distributions reported),numboot=desired number of bootstrap replicates,Atype=your_Atype,Syytype=your_Syytype,diagA=your_diagonal_of_A,parameter_signs=your_parameter_signs, and see other technical parameters in manual)', sep='')
    }
    if (evolmodel=="mvslouch"){
	chsedesc<-paste(chsedesc,'parametric.bootstrap(estimated.model=your_estimated_model,phyltree=your_phylogenetic_tree,values.to.bootstrap=c(vector of parameter names, e.g.: "evolutionary.regression", "optimal.regression"),regimes=your_regimes,root.regime=your_root_regime (does not need to be provided),M.error=your_measurement_error,predictors=your_predictors (if you want different than default (conditional on BM variables) conditional distributions reported),kY=number of response (i.e. OU) variables,numboot=desired number of bootstrap replicates,Atype=your_Atype,Syytype=your_Syytype,diagA=your_diagonal_of_A,parameter_signs=your_parameter_signs, and see other technical parameters in manual)', sep='')
    }
    chsedesc<-paste(chsedesc,"where  your_estimated_model_parameters are the estimated parameters 
    (can be found in the slot $BestModel$BestModel$ParamsInModel in the returned object), 
    your_phylogenetic_tree is your phylogeny in ape (phylo) format (passed to the parameter phyltree 
    when calling evol.model.est), your_measurement_error is your measurement error structure 
    (passed to the parameter M.error when calling evol.model.est, may be left NULL), 
    your_predictors are your indicated predictor variables (passed to the parameter predictors when calling 
    evol.model.est, may be left NULL)",sep="")    
    if ((evolmodel=="ouch") || (evolmodel=="mvslouch")){chsedesc<-paste(chsedesc,", 
    your_regimes is your regimes vector (passed to the parameter regimes when calling evol.model.est), 
    your_Atype is the class to which the A matrix belongs in the best found model 
    (can be found in the slot $BestModel$model$Atype in the returned object), 
    your_Syytype is the class to which the A matrix belongs in the best found model 
    (can be found in the slot $BestModel$model$Syytype in the returned object) and your_parameter_signs are the signs of parameters
    you assumed (if you did, read the manual files on this first!). ", sep="")
    }else{chsedesc<-paste(chsedesc,".",sep="")}
    chsedesc
}

.extract.model.key.properties<-function(evolmodel,estimated.model,k=NULL){
    lkey.properties<-list()    
    if (evolmodel=="bm"){lkey.properties$evolution.model<-"Brownian motion model"}
    if ((evolmodel=="ouch")||(evolmodel=="mvslouch")){
	if ((evolmodel=="ouch") && (k==1)){lkey.properties$evolution.model<-"Ornstein-Uhlenbeck model"}
	if ((evolmodel=="ouch") && (k>1)){lkey.properties$evolution.model<-"Ornstein-Uhlenbeck-Ornstein-Uhlenbeck model"}
	if (evolmodel=="mvslouch"){lkey.properties$evolution.model<-"Ornstein-Uhlenbeck-Brownian motion model"}
	lkey.properties$phylogenenetic.halflives<-list()
	lkey.properties$phylogenenetic.halflives$halflives<-estimated.model$ParamSummary$phyl.halflife$halflives["halflife",]
	lkey.properties$phylogenenetic.halflives$comment<-"Phylogenetic half-lives in the eigenvector (NOT trait) directions. "
	if (length(which(Re(estimated.model$ParamSummary$phyl.halflife$halflives["halflife",])>0))== length(estimated.model$ParamSummary$phyl.halflife$halflives["halflife",])){
	    lkey.properties$phylogenenetic.halflives$comment<-paste(lkey.properties$phylogenenetic.halflives$comment,"All halflives are positive so all traits are under stabilizing selection towards the optimum. See Bartoszek et. al. (2012) for details.",sep="")
	}else{
	    lkey.properties$phylogenenetic.halflives$comment<-paste(lkey.properties$phylogenenetic.halflives$comment,"There are negative halflives so all traits may be diverging. One needs to interpret the results very carefully. See Bartoszek et. al. (2012) for details.",sep="")
	}
	if (!is.null(estimated.model$ParamSummary$evolutionary.regression)){
	    lkey.properties$evolutionary.regression<-list()
	    lkey.properties$evolutionary.regression$regression.coefficents<-estimated.model$ParamSummary$evolutionary.regression
	    lkey.properties$evolutionary.regression$comment<-"The evolutionary regression between traits and the primary optimum for the suite of traits. This is the currently observed relationship. See Bartoszek et. al. (2012) and the help for the mvSLOUCH::mvslouchModel and mvSLOUCH::ouchModel functions for details."
	}
	if (!is.null(estimated.model$ParamSummary$optimal.regression)){
	    lkey.properties$optimal.regression<-list()
	    lkey.properties$optimal.regression$regression.coefficents<-estimated.model$ParamSummary$optimal.regression
	    lkey.properties$optimal.regression$comment<-"The optimal regression between traits and the primary optimum for the suite of traits. This is the currently optimal relationship that the traits want to attain. This however is not always possible, see the discussion in Bartoszek et. al.(2012)."
	}
	if (!is.null(estimated.model$ParamSummary$trait.regression)){
	    lkey.properties$trait.regression<-list()
	    lkey.properties$trait.regression$regression.coefficents<-estimated.model$ParamSummary$trait.regression
	    lkey.properties$trait.regression$comment<-"The regression of each trait on all of the others. See Bartoszek et. al. (2012) and the help for the mvSLOUCH::mvslouchModel and mvSLOUCH::ouchModel functions for details."
	}
	if (!is.null(estimated.model$ParamSummary$corr.matrix)){
	    lkey.properties$corr.matrix<-list()
	    lkey.properties$corr.matrix$correlation.matrix<-estimated.model$ParamSummary$corr.matrix
	    lkey.properties$corr.matrix$comment<-"The currently observed correlation between the traits. See Bartoszek et. al. (2012) and the help for the mvSLOUCH::mvslouchModel and mvSLOUCH::ouchModel functions for details."
	}
	if (!is.null(estimated.model$ParamSummary$stationary.corr.matrix)){
	    lkey.properties$stationary.corr.matrix<-list()
	    lkey.properties$stationary.corr.matrix$correlation.matrix<-estimated.model$ParamSummary$stationary.corr.matrix
	    lkey.properties$stationary.corr.matrix$comment<-"The limiting/stationary correlation between the traits. See Bartoszek et. al. (2012) and the help for the mvSLOUCH::mvslouchModel and mvSLOUCH::ouchModel functions for details."
	}
	lkey.properties$R2<-"Not implemented yet"
	lkey.properties$R2_conditional_on_predictors<-"Not implemented yet"
	if (!is.null(estimated.model$ParamSummary$RSS$R2)){lkey.properties$R2<-estimated.model$ParamSummary$RSS$R2}
	if (!is.null(estimated.model$ParamSummary$RSS$R2_conditional_on_predictors)){
	    lkey.properties$R2_conditional_on_predictors<-estimated.model$ParamSummary$RSS$R2_conditional_on_predictors;    
	    if (lkey.properties$R2=="Not implemented yet"){lkey.properties$R2<-"Not implemented yet"}
	}
	if (!is.null(estimated.model$ParamSummary$RSS$R2_non_phylogenetic_conditional_on_predictors)){
	    lkey.properties$R2_non_phylogenetic_conditional_on_predictors<-estimated.model$ParamSummary$RSS$R2_non_phylogenetic_conditional_on_predictors
	    if (lkey.properties$R2=="Not implemented yet"){lkey.properties$R2<-"Not implemented yet"}	    
	    if (lkey.properties$R2_conditional_on_predictors=="Not implemented yet"){lkey.properties$R2_conditional_on_predictors<-"Not implemented yet"}
	}	
    }
    lkey.properties$aic.c<- estimated.model$BestModel$ParamSummary$aic.c    
    lkey.properties$comment<-"See the slot $BestModel$BestModel$ParamSummary in the returned object for various other hopefully helpful summary statistics and information criteria"
    lkey.properties
}

.changeSigmatoSyy<-function(Sigma,Syytype,diagSyy,signsSyy,bcorrect0var=FALSE){
## function called in: estimBM.R, evolmodelest.R, modelparams.R, modelparamstransform.R, regimes.R
    if ((bcorrect0var) && (!.matrixcalc_is.diagonal.matrix(Sigma))){
    ## the rational is that if the variance is 0 the data is a constant, hence all covariances also should be 0
        v0var<-which(sapply(diag(Sigma),function(x){isTRUE(all.equal(x,0))},simplify=TRUE))
        bcorrect0var<-FALSE
        if (length(v0var)>0){## there are 0 variances
	    if (length(v0var)==nrow(Sigma)){Sigma<-matrix(0,nrow=nrow(Sigma),ncol=ncol(Sigma));v0var<-c()}
	    else{
	    	vNon0var<-sort(setdiff(1:nrow(Sigma),v0var))
		Sigma<-Sigma[-v0var,-v0var,drop=FALSE]
		bcorrect0var<-TRUE
	    }
        }
    }else{bcorrect0var<-FALSE}
    if (isTRUE(all.equal(sum(abs(Sigma)),0))){Syy<-matrix(0,nrow=nrow(Sigma),ncol=ncol(Sigma))}
    else{
        Syy<-Sigma
	if (.matrixcalc_is.diagonal.matrix(Sigma)){Syy<-diag(sqrt(diag(Sigma)),nrow=nrow(Sigma),ncol=nrow(Sigma))}
	else{   
    	    Syy=switch(Syytype,
        	SingleValueDiagonal={diag(sqrt(mean(diag(Sigma))),nrow=nrow(Sigma),ncol=ncol(Sigma))},
        	Diagonal={diag(sqrt(diag(Sigma)),nrow=nrow(Sigma),ncol=ncol(Sigma))},
        	Symmetric={
        	    eigS<-eigen(Sigma)
        	    eigS$vectors%*%diag(sqrt(eigS$values),nrow=nrow(Sigma),ncol=ncol(Sigma)) %*%t(eigS$vectors)
        	},
        	UpperTri={
		    ## factorization procedure taken from 
		    ## https://math.stackexchange.com/questions/2039477/cholesky-decompostion-upper-triangular-or-lower-triangular
		    k<-nrow(Sigma)
        	    P<-matrix(0,nrow=k,ncol=k);for (i in 1:k){P[k-i+1,i]<-1}## create permutation matrix with 1s on the anti-diagonal
        	    P%*%t(.my_chol(P%*%Sigma%*%P))%*%P 
        	},
        	LowerTri={
        	    t(.my_chol(Sigma))        	    
        	},
        	Any={Sigma},
        	.my_stop('Incorrect type for Syy provided! Admissable types are "SingleValueDiagonal", "Diagonal", "UpperTri", "LowerTri", "Symmetric", "Any".',TRUE)
        	)
    
	}    
	if (!is.null(diagSyy)){
    	    diag(Syy)=switch(diagSyy,
        	Positive={exp(diag(Syy))},
        	Negative={(-1)*exp(diag(Syy))}
    	    )
	}
	if (!is.null(signsSyy)){ 
    	    Syy[which(signsSyy=="-")]<- (-1)*exp(Syy[which(signsSyy=="-")])
    	    Syy[which(signsSyy=="+")]<- exp(Syy[which(signsSyy=="+")])
    	    signsSyy[which(signsSyy=="-")]<-NA
    	    signsSyy[which(signsSyy=="+")]<-NA
    	    if (lengths(which(!is.na(signsSyy)))>0){Syy[which(!is.na(signsSyy))]<-signsSyy[which(!is.na(signsSyy))]}
	}
    }    
    if (bcorrect0var){
	Syy<-cbind(Syy,matrix(0,nrow=nrow(Syy),ncol=length(v0var)))
	Syy<-Syy[,order(c(vNon0var,v0var)),drop=FALSE]
	Syy<-rbind(Syy,matrix(0,ncol=ncol(Syy),nrow=length(v0var)))
	Syy<-Syy[order(c(vNon0var,v0var)),,drop=FALSE]
    }
    Syy
}

.my_chol<-function(M){
## function called in estimBM.R evolmodelest.R matrixparametrizations.R
## M has to be symmetric--semi--positive--definte
## no check done for this done as the function is only called internally
## depending on platform two different types of calculations
## Debian with pivoting
## others without
## a different treatment takes place for Debian as the Debian flavor on CRAN
## raises an error to calls to chol(), noticed since 2019 XI 27
## other platforms do not do this
##    mchol<-matrix(NA,nrow=nrow(M),ncol=ncol(M))
##    platform<-R.Version()$platform 
##    if (grepl("deb", platform, ignore.case = TRUE)){## we are on Debian and some error seems to take place, so we pivot
##    	    org_warn<-options("warn")
##            options(warn=-1)
##            mchol<-chol(M,pivot=TRUE)
##            options(warn=org_warn$warn)
##    }else{## not Debian
##	mchol<-chol(M)
##    }
    mchol<-matrix(NA,nrow=nrow(M),ncol=ncol(M)) ## M has to be square
    tryCatch({mchol<-chol(M)},error=function(e){.my_message(paste("Caught:",e),FALSE);.my_message("\n",FALSE)})##,warning=function(e){.my_message(paste("Caught:",e),FALSE);.my_message("\n",FALSE)})
    if (is.na(mchol[1,1])){
    ## an error took place, most probably due to chol() considering M not to be of full rank
    ## therefore setting pivot=TRUE, to deal with such a matrix
    ## chol() throws a warning if M does not have full rank so we suppress it
	org_warn<-options("warn")
	options(warn=-1)
	tryCatch({mchol<-chol(M,pivot=TRUE)},error=function(e){.my_message(paste("Caught:",e),FALSE);.my_message("\n",FALSE)})##,warning=function(e){.my_message(paste("Caught:",e),FALSE);.my_message("\n",FALSE)})
	options(warn=org_warn$warn)    
        if (!is.null(attr(mchol,"pivot"))){attr(mchol,"pivot")<-NULL}
        if (!is.null(attr(mchol,"rank"))){attr(mchol,"rank")<-NULL}
    }
    mchol
}

.createStartPointsASyyB<-function(mData,tree_height,model_setup,kY,bDoB=FALSE){
	## Starting point motivated by results from 
	## K. Bartoszek and S. Sagitov. "Phylogenetic confidence intervals for the optimal trait value". Journal of Applied Probability 52.4 (2015), pp. 1115-1132.
    RatioForLambdas<-0.5
#    LambdaStart<-optim(1,function(x,tree_height){abs((1-exp(-2*x*tree_height))/(2*x)-1)},method="BFGS",tree_height=tree_height)$value
    LambdaStart<-optim(1,function(x,tree_height,RatioForLambdas){abs((1-exp(-2*x*tree_height)/(2*x))-RatioForLambdas)},method="BFGS",tree_height=tree_height,RatioForLambdas=RatioForLambdas)$value
    diagSyy<-NULL;signsSyy<-NULL
    if ((is.element("diagA",names(model_setup)))&&(!is.null(model_setup$diagA))&&(model_setup$Atype!="SymmetricPositiveDefinite")){
	        LambdaStart=switch(model_setup$diagA,
        	    Positive={exp(LambdaStart)},
        	    Negative={(-1)*exp(LambdaStart)},
        	    NULL=LambdaStart
    		)
    }
    
    if (is.element("diagSyy",names(model_setup))){diagSyy<-model_setup$diagSyy}
    if (is.element("signsSyy",names(model_setup$parameter_signs))){signsSyy<-model_setup$parameter_signs$signsSyy}		    
    Sigma<-(cov(mData,use="pairwise.complete.obs")[1:kY,1:kY])/RatioForLambdas
    vNASigma<-which(is.na(Sigma))
    if (length(vNASigma)>0){
        Sigma[vNASigma]<-runif(length(vNASigma))/10
    }
    
    LambdaStart<-diag(LambdaStart,ncol=kY,nrow=kY)
    if (is.element("signsA",names(model_setup$parameter_signs))){
	    signsA<-model_setup$parameter_signs$signsA
	    LambdaStart[which(signsA=="-")]<- (-1)*exp(LambdaStart[which(signsA=="-")])
	    LambdaStart[which(signsA=="+")]<- exp(LambdaStart[which(signsA=="+")])
    	    signsA[which(signsA=="-")]<-NA
    	    signsA[which(signsA=="+")]<-NA
    	    if (length(which(!is.na(signsA)))>0){LambdaStart[which(!is.na(signsA))]<-as.numeric(signsA[which(!is.na(signsA))])}	    
    }
    eigL<-eigen(LambdaStart)
    mP<-eigL$vectors
    maxTriesMP<-10
    while(isTRUE(all.equal(rcond(mP),0))&&(maxTriesMP>0)){mP<-jitter(mP);maxTriesMP<-maxTriesMP-1}
    if (isTRUE(all.equal(rcond(mP),0))){mP<-diag(1,nrow=nrow(mP),ncol=ncol(mP))}
    invmP<-solve(mP)
    vL<-eigL$values
    Sigma<-invmP%*%Sigma%*%t(invmP)
    mLtr<-matrix(0:(kY^2-1),nrow=kY,ncol=kY,byrow=TRUE)
    mLtr<-apply(mLtr,c(1,2),.CalcVlq,"vlambda"=vL,"t"=tree_height,"k"=kY)
    Sigma<-Sigma/mLtr
    Sigma<-mP%*%Sigma%*%t(mP)
    
	    
    ## resulting Sigma might not be sym-pos-def due to either "pairwise.complete.obs" if there are NA values, see ?cov
    ## or as a result of the performed transformation, but in this case is according to the formula for the covariance
    ## so should be fine, unless A taken as defective
    if ((!.matrixcalc_is.symmetric.matrix(Sigma)) || (!.matrixcalc_is.positive.definite(Sigma))){
	Sigma<-as.matrix(Matrix::nearPD(Sigma)$mat)
    }
    SyyStartPoint<-.changeSigmatoSyy(Sigma[1:kY,1:kY,drop=FALSE],model_setup$Syytype,diagSyy,signsSyy,FALSE)		    
    if (bDoB){
    	BStartPoint<-matrix(0,nrow=kY,ncol=ncol(mData)-kY)
	
	if (is.element("signsB",names(model_setup$parameter_signs))){
	    signsB<-model_setup$parameter_signs$signsB
    	    signsB[which(signsB=="-")]<-NA
    	    signsB[which(signsB=="+")]<-NA
    	    if (lengths(which(!is.na(signsB)))>0){BStartPoint[which(!is.na(signsB))]<-as.numeric(signsB[which(!is.na(signsB))])}	    
	}
	
	for (i in 1:kY){
	    vY<-mData[,i]-apply(mData[,(kY+1):ncol(mData),drop=FALSE],1,function(x,y){x%*%y},y=BStartPoint[i,])
	    vBols<-stats::lm(vY~mData[,(kY+1):ncol(mData)])$coefficients[-1]
	    names(vBols)<-NULL
	    if (is.element("signsB",names(model_setup$parameter_signs))){ BStartPoint[i,which(is.na(signsB[i,]))]<-vBols}
	    else{BStartPoint[i,]<-vBols}
	}
	
	if (is.element("signsB",names(model_setup$parameter_signs))){
	    signsB<-model_setup$parameter_signs$signsB
	    BStartPoint[which(signsB=="-")]<- sapply(BStartPoint[which(signsB=="-")],function(x){ifelse(x<0,x,0)})
	    BStartPoint[which(signsB=="+")]<- sapply(BStartPoint[which(signsB=="+")],function(x){ifelse(x>0,x,0)})
	}
	BStartPoint<-(-1)*LambdaStart%*%BStartPoint
##	lres<-list("A"=matrix(LambdaStart,kY,kY),"Syy"=SyyStartPoint,"B"=BStartPoint)
	lres<-list("A"=LambdaStart,"Syy"=SyyStartPoint,"B"=BStartPoint)
##    }else{lres<-list("A"=matrix(LambdaStart,kY,kY),"Syy"=SyyStartPoint)}
    }else{lres<-list("A"=LambdaStart,"Syy"=SyyStartPoint)}
    lres
}
