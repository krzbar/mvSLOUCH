## This file is part of mvSLOUCH

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this package or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .



BrownianMotionModel<-function(phyltree,mData,predictors=NULL,M.error=NULL,min_bl=0.0003){
    .internal_BrownianMotionModel(phyltree,mData,predictors,M.error,min_bl)
}


.internal_BrownianMotionModel<-function(phyltree,mData,predictors=NULL,M.error=NULL,min_bl=0.0003){
    EvolModel<-"bm"

    setres<-.set_sing_bl_threshold(min_bl)
    pcmbase_min_bl<-setres$pcmbase_min_bl
    pcmbase_skip<-setres$pcmbase_skip
##    phyltree$edge.length[which(phyltree$edge.length<min_bl)]<-min_bl

    res<-.PhyloSDEestim(phyltree,mData,kY=NULL,regimes=NULL,regimes.times=NULL,root.regime=NULL,params=list(EvolModel=EvolModel,EstimParams=list(calcCI=FALSE)),predictors=predictors,M.error=M.error)
    options(PCMBase.Threshold.Skip.Singular=pcmbase_min_bl)
    options(PCMBase.Skip.Singular=pcmbase_skip)
    res
}


.internal_ouchModel<-function(phyltree,mData,regimes=NULL,regimes.times=NULL,root.regime=NULL,predictors=NULL,M.error=NULL,Atype="Invertible",Syytype="UpperTri",diagA="Positive",estimate.root.state=FALSE,parameter_signs=NULL,lStartPoint=NULL,parscale=NULL,min_bl=0.0003,maxiter=c(10,100),b_warnInvertible=FALSE,b_warnAdiag=FALSE,b_warnAsymposdef=FALSE){
    if ((!is.vector(maxiter))|| (!is.numeric(maxiter)) || (length(maxiter)!=2)){
        maxiter<-c(10,100);
        .my_warning("WARNING: maxiter passed in a wrong way, setting it to default of c(10,100)",TRUE,FALSE)
    }

    setres<-.set_sing_bl_threshold(min_bl)
    pcmbase_min_bl<-setres$pcmbase_min_bl
    pcmbase_skip<-setres$pcmbase_skip
##    phyltree$edge.length[which(phyltree$edge.length<min_bl)]<-min_bl

    res<-.PhyloSDEestim(phyltree,mData,kY=NULL,regimes=regimes,regimes.times=regimes.times,root.regime=root.regime,params= list(EvolModel="ouch",EstimParams=c(list(Atype=Atype,Syytype=Syytype,diagA=diagA,diagSyy="Positive",calcCI=FALSE,lStartPoint=lStartPoint,optim_parscale=parscale),parameter_signs)),predictors=predictors,M.error=M.error,estimate.root.state=estimate.root.state,maxiter=c(maxiter[1],1,maxiter[2]))
    options(PCMBase.Threshold.Skip.Singular=pcmbase_min_bl)
    options(PCMBase.Skip.Singular=pcmbase_skip)
    .f_warningsA(b_warnInvertible=b_warnInvertible,b_warnAdiag=b_warnAdiag,Atype=Atype,diagA=diagA,b_warnAsymposdef=b_warnAsymposdef)
    res
}


ouchModel<-function(phyltree,mData,regimes=NULL,regimes.times=NULL,root.regime=NULL,predictors=NULL,M.error=NULL,Atype="Invertible",Syytype="UpperTri",diagA="Positive",estimate.root.state=FALSE,parameter_signs=NULL,start_point_for_optim=NULL,parscale=NULL,min_bl=0.0003,maxiter=c(10,100)){
    .internal_ouchModel(phyltree,mData,regimes,regimes.times,root.regime,predictors,M.error,Atype,Syytype,diagA,estimate.root.state,parameter_signs,start_point_for_optim,parscale,min_bl,maxiter=maxiter,b_warnInvertible=TRUE,b_warnAdiag=FALSE,b_warnAsymposdef=TRUE)
}


mvslouchModel<-function(phyltree,mData,kY,regimes=NULL,regimes.times=NULL,root.regime=NULL,predictors=NULL,M.error=NULL,Atype="Invertible",Syytype="UpperTri",diagA="Positive",estimate.root.state=FALSE,parameter_signs=NULL,start_point_for_optim=NULL,parscale=NULL,min_bl=0.0003,maxiter=c(10,50,100),estimateBmethod="ML"){
    .internal_mvslouchModel(phyltree,mData,kY,regimes,regimes.times,root.regime,predictors,M.error,Atype,Syytype,diagA,estimate.root.state,parameter_signs,start_point_for_optim,parscale,min_bl,maxiter=maxiter,estimateBmethod=estimateBmethod,b_warnInvertible=TRUE,b_warnAdiag=TRUE)
}


.internal_mvslouchModel<-function(phyltree,mData,kY,regimes=NULL,regimes.times=NULL,root.regime=NULL,predictors=NULL,M.error=NULL,Atype="Invertible",Syytype="UpperTri",diagA="Positive",estimate.root.state=FALSE,parameter_signs=NULL,lStartPoint=NULL,parscale=NULL,min_bl=0.0003,maxiter=c(10,50,100),estimateBmethod="ML",b_warnInvertible=FALSE,b_warnAdiag=FALSE){
    if (is.matrix(mData)){
	if (kY>=ncol(mData)){
	    .my_stop('Cannot have all variables as OU responses, i.e. kY>=ncol(mData). Set kY to a smaller (than number of columns in mData) number of OU responses.',TRUE)
	}
    }else{
	.my_stop('mData has to be a matrix!',TRUE)
    }

    if ((!is.vector(maxiter))|| (!is.numeric(maxiter)) || (length(maxiter)!=3)){
        maxiter<-c(10,50,100);
        .my_warning("WARNING: maxiter passed in a wrong way, setting it to default of c(10,50,100)",TRUE,FALSE)
    }

    EvolModel<-"mvslouch"

    setres<-.set_sing_bl_threshold(min_bl)
    pcmbase_min_bl<-setres$pcmbase_min_bl
    pcmbase_skip<-setres$pcmbase_skip
#    phyltree$edge.length[which(phyltree$edge.length<min_bl)]<-min_bl

    res<-.PhyloSDEestim(phyltree,mData,kY=kY,regimes=regimes,regimes.times=regimes.times,root.regime=root.regime,params=list(EvolModel=EvolModel,EstimParams=c(list(Atype=Atype,Syytype=Syytype,diagA=diagA,diagSyy="Positive",calcCI=FALSE,lStartPoint=lStartPoint,optim_parscale=parscale),parameter_signs)),predictors=predictors,M.error=M.error,estimate.root.state=estimate.root.state,maxiter=maxiter,cBestim_method=estimateBmethod)
    options(PCMBase.Threshold.Skip.Singular=pcmbase_min_bl)
    options(PCMBase.Skip.Singular=pcmbase_skip)
    .f_warningsA(b_warnInvertible=b_warnInvertible,b_warnAdiag=b_warnAdiag,Atype=Atype,diagA=diagA,b_warnAsymposdef=FALSE)
    res
}

.f_warningsA<-function(b_warnInvertible,b_warnAdiag,Atype,diagA,b_warnAsymposdef){
    if ((b_warnInvertible)&&(Atype=="Invertible")){
	.my_message('Atype is at the default "Invertible" setting. This is a highly inefficient and unrecommended setting. Please look into the possible constraints on A and choose the one best corresponding to the hypothesis on the relationship between the traits.')
    }
    if ((b_warnAdiag)&&(!is.null(diagA))&&(Atype%in%c("SymmetricPositiveDefinite", "DecomposablePositive", "DecomposableNegative", "DecomposableReal"))){ ##"Invertible"
	c_warnmessage<-""
	if (Atype=="SymmetricPositiveDefinite"){
	    c_warnmessage<-"A is defined as symmetric postive definite. Its diagonal is guaranteed to be positive. Using the diag A setting will result in the diagonal being exponentiated above this and could be in some cases detrimental to the optimization process. It is recommended to run with diagA=NULL."
	}else{
	    c_warnmessage<-"A is parametrized through its "
	    if (Atype%in%c("DecomposablePositive", "DecomposableNegative", "DecomposableReal")){c_warnmessage<-paste0(c_warnmessage,"eigen")}
	    ##if (Atype=="Invertible"){c_warnmessage<-paste0(c_warnmessage,"QR ")}
	    c_warnmessage<-paste0(c_warnmessage,"decomposition.	Using the diag A setting will result in the diagonal being exponentiated above this and could be in some cases (but performed simulations on this are not conclusive) slightly detrimental to the optimization process. It is recommended to experiment with diagA=NULL.")
	}
	.my_message(c_warnmessage)	
    }
    if ((!b_warnAdiag)&&(b_warnAsymposdef)&&(!is.null(diagA))&&(Atype=="SymmetricPositiveDefinite")){
	c_warnmessage<-"A is defined as symmetric postive definite. Its diagonal is guaranteed to be positive. Using the diag A setting will result in the diagonal being exponentiated above this and could be in some cases detrimental to the optimization process. It is recommended to run with diagA=NULL."
	.my_message(c_warnmessage)	
    }

    NULL
}

SummarizeBM<-function(phyltree,mData, modelParams,t=c(1),dof=NULL,M.error=NULL,predictors=NULL,min_bl=0.0003){
    bShouldPrint<-FALSE
    EvolModel<-"bm"

    setres<-.set_sing_bl_threshold(min_bl)
    pcmbase_min_bl<-setres$pcmbase_min_bl
    pcmbase_skip<-setres$pcmbase_skip
##    phyltree$edge.length[which(phyltree$edge.length<min_bl)]<-min_bl

    res<-.SummarizeFullPoint(NULL,mData=mData,PhylTree=phyltree,EvolModel=EvolModel,EstimationParams=list(calcCI=FALSE),regimes.times=NULL,regimes=NULL,modelParams=modelParams,t=t,dof=dof,bShouldPrint=bShouldPrint,LogLik=NULL,maxIter=NULL,tol=NULL,Merror=M.error,predictors=predictors)
    options(PCMBase.Threshold.Skip.Singular=pcmbase_min_bl)
    options(PCMBase.Skip.Singular=pcmbase_skip)
    res
}

SummarizeOUCH<-function(phyltree,mData,modelParams,regimes=NULL,regimes.times=NULL,t=c(1),dof=NULL,M.error=NULL,predictors=NULL,Atype="Invertible",Syytype="UpperTri",min_bl=0.0003){
## in summary there is no point in having a root regime as either it sits in the initial state or will have no effect

    bShouldPrint<-FALSE
    
    setres<-.set_sing_bl_threshold(min_bl)
    pcmbase_min_bl<-setres$pcmbase_min_bl
    pcmbase_skip<-setres$pcmbase_skip
##    phyltree$edge.length[which(phyltree$edge.length<min_bl)]<-min_bl

    res<-.SummarizeFullPoint(NULL,mData=mData,PhylTree=phyltree,EvolModel="ouch",EstimationParams=list(calcCI=FALSE,Atype=Atype,Syytype=Syytype),regimes.times=regimes.times,regimes=regimes,modelParams=modelParams,t=t,dof=dof,bShouldPrint=bShouldPrint,LogLik=NULL,maxIter=c(10,1,100),tol=c(0.01,0.01),Merror=M.error,predictors=predictors)
    options(PCMBase.Threshold.Skip.Singular=pcmbase_min_bl)
    options(PCMBase.Skip.Singular=pcmbase_skip)
    res
}

SummarizeMVSLOUCH<-function(phyltree,mData,modelParams,regimes=NULL,regimes.times=NULL,t=c(1),dof=NULL,M.error=NULL,predictors=NULL,Atype="Invertible",Syytype="UpperTri",min_bl=0.0003,maxiter=50){
## in summary there is no point in having a root regime as either it sits in the initial state or will have no effect

    bShouldPrint<-FALSE
    EvolModel<-"mvslouch"
    setres<-.set_sing_bl_threshold(min_bl)
    pcmbase_min_bl<-setres$pcmbase_min_bl
    pcmbase_skip<-setres$pcmbase_skip
##    phyltree$edge.length[which(phyltree$edge.length<min_bl)]<-min_bl

    res<-.SummarizeFullPoint(NULL,mData=mData,PhylTree=phyltree,EvolModel=EvolModel,EstimationParams=list(calcCI=FALSE,Atype=Atype,Syytype=Syytype),regimes.times=regimes.times,regimes=regimes,modelParams=modelParams,t=t,dof=dof,bShouldPrint=bShouldPrint,LogLik=NULL,maxIter=c(10,maxiter,100),tol=c(0.01,0.01),Merror=M.error,predictors=predictors)
    options(PCMBase.Threshold.Skip.Singular=pcmbase_min_bl)
    options(PCMBase.Skip.Singular=pcmbase_skip)
    res
}

.set_sing_bl_threshold<-function(min_bl=0.0003){
    pcmbase_min_bl<-PCMBase::PCMOptions()$PCMBase.Threshold.Skip.Singular
    pcmbase_skip<-PCMBase::PCMOptions()$PCMBase.Skip.Singular    
    new_min_bl_use<-max(c(pcmbase_min_bl,min_bl))
    options(PCMBase.Threshold.Skip.Singular=new_min_bl_use)
    options(PCMBase.Skip.Singular=TRUE)
    list(pcmbase_min_bl=pcmbase_min_bl,pcmbase_skip=pcmbase_skip)
}

simulBMProcPhylTree<-function(phyltree,X0=NULL,Sigma=NULL,dropInternal=TRUE,M.error=NULL,fullTrajectory=FALSE,jumpsetup=NULL,keep_tree=FALSE){
## if regimes is not null then X0 has to be a matrix 
    regimes<-NULL
    regimes.times<-NULL
    modelParams<-NULL
    evolmodel<-"bm"
    if (fullTrajectory){
	evolmodel<-"bmStep"
    }
    Simulparams<-NULL
    if (!is.null(jumpsetup)){Simulparams$jump<-jumpsetup}
    if (is.null(modelParams)){
	if (is.null(X0)||is.null(Sigma)){.my_stop("BM parameters are NULL",TRUE)}
	else{modelParams<-list(vX0=X0,Sxx=Sigma)}
    }
    res<-.simulVasicekProcPhylTree(phyltree,EvolModel=evolmodel,modelParams=modelParams,EstimationParams=NULL,regimes=regimes,regimes.times=regimes.times,Simulparams=Simulparams,dropInternal=dropInternal,M.error=M.error,bAllTrajectories=fullTrajectory,bKeep_tree=keep_tree)
    if (!fullTrajectory){res<-res$ExtantSample}
    res
}

simulOUCHProcPhylTree<-function(phyltree,modelParams,regimes=NULL,regimes.times=NULL,dropInternal=TRUE,M.error=NULL,fullTrajectory=FALSE,jumpsetup=NULL,keep_tree=FALSE){
    evolmodel<-"ouch"
    if (fullTrajectory){evolmodel<-"ouchStep"}
    Simulparams<-NULL
    if (!is.null(jumpsetup)){Simulparams$jump<-jumpsetup}
    res<-.simulVasicekProcPhylTree(phyltree,evolmodel,modelParams=modelParams,EstimationParams=NULL,regimes=regimes,regimes.times=regimes.times,Simulparams=Simulparams,dropInternal=dropInternal,M.error=M.error,bAllTrajectories=fullTrajectory,bKeep_tree=keep_tree)
    if (!fullTrajectory){res<-res$ExtantSample}
    res
}

simulMVSLOUCHProcPhylTree<-function(phyltree,modelParams,regimes=NULL,regimes.times=NULL,dropInternal=TRUE, M.error=NULL,fullTrajectory=FALSE,jumpsetup=NULL,keep_tree=FALSE){
    evolmodel<-"mvslouch"
    if (fullTrajectory){
	evolmodel<-"mvslouchStep"
    }
    Simulparams<-NULL
    if (!is.null(jumpsetup)){Simulparams$jump<-jumpsetup}
    res<-.simulVasicekProcPhylTree(phyltree,evolmodel,modelParams=modelParams,EstimationParams=NULL,regimes=regimes,regimes.times=regimes.times,Simulparams=Simulparams,dropInternal=dropInternal,M.error=M.error,bAllTrajectories=fullTrajectory,bKeep_tree=keep_tree)
    if (!fullTrajectory){res<-res$ExtantSample}
    res
}

