## This file is part of mvSLOUCH

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this package or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .

.Params.summary<-function(phyltree,modelParams,EvolModel,designToEstim,mData=NULL,t=1,LogLik=-Inf,npar0=0,RSS=NA,KnownParams=NULL,conf.level=0.95,vVars=NULL,conditional=FALSE,minLogLik=-Inf,bfullCI=FALSE){
       npar0<-.correct_npar0(npar0,EvolModel,designToEstim)
       tryCatch({
	modelParams$designToEstim<-designToEstim

	tree.height<-t
	if ((is.element("tree_height",names(phyltree)))&&(!is.na(phyltree$tree_height))){tree.height<-phyltree$tree_height}

	names(LogLik)<-c()
	lParamSummary=switch(EvolModel,
    		bm=.params.summary.bm(modelParams,mData,LogLik,RSS,n=phyltree$Ntips),
        	ouch=.params.summary.ouch(modelParams,mData,t,LogLik,phyltree$Ntips,npar0,RSS,tree.height,designToEstim),
#        	slouch=.params.summary.slouch(phyltree,modelParams,mData,LogLik,npar0,RSS,tree.height,designToEstim),
        	mvslouch=.params.summary.mvslouch(phyltree,modelParams,mData,t,LogLik,npar0,RSS,tree.height,designToEstim),
    		mvslouchtouchdetA0=.params.summary.mvslouch(phyltree,modelParams,mData,t,LogLik,npar0,RSS,tree.height,designToEstim)
    	    )

	if (is.element("regressCovar",names(modelParams))){
	## this will always be done, but the regression CIs are only calculated now
	    if (EvolModel=="mvslouchtoouchdetA0"){EvolModel<-"mvslouch"}
	    regressCovar<-modelParams$regressCovar
	    modelParams$paramPoint<-.cleanUpModelParams(modelParams) 
    	    tryCatch({
    		lParamSummary$confidence.interval<-.calcCI(modelParams,designToEstim,conf.level,regressCovar,bfullCI=FALSE)
	    },error=function(e){.my_message(paste("Error in confidence interval calculation",e),FALSE);.my_message("\n",FALSE)})
	}	
	lParamSummary
    },error=function(e){.my_message(paste("Error in parameter summary",e),FALSE);.my_message("\n",FALSE)})
}

.params.summary.bm<-function(modelParams,mData,LogLik,RSS,n){
    lParamsSummary<-vector("list",1)
    names(lParamsSummary)<-c("StS")
    lParamsSummary$StS<-modelParams$Sxx%*%t(modelParams$Sxx)
    numobs<- n*ncol(modelParams$Sxx)
    if (!is.null(mData)){numobs<- n*ncol(modelParams$Sxx)-length(which(is.na(c(mData))))}
    lParamsSummary$LogLik<-LogLik
    lParamsSummary$dof<- ncol(modelParams$Sxx)*(ncol(modelParams$Sxx)+1)/2+ncol(modelParams$Sxx)
    lParamsSummary$m2loglik<- -2*LogLik
    lParamsSummary$aic<- -2*LogLik+2*lParamsSummary$dof
    lParamsSummary$aic.c<- lParamsSummary$aic +2*lParamsSummary$dof*(lParamsSummary$dof+1)/(numobs-lParamsSummary$dof-1)
    lParamsSummary$sic<- lParamsSummary$m2loglik+log(numobs)*lParamsSummary$dof
    lParamsSummary$bic<-lParamsSummary$m2loglik+lParamsSummary$dof*log(numobs)
    lParamsSummary$RSS<-RSS
    lParamsSummary
}

.params.summary.ouch<-function(modelParams,mData=NULL,t=1,LogLik=-Inf,n=0,npar0=0,RSS=NA,tree.height=1,designToEstim=NULL){
    kY<-ncol(modelParams$A)
    lParamsSummary<-vector("list",18)   
    names(lParamsSummary)<-c("phyl.halflife","expmtA","mPsi.rotated","mPsi0.rotated","cov.matrix","corr.matrix","trait.regression","stationary.cov.matrix","stationary.corr.matrix","StS","LogLik","dof","m2loglik","aic","aic.c","sic","bic","RSS")
    lDecomps<-.decompEigenA.S(modelParams,NULL,NA,list(bCalcA=TRUE,bCovCalc=TRUE,dzetacalc=FALSE,lexptcalc=FALSE,kappacalc=FALSE,interceptcalc=FALSE),NULL)
    lParamsSummary$expmtA<-.calc.exptA(-t,lDecomps[[1]])
    lParamsSummary$mPsi.rotated<-apply(modelParams$mPsi,2,function(vPsi,expmtA){(diag(1,nrow(expmtA),ncol(expmtA))-expmtA)%*%vPsi},expmtA=lParamsSummary$expmtA)
    if (is.element("mPsi0",names(modelParams))&& !is.null(modelParams$mPsi0) && !is.na(modelParams$mPsi0[1])){
	lParamsSummary$mPsi0.rotated<-(diag(1,nrow(lParamsSummary$expmtA),ncol(lParamsSummary$expmtA))-lParamsSummary$expmtA)%*%modelParams$mPsi0
    }else{lParamsSummary$mPsi0.rotated<-NULL}
    lParamsSummary$cov.matrix<-.calc.cov.ouch.mv(t,lDecomps[[1]],lDecomps[[2]])
    lParamsSummary$corr.matrix<-.my_cov2cor(lParamsSummary$cov.matrix)
    if(kY>1){
	lParamsSummary$trait.regression<-NULL
	tryCatch({
	    lParamsSummary$trait.regression<-sapply(1:kY,function(i,mCov){mCov[i,-i,drop=FALSE]%*%solve(mCov[-i,-i,drop=FALSE])},mCov=lParamsSummary$cov.matrix,simplify=FALSE)
	},error=function(e){.my_message(paste("Error in trait regression calculation",e),FALSE);.my_message("\n",FALSE)})
    }else{lParamsSummary$trait.regression<-NULL}
    lParamsSummary$phyl.halflife<-.calc.phyl.halflife(modelParams$A,tree.height)
    k<-ncol(modelParams$A)
    if (length(which(Re(lDecomps[[1]]$eigA$values)<=0))==0){
	hadInvL1pL2<-apply(matrix(0:(k^2-1),k,k,byrow=TRUE),c(1,2),.CalcVlqStat,vlambda=lDecomps[[1]]$eigA$values,k=k)
	lParamsSummary$stationary.cov.matrix<-Re(lDecomps[[1]]$eigA$vectors%*%
					    (hadInvL1pL2*(
						lDecomps[[1]]$invP%*%(lDecomps[[2]]$S11)%*%t(lDecomps[[1]]$invP)))%*%
					    t(lDecomps[[1]]$eigA$vectors))
	lParamsSummary$stationary.corr.matrix<-.my_cov2cor(lParamsSummary$stationary.cov.matrix)
    }else{
        lParamsSummary$stationary.cov.matrix<-NULL
        lParamsSummary$stationary.cov.matrix.comment<-"A has negative eigenvalues, stationary covariance does not exist"
        lParamsSummary$stationary.corr.matrix<-NULL
        lParamsSummary$stationary.corr.matrix.comment<-"A has negative eigenvalues, stationary correlation does not exist"
    }
    lParamsSummary$StS<-lDecomps[[2]]$S11
    lParamsSummary$LogLik<-LogLik
    lParamsSummary$dof<-npar0
    if (modelParams$designToEstim$psi){
	lParamsSummary$dof<-lParamsSummary$dof+nrow(modelParams$A)*ncol(modelParams$mPsi)
	if (is.element("signsmPsi",names(designToEstim))){
	    vtoNA<-c(which(designToEstim$signsmPsi=="+"),which(designToEstim$signsmPsi=="-"))
	    if (length(vtoNA)){designToEstim$signsmPsi[vtoNA]<-NA}
	    numfixed<-length(which(!is.na(designToEstim$signsmPsi)))
	    lParamsSummary$dof<-lParamsSummary$dof-numfixed
	}
    }
    if (modelParams$designToEstim$psi0){
	lParamsSummary$dof<-lParamsSummary$dof+nrow(modelParams$A)
	if (is.element("signsmPsi0",names(designToEstim))){
	    vtoNA<-c(which(designToEstim$signsmPsi0=="+"),which(designToEstim$signsmPsi0=="-"))
	    if (length(vtoNA)){designToEstim$signsmPsi0[vtoNA]<-NA}
	    numfixed<-length(which(!is.na(designToEstim$signsmPsi0)))
	    lParamsSummary$dof<-lParamsSummary$dof-numfixed
	}
    }
    if (!modelParams$designToEstim$y0AncState && modelParams$designToEstim$y0){
	lParamsSummary$dof<-lParamsSummary$dof+nrow(modelParams$A)
	if (is.element("signsvY0",names(designToEstim))){
	    vtoNA<-c(which(designToEstim$signsvY0=="+"),which(designToEstim$signsvY0=="-"))
	    if (length(vtoNA)){designToEstim$signsvY0[vtoNA]<-NA}
	    numfixed<-length(which(!is.na(designToEstim$signsvY0)))
	    lParamsSummary$dof<-lParamsSummary$dof-numfixed
	}	
    }
    lParamsSummary$m2loglik<- -2*LogLik
    numNAdata<-0
    if (!is.null(mData)){numNAdata<-length(which(is.na(mData)))}
    lParamsSummary$aic<-lParamsSummary$m2loglik+2*lParamsSummary$dof
    lParamsSummary$aic.c<-lParamsSummary$aic+2*lParamsSummary$dof*(lParamsSummary$dof+1)/(nrow(modelParams$A)*n-numNAdata-lParamsSummary$dof-1)
    lParamsSummary$sic<-lParamsSummary$m2loglik+log(nrow(modelParams$A)*n-numNAdata)*lParamsSummary$dof
    lParamsSummary$bic<-lParamsSummary$m2loglik+lParamsSummary$dof*log(nrow(modelParams$A)*n-numNAdata)
    lParamsSummary$RSS<-RSS
    lParamsSummary
}

.params.summary.mvslouch<-function(phyltree,modelParams,mData=NULL,t=1,LogLik=-Inf,npar0=0,RSS=NA,tree.height=1,designToEstim=NULL){
    n<-phyltree$Ntips

    lParamsSummary<-vector("list",26)   
    names(lParamsSummary)<-c("phyl.halflife","expmtA","optimal.regression","mPsi.rotated","mPsi0.rotated","cov.matrix","corr.matrix","conditional.cov.matrix","conditional.corr.matrix","stationary.cov.matrix","stationary.corr.matrix","optima.cov.matrix","optima.corr.matrix","cov.with.optima","corr.with.optima","evolutionary.regression","trait.regression","StS","LogLik","dof","m2loglik","aic","aic.c","sic","bic","RSS")

    lDecomps<-.decompEigenA.S(modelParams,NULL,NA,list(bCalcA=TRUE,bCovCalc=TRUE,dzetacalc=FALSE,lexptcalc=FALSE,kappacalc=FALSE,interceptcalc=FALSE),NULL)

    lParamsSummary$expmtA<-.calc.exptA(-t,lDecomps[[1]])
    lParamsSummary$mPsi.rotated<-apply(modelParams$mPsi,2,function(vPsi,expmtA){(diag(1,nrow(expmtA),ncol(expmtA))-expmtA)%*%vPsi},expmtA=lParamsSummary$expmtA)
    if (is.element("mPsi0",names(modelParams))&& !is.null(modelParams$mPsi0) && !is.na(modelParams$mPsi0[1])){
	lParamsSummary$mPsi0.rotated<-(diag(1,nrow(lParamsSummary$expmtA),ncol(lParamsSummary$expmtA))-lParamsSummary$expmtA)%*%modelParams$mPsi0
    }else{lParamsSummary$mPsi0.rotated<-NULL}
    lParamsSummary$cov.matrix<-.calc.cov.slouch.mv(t,lDecomps[[1]],lDecomps[[2]])
    lParamsSummary$corr.matrix<-.my_cov2cor(lParamsSummary$cov.matrix)
    lParamsSummary$phyl.halflife<-.calc.phyl.halflife(modelParams$A,tree.height)
    kY<-ncol(modelParams$A)
    kX<-ncol(modelParams$B)
    lParamsSummary$evolutionary.regression<-lParamsSummary$cov.matrix[1:kY,(kY+1):(kY+kX)]%*%solve(lParamsSummary$cov.matrix[(kY+1):(kY+kX),(kY+1):(kY+kX)])
    lParamsSummary$trait.regression<-NULL
    tryCatch({
	lParamsSummary$trait.regression<-sapply(1:kY,function(i,mCov){mCov[i,-i,drop=FALSE]%*%solve(mCov[-i,-i,drop=FALSE])},mCov=lParamsSummary$cov.matrix,simplify=FALSE)
    },error=function(e){.my_message(paste("Error in trait regression calculation",e),FALSE);.my_message("\n",FALSE)})
    lParamsSummary$conditional.cov.matrix<-lParamsSummary$cov.matrix[1:kY,1:kY]-lParamsSummary$cov.matrix[1:kY,(kY+1):(kY+kX)]%*%solve(lParamsSummary$cov.matrix[(kY+1):(kY+kX),(kY+1):(kY+kX)])%*%lParamsSummary$cov.matrix[(kY+1):(kY+kX),1:kY]
    lParamsSummary$conditional.corr.matrix<-.my_cov2cor(lParamsSummary$conditional.cov.matrix)
    if (!is.na(lDecomps[[1]]$A1B[1])&&(length(which(Re(lDecomps[[1]]$eigA$values)<=0))==0)){
	hadInvL1pL2<-apply(matrix(0:(kY^2-1),kY,kY,byrow=TRUE),c(1,2),.CalcVlqStat,vlambda=lDecomps[[1]]$eigA$values,k=kY)
	lParamsSummary$stationary.cov.matrix<-Re(lDecomps[[1]]$eigA$vectors%*%
					    (hadInvL1pL2*(
						lDecomps[[1]]$invP%*%(
						    lDecomps[[2]]$S11+
						    lDecomps[[2]]$S12%*%t(lDecomps[[1]]$A1B)+
						    lDecomps[[1]]$A1B%*%lDecomps[[2]]$S21+
				    		    lDecomps[[1]]$A1B%*%lDecomps[[2]]$S22%*%t(lDecomps[[1]]$A1B)
						)%*%t(lDecomps[[1]]$invP)))%*%t(lDecomps[[1]]$eigA$vectors))
        lParamsSummary$stationary.corr.matrix<-.my_cov2cor(lParamsSummary$stationary.cov.matrix)
        lParamsSummary$optimal.regression<- (-1)*lDecomps[[1]]$A1B
    }else{
        lParamsSummary$stationary.cov.matrix<-NULL
        lParamsSummary$stationary.cov.matrix.comment<-"A has negative or zero eigenvalues, stationary covariance does not exist."
        lParamsSummary$stationary.corr.matrix<-NULL
        lParamsSummary$stationary.corr.matrix.comment<-"A has negative or zero eigenvalues, stationary correlation does not exist."
        lParamsSummary$optimal.regression<- NULL
        lParamsSummary$optimal.regression.comment<- "A has negative eigenvalues, optimal regression does not exist."
        lParamsSummary$A1B<-lDecomps[[1]]$A1B
    }
    if(!is.na(lDecomps[[1]]$A1B[1])){
	lParamsSummary$optima.cov.matrix<-t*lDecomps[[1]]$A1B%*%lDecomps[[2]]$S22%*%t(lDecomps[[1]]$A1B)
	lParamsSummary$optima.corr.matrix<-.my_cov2cor(lParamsSummary$optima.cov.matrix)
	lParamsSummary$cov.with.optima<-lParamsSummary$cov.matrix[1:kY,(kY+1):(kY+kX)]%*%t(lDecomps[[1]]$A1B)*(-1)
        lParamsSummary$corr.with.optima<-apply(matrix(0:((kY)^2-1),kY,kY,byrow=TRUE),c(1,2),function(ij,kY,mcov.traits,mcov.optima,mcov.with){i<-ij%/%kY+1;j<-ij%%kY+1;mcov.with[i,j]/(sqrt(mcov.traits[i,i]*mcov.optima[j,j]))},kY=kY,mcov.traits=lParamsSummary$cov.matrix,mcov.optima=lParamsSummary$optima.cov.matrix,mcov.with=lParamsSummary$cov.with.optima)
    }
    else{
	lParamsSummary$optima.cov.matrix<-"A has 0 eigenvalues, optima covariance matrix cannot be calculated in the current implementation."
	lParamsSummary$optima.corr.matrix<-"A has 0 eigenvalues, optima correlation matrix cannot be calculated in the current implementation."
        lParamsSummary$cov.with.optima<-"A has 0 eigenvalues, covariance with optima matrix cannot be calculated in the current implementation."
        lParamsSummary$corr.with.optima<-"A has 0 eigenvalues, correlation with optima matrix cannot be calculated in the current implementation."
    }

    lParamsSummary$StS<-rbind(cbind(lDecomps[[2]]$S11,lDecomps[[2]]$S12),cbind(lDecomps[[2]]$S21,lDecomps[[2]]$S22))
    lParamsSummary$LogLik<-LogLik
    ## no correction for signs Sxx
    lParamsSummary$dof<-npar0 + nrow(modelParams$Sxx)+nrow(modelParams$Sxx)*(1+nrow(modelParams$Sxx))/2 ## Sxx is for now estimated directly
    if (modelParams$designToEstim$B){
	lParamsSummary$dof<-lParamsSummary$dof+nrow(modelParams$B)*ncol(modelParams$B)
	if (is.element("signsB",names(designToEstim))){
	    vtoNA<-c(which(designToEstim$signsB=="+"),which(designToEstim$signsB=="-"))
	    if (length(vtoNA)){designToEstim$signsB[vtoNA]<-NA}
	    numfixed<-length(which(!is.na(designToEstim$signsB)))
	    lParamsSummary$dof<-lParamsSummary$dof-numfixed
	}	
    }    
    if (modelParams$designToEstim$psi){
	lParamsSummary$dof<-lParamsSummary$dof+nrow(modelParams$A)*ncol(modelParams$mPsi)
    	if (is.element("signsmPsi",names(designToEstim))){
	    vtoNA<-c(which(designToEstim$signsmPsi=="+"),which(designToEstim$signsmPsi=="-"))
	    if (length(vtoNA)){designToEstim$signsmPsi[vtoNA]<-NA}
	    numfixed<-length(which(!is.na(designToEstim$signsmPsi)))
	    lParamsSummary$dof<-lParamsSummary$dof-numfixed
	}
    }    
    if (modelParams$designToEstim$psi0){
	lParamsSummary$dof<-lParamsSummary$dof+nrow(modelParams$A)
    	if (is.element("signsmPsi0",names(designToEstim))){
	    vtoNA<-c(which(designToEstim$signsmPsi0=="+"),which(designToEstim$signsmPsi0=="-"))
	    if (length(vtoNA)){designToEstim$signsmPsi0[vtoNA]<-NA}
	    numfixed<-length(which(!is.na(designToEstim$signsmPsi0)))
	    lParamsSummary$dof<-lParamsSummary$dof-numfixed
	}
    }    
    if (!modelParams$designToEstim$y0AncState && modelParams$designToEstim$y0){
	lParamsSummary$dof<- lParamsSummary$dof+nrow(modelParams$A)
    	if (is.element("signsvY0",names(designToEstim))){
	    vtoNA<-c(which(designToEstim$signsvY0=="+"),which(designToEstim$signsvY0=="-"))
	    if (length(vtoNA)){designToEstim$signsvY0[vtoNA]<-NA}
	    numfixed<-length(which(!is.na(designToEstim$signsvY0)))
	    lParamsSummary$dof<-lParamsSummary$dof-numfixed
	}	
    }
    lParamsSummary$m2loglik<- -2*LogLik
    lParamsSummary$aic<-lParamsSummary$m2loglik+2*lParamsSummary$dof
    numNAdata<-0
    if (!is.null(mData)){numNAdata<-length(which(is.na(mData)))}
    lParamsSummary$aic.c<-lParamsSummary$aic+2*lParamsSummary$dof*(lParamsSummary$dof+1)/((nrow(modelParams$A)+nrow(modelParams$Sxx))*n-numNAdata-lParamsSummary$dof-1)
    lParamsSummary$sic<-lParamsSummary$m2loglik+log((nrow(modelParams$A)+nrow(modelParams$Sxx))*n-numNAdata)*lParamsSummary$dof
    lParamsSummary$bic<-lParamsSummary$m2loglik+lParamsSummary$dof*log((nrow(modelParams$A)+nrow(modelParams$Sxx))*n-numNAdata)
    lParamsSummary$RSS<-RSS
    lParamsSummary
}

.norm.max<-function(M){max(abs(M))}

.calc.phyl.halflife<-function(A,tree.height){
    mPhylHalfLife<-matrix(NA,nrow=3,ncol=nrow(A))
    rownames(mPhylHalfLife)<-c("eigenvalues","halflife","%treeheight")
    eigA<-eigen(A)
    mPhylHalfLife[1,]<-eigA$values
    mPhylHalfLife[2,]<- log(2)/(Re(eigA$values)) 
    mPhylHalfLife[3,]<- 100*(mPhylHalfLife[2,]/tree.height)
    ## Idea for bound taken from von Lohan, Sensitivity of the matrix exponential remove det part
    list("directions"=eigA$vectors,"halflives"=mPhylHalfLife,"halflifeLowerbounds"=c(Re(log(2)/sum(eigA$values))))
}

.SummarizeFullPoint<-function(vPoint,mData,PhylTree,EvolModel,regimes.times,regimes,EstimationParams=NULL,modelParams=NULL,t=c(1),dof=NULL,calcCI=NULL,sigmaRule=NULL,tol=c(0.01,0.01),maxIter=c(10,50,100),bShouldPrint=FALSE,Merror=NULL,predictors=NULL,minLogLik=-Inf,LogLik=NULL){
## called in wrappers.R
## there is never a call with mData==NULL
## vPoint and modelParams CANNOT be NULL at the same time
    lEvalPoint<-NA

    

## =========================================================================================================================
    ## checking if mData is a matrix, has to be for PCMBase
    if (!inherits(PhylTree,"phylo")){.my_stop("The phylogenetic tree has to be of the phylo format (i.e. ape).",TRUE)}
    mData<-.check_input_trait_data(mData,NULL,PhylTree$tip.label)
    n<-nrow(mData)
    kYX<-ncol(mData)
## =========================================================================================================================
    
    if (is.null(vPoint) && is.null(modelParams)){.my_stop("Error in .SummarizeFullPoint call: both vPoint and modelParams are NULL. Cannot summarize the observed point. Provide either the model parameters or their transformation for the optimizer.",TRUE)}
    
    if (EvolModel=="slouch"){EvolModel<-"mvslouch";.my_warning("WARNING: Call to slouch not yet implemented, using mvslouch model. You might need to run again with correct input structure. \n",TRUE,FALSE)}  

    if (is.null(EstimationParams)){EstimationParams<-list()}

    if (!is.null(calcCI)){
	EstimationParams$calcCI<-calcCI
	if (!is.null(sigmaRule)){
	    if (!is.element("designToEstim",names(EstimationParams))){EstimationParams$desginToEstim<-list()}
    	    EstimationParams$desginToEstim$sigmaRule<-sigmaRule
        }
    }
    
## =========================================================================================================================
    ## checking if M.error is of the correct type and making it appropriate for further use
    Merror<-.createMeasurementError(Merror,n,kYX)
## =========================================================================================================================

    if (!is.element("kY",names(EstimationParams))){
        if ((EvolModel=="mvslouch") && (is.null(modelParams))){
    	    if (is.element("A",names(vPoint))){EstimationParams$kY<-1}
	    else{EstimationParams$kY<-which(names(vPoint)=="Astart")-which(names(vPoint)=="Aend")+1}
	}
	else{
	    EstimationParams$kY=switch(EvolModel,
    		bm=0,
        	ouch=kYX,
        	slouch=1,
        	mvslouch=nrow(modelParams$A)
    	    )
	}
    }

    if (!is.element("kX",names(EstimationParams))){
	    EstimationParams$kX=switch(EvolModel,
    	    	    bm=kYX,
    		    ouch=0,
        	    slouch=kYX-1,
        	    mvslouch=kYX-EstimationParams$kY
		)    
    }

    root.regime<-NULL ## in summary there is no point in having a root regime as either it sits in the initial state or will have no effect
    tmpkX<-EstimationParams$kX; if (tmpkX==0){tmpkX<-ncol(mData)}
    regimesList<-.InitialRegimeSetup(phyltree=PhylTree,regimes=regimes,regimes.times=regimes.times,mData=mData,kX=tmpkX,kYX=ncol(mData),root.regime=root.regime,M.error=Merror)
    bOKregimes<-regimesList$bOKregimes
    if (!bOKregimes){.my_stop("The regimes, regimes.times and phyltree variables are not consistant. Cannot perform estimation. Please see manual on their format.",TRUE)}
    regimes<-regimesList$regimes
    regimes.times<-regimesList$regimes.times
    regimes.types<-regimesList$regimes.types
    root.regime<-regimesList$root.regime
    regimes.types.orig<-regimesList$regimes.types.orig
    PhylTree<-regimesList$phyltree
    pcmbase_model_box<-regimesList$pcmbase_model_box ## this is the model object in which parameters for all regimes sit in for PCMBase
    mData<-.check_input_trait_data(mData,n=PhylTree$Ntips,vSpeciesLabels=PhylTree$tip.label)
    ## mData has to be a matrix of PCMBase
    ## order of rows in mData has to be the same as the order of tips in the phylogenetic tree

    ## This is essentially a dummy object that might become useful in the future when B for the mvSLOUCH model is estimated by GLS
    ## now it is used to pass the times of species
    lPrecalculates<-list(vSpecies_times=NA)
    lPrecalculates$vSpecies_times<-PhylTree$time.of.nodes[PhylTree$tip_species_index] 

    tree_height<-NA
    if ((is.element("tree_height",names(PhylTree)))&&(!is.na(PhylTree$tree_height))){tree_height<-PhylTree$tree_height}
    if (!is.na(tree_height)){
	if(!any(sapply(t,function(t1,t2){isTRUE(all.equal(t1,t2))},t2=tree_height,simplify=TRUE))){
	    t<-c(t,tree_height)
	}
    }
    
    ## we do not need estimate.root.state here 
    ## the dof is to be provided by the user in the wrapper function that calls this function	
    ## in the line below Merror=Merror was removed as this now sits in the pcmbase_model_box object
    EstimationParams<-.set.estimparams(params=list(pcmbase_model_box=pcmbase_model_box,method="glsgc",EvolModel=EvolModel,EstimParams=EstimationParams,tol=tol,maxIter=maxIter,bShouldPrint=bShouldPrint),kY=EstimationParams$kY,kX=EstimationParams$kX,numregs=length(regimes.types))

    if (!is.null(predictors)){
        EstimationParams$predictors<-predictors
        predictors<-colnames(mData)[predictors]
    }

    if (EstimationParams$designToEstim$y0AncState){
        if (!is.null(root.regime)){EstimationParams$designToEstim$y0Regime<-which(regimes.types.orig==root.regime)}
        else{
	    if (is.element("y0Regime",names(EstimationParams$designToEstim))){EstimationParams$designToEstim$y0Regime<-which(regimes.types.orig==EstimationParams$designToEstim$y0Regime)}       
    	    else{EstimationParams$designToEstim$y0Regime<-regimes[[1]][1]} ## need to choose something anyway ...    
    	}
    }
    
    if (is.null(vPoint)){
        if (EvolModel=="mvslouch"){
            if(is.element("Fixed",names(EstimationParams))&& is.element("Sxx",names(EstimationParams$Fixed))){Sxx<-EstimationParams$Fixed$Sxx}else{if(is.element("Sxx",names(modelParams))){Sxx<-modelParams$Sxx}}
            if(is.element("Fixed",names(EstimationParams))&& is.element("vX0",names(EstimationParams$Fixed))){vX0<-EstimationParams$Fixed$vX0}else{if(is.element("vX0",names(modelParams))){vX0<-modelParams$vX0}}
            if(is.element("Fixed",names(EstimationParams))&& is.element("B",names(EstimationParams$Fixed))){B<-EstimationParams$Fixed$B}else{if(is.element("B",names(modelParams))){B<-modelParams$B}}
    	}
    	if ((EvolModel=="ouch")||(EvolModel=="mvslouch")){
    	    if(is.element("Fixed",names(EstimationParams))&& is.element("vY0",names(EstimationParams$Fixed))){vY0<-EstimationParams$Fixed$vY0}else{if(is.element("vY0",names(modelParams))){vY0<-modelParams$vY0}}
    	    if(is.element("Fixed",names(EstimationParams))&& is.element("mPsi",names(EstimationParams$Fixed))){mPsi<-EstimationParams$Fixed$mPsi}else{if(is.element("mPsi",names(modelParams))){mPsi<-modelParams$mPsi}}
    	    mPsi0<-matrix(0,ncol=1,nrow=EstimationParams$kY)
    	    if(is.element("Fixed",names(EstimationParams))&& is.element("mPsi0",names(EstimationParams$Fixed))){mPsi0<-EstimationParams$Fixed$mPsi0}else{if(is.element("mPsi0",names(modelParams))&&!is.null(modelParams$mPsi0) &&!is.na(modelParams$mPsi0[1])){mPsi0<-modelParams$mPsi0}}
	}
    }
    EstimationParams<-.beginEstimationParams(EvolModel,EstimationParams,mData,PhylTree)
    if (is.null(vPoint)){
        if (!is.element("Fixed",names(EstimationParams))){EstimationParams$Fixed<-list()}
    	if (EvolModel=="mvslouch"){
    	    if (!is.element("Sxx",names(EstimationParams$Fixed))){EstimationParams$Fixed$Sxx<-Sxx}else{EstimationParams$Fixed$Sxx<-modelParams$Sxx}
	    if (!is.element("vX0",names(EstimationParams$Fixed))){EstimationParams$Fixed$vX0<-vX0}else{EstimationParams$Fixed$vX0<-modelParams$vX0}
    	    if (!is.element("B",names(EstimationParams$Fixed))){EstimationParams$Fixed$B<-B}else{EstimationParams$Fixed$B<-modelParams$B}
    	}
    	if ((EvolModel=="ouch")||(EvolModel=="mvslouch")){
    	    if (!is.element("vY0",names(EstimationParams$Fixed))){EstimationParams$Fixed$vY0<-vY0}else{EstimationParams$Fixed$vY0<-modelParams$vY0}
    	    if (!is.element("mPsi",names(EstimationParams$Fixed))){EstimationParams$Fixed$mPsi<-mPsi}else{EstimationParams$Fixed$mPsi<-modelParams$mPsi}
	    if (!is.element("mPsi0",names(EstimationParams$Fixed))){EstimationParams$Fixed$mPsi0<-mPsi0}else{if(is.element("mPsi0",names(modelParams))&&!is.null(modelParams$mPsi0) &&!is.na(modelParams$mPsi0[1])){EstimationParams$Fixed$mPsi0<-modelParams$mPsi0}else{EstimationParams$Fixed$mPsi0<-matrix(0,ncol=1,nrow=EstimationParams$kY)}}
        }
    }

    if (is.null(modelParams)){
        modelParams<-.par_transform_withPCMBase(vPoint,EstimationParams,EvolModel)
        modelParams$vPoint<-vPoint
	bFull<-TRUE
    }else{
	modelParams$pcmbase_model_box<-EstimationParams$pcmbase_model_box
	bFull<-FALSE
    }
        
    if (!is.element("regimeTimes",names(modelParams))){modelParams$regimeTimes<-regimes.times}
    if (!is.element("regimes",names(modelParams))){modelParams$regimes<-regimes}
    if (!is.element("regimeTypes",names(modelParams))){modelParams$regimeTypes<-regimes.types}
    if (!is.element("regimes.types.orig",names(modelParams))){modelParams$regimes.types.orig<-regimes.types.orig} ## this has to be passed as PCMBase will use original regime names while in GLS we just use the regime indices
    
    modelParams$pcmbase_model_box<-.update_pcmbase_box_params(modelParams,EvolModel)
    
    if (!is.element("M_error",names(modelParams))){modelParams$M_error<-Merror}
    lEvalPoint<-vector("list",length(t))

    ## mData is assumed to be in the following format
    ## number of rows is the number of species, first columns 1:kY response, columns (kY+1):(kY+kX) predictor (for mvslouch)
    if (!is.null(mData)){
        mY<-(mData[,1:EstimationParams$kY,drop=FALSE])            
    }else{mY<-NA}

    if (is.null(vPoint)){
        if (is.element("A",names(modelParams))){EstimationParams$Fixed$A<-modelParams$A}
        if (is.element("Syy",names(modelParams))){EstimationParams$Fixed$Syy<-modelParams$Syy}
        
        if (is.null(LogLik)){lEvalPoint<-sapply(t,function(tcalc,EvolModel,PhylTree,mData,mY,modelParams,lPrecalculates,EstimationParams,tol,maxIter,bShouldPrint,minLogLik){
    		lres<-list(t=NA)
    		tryCatch({lres<-.EvaluatePoint(EvolModel,PhylTree,mData,mY,modelParams,lPrecalculates,EstimationParams,tol[2],maxIter[2],bShouldPrint,FALSE,list(A=EstimationParams$Fixed$A,Syy=EstimationParams$Fixed$Syy),TRUE,TRUE,minLogLik=minLogLik,t=tcalc)},error=function(e){.my_message("Error in evaluating point to summarize: \n",FALSE);.my_message(e,FALSE);.my_message("\n",FALSE)})
    		lres
    	    },EvolModel=EvolModel,PhylTree=PhylTree,mData=mData,mY=mY,modelParams=modelParams,lPrecalculates=lPrecalculates,EstimationParams=EstimationParams,tol=tol,maxIter=maxIter,bShouldPrint=bShouldPrint,minLogLik=minLogLik,simplify=FALSE)}
        else{lEvalPoint<-sapply(t,function(tcalc,EvolModel,PhylTree,mData,mY,modelParams,lPrecalculates,EstimationParams,tol,maxIter,bShouldPrint,minLogLik){
    		lres<-list(t=NA)
    		tryCatch({lres<-.EvaluatePoint(EvolModel,PhylTree,mData,mY,modelParams,lPrecalculates,EstimationParams,tol[2],maxIter[2],bShouldPrint,FALSE,list(A=EstimationParams$Fixed$A,Syy=EstimationParams$Fixed$Syy),FALSE,TRUE,minLogLik=LogLik,t=tcalc)},error=function(e){.my_message("Error in evaluating point to summarize: \n",FALSE);.my_message(e,FALSE);.my_message("\n",FALSE)})
		lres    	    
    	    },EvolModel=EvolModel,PhylTree=PhylTree,mData=mData,mY=mY,modelParams=modelParams,lPrecalculates=lPrecalculates,EstimationParams=EstimationParams,tol=tol,maxIter=maxIter,bShouldPrint=bShouldPrint,minLogLik=minLogLik,simplify=FALSE)}
    }
    else{
        if (is.null(LogLik)){lEvalPoint<-sapply(t,function(tcalc){
    		lres<-list(t=NA)
    		tryCatch({lres<-.EvaluatePoint(EvolModel,PhylTree,mData,mY,modelParams,lPrecalculates,EstimationParams,tol[2],maxIter[2],bShouldPrint,bFull,NULL,TRUE,TRUE,minLogLik=minLogLik,t=tcalc)},error=function(e){.my_message("Error in evaluating point to summarize: \n",FALSE);.my_message(e,FALSE);.my_message("\n",FALSE)})
    		lres
    	    },simplify=FALSE)}
        else{lEvalPoint<-sapply(t,function(tcalc){
    		lres<-list(t=NA)
    		tryCatch({lres<-.EvaluatePoint(EvolModel,PhylTree,mData,mY,modelParams,lPrecalculates,EstimationParams,tol[2],maxIter[2],bShouldPrint,bFull,NULL,FALSE,TRUE,minLogLik=LogLik,t=tcalc)},error=function(e){.my_message("Error in evaluating point to summarize: \n",FALSE);.my_message(e,FALSE);.my_message("\n",FALSE)})
    		lres
    	    },simplify=FALSE)}
    }
    
    ## here we have names of list entries as numbers casted to a string
##    names(lEvalPoint)<-as.character(t)
    names(lEvalPoint)<-as.character(paste0("t_",t))
    for (i in 1:length(t)){
        lEvalPoint[[i]]$modelParams<-.cleanUpModelParams(lEvalPoint[[i]]$modelParams)
        lEvalPoint[[i]]<-.correct.names(lEvalPoint[[i]],regimes.types.orig,if(EvolModel!="bm"){colnames(mData)[1:EstimationParams$kY]}else{NULL},if ((EvolModel=="mvslouch")||(EvolModel=="bm")||(EvolModel=="slouch")){colnames(mData)[(EstimationParams$kY+1):(EstimationParams$kY+EstimationParams$kX)]}else{NULL},predictors,EvolModel)
    }
    if (!is.null(dof)){
        numNAdata<-length(which(is.na(mData)))
	for (i in 1:length(t)){
	    lEvalPoint[[i]]$PointSummary$dof<-dof
	    lEvalPoint[[i]]$PointSummary$aic<-lEvalPoint[[i]]$PointSummary$m2loglik+2*dof
	    lEvalPoint[[i]]$PointSummary$aic.c<-lEvalPoint[[i]]$PointSummary$aic+2*dof*(dof+1)/(kYX*n-numNAdata-dof-1)
	    lEvalPoint[[i]]$PointSummary$sic<-lEvalPoint[[i]]$PointSummary$m2loglik+log(kYX*n-numNAdata)*dof
	    lEvalPoint[[i]]$PointSummary$bic<-lEvalPoint[[i]]$PointSummary$sic
	}
    }
    lEvalPoint
}

.cleanUpModelParams<-function(modelParams){
    if (is.element("eigenSignsA",names(modelParams))){modelParams$eigenSignsA<-NULL}
    if (is.element("GivensQCsignsA",names(modelParams))){modelParams$GivensQCsignsA<-NULL}
    if (is.element("regimes",names(modelParams))){modelParams$regimes<-NULL}
    if (is.element("regimeTypes",names(modelParams))){modelParams$regimeTypes<-NULL}
    if (is.element("regimeTimes",names(modelParams))){modelParams$regimeTimes<-NULL}
    if (is.element("regimes.types",names(modelParams))){modelParams$regimes.types<-NULL}
    if (is.element("regimes.times",names(modelParams))){modelParams$regimes.times<-NULL}
    if (is.element("mCovPhyl",names(modelParams))){modelParams$mCovPhyl<-NULL}
    if (is.element("invSXX",names(modelParams))){modelParams$invSXX<-NULL}
    if (is.element("intercept",names(modelParams))){modelParams$intercept<-NULL}
    if (is.element("precalcMatrices",names(modelParams))){modelParams$precalcMatrices<-NULL}
    if (is.element("paramPoint",names(modelParams))){modelParams$paramPoint<-NULL}    
    if (is.element("regressCovar",names(modelParams))){modelParams$regressCovar<-NULL}    
    if (is.element("EstimParams",names(modelParams))){modelParams$EstimParams<-NULL}    
    if (is.element("kY",names(modelParams))){modelParams$kY<-NULL}    
    if (is.element("kX",names(modelParams))){modelParams$kX<-NULL}    
    if (is.element("method",names(modelParams))){modelParams$method<-NULL}        
    if (is.element("tol",names(modelParams))){modelParams$tol<-NULL}    
    if (is.element("maxIter",names(modelParams))){modelParams$maxIter<-NULL}    
    if (is.element("bShouldPrint",names(modelParams))){modelParams$bShouldPrint<-NULL}    
    if (is.element("EvolModel",names(modelParams))){modelParams$EvolModel<-NULL}    
    if (is.element("maxTries",names(modelParams))){modelParams$maxTries<-NULL}    
    if (is.element("process",names(modelParams))){modelParams$process<-NULL}    
    if (is.element("procparams",names(modelParams))){modelParams$procparams<-NULL}    
    if (is.element("minLogLik",names(modelParams))){modelParams$minLogLik<-NULL}    
    if (is.element("designToEstim",names(modelParams))){modelParams$designToEstim<-NULL}    
    if (is.element("Merror",names(modelParams))){modelParams$Merror<-NULL}    
    
    if (is.element("regimes.types.orig",names(modelParams))){modelParams$regimes.types.orig<-NULL}
    if (is.element("pcmbase_model_regimes",names(modelParams))){modelParams$pcmbase_model_regimes<-NULL}
        
    modelParams
}

.correct_npar0<-function(npar0,evolmodel,designToEstim=NULL){
    if (evolmodel=="mvslouchtoouchdetA0"){evolmodel<-"mvslouch"}
    res<-npar0
    if ((!is.null(designToEstim))&&((evolmodel=="ouch")||(evolmodel=="mvslouch"))){
    ## we do correction here for non-GLS estimated parameters
	if ((is.element("signsAnonGLS",names(designToEstim)))&&(designToEstim$Atype!="Any")){
	    vtoNA<-c(which(designToEstim$signAnonGLS=="+"),which(designToEstim$signsAnonGLS=="-"))
	    if (length(vtoNA)){designToEstim$signsAnonGLS[vtoNA]<-NA}
	    numfixed<-length(which(!is.na(designToEstim$signsAnonGLS)))	    
	    ##npar0<-npar0-numfixed
	    res<-npar0-numfixed
	}
	if ((is.element("signsBnonGLS",names(designToEstim)))&&(!designToEstim$B)&&(designToEstim$Btype!="Any")){
	    vtoNA<-c(which(designToEstim$signsBnonGLS=="+"),which(designToEstim$signsBnonGLS=="-"))
	    if (length(vtoNA)){designToEstim$signsBnonGLS[vtoNA]<-NA}
	    numfixed<-length(which(!is.na(designToEstim$signsBnonGLS)))	    
	    ##npar0<-npar0-numfixed
	    res<-npar0-numfixed
	}
	if ((is.element("signsSyynonGLS",names(designToEstim)))&&(designToEstim$Syytype!="Any")){
	    vtoNA<-c(which(designToEstim$signSyynonGLS=="+"),which(designToEstim$signsSyynonGLS=="-"))
	    if (length(vtoNA)){designToEstim$signsSyynonGLS[vtoNA]<-NA}
	    numfixed<-length(which(!is.na(designToEstim$signsSyynonGLS)))	    
	    ##npar0<-npar0-numfixed
	    res<-npar0-numfixed
	}	
    }    
    res
}
