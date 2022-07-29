## This file is part of mvSLOUCH

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this package or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .

.calc_species_mean<-function(evolmodel,expmtA,exptjA,regimes,vX0=NULL,vY0=NULL,mPsi=NULL,A1B=NULL,mPsi0=NULL){
## called in loglik.R in order to get RSS
    vMean=switch(evolmodel,
	bm=vX0, ## BM with drift here also
	ouch=.calc.mean.ouch.mv(expmtA,vY0,mPsi,mPsi0,exptjA,regimes),
#	slouch=(),
	mvslouch=.calc.mean.slouch.mv(expmtA,A1B,vY0,vX0,mPsi,mPsi0,exptjA,regimes)
    )
    vMean
}

.calc_species_covariance<-function(evolmodel,species_time,lAcalcs,lScalcs,method="minus.v",tol=1e-10,stationary=FALSE,b_correct_forPD=TRUE){
    mCov=switch(evolmodel, 
	bm=species_time*lScalcs$S22, ## BM with drift here also
	ouch=.calc.cov.ouch.mv(species_time,lAcalcs,lScalcs,stationary=stationary,method=method,b_correct_forPD=b_correct_forPD),
#	slouch=(),
	mvslouch=.calc.cov.slouch.mv(species_time,lAcalcs,lScalcs,tol=tol,method=method,b_correct_forPD=b_correct_forPD)
    )
    if (b_correct_forPD&&((!.matrixcalc_is.symmetric.matrix(mCov)) || (!.matrixcalc_is.positive.definite(mCov)))){
	tryCatch({mCov<-as.matrix(Matrix::nearPD(mCov)$mat)},error=function(e){.my_message(e,FALSE)})    
    }
    mCov
}

.calc.mean.bm.mv<-function(vX0){vX0}
.calc.cov.bm.mv<-function(species_time,StS,b_correct_forPD=TRUE){
    mCov<-species_time*StS
    if (b_correct_forPD && ((!.matrixcalc_is.symmetric.matrix(mCov)) || (!.matrixcalc_is.positive.definite(mCov)))){
        tryCatch({mCov<-as.matrix(Matrix::nearPD(mCov)$mat)},error=function(e){.my_message(e,FALSE)})
    }
    mCov
}
    

.calc.mean.ouch.mv<-function(expmtA,vY0,mPsi,mPsi0=NULL,exptjA=NULL,regimes=NULL){
## called in simulVasicekproc.R, simulVasicekprocphyl.R
## t is not needed to be passed, it is in the expmtA not explicitely
## vPsia has to be given as a vector
    kY<-length(vY0)
    vMean<-expmtA%*%vY0
    if (!(is.null(exptjA))){
	## we have for length - 1 as we only take differences
	vPsiSum<-apply(matrix(sapply(1:(length(exptjA)-1),function(x){(exptjA[[x+1]]-exptjA[[x]])%*%mPsi[,regimes[x]]},simplify=TRUE),nrow=length(vY0)),1,sum) 
	vMean<-vMean+vPsiSum ## a trick here to reuse the calculations from vMean
    }else{vMean<-vMean+(diag(1,kY,kY)-expmtA)%*%mPsi[,1]} ## there is only one optimum ie one Psi
    if (!(is.null(mPsi0))){vMean<-vMean+(diag(1,kY,kY)-expmtA)%*%mPsi0}
    vMean[which(abs(vMean)<1e-15)]<-0
    names(vMean)<-NULL
    vMean
}

.calc.cov.ouch.mv<-function(t,lAcalcs,lScalcs,stationary=FALSE,method="minus.v",b_correct_forPD=TRUE){
## called in simulVasicekproc.R, simulVasicekprocphyl.R
    invA<-lAcalcs$invA
    kY<-nrow(lAcalcs$A)
    mCov<-matrix(NA,kY,kY)
    if (method=="minus.v"){
        Int11<-.calc.integral.evAStevA(-t,lScalcs$S11,lAcalcs)
	mCov<-Int11
    }
    if (method=="plus.v"){
        expmtA<-.calc.exptA(-t,lAcalcs)
## ----------------------------------------------------------------------------------
	Int11<-.calc.integral.evAStevA(t,lScalcs$S11,lAcalcs)
## ----------------------------------------------------------------------------------    
	if (!stationary){ mCov<-expmtA%*%Int11%*%t(expmtA)}
	else{mCov<-Int11}    
    }
    mCov<-(mCov+t(mCov))/2            
    mCov[which(abs(mCov)<1e-15)]<-0
    colnames(mCov)<-NULL
    rownames(mCov)<-NULL
    if (b_correct_forPD && ((!.matrixcalc_is.symmetric.matrix(mCov)) || (!.matrixcalc_is.positive.definite(mCov)))){
        tryCatch({mCov<-as.matrix(Matrix::nearPD(mCov)$mat)},error=function(e){.my_message(e,FALSE)})
    }
    mCov
}

.calc.mean.slouch.mv<-function(expmtA,A1B,vY0,vX0,mPsi,mPsi0=NULL,exptjA=NULL,regimes=NULL){
## called in simulVasicekproc.R, simulVasicekprocphyl.R
## t is not needed to be passed, it is in the expmtA not explicitely
## vPsia has to be given as a vector
    kY<-length(vY0)
    vMean<-rbind(-A1B%*%vX0,as.matrix(vX0,ncol=1))
    if (!(is.null(exptjA))){
	## we have for length - 1 as we only take differences
	vPsiSum<-apply(matrix(sapply(1:(length(exptjA)-1),function(x){(exptjA[[x+1]]-exptjA[[x]])%*%mPsi[,regimes[x]]},simplify=TRUE),nrow=length(vY0)),1,sum) 
	vMean[1:kY,]<-vMean[1:kY,]+expmtA%*%(vY0-vMean[1:kY,])+vPsiSum
    }else{vMean[1:kY,]<-vMean[1:kY,]+expmtA%*%(vY0-mPsi[,1]-vMean[1:kY,])+mPsi[,1]} ## there is only one optimum ie one Psi
    if (!(is.null(mPsi0))){vMean[1:kY,]<-vMean[1:kY,]+(diag(1,kY,kY)-expmtA)%*%mPsi0}
    vMean[which(abs(vMean)<1e-15)]<-0
    names(vMean)<-NULL
    vMean
}


.calc.cov.slouch.mv<-function(t,lAcalcs,lScalcs,tol=1e-10,method="minus.v",b_correct_forPD=TRUE){
## called in simulVasicekproc.R, simulVasicekprocphyl.R
## once again invP can be calculated from eigA but more effective to do this once
## same with invA
## the idea of the function is just to glue everything together and only calculate time depedent bits
    A1B<-lAcalcs$A1B
    invA<-lAcalcs$invA
    kY<-nrow(A1B)
    kX<-ncol(A1B)
    tA1B<-t(A1B) 
    exptA<-.calc.exptA(t,lAcalcs)
    expmtA<-.calc.exptA(-t,lAcalcs)
    A1BS22tA1B<-A1B%*%lScalcs$S22%*%tA1B
## ----------------------------------------------------------------------------------
    if (method=="plus.v"){
	Int11<-.calc.integral.evAStevA(t,lScalcs$S11,lAcalcs)
	Int21<-.calc.integral.evAStevA(t,A1B%*%lScalcs$S21,lAcalcs)
	Int12<-.calc.integral.evAStevA(t,lScalcs$S12%*%tA1B,lAcalcs)
	Int22<-.calc.integral.evAStevA(t,A1BS22tA1B,lAcalcs) 
    }
    if (method=="minus.v"){
	Int11<-.calc.integral.evAStevA(-t,lScalcs$S11,lAcalcs)
	Int21<-.calc.integral.evAStevA(-t,A1B%*%lScalcs$S21,lAcalcs)
	Int12<-.calc.integral.evAStevA(-t,lScalcs$S12%*%tA1B,lAcalcs)
	Int22<-.calc.integral.evAStevA(-t,A1BS22tA1B,lAcalcs) 
    }
## ----------------------------------------------------------------------------------    

    CovXX<-t*lScalcs$S22
    CovYX<-(diag(1,nrow=kY,ncol=kY)-expmtA)%*%invA%*%(lScalcs$S12+A1B%*%lScalcs$S22)-t*A1B%*%lScalcs$S22
    CovXY<-t(CovYX)
    CovYY<-matrix(NA,nrow=nrow(A1B),ncol=nrow(A1B))
    if (method=="plus.v"){
	CovYY<-expmtA%*%(Int11+Int12+Int21+Int22)%*%t(expmtA)  - A1B%*%CovXY - CovYX%*%tA1B-t*A1BS22tA1B
    }
    if (method=="minus.v"){
	CovYY<-(Int11+Int12+Int21+Int22)  - A1B%*%CovXY - CovYX%*%tA1B-t*A1BS22tA1B
    }
    mCov<-rbind(cbind(CovYY,CovYX),cbind(CovXY,CovXX))

    mCov[which(abs(mCov)<1e-15)]<-0
    mCov<-(mCov+t(mCov))/2            
    colnames(mCov)<-NULL
    rownames(mCov)<-NULL
    if (b_correct_forPD&&((!.matrixcalc_is.symmetric.matrix(mCov)) || (!.matrixcalc_is.positive.definite(mCov)))){
        tryCatch({mCov<-as.matrix(Matrix::nearPD(mCov)$mat)},error=function(e){.my_message(e,FALSE)})
    }

    mCov
}



