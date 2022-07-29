## This file is part of mvSLOUCH

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this package or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .

.calc.phyl.cov<-function(EvolModel,modelParams,mTreeDist=NULL,vSpeciesTime=NULL,vSpeciesPairs=NULL,mAncestorTimes=NULL){
    mPhylCov=switch(EvolModel,
	bm=.bm.phyl.cov(mAncestorTimes,S=modelParams$Sxx), 
	ouch=.ouch.phyl.cov(mTreeDist,vSpeciesTime,modelParams,vSpeciesPairs), 
	slouch=.mvslouch.phyl.cov(mTreeDist,vSpeciesTime,modelParams,vSpeciesPairs),
	mvslouch=.mvslouch.phyl.cov(mTreeDist,vSpeciesTime,modelParams,vSpeciesPairs)
    )
    if ((is.element("Merror",names(modelParams)))&&(!(is.null(modelParams$Merror)))&&(!all(is.na(modelParams$Merror)))){mPhylCov<-mPhylCov+modelParams$Merror}
    mPhylCov[which(is.na(mPhylCov))]<-0
    mPhylCov[which(abs(mPhylCov)<1e-15)]<-0
    if ((!.matrixcalc_is.symmetric.matrix(mPhylCov)) || (!.matrixcalc_is.positive.definite(mPhylCov))){
        mPhylCov<-as.matrix(Matrix::nearPD(mPhylCov)$mat)
    }

    mPhylCov
}

.bm.phyl.cov<-function(mAncestorTimes,StS=NULL,S=NULL){
    if (is.null(StS)){StS<-S%*%t(S)}
    mPhylCov<-mAncestorTimes%x%StS
    if ((!.matrixcalc_is.symmetric.matrix(mPhylCov)) || (!.matrixcalc_is.positive.definite(mPhylCov))){
        mPhylCov<-as.matrix(Matrix::nearPD(mPhylCov)$mat)
    }
    mPhylCov<-(mPhylCov+t(mPhylCov))/2
    mPhylCov
}

.ouch.phyl.cov<-function(mTreeDist,vSpeciesDist,modelParams,vSpeciesPairs,ultrametric=FALSE){
    n<-nrow(mTreeDist)
    kY<-ncol(modelParams$A)
    mPhylCov<-matrix(NA,nrow=n*kY,ncol=n*kY)
    lCovMats<-sapply(vSpeciesPairs,function(spPair){
	s1<-(spPair-1)%/%n+1
	s2<-(spPair-1)%%n+1	
	mexpmt1A<-.calc.exptA((-1)*mTreeDist[s1,s2],modelParams$precalcMatrices[[1]])
	if(ultrametric){mexpmt2AT<-t(mexpmt1A)} ## if the tree is ultrametric then this will be the same
	else{mexpmt2AT<-t(.calc.exptA((-1)*mTreeDist[s2,s1],modelParams$precalcMatrices[[1]]))} 	
	mG<-.calc.cov.ouch.mv(vSpeciesDist[s1]-mTreeDist[s1,s2],modelParams$precalcMatrices[[1]],modelParams$precalcMatrices[[2]],method="minus.v")
	mexpmt1A%*%mG%*%mexpmt2AT		
    },simplify=FALSE)
    for (i in 1:length(vSpeciesPairs)){
	spPair<-vSpeciesPairs[i]
	s1<-(spPair-1)%/%n+1 
	s2<-(spPair-1)%%n+1	
	mPhylCov[((s1-1)*kY+1):(s1*kY),((s2-1)*kY+1):(s2*kY)]<-lCovMats[[i]]
	if(s1!=s2){mPhylCov[((s2-1)*kY+1):(s2*kY),((s1-1)*kY+1):(s1*kY)]<-t(lCovMats[[i]])}	
    }
    mPhylCov<-(mPhylCov+t(mPhylCov))/2
    if (min(Re(eigen(mPhylCov)$values))<0){
	mPhylCov<-matrix(NA,nrow=n*kY,ncol=n*kY)
	lCovMats<-sapply(vSpeciesPairs,function(spPair){
	    s1<-(spPair-1)%/%n+1
	    s2<-(spPair-1)%%n+1	
	    mexpmt1A<-.calc.exptA((-1)*mTreeDist[s1,s2],modelParams$precalcMatrices[[1]])
	    if(ultrametric){mexpmt2AT<-t(mexpmt1A)} ## if the tree is ultrametric then this will be the same
	    else{mexpmt2AT<-t(.calc.exptA((-1)*mTreeDist[s2,s1],modelParams$precalcMatrices[[1]]))} 	
	    mG<-.calc.cov.ouch.mv(vSpeciesDist[s1]-mTreeDist[s1,s2],modelParams$precalcMatrices[[1]],modelParams$precalcMatrices[[2]],method="plus.v")
	    mexpmt1A%*%mG%*%mexpmt2AT		
        },simplify=FALSE)
	for (i in 1:length(vSpeciesPairs)){
	    spPair<-vSpeciesPairs[i]
	    s1<-(spPair-1)%/%n+1 
	    s2<-(spPair-1)%%n+1	
	    mPhylCov[((s1-1)*kY+1):(s1*kY),((s2-1)*kY+1):(s2*kY)]<-lCovMats[[i]]
	    if(s1!=s2){mPhylCov[((s2-1)*kY+1):(s2*kY),((s1-1)*kY+1):(s1*kY)]<-t(lCovMats[[i]])}	
	}
	mPhylCov<-(mPhylCov+t(mPhylCov))/2
    }
    if ((!.matrixcalc_is.symmetric.matrix(mPhylCov)) || (!.matrixcalc_is.positive.definite(mPhylCov))){
        mPhylCov<-as.matrix(Matrix::nearPD(mPhylCov)$mat)
    }
    mPhylCov
}

.mvslouch.phyl.cov<-function(mTreeDist,vSpeciesDist,modelParams,vSpeciesPairs,ultrametric=FALSE){
    ## the mTreeDist matrix is not symmetrical so we have to have a convention for it (tree is not ultrametric)
    ## perhaps do some simplification in case of ultrametricity
    ## if not then structure M[i,j] the time from species i to the last common ancestor of species j
    ## vSpeciesPairs is just the set of indexes but as it will always be used it is precalculated
    n<-nrow(mTreeDist)
    kY<-ncol(modelParams$A)
    kX<-ncol(modelParams$B)
    lCovMats<-sapply(vSpeciesPairs,function(spPair){
	s1<-(spPair-1)%/%n+1
	s2<-(spPair-1)%%n+1	
	mCovta<-.calc.cov.slouch.mv(vSpeciesDist[s1]-mTreeDist[s1,s2],modelParams$precalcMatrices[[1]],modelParams$precalcMatrices[[2]],method="minus.v")
	mexpmt1A<-.calc.exptA((-1)*mTreeDist[s1,s2],modelParams$precalcMatrices[[1]]) 
	if(ultrametric){mexpmt2AT<-t(mexpmt1A)} ## if the tree is ultrametric then this will be the same
	else{mexpmt2AT<-t(.calc.exptA((-1)*mTreeDist[s2,s1],modelParams$precalcMatrices[[1]]))} 

	mCovY1Y2<-mexpmt1A%*%mCovta[1:kY,1:kY]%*%mexpmt2AT + mexpmt1A%*%mCovta[1:kY,(kY+1):(kY+kX)]%*%t(modelParams$precalcMatrices[[1]]$A1B)%*%(mexpmt2AT-diag(x=1,nrow=kY,ncol=kY)) + (mexpmt1A-diag(x=1,nrow=kY,ncol=kY))%*%modelParams$precalcMatrices[[1]]$A1B%*%mCovta[(kY+1):(kY+kX),1:kY]%*%mexpmt2AT + (mexpmt1A-diag(x=1,nrow=kY,ncol=kY))%*%modelParams$precalcMatrices[[1]]$A1B%*%mCovta[(kY+1):(kY+kX),(kY+1):(kY+kX)]%*%t(modelParams$precalcMatrices[[1]]$A1B)%*%(mexpmt2AT-diag(x=1,nrow=kY,ncol=kY))
	mCovY1T2<-mexpmt1A%*%mCovta[1:kY,(kY+1):(kY+kX)] + (mexpmt1A-diag(x=1,nrow=kY,ncol=kY))%*%modelParams$precalcMatrices[[1]]$A1B%*%mCovta[(kY+1):(kY+kX),(kY+1):(kY+kX)]
	mCovT1Y2<-mCovta[(kY+1):(kY+kX),1:kY]%*%mexpmt2AT + mCovta[(kY+1):(kY+kX),(kY+1):(kY+kX)]%*%t(modelParams$precalcMatrices[[1]]$A1B)%*%(mexpmt2AT-diag(x=1,nrow=kY,ncol=kY))
	mCovT1T2<-mCovta[(kY+1):(kY+kX),(kY+1):(kY+kX)]
	rbind(cbind(mCovY1Y2,mCovY1T2),cbind(mCovT1Y2,mCovT1T2))
    },simplify =FALSE)
    ## Now we have a list with all of the matrix components needed -> we need to make a matrix out of it ....
    mPhylCov<-matrix(0,nrow=n*(kY+kX),ncol=n*(kY+kX))    
    for (i in 1:length(vSpeciesPairs)){
	spPair<-vSpeciesPairs[i]
	s1<-(spPair-1)%/%n+1 
	s2<-(spPair-1)%%n+1	
	mPhylCov[((s1-1)*(kY+kX)+1):(s1*(kY+kX)),((s2-1)*(kY+kX)+1):(s2*(kY+kX))]<-lCovMats[[i]]
	if(s1!=s2){mPhylCov[((s2-1)*(kY+kX)+1):(s2*(kY+kX)),((s1-1)*(kY+kX)+1):(s1*(kY+kX))]<-t(lCovMats[[i]])}
    }
    mPhylCov<-(mPhylCov+t(mPhylCov))/2
    if (min(Re(eigen(mPhylCov)$values))<0){
	lCovMats<-sapply(vSpeciesPairs,function(spPair){
	    s1<-(spPair-1)%/%n+1
	    s2<-(spPair-1)%%n+1	
	    mCovta<-.calc.cov.slouch.mv(vSpeciesDist[s1]-mTreeDist[s1,s2],modelParams$precalcMatrices[[1]],modelParams$precalcMatrices[[2]],method="plus.v")
	    mexpmt1A<-.calc.exptA((-1)*mTreeDist[s1,s2],modelParams$precalcMatrices[[1]]) 
	    if(ultrametric){mexpmt2AT<-t(mexpmt1A)} ## if the tree is ultrametric then this will be the same
	    else{mexpmt2AT<-t(.calc.exptA((-1)*mTreeDist[s2,s1],modelParams$precalcMatrices[[1]]))} 

	    mCovY1Y2<-mexpmt1A%*%mCovta[1:kY,1:kY]%*%mexpmt2AT + mexpmt1A%*%mCovta[1:kY,(kY+1):(kY+kX)]%*%t(modelParams$precalcMatrices[[1]]$A1B)%*%(mexpmt2AT-diag(x=1,nrow=kY,ncol=kY)) + (mexpmt1A-diag(x=1,nrow=kY,ncol=kY))%*%modelParams$precalcMatrices[[1]]$A1B%*%mCovta[(kY+1):(kY+kX),1:kY]%*%mexpmt2AT + (mexpmt1A-diag(x=1,nrow=kY,ncol=kY))%*%modelParams$precalcMatrices[[1]]$A1B%*%mCovta[(kY+1):(kY+kX),(kY+1):(kY+kX)]%*%t(modelParams$precalcMatrices[[1]]$A1B)%*%(mexpmt2AT-diag(x=1,nrow=kY,ncol=kY))
	    mCovY1T2<-mexpmt1A%*%mCovta[1:kY,(kY+1):(kY+kX)] + (mexpmt1A-diag(x=1,nrow=kY,ncol=kY))%*%modelParams$precalcMatrices[[1]]$A1B%*%mCovta[(kY+1):(kY+kX),(kY+1):(kY+kX)]
	    mCovT1Y2<-mCovta[(kY+1):(kY+kX),1:kY]%*%mexpmt2AT + mCovta[(kY+1):(kY+kX),(kY+1):(kY+kX)]%*%t(modelParams$precalcMatrices[[1]]$A1B)%*%(mexpmt2AT-diag(x=1,nrow=kY,ncol=kY))
	    mCovT1T2<-mCovta[(kY+1):(kY+kX),(kY+1):(kY+kX)]
	    rbind(cbind(mCovY1Y2,mCovY1T2),cbind(mCovT1Y2,mCovT1T2))
	},simplify =FALSE)
	## Now we have a list with all of the matrix components needed -> we need to make a matrix out of it ....
	mPhylCov<-matrix(0,nrow=n*(kY+kX),ncol=n*(kY+kX))    
	for (i in 1:length(vSpeciesPairs)){
	    spPair<-vSpeciesPairs[i]
	    s1<-(spPair-1)%/%n+1 
	    s2<-(spPair-1)%%n+1	
	    mPhylCov[((s1-1)*(kY+kX)+1):(s1*(kY+kX)),((s2-1)*(kY+kX)+1):(s2*(kY+kX))]<-lCovMats[[i]]
	    if(s1!=s2){mPhylCov[((s2-1)*(kY+kX)+1):(s2*(kY+kX)),((s1-1)*(kY+kX)+1):(s1*(kY+kX))]<-t(lCovMats[[i]])}
	}
	mPhylCov<-(mPhylCov+t(mPhylCov))/2
    }
    if ((!.matrixcalc_is.symmetric.matrix(mPhylCov)) || (!.matrixcalc_is.positive.definite(mPhylCov))){
        mPhylCov<-as.matrix(Matrix::nearPD(mPhylCov)$mat)
    }
    mPhylCov
}
