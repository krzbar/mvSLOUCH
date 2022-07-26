## This file is part of mvSLOUCH

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this package or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .


.decompEigenA.S<-function(modelParams,lPrecalculates,designToEstim,toCalc=list(bCalcA=TRUE,bCovCalc=FALSE,dzetacalc=FALSE,lexptcalc=FALSE,kappacalc=FALSE,interceptcalc=FALSE),mXmX0=NULL){
## function precalculates all the matrices needed to calculate the covariance
## ie all the bits that are not time dependent
## the function does not consider the situation A=0 ie degenerate to BM
## called in estimGLSGC.R, getESS.R, loglik.R, matrixexps.R (not lPrecalculates), modelparams.R, OUphylregression.R
## modelparamsummary.R (not lPrecalculates), precalcs.R, simulVasicekproc.R, simulVasicekprocphyl.R

##---------------------------Prepare all the variables ----------------------------------
    regimes<-modelParams$regimes
    regimeTimes<-modelParams$regimeTimes
    kY<-NA;kX<-NA
    if (is.element("A",names(modelParams))){A<-modelParams$A;if((!any(is.na(A)))&&(is.matrix(A))){kY<-nrow(A)}}else{A<-NULL}
    if (is.element("B",names(modelParams))){B<-modelParams$B;if((!any(is.na(B)))&&(is.matrix(B))){kX<-ncol(B)}}else{B<-NULL}
    if (is.element("Syy",names(modelParams))){Syy<-modelParams$Syy}else{Syy<-NULL}
    if (is.element("Syx",names(modelParams))){Syx<-modelParams$Syx}else{Syx<-NULL}
    if (is.element("Sxy",names(modelParams))){Sxy<-modelParams$Sxy}else{Sxy<-NULL}
    if (is.element("Sxx",names(modelParams))){Sxx<-modelParams$Sxx}else{Sxx<-NULL}
    if (is.element("mPsi",names(modelParams))){mPsi<-modelParams$mPsi}else{mPsi<-NULL}
    if (is.element("mPsi0",names(modelParams))){mPsi0<-modelParams$mPsi0}else{mPsi0<-NULL}
    if (is.element("vY0",names(modelParams))){vY0<-modelParams$vY0}else{vY0<-NULL}
    if (is.element("vX0",names(modelParams))){vX0<-modelParams$vX0}else{vX0<-NULL}
    if (!is.null(lPrecalculates)){
	if (is.element("vSpecies_times",names(lPrecalculates))){vSpecies_times<-lPrecalculates$vSpecies_times}else{vSpecies_times<-NULL}
	if (!is.element("vSpecies_times",names(lPrecalculates))){
	    if (is.element("mSpecDist",names(lPrecalculates))){vSpecies_times<-lPrecalculates$mSpecDist[nrow(lPrecalculates$mSpecDist),]}else{vSpecies_times<-NULL}
	}
	if (is.element("mTreeDist",names(lPrecalculates))){mTreeDist<-lPrecalculates$mTreeDist}else{mTreeDist<-NULL}
	if (is.element("invmAncestorTimes",names(lPrecalculates))){invmAncestorTimes<-lPrecalculates$invmAncestorTimes}else{invmAncestorTimes<-NULL}
	if (is.element("mAncestorTimes",names(lPrecalculates))){mAncestorTimes<-lPrecalculates$mAncestorTimes}else{mAncestorTimes<-NULL}
	if (is.element("vSpeciesPairs",names(lPrecalculates))){vSpeciesPairs<-lPrecalculates$vSpeciesPairs}else{vSpeciesPairs<-NULL}
    }
    bCovCalc<-toCalc$bCovCalc
    dzetacalc<-toCalc$dzetacalc    
    lexptcalc<-toCalc$lexptcalc
    kappacalc<-toCalc$kappacalc    
    bcalcA<-toCalc$bCalcA 
    interceptcalc<-toCalc$interceptcalc
    if (!is.element("bCalcG",names(toCalc))){bCalcG<-FALSE}else{bCalcG<-toCalc$bCalcG}
## --------------------------------------------------------------------------------    
    lReturn<-vector("list",6)
    lReturn[[6]]<-list(NA,NA)
    if (is.element("precalcMatrices",names(modelParams))){lReturn<-modelParams$precalcMatrices;ldecompEigenA.precalc<-lReturn[[1]];mKappa<-lReturn[[4]]$mKappa;S12<-lReturn[[2]]$S12;invS22<-lReturn[[2]]$invS22}
    else{modelParams$precalcMatrices<-list();ldecompEigenA.precalc=NULL}
    
    if ((bcalcA)&&(!is.null(A))){
	ldecompEigenA.precalc<-vector("list",7)
	names(ldecompEigenA.precalc)<-c("A","invA","A1B","eigA","invP","decomp","TwoByTwo")
	ldecompEigenA.precalc$A<-A
	eigA<-eigen(A)
	if (!isTRUE(all.equal(rcond(A),0))){
	    ldecompEigenA.precalc$invA<-solve(A)
	    if((!(is.null(B)))&&(!any(is.na(B)))){ldecompEigenA.precalc$A1B<-ldecompEigenA.precalc$invA%*%B}
	}else{ldecompEigenA.precalc$invA<-matrix(NA,nrow(A),nrow(A));if((!(is.null(B)))&&(!any(is.na(B)))){ldecompEigenA.precalc$A1B<-matrix(NA,nrow(A),ncol(B))}}
	if (!isTRUE(all.equal(rcond(eigA$vectors),0))){ldecompEigenA.precalc$invP<-solve(eigA$vectors)}else{.my_stop("ERROR: the matrix A does not have an eigendecomposition.",FALSE)}
	if (any(class(ldecompEigenA.precalc$invP)=="matrix")){## the matrix of eigenvectors is invertible	
	## if there was an error we should have an object of class try-error
	    ldecompEigenA.precalc$eigA<-eigA
	    ldecompEigenA.precalc$decomp<-TRUE    
        }else{
    	    .my_message("Warning: the matrix A does not have an eigendecomposition. \n",FALSE)
    	    ldecompEigenA.precalc$decomp<-FALSE;if (nrow(A)==2){ldecompEigenA.precalc$TwoByTwo<-TRUE}else{ldecompEigenA.precalc$TwoByTwo<-FALSE}
	}
	lReturn[[1]]<-ldecompEigenA.precalc
	if (is.element("precalcMatrices",names(modelParams))){modelParams$precalcMatrices[[1]]<-lReturn[[1]]}
    }

    if (bCovCalc){
	lSs<-vector("list",5)
	names(lSs)<-c("S11","S12","S21","S22","invS22")
	if (!(is.null(Syy)||(any(is.na(Syy))))){lSs$S11<-Syy%*%t(Syy)}
        if (!((is.null(lSs$S11))||(is.null(Syx))||(is.null(Sxy))||(any(is.na(Sxy)))||(any(is.na(Sxy))))){lSs$S11<-lSs$S11+Syx%*%t(Syx)}
        if (!((is.null(Syy))||(is.null(Syx))||(is.null(Sxy))||(is.null(Sxx))||(any(is.na(Syy)))||(any(is.na(Sxy)))||(any(is.na(Sxy)))||(any(is.na(Sxx))))){lSs$S12<-Syy%*%t(Sxy)+Syx%*%t(Sxx)}
        if (!((is.null(Syy))||(is.null(Syx))||(is.null(Sxy))||(is.null(Sxx))||(any(is.na(Syy)))||(any(is.na(Sxy)))||(any(is.na(Sxy)))||(any(is.na(Sxx))))){lSs$S21<-Sxy%*%t(Syy)+Sxx%*%t(Syx)}
        if (!((is.null(Sxx))||(any(is.na(Sxx))))){
    	    lSs$S22<-Sxx%*%t(Sxx)
	    if (!((is.null(Syx))||(is.null(Sxy))||(any(is.na(Sxy)))||(any(is.na(Sxy))))){lSs$S22<- lSs$S22+Sxy%*%t(Sxy)}
	    lSs$invS22<-solve(lSs$S22)
	    lSs$invS22[which(abs(lSs$invS22)<1e-15)]<-0
	}
	invS22<-lSs$invS22
	S12<-lSs$S12	    
	
	lReturn[[2]]<-lSs
	if (is.element("precalcMatrices",names(modelParams))){modelParams$precalcMatrices[[2]]<-lReturn[[2]]}
    }

    if (lexptcalc){    
        lexpmtA<-sapply(vSpecies_times,function(t){.calc.exptA(-t,ldecompEigenA.precalc)},simplify=FALSE) 
        if (!(is.null(regimeTimes))){
    	    lexptjA<-sapply(1:length(regimeTimes),function(i,regimeTimes,vSpecDist){tjs<-regimeTimes[[i]];specT<-vSpecDist[i];sapply(tjs,function(t,specT){.calc.exptA(t-specT,ldecompEigenA.precalc)},specT=specT,simplify=FALSE)},regimeTimes=regimeTimes,vSpecDist=vSpecies_times,simplify=FALSE)}
        else{ 
    	    lexptjA<-sapply(1:length(vSpecies_times),function(i,k){list(lexpmtA[[i]],diag(1,nrow=k,ncol=k))},k=nrow(A),simplify=FALSE)
        }   
	lReturn[[3]]<-list(lexpmtA=lexpmtA,lexptjA=lexptjA)
	if (is.element("precalcMatrices",names(modelParams))){modelParams$precalcMatrices[[3]]<-lReturn[[3]]}
    }
    lReturn[[4]]<-list(lDzeta=NA,lDzetaKappa=NA,mKappa=NA)

    
    if ((dzetacalc)&&(!designToEstim$YnonCondX)&&(designToEstim$B)){   
       ## precalc all dzeta matrix and DzetaKappa matrix
       ## here A HAS to be invertible for this to work
       bSimpReg<-designToEstim$SimpReg
       n<-length(vSpecies_times)
       kY<-nrow(A)
       kX<-nrow(mXmX0)
       lDzeta<-sapply(vSpecies_times,.dzeta.matrix.t,"precalcs"=ldecompEigenA.precalc,"A"=A,"invA"=ldecompEigenA.precalc$invA,simplify=FALSE)
       lDzetaIJ<-NULL
       if (!bSimpReg){  
    	    vNAX<-which(is.na(c(mXmX0)))
    	    if ((length(vNAX)>0)&&(designToEstim$BFullXNA)){## if there are missing values this is more complicated		
    		mKappa<-.Kappa.matrix.NA(vNAX,mAncestorTimes,lSs$S22,lSs$S22,n,kX,mXmX0)
    	    }else{mKappa<-.Kappa.matrix(invmAncestorTimes,mXmX0,n,diag(1,kX,kX))}
	    if (kX==1){mKappa<-matrix(mKappa,ncol=n,nrow=1)}
    	    if (is.null(B)||any(is.na(B))){
    		lDzetaKappa<-sapply(1:n,function(i,mTreeDist,mAncestorTimes,precalcs,invA,mKappa,n,kY,kX){
                	    .dzetaKappa.matrix(mTreeDist[i,],mAncestorTimes[i,],precalcs,invA,mKappa,n,kY,kX)},
                	mTreeDist=mTreeDist,mAncestorTimes=mAncestorTimes,precalcs=ldecompEigenA.precalc,invA=ldecompEigenA.precalc$invA,mKappa=mKappa,n=n,kY=kY,kX=kX,simplify=FALSE)
    		if (bCalcG){
    	    	    lDzetaIJ<-sapply(1:n,function(i,mTreeDist,mAncestorTimes,precalcs,invA,n,kY,kX){
    	        	    .dzetaIJ.matrix(mTreeDist[i,],mAncestorTimes[i,],precalcs,invA,n,kY,kX)},
    	                mTreeDist=mTreeDist,mAncestorTimes=mAncestorTimes,precalcs=ldecompEigenA.precalc,invA=ldecompEigenA.precalc$invA,n=n,kY=kY,kX=kX,simplify=FALSE)
                }    	                                                                                                    
    	    }else{lDzetaKappa<-NA}
       }else{
    	   mKappa<-NA
	   if (is.null(B)||any(is.na(B))){	
    		lDzetaKappa<-sapply(vSpecies_times,.dzetaKappa.matrix.simp.t,"precalcs"=ldecompEigenA.precalc,"A"=A,"invA"=ldecompEigenA.precalc$invA,simplify=FALSE)    
	   }else{lDzetaKappa<-NA}
       }
       lReturn[[4]]<-list(lDzeta=lDzeta,lDzetaKappa=lDzetaKappa,mKappa=mKappa,lDzetaIJ=lDzetaIJ)
       if (is.element("precalcMatrices",names(modelParams))){modelParams$precalcMatrices[[4]]<-lReturn[[4]]}
    }
    intercept<-NA
    mKappaIntercept<-NA
    if (interceptcalc){
    ## intercept is calculated as a vector, i.e. the response is row stacked
	kY<-nrow(lSs$S11)      
	if (!is.null(mPsi) && (!is.na(mPsi[1]))){vAncPsi<-matrix(mPsi[,designToEstim$y0Regime],ncol=1,nrow=kY)} ## done here so no needless passing of designToEstim
	else {vAncPsi<-matrix(NA,ncol=1,nrow=kY)}
	if (!is.null(mPsi0) &&(!is.na(mPsi0[1]))){vAncPsi<-vAncPsi+mPsi0}
        bSimpReg<-designToEstim$SimpReg
        n<-length(lPrecalculates$vSpecies_times)
        
	intercept<-rep(0,n*kY)
	if (!designToEstim$y0){intercept<-intercept+c(sapply(modelParams$precalcMatrices[[3]]$lexpmtA,function(expmtA,vY0){expmtA%*%vY0},vY0=vY0,simplify=TRUE))}
	else{
	    if (!designToEstim$psi && designToEstim$y0AncState){intercept<-intercept+c(sapply(modelParams$precalcMatrices[[3]]$lexpmtA,function(expmtA,AncPsi){expmtA%*%AncPsi},AncPsi=vAncPsi,simplify=TRUE))}
	    if (!is.null(mPsi0) && designToEstim$psi0 && !designToEstim$psi && designToEstim$y0AncState){intercept<-intercept+c(sapply(modelParams$precalcMatrices[[3]]$lexpmtA,function(expmtA,mPsi0){expmtA%*%mPsi0},mPsi0=mPsi0,simplify=TRUE))}	
	}
	
	if (!designToEstim$psi){
    	    intercept<-intercept+c(sapply(1:n,function(i,mPsi,lexptjA){
    		vRegs<-regimes[[i]]
		Reduce('+',c(sapply(1:(length(vRegs)),function(reg,mPsi,mexptjA){(mexptjA[[reg+1]]-mexptjA[[reg]])%*%mPsi[,vRegs[reg]]},mPsi=mPsi,mexptjA=lexptjA[[i]],simplify=TRUE)))    		
    	    },mPsi=mPsi,lexptjA=modelParams$precalcMatrices[[3]]$lexptjA,simplify=TRUE))
    	}
	if (!is.null(mPsi0) && (!designToEstim$psi0)){
	    intercept<-intercept+c(sapply(modelParams$precalcMatrices[[3]]$lexpmtA,function(expmtA,mPsi0){(diag(1,ncol(expmtA),nrow(expmtA))-expmtA)%*%mPsi0},mPsi0=mPsi0,simplify=TRUE))	
    	}   
    	 	
    	if (!is.null(mXmX0)){ ## this is checking whether we have have mvslouch or ouch
	    if(designToEstim$B){ ## we do not know B matrix
    		if (!bSimpReg){
    		    vNAX<-which(is.na(c(mXmX0)))
    		    if ((length(vNAX)>0)&&(designToEstim$BFullXNA)){## if there are missing values this is more complicated		
    			mKappaIntercept<-.Kappa.matrix.NA(vNAX,mAncestorTimes,lSs$S22,lSs$S12,n,kX,mXmX0)
    		    }else{mKappaIntercept<-.Kappa.matrix(invmAncestorTimes,mXmX0,n,lSs$S12%*%lSs$invS22)}
    		    if (kY==1){mKappaIntercept<-matrix(mKappaIntercept,ncol=n,nrow=1)}
    		    intercept<-intercept+c(sapply(1:n,function(i,mTreeDist,mAncestorTimes,precalcs,invA,mKappa,n,kY){
					.dzetaKappa.intercept(mTreeDist[i,],mAncestorTimes[i,],precalcs,invA,mKappa,n,kY)
				},mTreeDist=mTreeDist,mAncestorTimes=mAncestorTimes,precalcs=ldecompEigenA.precalc,invA=ldecompEigenA.precalc$invA,mKappa=mKappaIntercept,n=n,kY=kY,simplify=TRUE))
    		}else{
    		    intercept<-intercept+c(sapply(1:n,function(i,precalcs,invA,S12,invS22,mXmX0,kY){
				    .dzetaKappa.simp.intercept(vSpecies_times[i],
				"precalcs"=precalcs,"invA"=invA,S12=S12,invS22=invS22,mXmX0i=mXmX0[,i],kY=kY,bFullNA=designToEstim$BFullXNA)},
    		     "precalcs"=ldecompEigenA.precalc,"invA"=ldecompEigenA.precalc$invA,S12=S12,invS22=invS22,mXmX0=mXmX0,kY=kY,simplify=TRUE))      
    		}
	    }else{## we know B matrix
		if (designToEstim$UseX0){intercept<-intercept+c(sapply(modelParams$precalcMatrices[[3]]$lexpmtA,function(mexpmtA,kY,A1B,vX0){(mexpmtA-diag(1,kY,kY))%*%A1B%*%vX0},kY=kY,A1B=modelParams$precalcMatrices[[1]]$A1B,vX0=vX0,simplify=TRUE))}
		if (designToEstim$YnonCondX){## we do the GLS non-conditional on the BM predictors. This can only have a positive impact if B is not in the GLS		    		    
		    if (designToEstim$y0 && designToEstim$UseX0 && designToEstim$y0AncState){intercept<-intercept+c(sapply(modelParams$precalcMatrices[[3]]$lexpmtA,function(expmtA,A1B,vX0){(-1)*expmtA%*%A1B%*%vX0},A1B=modelParams$precalcMatrices[[1]]$A1B,vX0=vX0,simplify=TRUE))}
		}else{## we do the GLS conditional on BM predictors		    
    		    if (designToEstim$UseX0 && !designToEstim$y0OnlyFixed){intercept<-intercept+c(sapply(modelParams$precalcMatrices[[3]]$lexpmtA,function(mexpmtA,A1B,vX0){mexpmtA%*%A1B%*%vX0},A1B=modelParams$precalcMatrices[[1]]$A1B,vX0=vX0,simplify=TRUE))}
    		    if (!bSimpReg){
    		        mCovPhyl<-.calc.phyl.cov(mTreeDist,vSpecies_times,NULL,vSpeciesPairs,"mvslouch",modelParams)
			lReturn[[6]][[1]]<-mCovPhyl
    			SXX<-mCovPhyl[c(sapply(((1:n)-1)*(kY+kX),function(x,v){x+v},"v"=(kY+1):(kY+kX))),c(sapply(((1:n)-1)*(kY+kX),function(x,v){x+v},"v"=(kY+1):(kY+kX)))]
			SYX<-mCovPhyl[c(sapply(((1:n)-1)*(kY+kX),function(x,v){x+v},"v"=1:kY)),c(sapply(((1:n)-1)*(kY+kX),function(x,v){x+v},"v"=(kY+1):(kY+kX)))]
			vNAX<-which(is.na(c(mXmX0)))
			if (length(vNAX)>0){SYX<-SYX[,-vNAX,drop=FALSE];SXX<-SXX[-vNAX,-vNAX,drop=FALSE];vXmX0<-c(mXmX0)[-vNAX]}
			else{vXmX0<-c(mXmX0)}
			invSXX<-solve(SXX)		
			lReturn[[6]][[2]]<-invSXX
			intercept<-intercept+SYX%*%invSXX%*%vXmX0 
		    }else{
			test_t<-sum(abs(vSpecies_times-vSpecies_times[1]))
			phyl_ultra<-FALSE
			if (isTRUE(all.equal(test_t,0))){phyl_ultra<-TRUE}
			if(phyl_ultra){vtmpSpecies_times<-vSpecies_times[1]}else{vtmpSpecies_times<-vSpecies_times}
			lScovs<-sapply(vtmpSpecies_times,function(spec_time,modelParams){
    				mS<-.calc.cov.slouch.mv(spec_time,modelParams$precalcMatrices[[1]],modelParams$precalcMatrices[[2]])    	    
				syt<-mS[1:kY,(kY+1):(kY+kX)]
				invsxx<-solve(mS[(kY+1):(kY+kX),(kY+1):(kY+kX)])
				list(syt=syt,invsxx=invsxx,mS=mS)
    			    },modelParams=modelParams,simplify=FALSE)
    			lReturn[[6]][[1]]<-lScovs
			if (phyl_ultra){syt<-lScovs[[1]]$syt;invsxx<-lScovs[[1]]$invsxx}
			for (i in 1:n){
			    if (!phyl_ultra){
				syt<-lScovs[[i]]$syt
				invsxx<-lScovs[[i]]$invsxx
			    }
			    vX<-c(mXmX0)[((i-1)*kX+1):(i*kX)] 
			    vNAX<-vX
			    if (length(vNAX)>0){
				if (length(vNAX)<kX){
				    vX<-vX[-vNAX];
				    syt<-syt[,-vNAX,drop=FALSE]
				    invsxx<-invsxx[-vNAX,-vNAX,drop=FALSE]
				}
				else{vX<-rep(0,kX)}
			    }
			    intercept[((i-1)*kY+1):(i*kY)]<-intercept[((i-1)*kY+1):(i*kY)]+syt%*%invsxx%*%vX
			}
		    }
    		}
    	    }
    	}
	intercept[which(abs(intercept)<1e-15)]<-0
    }

    lReturn[[5]]<-list(intercept=intercept,mKappaIntercept=mKappaIntercept)
    lReturn
}

.dzeta.matrix.t<-function(t,precalcs,invA,A=NULL){
    ## remove after working and copied to FilesForGLSforB/
    ExpmtA<-.calc.exptA(t=-t,precalcs)
    M<-(ExpmtA-diag(1,nrow=nrow(ExpmtA),ncol=nrow(ExpmtA)))%*%invA
    M[which(abs(M)<1e-15)]<-0
    M
}


.dzeta.matrix.tv2<-function(t,precalcs,invA,A=NULL){
    ## remove after working and copied to FilesForGLSforB/
    ExpmtA<-.calc.exptA(t=-t,precalcs)
    M<-ExpmtA%*%invA
    M[which(abs(M)<1e-15)]<-0
    M
}

.dzetaKappa.matrix<-function(tiij,taij,precalcs,invA,mKappa,n,kY,kX){
    ## remove after working and copied to FilesForGLSforB/
## at the moment we do not allow precalcs to be NULL
    sumDzetaKappa<-sapply(1:n,function(j,tiij,taij,precalcs,invA,mKappa,kY){
	((.calc.exptA(t=-tiij[j],precalcs)%*%.calc.exptA(t=-taij[j],precalcs)%*%invA%*%
	(.calc.exptA(t=taij[j],precalcs)-diag(1,kY,kY))-taij[j]*diag(1,kY,kY))%*%invA)%x%matrix(mKappa[,j],nrow=1)},
	tiij=tiij,taij=taij,precalcs=precalcs,invA=invA,mKappa=mKappa,kY=kY,simplify=FALSE)
    M<-Reduce('+',sumDzetaKappa,matrix(0,ncol=kY*kX,nrow=kY))    
    M[which(abs(M)<1e-15)]<-0
    M
}

.dzetaIJ.matrix<-function(tiij,taij,precalcs,invA,n,kY,kX){
    ## remove after working and copied to FilesForGLSforB/
## at the moment we do not allow precalcs to be NULL
    DzetaJ<-sapply(1:n,function(j,tiij,taij,precalcs,invA,kY){
	((.calc.exptA(t=-tiij[j],precalcs)%*%.calc.exptA(t=-taij[j],precalcs)%*%invA%*%
	(.calc.exptA(t=taij[j],precalcs)-diag(1,kY,kY))-taij[j]*diag(1,kY,kY))%*%invA)},
	tiij=tiij,taij=taij,precalcs=precalcs,invA=invA,kY=kY,simplify=FALSE)
    DzetaJ
}

.dzetaKappa.matrix.simp.t<-function(t,precalcs,invA,A=NULL){
    ## remove after working and copied to FilesForGLSforB/
    ExpmtA<-.calc.exptA(t=-t,precalcs)
    ExptA<-.calc.exptA(t=t,precalcs)
    M<-ExpmtA%*%invA%*%(ExptA-diag(1,nrow=nrow(ExptA),ncol=nrow(ExptA)))%*%invA*(1/t)-invA
    M[which(abs(M)<1e-15)]<-0
    M
}


.dzetaKappa.simp.intercept<-function(t,precalcs,invA,S12,invS22,mXmX0i,kY,bFullNA){
    ## remove after working and copied to FilesForGLSforB/
    corrNA<-matrix(0,nrow=nrow(S12),ncol=ncol(S12))
    vNAX<-which(is.na(mXmX0i))
    mXmX0i<-matrix(mXmX0i,ncol=1,nrow=length(mXmX0i))
    if ((length(vNAX)>0)&&bFullNA){
	    if (length(vNAX)<length(mXmX0i)){mXmX0i<-mXmX0i[vNAX];
	    S12<-S12[,-vNAX,drop=FALSE]
	    ##if(nrow(S12)>1){S12<-S12[,-vNAX]}else{S12<-matrix(S12[,-vNAX],nrow=1)};
	    corrNA<-S12%*%invS22[-vNAX,-vNAX,drop=FALSE]%*%mXmX0i}
    }else{if (length(vNAX)==0){corrNA<-S12%*%invS22%*%mXmX0i}}
    M<-(1/t)*.calc.exptA(t=-t,precalcs)%*%invA%*%(.calc.exptA(t=t,precalcs)-diag(1,kY,kY))%*%corrNA
    M[which(abs(M)<1e-15)]<-0
    M
}

.dzetaKappa.intercept<-function(tiij,taij,precalcs,invA,mKappa,n,kY,invm){
    ## remove after working and copied to FilesForGLSforB/
## at the moment we do not allow precalcs to be NULL
    sumDzetaKappaInterc<-sapply(1:n,function(j,tiij,taij,precalcs,mKappa,kY){
	(.calc.exptA(t=-tiij[j],precalcs)%*%.calc.exptA(t=-taij[j],precalcs)%*%invA%*%(.calc.exptA(t=taij[j],precalcs)-diag(1,kY,kY)))%*%mKappa[,j]},
	tiij=tiij,taij=taij,precalcs=precalcs,mKappa=mKappa,kY=kY,simplify=FALSE)
    M<-Reduce('+',sumDzetaKappaInterc,rep(0,kY))    
    M[which(abs(M)<1e-15)]<-0
    M
}

.dzetaKappa.matrix2<-function(tiij,taij,precalcs,invA,mKappa,n,kY,kX){
    ## remove after working and copied to FilesForGLSforB/
## at the moment we do not allow precalcs to be NULL
    M<-sapply(1:n,function(j,tiij,taij,precalcs,invA,mKappa,kY){
	((.calc.exptA(t=-tiij[j],precalcs)%*%invA%*%
	(.calc.exptA(t=taij[j],precalcs)-diag(1,kY,kY))-diag(taij[j],kY,kY))%*%invA)},
	tiij=tiij,taij=taij,precalcs=precalcs,invA=invA,mKappa=mKappa,kY=kY,simplify=FALSE)
    M[which(abs(M)<1e-15)]<-0
    M
}

.dzeta.matrix.t1t2<-function(t2,t1,precalcs,invA,A=NULL){
    ## remove after working and copied to FilesForGLSforB/
    M<-.calc.exptA(t=-t1,precalcs)%*%invA%*%.dzeta.matrix.t(-t2,precalcs,invA,A)
    M[which(abs(M)<1e-15)]<-0
    M
}

.dzeta.matrix.t1t2v2<-function(t2,t1,precalcs,invA,A=NULL){
    ## remove after working and copied to FilesForGLSforB/
    M<-.calc.exptA(t=-t1,precalcs)%*%invA%*%(.calc.exptA(t=t2,precalcs)-diag(1,nrow=nrow(invA),ncol=ncol(invA)))%*%invA
    M[which(abs(M)<1e-15)]<-0
    M
}

.Kappa.matrix<-function(invT,mXmX0,n,S1S2){
    ## remove after working and copied to FilesForGLSforB/
    mXmX0<-apply(mXmX0,2,function(x,S1S2){S1S2%*%x},S1S2=S1S2)
    if (nrow(S1S2)==1){mXmX0<-matrix(mXmX0,nrow=1,ncol=n)}
    M<-apply(invT,1,function(tj,mXmX0,n){
	TM<-matrix(tj,ncol=n,nrow=nrow(mXmX0),byrow=T)*mXmX0
	TM<-apply(TM,2,function(x){if (length(which(is.na(x)))>0){rep(0,length(x))}else{x}})
	if (nrow(mXmX0)==1){TM<-matrix(TM,nrow=1,ncol=ncol(mXmX0))}
	apply(TM,1,sum)
    },mXmX0=mXmX0,n=n)
    M[which(abs(M)<1e-15)]<-0
    M
}

.Kappa.matrix.NA<-function(vNAX,mAncestorTimes,S22,mS,n,kX,mXmX0){
## remove after working and copied to FilesForGLSforB/
    if (length(vNAX)>0){
	invCovXX<-solve((mAncestorTimes%x%S22)[-vNAX,-vNAX,drop=FALSE])
	mXmX0[-vNAX]<-invCovXX%*%(c(mXmX0)[-vNAX]) 
    }else{
        invCovXX<-solve(mAncestorTimes)%x%solve(S22)
	mXmX0<-matrix(invCovXX%*%(c(mXmX0)),ncol=n,nrow=kX)
    }
    res<-sapply(1:n,function(j,mS,mXmX0,n,kX){
        vX<-mXmX0[,j]
        M<-rep(0,nrow(mS))
    	vNAXr<-which(is.na(vX))
    	if (length(vNAXr)>0){
    	    if (length(vNAXr)<kX){
    		vX<-vX[-vNAXr];
    		mS<-mS[,-vNAXr,drop=FALSE]
    		#if(nrow(mS)>1){mS<-mS[,-vNAXr]}else{mS<-matrix(mS[,-vNAXr],nrow=1)};
    		M<-mS%*%vX
    	    }
    	}
    	else{M<-mS%*%vX}
    	M[which(abs(M)<1e-15)]<-0
	M
    },mS=mS,mXmX0=mXmX0,n=n,kX=kX,simplify=TRUE)## glues by column
    res[which(abs(res)<1e-15)]<-0
    if (nrow(mS)==1){res<-matrix(res,nrow=1,ncol=n)}
    res
}

.CalcVlqStat<-function(lq,vlambda,k){
    l<-lq%/%k+1
    q<-lq%%k+1
    sumllq<-vlambda[l]+vlambda[q]
    1/sumllq
}
                
.CalcVlq<-function(lq,vlambda,t,k){
    l<-lq%/%k+1
    q<-lq%%k+1
    sumllq<-vlambda[l]+vlambda[q]
    if (sumllq==0){t}else{(1-exp(-1*sumllq*t))/sumllq}
}
                                
.CalcVlq2<-function(lq,vlambda,t,k){
    ## t > 0 assumed
    l<-lq%/%k+1
    q<-lq%%k+1
    sumllq<-vlambda[l]+vlambda[q]
    if (sumllq==0){t}else{(exp(sumllq*t)-1)/sumllq}
}
                                                    