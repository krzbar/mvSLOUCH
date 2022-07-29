## This file is part of mvSLOUCH

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this package or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .

.calcCI<-function(params,designToEstim,conf.level=0.95,regressCovar=NULL,bfullCI=FALSE){
## called in: .ParamsSummary(modelparamssummary.R)
## mData is a dummy parameter at the moment

    lCI<-NULL
    if (bfullCI){
        .my_message("This confidence interval option has been removed due to numeric instabilities and long running times. Please use the parametric bootstrap option.\n",TRUE)
    }
## -------------------------------------------------------------------------------
    if ((!is.null(regressCovar))&&((nrow(regressCovar)==0)||(ncol(regressCovar)==0)||is.na(regressCovar[1]))){regressCovar<-NULL}
    if (!is.null(regressCovar)){
        v_names_regressCovar<-rep(NA,nrow(regressCovar))
        if (is.null(lCI)){lCI<-vector("list",1);names(lCI)<-"regression.summary"}
	lCI$regression.summary<-list()
	
	vRegCIs<-qnorm(1-(1-conf.level)/2)*sqrt(diag(regressCovar))
	currVar<-length(vRegCIs) ## we go from the end as at the beginning we could have some fixed effects
	if (is.element("BX0",names(designToEstim)) && designToEstim$BX0){
	## NOT CORRECTED HERE IN ANY WAY FOR signsB or signsBX0 or signsvX0
	    trueBX0<-params$B%*%params$vX0
	    lCI$regression.summary$BX0.regression.confidence.interval<-matrix(NA,ncol=3,nrow=length(trueBX0))
	    colnames(lCI$regression.summary$BX0.regression.confidence.interval)<-c("Lower.end","Estimated.Point","Upper.end")
	    lCI$regression.summary$BX0.regression.confidence.interval[,"Estimated.Point"]<-trueBX0
	    lCI$regression.summary$BX0.regression.confidence.interval[,"Lower.end"]<-trueBX0-vRegCIs[(currVar-length(trueBX0)+1):currVar]
	    lCI$regression.summary$BX0.regression.confidence.interval[,"Upper.end"]<-trueBX0+vRegCIs[(currVar-length(trueBX0)+1):currVar]
	    currVar<-currVar-length(trueBX0)
	    v_names_regressCovar[(currVar-length(trueBX0)+1):currVar]<-paste0("BX0_",1:length(trueBX0))
	}
	if (is.element("X0",names(designToEstim)) && designToEstim$X0){
	    lenGLSvX0<-length(params$vX0)
	    if ((is.element("signsvX0",names(designToEstim)))&&(length(which(!is.na(designToEstim$signsvX0))>0))){
		lenGLSvX0<-lenGLSvX0-length(which(!is.na(designToEstim$signsvX0)))
	    }
	    if (lenGLSvX0>0){
		lCI$regression.summary$X0.regression.confidence.interval<-matrix(NA,ncol=3,nrow=lenGLSvX0)
		colnames(lCI$regression.summary$X0.regression.confidence.interval)<-c("Lower.end","Estimated.Point","Upper.end")
		lCI$regression.summary$X0.regression.confidence.interval[,"Estimated.Point"]<-params$vX0
		lCI$regression.summary$X0.regression.confidence.interval[,"Lower.end"]<-params$vX0
		lCI$regression.summary$X0.regression.confidence.interval[,"Upper.end"]<-params$vX0		
		if ((is.element("signsvX0",names(designToEstim)))&&(length(which(is.na(designToEstim$signsvX0))>0))){
		    lCI$regression.summary$X0.regression.confidence.interval[which(is.na(designToEstim$signsvX0)),"Lower.end"]<-lCI$regression.summary$X0.regression.confidence.interval[which(is.na(designToEstim$signsvX0)),"Lower.end"]-vRegCIs[(currVar-lenGLSvX0+1):currVar]
		    lCI$regression.summary$X0.regression.confidence.interval[which(is.na(designToEstim$signsvX0)),"Upper.end"]<-lCI$regression.summary$X0.regression.confidence.interval[which(is.na(designToEstim$signsvX0)),"Upper.end"]+vRegCIs[(currVar-lenGLSvX0+1):currVar]
		}else{
		    lCI$regression.summary$X0.regression.confidence.interval[,"Lower.end"]<-lCI$regression.summary$X0.regression.confidence.interval[,"Lower.end"]-vRegCIs[(currVar-lenGLSvX0+1):currVar]
		    lCI$regression.summary$X0.regression.confidence.interval[,"Upper.end"]<-lCI$regression.summary$X0.regression.confidence.interval[,"Upper.end"]+vRegCIs[(currVar-lenGLSvX0+1):currVar]
		}
		v_names_regressCovar[(currVar-lenGLSvX0+1):currVar]<-paste0("X0_",1:lenGLSvX0)
		currVar<-currVar-lenGLSvX0
	    }
	}
	if (is.element("B",names(designToEstim)) && designToEstim$B){
	    lenGLSB<-length(c(params$B))
	    if ((is.element("signsB",names(designToEstim)))&&(length(which(!is.na(designToEstim$signsB))>0))){
		lenGLSB<-lenGLSB-length(which(!is.na(designToEstim$signsB)))
	    }	    	    
	    if (lenGLSB>0){
		if (ncol(params$B)==1){
	    	    lCI$regression.summary$B.regression.confidence.interval<-matrix(NA,ncol=3,nrow=nrow(params$B))
		    colnames(lCI$regression.summary$B.regression.confidence.interval)<-c("Lower.end","Estimated.Point","Upper.end")
		    lCI$regression.summary$B.regression.confidence.interval[,"Estimated.Point"]<-params$B[,1]
		    lCI$regression.summary$B.regression.confidence.interval[,"Lower.end"]<-params$B[,1]
		    lCI$regression.summary$B.regression.confidence.interval[,"Upper.end"]<-params$B[,1]
		    if ((is.element("signsB",names(designToEstim)))&&(length(which(is.na(designToEstim$signsB))>0))){
			lCI$regression.summary$B.regression.confidence.interval[which(is.na(designToEstim$signsB)),"Lower.end"]<-lCI$regression.summary$B.regression.confidence.interval[which(is.na(designToEstim$signsB)),"Lower.end"]-vRegCIs[(currVar-lenGLSB+1):currVar]
			lCI$regression.summary$B.regression.confidence.interval[which(is.na(designToEstim$signsB)),"Upper.end"]<-lCI$regression.summary$B.regression.confidence.interval[which(is.na(designToEstim$signsB)),"Upper.end"]+vRegCIs[(currVar-lenGLSB+1):currVar]
		    }else{
			lCI$regression.summary$B.regression.confidence.interval[,"Lower.end"]<-lCI$regression.summary$B.regression.confidence.interval[,"Lower.end"]-vRegCIs[(currVar-lenGLSB+1):currVar]
			lCI$regression.summary$B.regression.confidence.interval[,"Upper.end"]<-lCI$regression.summary$B.regression.confidence.interval[,"Upper.end"]+vRegCIs[(currVar-lenGLSB+1):currVar]
		    }
		}else{
	    	    lCI$regression.summary$B.regression.confidence.interval<-vector("list",3)
		    names(lCI$regression.summary$B.regression.confidence.interval)<-c("Lower.end","Estimated.Point","Upper.end")
	    	    lCI$regression.summary$B.regression.confidence.interval$Estimated.Point<-params$B
	    	    lCI$regression.summary$B.regression.confidence.interval$Lower.end<-params$B	    	    
		    lCI$regression.summary$B.regression.confidence.interval$Upper.end<-params$B
		    if ((is.element("signsB",names(designToEstim)))&&(length(which(is.na(designToEstim$signsB))>0))){
			## CHECK IF THIS IS CORRECT WHEN B IS BY GLS
			mBCI<-matrix(NA,nrow=nrow(params$B),ncol=ncol(params$B))
			mtmpBCI<-t(mBCI)
			mtmpBCI[which(is.na(designToEstim$signsB))]<-vRegCIs[(currVar-lenGLSB+1):currVar]
			mBCI<-t(mtmpBCI)
			lCI$regression.summary$B.regression.confidence.interval$Lower.end<-lCI$regression.summary$B.regression.confidence.interval$Lower.end-mBCI
		    	lCI$regression.summary$B.regression.confidence.interval$Upper.end<-lCI$regression.summary$B.regression.confidence.interval$Upper.end+mBCI
		    }else{
			lCI$regression.summary$B.regression.confidence.interval$Lower.end<-lCI$regression.summary$B.regression.confidence.interval$Lower.end-matrix(vRegCIs[(currVar-lenGLSB+1):currVar],nrow=nrow(params$B),ncol=ncol(params$B),byrow=TRUE)
			lCI$regression.summary$B.regression.confidence.interval$Upper.end<-lCI$regression.summary$B.regression.confidence.interval$Upper.end+matrix(vRegCIs[(currVar-lenGLSB+1):currVar],nrow=nrow(params$B),ncol=ncol(params$B),byrow=TRUE)
		    }
		}
	    	v_names_regressCovar[(currVar-lenGLSB+1):currVar]<-paste0("B_",1:lenGLSB)
	    	currVar<-currVar-lenGLSB	    	
	    }
	}
	if (is.element("psi",names(designToEstim)) && designToEstim$psi){
	    lenGLSmPsi<-length(c(params$mPsi))
	    if ((is.element("signsmPsi",names(designToEstim)))&&(length(which(!is.na(designToEstim$signsmPsi))>0))){
		lenGLSmPsi<-lenGLSmPsi-length(which(!is.na(designToEstim$signsmPsi)))
	    }	    	    
	    if (lenGLSmPsi>0){
		if (ncol(params$mPsi)==1){
	    	    lCI$regression.summary$mPsi.regression.confidence.interval<-matrix(NA,ncol=3,nrow=nrow(params$mPsi))
		    colnames(lCI$regression.summary$mPsi.regression.confidence.interval)<-c("Lower.end","Estimated.Point","Upper.end")
		    lCI$regression.summary$mPsi.regression.confidence.interval[,"Estimated.Point"]<-params$mPsi[,1]
		    lCI$regression.summary$mPsi.regression.confidence.interval[,"Lower.end"]<-params$mPsi[,1]
		    lCI$regression.summary$mPsi.regression.confidence.interval[,"Upper.end"]<-params$mPsi[,1]
		    if ((is.element("signsmPsi",names(designToEstim)))&&(length(which(is.na(designToEstim$signsmPsi))>0))){
		    	lCI$regression.summary$mPsi.regression.confidence.interval[which(is.na(designToEstim$signsmPsi)),"Lower.end"]<-lCI$regression.summary$mPsi.regression.confidence.interval[which(is.na(designToEstim$signsmPsi)),"Lower.end"]-vRegCIs[(currVar-lenGLSmPsi+1):currVar]
			lCI$regression.summary$mPsi.regression.confidence.interval[which(is.na(designToEstim$signsmPsi)),"Upper.end"]<-lCI$regression.summary$mPsi.regression.confidence.interval[which(is.na(designToEstim$signsmPsi)),"Upper.end"]+vRegCIs[(currVar-lenGLSmPsi+1):currVar]
		    }else{
			lCI$regression.summary$mPsi.regression.confidence.interval[,"Lower.end"]<-lCI$regression.summary$mPsi.regression.confidence.interval[,"Lower.end"]-vRegCIs[(currVar-lenGLSmPsi+1):currVar]
			lCI$regression.summary$mPsi.regression.confidence.interval[,"Upper.end"]<-lCI$regression.summary$mPsi.regression.confidence.interval[,"Upper.end"]+vRegCIs[(currVar-lenGLSmPsi+1):currVar]
		    }
		}else{
	    	    lCI$regression.summary$mPsi.regression.confidence.interval<-vector("list",3)
		    names(lCI$regression.summary$mPsi.regression.confidence.interval)<-c("Lower.end","Estimated.Point","Upper.end")
		    lCI$regression.summary$mPsi.regression.confidence.interval$Estimated.Point<-params$mPsi
		    lCI$regression.summary$mPsi.regression.confidence.interval$Lower.end<-params$mPsi
		    lCI$regression.summary$mPsi.regression.confidence.interval$Upper.end<-params$mPsi
		    if ((is.element("signsmPsi",names(designToEstim)))&&(length(which(is.na(designToEstim$signsmPsi))>0))){	
			mPsiCI<-matrix(NA,ncol=ncol(params$mPsi),nrow=nrow(params$mPsi))
			mPsiCI[which(is.na(designToEstim$signsmPsi))]<-vRegCIs[(currVar-lenGLSmPsi+1):currVar]
			lCI$regression.summary$mPsi.regression.confidence.interval$Lower.end<-lCI$regression.summary$mPsi.regression.confidence.interval$Lower.end-mPsiCI
			lCI$regression.summary$mPsi.regression.confidence.interval$Upper.end<-lCI$regression.summary$mPsi.regression.confidence.interval$Upper.end+mPsiCI
		    }else{
			lCI$regression.summary$mPsi.regression.confidence.interval$Lower.end<-lCI$regression.summary$mPsi.regression.confidence.interval$Lower.end-matrix(vRegCIs[(currVar-lenGLSmPsi+1):currVar],nrow=nrow(params$mPsi),ncol=ncol(params$mPsi),byrow=FALSE)
			lCI$regression.summary$mPsi.regression.confidence.interval$Upper.end<-lCI$regression.summary$mPsi.regression.confidence.interval$Upper.end+matrix(vRegCIs[(currVar-lenGLSmPsi+1):currVar],nrow=nrow(params$mPsi),ncol=ncol(params$mPsi),byrow=FALSE)
		    }
		
		}
		v_names_regressCovar[(currVar-lenGLSmPsi+1):currVar]<-paste0("mPsi_",1:lenGLSmPsi)
		currVar<-currVar-lenGLSmPsi
	    }
	}
	if (is.element("psi0",names(designToEstim)) && designToEstim$psi0){
	    lenGLSmPsi0<-length(params$mPsi0)
	    if ((is.element("signsmPsi0",names(designToEstim)))&&(length(which(!is.na(designToEstim$signsmPsi0))>0))){
		lenGLSmPsi0<-lenGLSmPsi0-length(which(!is.na(designToEstim$signsmPsi0)))
	    }
	    if (lenGLSmPsi0>0){
		lCI$regression.summary$mPsi0.regression.confidence.interval<-matrix(NA,ncol=3,nrow=nrow(params$mPsi0))
		colnames(lCI$regression.summary$mPsi0.regression.confidence.interval)<-c("Lower.end","Estimated.Point","Upper.end")
		lCI$regression.summary$mPsi0.regression.confidence.interval[,"Estimated.Point"]<-params$mPsi0
		lCI$regression.summary$mPsi0.regression.confidence.interval[,"Lower.end"]<-params$mPsi0
		lCI$regression.summary$mPsi0.regression.confidence.interval[,"Upper.end"]<-params$mPsi0
		if ((is.element("signsmPsi0",names(designToEstim)))&&(length(which(is.na(designToEstim$signsmPsi0))>0))){	
		    lCI$regression.summary$mPsi0.regression.confidence.interval[which(is.na(designToEstim$signsmPsi0)),"Lower.end"]<-lCI$regression.summary$mPsi0.regression.confidence.interval[which(is.na(designToEstim$signsmPsi0)),"Lower.end"]-vRegCIs[(currVar-lenGLSmPsi0+1):currVar]
		    lCI$regression.summary$mPsi0.regression.confidence.interval[which(is.na(designToEstim$signsmPsi0)),"Upper.end"]<-lCI$regression.summary$mPsi0.regression.confidence.interval[which(is.na(designToEstim$signsmPsi0)),"Upper.end"]+vRegCIs[(currVar-lenGLSmPsi0+1):currVar]
		}else{
		    lCI$regression.summary$mPsi0.regression.confidence.interval[,"Lower.end"]<-lCI$regression.summary$mPsi0.regression.confidence.interval[,"Lower.end"]-vRegCIs[(currVar-lenGLSmPsi0+1):currVar]
		    lCI$regression.summary$mPsi0.regression.confidence.interval[,"Upper.end"]<-lCI$regression.summary$mPsi0.regression.confidence.interval[,"Upper.end"]+vRegCIs[(currVar-lenGLSmPsi0+1):currVar]
		}
		v_names_regressCovar[(currVar-lenGLSmPsi0+1):currVar]<-paste0("mPsi0_",1:lenGLSmPsi0)
		currVar<-currVar-lenGLSmPsi0		
	    }
	}
	if (is.element("y0",names(designToEstim)) && designToEstim$y0 && !designToEstim$y0AncState ){
	    lenGLSvY0<-length(params$vY0)
	    if ((is.element("signsmvY0",names(designToEstim)))&&(length(which(!is.na(designToEstim$signsvY0))>0))){
		lenGLSmPsi0<-lenGLSvY0-length(which(!is.na(designToEstim$signsvY0)))
	    }
	    if (lenGLSvY0>0){
		lCI$regression.summary$Y0.regression.confidence.interval<-matrix(NA,ncol=3,nrow=length(params$vY0))
		colnames(lCI$regression.summary$Y0.regression.confidence.interval)<-c("Lower.end","Estimated.Point","Upper.end")
		lCI$regression.summary$Y0.regression.confidence.interval[,"Estimated.Point"]<-params$vY0
		lCI$regression.summary$Y0.regression.confidence.interval[,"Lower.end"]<-params$vY0
	        lCI$regression.summary$Y0.regression.confidence.interval[,"Upper.end"]<-params$vY0
		if ((is.element("signsmPsi0",names(designToEstim)))&&(length(which(is.na(designToEstim$signsmPsi0))>0))){	
		    lCI$regression.summary$Y0.regression.confidence.interval[which(is.na(designToEstim$signsmPsi0)),"Lower.end"]<-lCI$regression.summary$Y0.regression.confidence.interval[which(is.na(designToEstim$signsmPsi0)),"Lower.end"]-vRegCIs[(currVar-lenGLSvY0+1):currVar]
		    lCI$regression.summary$Y0.regression.confidence.interval[which(is.na(designToEstim$signsmPsi0)),"Upper.end"]<-lCI$regression.summary$Y0.regression.confidence.interval[which(is.na(designToEstim$signsmPsi0)),"Upper.end"]+vRegCIs[(currVar-lenGLSvY0+1):currVar]
		}else{
		    lCI$regression.summary$Y0.regression.confidence.interval[,"Lower.end"]<-lCI$regression.summary$Y0.regression.confidence.interval[,"Lower.end"]-vRegCIs[(currVar-lenGLSvY0+1):currVar]
		    lCI$regression.summary$Y0.regression.confidence.interval[,"Upper.end"]<-lCI$regression.summary$Y0.regression.confidence.interval[,"Upper.end"]+vRegCIs[(currVar-lenGLSvY0+1):currVar]
		}
	    	v_names_regressCovar[(currVar-lenGLSvY0+1):currVar]<-paste0("Y0_",1:lenGLSvY0)
		currVar<-currVar-lenGLSvY0	    	
	    }
	}
	colnames(regressCovar)<-v_names_regressCovar
	rownames(regressCovar)<-v_names_regressCovar
	lCI$regression.summary$regression.covariance.matrix<-regressCovar
	lCI$regression.summary$regression.confidence.interval.comment<-"These are confidence intervals for parameters estimated by a GLS conditional on the A and diffusion matrix parameters. In the full covariance matrix of the regression estimators the matrix parameters are entered column wise for the deterministic optimum and row wise for the B matrix. Be careful if some of the parameters were set at fixed values in the user's input, i.e. check if the variances correctly correspond to the presented 95% CIs. These are calculated  as estimate =/- (qnorm(0.975), i.e. 0.975 quantile of N(0,1))*(square root of appropriate diagonal element of the covariance matrix), i.e. standard deviation."
    }
    lCI
}

