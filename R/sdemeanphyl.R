.calc.phyl.mean<-function(vSpecDist,EvolModel,modelParams){
## called in OUphylregression.R
    vMean=switch(EvolModel,
	bm=.bm.phyl.mean(modelParams$vX0,length(vSpecDist)), 
	ouch=.ouch.phyl.mean(vSpecDist,modelParams), 
	slouch=.mvslouch.phyl.mean(vSpecDist,modelParams),
	mvslouch=.mvslouch.phyl.mean(vSpecDist,modelParams)
    )
    vMean[which(abs(vMean)<1e-15)]<-0
    vMean
}

.bm.phyl.mean<-function(vY0,n){rep(vY0,n)} 

.ouch.phyl.mean<-function(vSpecDist,modelParams){
    if (is.null(modelParams$precalcMatrices[[3]])){
	lexpmtA<-sapply(vSpecDist,function(t){.calc.exptA(-t,modelParams$precalcMatrices[[1]])},simplify=FALSE) ## vSpecDist  will be the times of current species
        if (!(is.null(modelParams$regimeTimes))){lexptjA<-
    	    sapply(1:length(modelParams$regimeTimes),function(i,regimeTimes,vSpecDist){tjs<-regimeTimes[[i]];specT<-vSpecDist[i];sapply(tjs,function(t,specT){.calc.exptA(t-specT,modelParams$precalcMatrices[[1]])},specT=specT,simplify=FALSE)},regimeTimes=modelParams$regimeTimes,vSpecDist=vSpecDist,simplify=FALSE)}
        else{ lexptjA<-sapply(1:length(vSpecDist),function(i,k){list(lexpmtA[[i]],diag(1,nrow=k,ncol=k))},k=length(modelParams$vY0),simplify=FALSE) }
    }else{
	lexpmtA<-modelParams$precalcMatrices[[3]]$lexpmtA
	lexptjA<-modelParams$precalcMatrices[[3]]$lexptjA
    }
    if ((ncol(modelParams$mPsi)>1)&&(!is.null(colnames(modelParams$mPsi)))){modelParams$mPsi<-modelParams$mPsi[,order(colnames(modelParams$mPsi)),drop=FALSE]}
    c(sapply(1:length(vSpecDist),function(s){
	.calc.mean.ouch.mv(lexpmtA[[s]],modelParams$vY0,modelParams$mPsi,modelParams$mPsi0,
	{if(is.null(modelParams$regimeTimes)){NULL}else{lexptjA[[s]]}},{if(is.null(modelParams$regimeTimes)){NULL}else{modelParams$regimes[[s]]}})
    })) 
}



.mvslouch.phyl.mean<-function(vSpecDist,modelParams){
    if (is.null(modelParams$precalcMatrices[[3]])){
	lexpmtA<-sapply(vSpecDist,function(t){.calc.exptA(-t,modelParams$precalcMatrices[[1]])},simplify=FALSE) ##  vSpecDist will be the times of current species
        if (!(is.null(modelParams$regimeTimes))){lexptjA<-
    		    sapply(1:length(modelParams$regimeTimes),function(i,regimeTimes,vSpecDist){tjs<-regimeTimes[[i]];specT<-vSpecDist[i];sapply(tjs,function(t,specT){.calc.exptA(t-specT,modelParams$precalcMatrices[[1]])},specT=specT,simplify=FALSE)},regimeTimes=modelParams$regimeTimes,vSpecDist=vSpecDist,simplify=FALSE)}
        else{ lexptjA<-sapply(1:length(vSpecDist),function(i,k){list(lexpmtA[[i]],diag(1,nrow=k,ncol=k))},k=length(modelParams$vY0),simplify=FALSE)}        
    }else{
    	lexpmtA<-modelParams$precalcMatrices[[3]]$lexpmtA
	lexptjA<-modelParams$precalcMatrices[[3]]$lexptjA
    }    
    if ((ncol(modelParams$mPsi)>1)&&(!is.null(colnames(modelParams$mPsi)))){modelParams$mPsi<-modelParams$mPsi[,order(colnames(modelParams$mPsi)),drop=FALSE]}
    c(sapply(1:length(vSpecDist),function(s){
	.calc.mean.slouch.mv(lexpmtA[[s]],modelParams$precalcMatrices[[1]]$A1B,modelParams$vY0,modelParams$vX0,modelParams$mPsi,modelParams$mPsi0,
	{if(is.null(modelParams$regimeTimes)){NULL}else{lexptjA[[s]]}},{if(is.null(modelParams$regimeTimes)){NULL}else{modelParams$regimes[[s]]}})
    })) 
}


