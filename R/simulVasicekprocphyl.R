## This file is part of mvSLOUCH

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this package or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .


.simulVasicekProcPhylTree<-function(phyltree,EvolModel,modelParams,EstimationParams=NULL,regimes=NULL,regimes.times=NULL,Simulparams=NULL,dropInternal=TRUE,bAllTrajectories=FALSE,M.error=NULL,bKeep_tree=FALSE){
## This is a dummy object that might become useful in the future when B for the mvSLOUCH model is estimated by GLS
    lPrecalculates<-NULL 

    bJumpAtNode<-FALSE
    if (is.element("jump",names(Simulparams))){bJumpAtNode<-TRUE}
## -----------------------------------------------------------------------    

    if ((EvolModel=="bm")||(EvolModel=="bmStep")){
	regimes<-NULL
	regimes.times<-NULL
    }
    if ((EvolModel=="mvslouch") && .is_det0(modelParams$A,NULL)){
	EvolModel<-"ouch"
	modelParams<-.mvslouch_to_ouch_model(modelParams)
	.my_warning("WARNING: A is a singular matrix, simulating under the OUOU model and not OUBM model so that distribution of output agrees with provided model parameters.",TRUE,FALSE)
    }

    regimesList<-.InitialRegimeSetup(phyltree=phyltree,regimes=regimes,regimes.times=regimes.times,mData=NULL,kX=NULL,kYX=NULL,root.regime=NULL,M.error=NULL)
    bOKregimes<-regimesList$bOKregimes
    if (!bOKregimes){.my_stop("The regimes, regimes.times and phyltree variables are not consistant. Cannot perform estimation. Please see manual on their format.",TRUE)}
    regimes<-regimesList$regimes
    regimes.times<-regimesList$regimes.times
    regimes.types<-regimesList$regimes.types
    root.regime<-regimesList$root.regime
    regimes.types.orig<-regimesList$regimes.types.orig
    phyltree<-regimesList$phyltree

    if (bAllTrajectories && (is.null(Simulparams) || !is.element("step",names(Simulparams)))){
	step<-min(c(0.001,phyltree$tree_height/2000))
	if (is.null(Simulparams)){Simulparams<-list(step=step)}
	else{Simulparams$step<-step}
    }
    modelParams$regimes<-regimes
    modelParams$regimeTimes<-regimes.times

    inumber_species<-phyltree$Ntips
    inumber_nodes<-phyltree$Ntips+phyltree$Nnode
    inumber_internal_nodes<-phyltree$Nnode
    inumber_edges<-nrow(phyltree$edge)
    
     params<-list()
     params$EvolModel<-EvolModel
     if (params$EvolModel=="slouch"){params$EvolModel<-"mvslouch";.my_warning("Call to slouch not yet implemented, using mvslouch model. You might need to run again with correct input structure. \n",TRUE,FALSE)}  
     if (!is.element("method",names(params))){params$method<-"glsgc"}
     if (params$EvolModel=="bm"){params$method<-"maxlik"}
     if (!is.element("bShouldPrint",names(params))){params$bShouldPrint<-TRUE}
     if (!is.element("tol",names(params))){params$tol<-.set.tol(params$method)}
     if (!is.element("maxIter",names(params))){params$maxIter<-.set.maxiter(params$method)}
     if (!is.element("maxTries",names(params))){params$maxTries<-10}
     if (!is.element("minLogLik",names(params))){params$minLogLik<- -Inf}                                                           

    
    lPrecalculates<-list(vSpecies_times=NA)
    lPrecalculates$vSpecies_times<-phyltree$time.of.nodes[phyltree$tip_species_index] 
    
    if ((EvolModel=="bm")||(EvolModel=="bmStep")){
	x0<-modelParams$vX0;kX<-nrow(modelParams$Sxx)
	EstimationParams<-.set.estimparams(params,0,kX,length(regimes.types))
    }
    if ((EvolModel=="ouch")||(EvolModel=="ouchStep")){
	kY<-nrow(modelParams$A)
	x0<-modelParams$vY0
	EstimationParams<-.set.estimparams(params,kY,0,length(regimes.types))
        modelParams$precalcMatrices<-.decompEigenA.S(modelParams,lPrecalculates,NA,list(bCalcA=TRUE,bCovCalc=TRUE,dzetacalc=FALSE,lexptcalc=TRUE,kappacalc=FALSE,interceptcalc=FALSE),NULL)
    }
    if ((EvolModel=="mvslouch")||(EvolModel=="mvslouchStep")){
	kY<-nrow(modelParams$A)
	kX<-ncol(modelParams$B)
    	x0<-c(modelParams$vY0,modelParams$vX0)
    	EstimationParams<-.set.estimparams(params,kY,kX,length(regimes.types))
        modelParams$precalcMatrices<-.decompEigenA.S(modelParams,lPrecalculates,EstimationParams$designToEstim,list(bCalcA=TRUE,bCovCalc=TRUE,dzetacalc=FALSE,lexptcalc=TRUE,kappacalc=FALSE,interceptcalc=FALSE),NULL)
    }
    
    RootNode<-phyltree$root_index

    mTreeTraject<-matrix(NA,ncol=length(x0),nrow=inumber_nodes)
    vnodelabels<-paste("node_",1:inumber_nodes,sep="")
    vnodelabels[phyltree$tip_species_index]<-phyltree$tip.label
    rownames(mTreeTraject)<-vnodelabels
    colnames(mTreeTraject)<-names(x0)
    
    if (!is.na(x0[1])){ mTreeTraject[RootNode,]<-x0} ## the state of the root is x0    
    if (bJumpAtNode){vJumpNode<-rep(NA,nrow(mTreeTraject))}
    if (bAllTrajectories){
    	if (is.element(EvolModel,c("bm","bmStep","ouch","ouchStep","mvslouch","mvslouchStep"))){
	    ## we do the simulation on each edge
	    lFullTraject<-vector("list",inumber_edges)
	    if (!substr(EvolModel,start=nchar(EvolModel)-3,stop=nchar(EvolModel))=="Step"){
	    	orgEvolModel<-EvolModel
		EvolModel<-paste(EvolModel,"Step",sep="")
	    }
	}else{
	    .my_warning("Cannot yet simulate full trajectory for the chosen model",TRUE,TRUE)
	    bAllTrajectories<-FALSE
	}
    }else{lFullTraject<-NA}

    for (Term in phyltree$tip_species_index){## for each terminal node -- defines a lineage
	vTermLineage<-phyltree$path.from.root[[Term]]$nodes ## we want it to start from start of tree, but DO NOT need to reverse as due to setup of phyltree
	i<-1
	while(!is.na(mTreeTraject[vTermLineage[i],1])){i<-i+1}
	if(i<length(vTermLineage)+1){## just to check if we haven't run over
	    for (j in i:length(vTermLineage)){
    		branch_index<-intersect(which(phyltree$edge[,1]==vTermLineage[j-1]),which(phyltree$edge[,2]==vTermLineage[j]))
    		currbranch<-branch_index
    		if (bAllTrajectories){
    		    lFullTraject[[currbranch]]<-vector("list",4)
    		    names(lFullTraject[[currbranch]])<-c("branch","nodenames","trajectory","branch_index")
    		    lFullTraject[[currbranch]]$branch<-c(vTermLineage[j-1],vTermLineage[j])
    		    lFullTraject[[currbranch]]$nodenames<-vnodelabels[c(vTermLineage[j-1],vTermLineage[j])]
    		    lFullTraject[[currbranch]]$branch_index<-branch_index
    		}
    		Xprev<-mTreeTraject[vTermLineage[j-1],]
		timeDiff<-phyltree$edge.length[branch_index]

		if (bJumpAtNode){		    
		    bDoJump<-FALSE
		    if (Simulparams$jump$jumptype=="ForBoth"){
                        if (is.na(vJumpNode[vTermLineage[j-1]])){
                            if (runif(1)<Simulparams$jump$jumpprob){bDoJump<-TRUE}
                            else{bDoJump<-FALSE;vJumpNode[vTermLineage[j-1]]<-0}
                        }
                    }
		    if (Simulparams$jump$jumptype=="RandomLineage"){
			if (is.na(vJumpNode[vTermLineage[j-1]])){vJumpNode[vTermLineage[j-1]]<-sample(c(1,2),1)}
			if ((vJumpNode[vTermLineage[j-1]]==1)||(vJumpNode[vTermLineage[j-1]]==3)){bDoJump<-TRUE}
			if (vJumpNode[vTermLineage[j-1]]==2){bDoJump<-FALSE;vJumpNode[vTermLineage[j-1]]<-3}
		    }
		    if (Simulparams$jump$jumptype=="JumpWithProb"){
		    	if (is.na(vJumpNode[vTermLineage[j-1]])||(vJumpNode[vTermLineage[j-1]]==1)){if (runif(1)<Simulparams$jump$jumpprob){bDoJump<-TRUE}}
			if ((!is.na(vJumpNode[vTermLineage[j-1]]))&&(vJumpNode[vTermLineage[j-1]]==2)){bDoJump<-FALSE;vJumpNode[vTermLineage[j-1]]<-0}
		    }
		    if (bDoJump){		    	
		    	## add according to distribution
		    	if (Simulparams$jump$jumpdistrib=="Normal"){Xprev<-Xprev+mvtnorm::rmvnorm(1,mean=Simulparams$jump$vMean,sigma=Simulparams$jump$mCov)}
			
			## any post-add trajectory corrections
			if (Simulparams$jump$jumptype=="ForBoth"){mTreeTraject[vTermLineage[j-1],]<-Xprev;vJumpNode[vTermLineage[j-1]]<-0}
			if (Simulparams$jump$jumptype=="RandomLineage"){vJumpNode[vTermLineage[j-1]]<-0}
			if (Simulparams$jump$jumptype=="JumpWithProb"){
			    if (is.na(vJumpNode[vTermLineage[j-1]])){vJumpNode[vTermLineage[j-1]]<-1}
			    else{vJumpNode[vTermLineage[j-1]]<-vJumpNode[vTermLineage[j-1]]+1}
			}
		    }
		}
## -----------------------------------------------------------------------    
		if (EvolModel=="slouch"){}## Empty at the moment
		if (EvolModel=="bm"){
		    vMean<-Xprev[1:kX]
		    mCov<-timeDiff*modelParams$Sxx%*%t(modelParams$Sxx)
		    Xdrawn<-mvtnorm::rmvnorm(n=1,mean=vMean,sigma=mCov) 
		    Xdrawn<-Xdrawn[nrow(Xdrawn),]
		}
		if (EvolModel=="bmStep"){
			## in principle (by phylo object standard) this should be Term
			itermNum<-which(phyltree$tip_species_index==Term)
			if (itermNum!=Term){.my_warning("WARNING: non-standard (for phylo object) numbering of tip species!",TRUE,FALSE)}
	    		vX0<-Xprev[1:kX]
	    		tmpmodelparams<-modelParams
	    		tmpmodelparams$vX0<-matrix(vX0,ncol=1,nrow=kX)
			Xdrawn<-.bm.simulate(step=ifelse(is.null(Simulparams),timeDiff,Simulparams$step),duration=timeDiff,modelParams=tmpmodelparams,regimes=NULL,regimes.times=NULL,mCov=NULL)
			if (bAllTrajectories){
			    Xdrawn[,1]<-Xdrawn[,1]+phyltree$time.of.nodes[vTermLineage[j-1]]
			    lFullTraject[[currbranch]]$trajectory<-Xdrawn[1:(nrow(Xdrawn)-1),1:ncol(Xdrawn),drop=FALSE]
			}
			Xdrawn<-Xdrawn[nrow(Xdrawn)-1,2:ncol(Xdrawn)]
		}
		if (EvolModel=="ouchStep"){
			## in principle (by phylo object standard) this should be Term
			itermNum<-which(phyltree$tip_species_index==Term)
			if (itermNum!=Term){.my_warning("WARNING: non-standard (for phylo object) numbering of tip species!",TRUE,FALSE)}
			vWhichTimes<-intersect(which(regimes.times[[itermNum]]>=phyltree$time.of.nodes[vTermLineage[j-1]]),which(regimes.times[[itermNum]]<=phyltree$time.of.nodes[vTermLineage[j]]))
			regimesCurrTimes<-regimes.times[[itermNum]][vWhichTimes[1:(length(vWhichTimes))]]		
			regimesCurr<-regimes[[itermNum]][vWhichTimes[1:(length(vWhichTimes)-1)]]		
	    		vY0<-Xprev[1:kY]
	    		tmpmodelparams<-modelParams
	    		tmpmodelparams$vY0<-matrix(vY0,ncol=1,nrow=kY)
	    		regimesCurrTimes<-regimesCurrTimes-phyltree$time.of.nodes[vTermLineage[j-1]]
			Xdrawn<-.ouch.simulate(step=ifelse(is.null(Simulparams),timeDiff,Simulparams$step),duration=timeDiff,modelParams=tmpmodelparams,regimes=regimesCurr,regimes.times=regimesCurrTimes,mCov=NULL)			
			if (bAllTrajectories){
			    Xdrawn[,1]<-Xdrawn[,1]+phyltree$time.of.nodes[vTermLineage[j-1]]
			    lFullTraject[[currbranch]]$trajectory<-Xdrawn[1:(nrow(Xdrawn)-1),1:ncol(Xdrawn),drop=FALSE]
			}
			Xdrawn<-Xdrawn[nrow(Xdrawn)-1,2:ncol(Xdrawn)]
		}
		if (EvolModel=="mvslouchStep"){
			## in principle (by phylo object standard) this should be Term
			itermNum<-which(phyltree$tip_species_index==Term)
			if (itermNum!=Term){.my_warning("WARNING: non-standard (for phylo object) numbering of tip species!",TRUE,FALSE)}
			vWhichTimes<-intersect(which(regimes.times[[itermNum]]>=phyltree$time.of.nodes[vTermLineage[j-1]]),which(regimes.times[[itermNum]]<=phyltree$time.of.nodes[vTermLineage[j]]))
			regimesCurrTimes<-regimes.times[[itermNum]][vWhichTimes[1:(length(vWhichTimes))]]		
			regimesCurr<-regimes[[itermNum]][vWhichTimes[1:(length(vWhichTimes)-1)]]		
	    		vY0<-Xprev[1:kY]
	    		vX0<-Xprev[(kY+1):(kY+kX)]
	    		tmpmodelparams<-modelParams
	    		tmpmodelparams$vY0<-matrix(vY0,ncol=1,nrow=kY)
	    		tmpmodelparams$vX0<-matrix(vX0,ncol=1,nrow=kX)
	    		regimesCurrTimes<-regimesCurrTimes-phyltree$time.of.nodes[vTermLineage[j-1]]
			Xdrawn<-.mvslouch.simulate(step=ifelse(is.null(Simulparams),timeDiff,Simulparams$step),duration=timeDiff,modelParams=tmpmodelparams,regimes=regimesCurr,regimes.times=regimesCurrTimes,mCov=NULL)
			if (bAllTrajectories){
			    Xdrawn[,1]<-Xdrawn[,1]+phyltree$time.of.nodes[vTermLineage[j-1]]
			    lFullTraject[[currbranch]]$trajectory<-Xdrawn[1:(nrow(Xdrawn)-1),1:ncol(Xdrawn),drop=FALSE]
			}
			Xdrawn<-Xdrawn[nrow(Xdrawn)-1,2:ncol(Xdrawn)]
		}

		if ((EvolModel=="ouch")||(EvolModel=="mvslouch")){
		    ## in principle (by phylo object standard) this should be Term
		    itermNum<-which(phyltree$tip_species_index==Term)
		    if (itermNum!=Term){.my_warning("WARNING: non-standard (for phylo object) numbering of tip species!",TRUE,FALSE)}
		    vWhichTimes<-intersect(which(regimes.times[[itermNum]]>=phyltree$time.of.nodes[vTermLineage[j-1]]),which(regimes.times[[itermNum]]<=phyltree$time.of.nodes[vTermLineage[j]]))
	    	    vY0<-Xprev[1:kY]
	    	    mPsi<-modelParams$mPsi
		    if(ncol(mPsi)>1){mPsi<-mPsi[,order(colnames(mPsi)),drop=FALSE]}
	    	    mPsi0<-modelParams$mPsi0

		    expmtAcurr<-.calc.exptA(t=-timeDiff,modelParams$precalcMatrices[[1]])   

		    exptjA<-sapply(c(regimes.times[[itermNum]][vWhichTimes]-regimes.times[[itermNum]][vWhichTimes[1]]-timeDiff),function(t,precalc){.calc.exptA(t=t,precalc)},precalc=modelParams$precalcMatrices[[1]],simplify=FALSE)		    
	    	    regimesCurr<-regimes[[itermNum]][vWhichTimes[1:(length(vWhichTimes)-1)]]		
	    	    names(regimesCurr)<-NULL
		    if (EvolModel=="ouch"){
		        vMean<-.calc.mean.ouch.mv(expmtAcurr,vY0,mPsi,mPsi0,exptjA,regimesCurr)
		        mCov<-.calc.cov.ouch.mv(timeDiff,modelParams$precalcMatrices[[1]],modelParams$precalcMatrices[[2]])
	    	    }		
		    if (EvolModel=="mvslouch"){## we do not simulate the trajectory -> just draw from appropriate distribution
			vX0<-Xprev[(kY+1):(length(Xprev))]
			A1B<-modelParams$precalcMatrices[[1]]$A1B		
    			vMean<-.calc.mean.slouch.mv(expmtAcurr,A1B,vY0,vX0,mPsi,mPsi0,exptjA,regimesCurr)
			mCov<-.calc.cov.slouch.mv(timeDiff,modelParams$precalcMatrices[[1]],modelParams$precalcMatrices[[2]])		
		    }
		    
		    Xdrawn<-mvtnorm::rmvnorm(n=1,mean=vMean,sigma=mCov) 
		    Xdrawn<-Xdrawn[nrow(Xdrawn),]
		}
		mTreeTraject[vTermLineage[j],]<-Xdrawn
	    }
	}else{.my_stop("Something is wrong with tree---seems to be a network",TRUE)}	
    }    

    
    if (EvolModel=="bm"){if (!is.null(colnames(modelParams$Sxx))){colnames(mTreeTraject)<-colnames(modelParams$Sxx)}}
    if (EvolModel=="ouch"){if (!is.null(colnames(modelParams$A))){colnames(mTreeTraject)<-colnames(modelParams$A)}}
    if (EvolModel=="mvslouch"){if ((!is.null(colnames(modelParams$A)))&&(!is.null(colnames(modelParams$Sxx)))){colnames(mTreeTraject)<-c(colnames(modelParams$A),colnames(modelParams$Sxx))}}

    if (!is.null(M.error)){
	M.error<-.createMeasurementError(M.error,phyltree$Ntips,ncol(mTreeTraject))
    	for (i in phyltree$tip_species_index){
    	    Term<-phyltree$tip_species_index[i]
    	    mErrorCov<-matrix(NA,ncol(mTreeTraject),ncol(mTreeTraject))
	    if (is.list(M.error)){mErrorCov<-M.error[[i]]}else{mErrorCov<-M.error}
    	    m.errors<-mvtnorm::rmvnorm(1,mean=rep(0,ncol(mErrorCov)),sigma=mErrorCov)
	    mTreeTraject[Term,]<-mTreeTraject[Term,]+m.errors
	    if (bAllTrajectories){
		branch_index<-which(phyltree$edge[,2]==Term)
		lFullTraject[[branch_index]]$trajectory[nrow(lFullTraject[[branch_index]]$trajectory),2:ncol(lFullTraject[[branch_index]]$trajectory)]<-lFullTraject[[branch_index]]$trajectory[nrow(lFullTraject[[branch_index]]$trajectory),2:ncol(lFullTraject[[branch_index]]$trajectory)]+m.errors
	    }
	}
    }

    if (dropInternal){mTreeTraject<-mTreeTraject[phyltree$tip_species_index,,drop=FALSE]}

    if (is.null(colnames(mTreeTraject))){colnames(mTreeTraject)<-paste("trait_",1:ncol(mTreeTraject),sep="")}
    if (bAllTrajectories){
	lFullTraject<-sapply(lFullTraject,function(ltraj,vtraitnames){colnames(ltraj$trajectory)<-c("time",vtraitnames);ltraj},vtraitnames=colnames(mTreeTraject),simplify=FALSE)}
    simulReturn<-vector("list",3)
    names(simulReturn)<-c("Tree","ExtantSample","FullTrajectory")
    if (bKeep_tree){simulReturn$Tree<-phyltree}else{simulReturn$Tree<-NULL}
    simulReturn$ExtantSample<-mTreeTraject
    simulReturn$FullTrajectory<-lFullTraject
    simulReturn
}


drawPhylProcess<-function(PhylTraitProcess,vTraitsToPlot=NULL,vColours="black",plotlayout=c(1,1),additionalfigs=FALSE,modelParams=NULL,EvolModel=NULL,xlimits=NULL,ylimits=NULL){
    ## prepare data matrix to plot
    mData<-NA
    if (is.list(PhylTraitProcess)){
	if (is.element("FullTrajectory",names(PhylTraitProcess))){
	    lPhylTraject<-PhylTraitProcess$FullTrajectory
	    mData<-lPhylTraject[[1]]$trajectory
	    if (length(lPhylTraject)>1){for (i in 2:length(lPhylTraject)){mData<-rbind(mData,lPhylTraject[[i]]$trajectory)}}
	}
    }else{if (is.matrix(PhylTraitProcess)){mData<-PhylTraitProcess}}
    
    if (!is.na(mData[1])){    
	ntraits<-ncol(mData)-1
	if (is.null(vTraitsToPlot)){vTraitsToPlot<-2:(ntraits+1)}
	if (!is.numeric(vTraitsToPlot)){.my_stop("The trait column numbers in vTraitsToPlot have to be numeric!",TRUE)}
	if (length(setdiff(vTraitsToPlot,2:(ntraits+1)))>0){.my_stop(paste("Wrong trait column numbers provided in vTraitsToPlot object. The trait column numbers has to be a vector of integers in the range 2 to ",ntraits+1,".",sep=""),TRUE)}
	if (length(vTraitsToPlot)>length(unique(vTraitsToPlot))){.my_warning("There are traits that will be plotted multiple times!",TRUE,FALSE)}
	## check if we have the colors
	if ((length(vColours)==0)||(is.na(vColours[1]))||(is.null(vColours))){vColours<-"black"}
	if (length(vColours)<ntraits){vColours<-rep(vColours,length.out=ntraits)}
    
	## prepare plot area
	if (prod(plotlayout)<ntraits){plotlayout<-c(1,length(vTraitsToPlot));.my_warning("WARNING : Possibly ugly plot layout, change parameter plotlayout! \n",TRUE,FALSE)}
    
	## plot
	par(mfrow=(plotlayout))
	for (i in vTraitsToPlot){
	    minmaxx<-NA
	    if (!is.null(xlimits)){
		if ((is.vector(xlimits))&&(length(xlimits)==2)){minmaxx<-xlimits}
		if ((is.list(xlimits))&&(length(xlimits)==1)&&(is.vector(xlimits[[1]]))&&(length(xlimits[[1]])==2)){minmaxx<-xlimits[[1]]}
		if ((is.list(xlimits))&&(length(xlimits)==ntraits)&&(is.vector(xlimits[[i-1]]))&&(length(xlimits[[i-1]])==2)){minmaxx<-xlimits[[i-1]]}
		if ((is.matrix(xlimits))&&(ncol(xlimits)==2)&&(nrow(xlimits)==1)){minmaxx<-xlimits[1,]}
		if ((is.matrix(xlimits))&&(ncol(xlimits)==2)&&(nrow(xlimits)==ntraits)){minmaxx<-xlimits[i-1,]}		
	    }
	    if (any(is.na(minmaxx))){
		minmaxx<-c(min(mData[,i]),max(mData[,i]))
	    }	    
	    minmaxy<-NA
	    if (!is.null(ylimits)){
		if ((is.vector(ylimits))&&(length(ylimits)==2)){minmaxy<-ylimits}
		if ((is.list(ylimits))&&(length(ylimits)==1)&&(is.vector(ylimits[[1]]))&&(length(ylimits[[1]])==2)){minmaxy<-ylimits[[1]]}
		if ((is.list(ylimits))&&(length(ylimits)==ntraits)&&(is.vector(ylimits[[i-1]]))&&(length(ylimits[[i-1]])==2)){minmaxy<-ylimits[[i-1]]}
		if ((is.matrix(ylimits))&&(ncol(ylimits)==2)&&(nrow(ylimits)==1)){minmaxy<-ylimits[1,]}
		if ((is.matrix(ylimits))&&(ncol(ylimits)==2)&&(nrow(ylimits)==ntraits)){minmaxy<-ylimits[i-1,]}		
	    }
	    if (any(is.na(minmaxy))){
	    	minmaxy<-c(min(mData[,1]),max(mData[,1]))
	    }	    

	    plot(mData[,i],mData[,1],col=vColours[i-1],pch=19,cex=0.2,main="",xlab="",ylab="",frame.plot=FALSE,axes=FALSE,ylim=rev(minmaxy),xlim=minmaxx)
	    if (additionalfigs && !is.null(modelParams) && !is.null(EvolModel)){
		if (EvolModel=="bm"){abline(v=modelParams$vX0[i-1,1],lty=2,lwd=1.5);mtext(side=3,at=modelParams$vX0[i-1,1],text=expression(X[0]),cex=2)}
		if (EvolModel=="ouou"){abline(v=modelParams$vY0[i-1,1],lty=2,lwd=1.5);abline(v=modelParams$mPsi[i-1,1],lty=2,lwd=1.5);mtext(side=3,at=modelParams$vY0[i-1,1],text=expression(X[0]),cex=2);mtext(side=3,at=modelParams$mPsi[i-1,1],text=expression(theta),cex=2)}
		if (EvolModel=="mvslouch"){
		    kY<-nrow(modelParams$A);kX<-nrow(modelParams$Sxx)
		    if (i-1<=kY){abline(v=modelParams$vY0[i-1,1],lty=2,lwd=1.5);abline(v=modelParams$mPsi[i-1,1],lty=2,lwd=1.5);mtext(side=3,at=modelParams$vY0[i-1,1],text=expression(Y[0]),cex=2);mtext(side=3,at=modelParams$mPsi[i-1,1],text=expression(psi),cex=2)}
		    else{abline(v=modelParams$vX0[i-1-kY,1],lty=2,lwd=1.5);mtext(side=3,at=modelParams$vX0[i-1-kY,1],text=expression(X[0]),cex=2)}
		}
	    }
	}
    }else{.my_stop("Error: wrong data provided to plotting function",TRUE)}
    NA        
}

