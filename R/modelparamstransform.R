## This file is part of mvSLOUCH

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this package or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .


.par_transform_withPCMBase<-function(params,EstimationParams,evolmodel,model_params=NULL){
    model_params<-.par.transform(params,EstimationParams,model_params)
    model_params$pcmbase_model_box<-.update_pcmbase_box_params(model_params,evolmodel)
    model_params
}

.update_pcmbase_box_params<-function(model_params,evolmodel,vDo=NULL){
## callled in OUphylregression.R
    pcmbase_model_box=switch(evolmodel,
                bm=.update_pcmbase_box_params_bm(model_params,vDo=NULL),
                ouch=.update_pcmbase_box_params_ouch(model_params,vDo=NULL),
#               slouch=.update_pcmbase_box_params_mvslouch(model_params,vDo=NULL),
                mvslouch=.update_pcmbase_box_params_mvslouch(model_params,vDo=NULL)
            )
    pcmbase_model_box
}

.update_pcmbase_box_params_bm<-function(model_params,vDo=NULL){
## at the moment in the BM case we do not use this functionality
    if (is.null(vDo)){vDo<-c("H"=FALSE,"Theta"=FALSE,"Sigma_x"=TRUE,"X0"=TRUE)}

    if (is.element("H",names(model_params$pcmbase_model_box))){model_params$pcmbase_model_box[which(names(model_params$pcmbase_model_box)=="H")]<-NULL}
    if (is.element("Theta",names(model_params$pcmbase_model_box))){model_params$pcmbase_model_box[which(names(model_params$pcmbase_model_box)=="Theta")]<-NULL}

    if (vDo["Sigma_x"]){
	for (cpcmbase_reg in names(model_params$pcmbase_model_box$Sigma_x[1,1,])){
	    model_params$pcmbase_model_box$Sigma_x[,,which(names(model_params$pcmbase_model_box$Sigma_x[1,1,])==cpcmbase_reg)]<- .changeSigmatoSyy(model_params$Sxx%*%t(model_params$Sxx),"UpperTri",NULL,NULL,FALSE)
	}
    }
    if (vDo["X0"]){
	model_params$pcmbase_model_box<-.set_pcmbase_model_box_X0(model_params$pcmbase_model_box, model_params$vX0[,1])
    }
    class(model_params$pcmbase_model_box)<-"BM"
    model_params$pcmbase_model_box
}


.update_pcmbase_box_params_ouch<-function(model_params,vDo=NULL){
    if (is.null(vDo)){vDo<-c("H"=TRUE,"Theta"=TRUE,"Sigma_x"=TRUE,"X0"=TRUE)}
    if ((length(model_params$mPsi)==1)&&(is.na(model_params$mPsi[1]))){vDo["Theta"]<-FALSE}
    if ((length(model_params$vY0)==1)&&(is.na(model_params$vY0[1]))){vDo["X0"]<-FALSE}

    if (vDo["H"]){
	for (cpcmbase_reg in names(model_params$pcmbase_model_box$H[1,1,])){
	    model_params$pcmbase_model_box$H[,,which(names(model_params$pcmbase_model_box$H[1,1,])==cpcmbase_reg)]<- model_params$A
	}
    }

    if (vDo["Sigma_x"]){
	for (cpcmbase_reg in names(model_params$pcmbase_model_box$Sigma_x[1,1,])){
	    model_params$pcmbase_model_box$Sigma_x[,,which(names(model_params$pcmbase_model_box$Sigma_x[1,1,])==cpcmbase_reg)]<- .changeSigmatoSyy(model_params$Syy%*%t(model_params$Syy),"UpperTri",NULL,NULL,FALSE)
	}
    }

    if (vDo["Theta"]){
    	for (cpcmbase_reg in names(model_params$pcmbase_model_box$Theta[1,])){
    	    ## strsplit is in base
    	    reg_name<-strsplit(cpcmbase_reg,split="_merrorregime_node")[[1]][1]
	    index_Psi_reg<-NA
	    ## PCMBase has regimes called by name, in mvSLOUCH mPsi has them by index of position in $regimes.types.orig
	    if (is.element("regimes.types.orig",names(model_params))){
		index_Psi_reg<-which(model_params$regimes.types.orig==reg_name)
	    }else{
		if (is.element(reg_name,colnames(model_params$mPsi))){
		    index_Psi_reg<-which(colnames(model_params$mPsi)==reg_name)
		}else{
		    if (length(names(model_params$pcmbase_model_box$Theta[1,]))==ncol(model_params$mPsi)){
			index_Psi_reg<-which(names(model_params$pcmbase_model_box$Theta[1,])==cpcmbase_reg)		    
		    }else{
			.my_stop("Cannot identify regimes in pcmabase_model_box_update")
		    }
		}
	    }
	    model_params$pcmbase_model_box$Theta[,which(names(model_params$pcmbase_model_box$Theta[1,])==cpcmbase_reg)] <- model_params$mPsi[,index_Psi_reg]
    	}
    }
    if (vDo["X0"]){
	model_params$pcmbase_model_box<-.set_pcmbase_model_box_X0(model_params$pcmbase_model_box,model_params$vY0[,1])
    }
    model_params$pcmbase_model_box
}

.update_pcmbase_box_params_mvslouch<-function(model_params,vDo=NULL){
    if (is.null(vDo)){vDo<-c("H"=TRUE,"Theta"=TRUE,"Sigma_x"=TRUE,"X0"=TRUE)}
    if ((length(model_params$mPsi)==1)&&(is.na(model_params$mPsi[1]))){vDo["Theta"]<-FALSE}
    if ((length(model_params$vY0)==1)&&(is.na(model_params$vY0[1]))){vDo["X0"]<-FALSE}
    if ((length(model_params$vX0)==1)&&(is.na(model_params$vX0[1]))){vDo["X0"]<-FALSE}

    kY<-ncol(model_params$A)
    kX<-ncol(model_params$Sxx)
    if (vDo["H"]){
	if(is.na(model_params$B[1])&&(length(model_params$B)==1)){model_params$B<-matrix(0,nrow=kY,ncol=kX)}
	mAB<-rbind(cbind(model_params$A,model_params$B),cbind(matrix(0,nrow=kX,ncol=kY),matrix(0,nrow=kX,ncol=kX)))
	for (cpcmbase_reg in names(model_params$pcmbase_model_box$H[1,1,])){
	    model_params$pcmbase_model_box$H[,,which(names(model_params$pcmbase_model_box$H[1,1,])==cpcmbase_reg)]<- mAB
	}
    }

    if (vDo["Sigma_x"]){
	mV<-rbind(cbind(model_params$Syy%*%t(model_params$Syy)+model_params$Syx%*%t(model_params$Syx),model_params$Syy%*%t(model_params$Sxy)+model_params$Syx%*%t(model_params$Sxx)),cbind(model_params$Sxy%*%t(model_params$Syy)+model_params$Sxx%*%t(model_params$Syx),model_params$Sxy%*%t(model_params$Sxy)+model_params$Sxx%*%t(model_params$Sxx)))
	for (cpcmbase_reg in names(model_params$pcmbase_model_box$Sigma_x[1,1,])){	        
	    model_params$pcmbase_model_box$Sigma_x[,,which(names(model_params$pcmbase_model_box$Sigma_x[1,1,])==cpcmbase_reg)]<- .changeSigmatoSyy(mV,"UpperTri",NULL,NULL,FALSE)
	}
    }

    if (vDo["Theta"]){
    	for (cpcmbase_reg in names(model_params$pcmbase_model_box$Theta[1,])){
    	    ## strsplit is in base
    	    reg_name<-strsplit(cpcmbase_reg,split="_merrorregime_node")[[1]][1]
	    index_Psi_reg<-NA
	    ## PCMBase has regimes called by name, in mvSLOUCH mPsi has them by index of position in $regimes.types.orig
	    if (is.element("regimes.types.orig",names(model_params))){
		index_Psi_reg<-which(model_params$regimes.types.orig==reg_name)
	    }else{
		if (is.element(reg_name,colnames(model_params$mPsi))){
		    index_Psi_reg<-which(colnames(model_params$mPsi)==reg_name)
		}else{
		    if (length(names(model_params$pcmbase_model_box$Theta[1,]))==ncol(model_params$mPsi)){
			index_Psi_reg<-which(names(model_params$pcmbase_model_box$Theta[1,])==cpcmbase_reg)		    
		    }else{
			.my_stop("Cannot identify regimes in pcmabase_model_box_update")
		    }
		}	
	    }
	    model_params$pcmbase_model_box$Theta[,which(names(model_params$pcmbase_model_box$Theta[1,])==cpcmbase_reg)]<- c(model_params$mPsi[,index_Psi_reg],rep(0,kX))
    	}
    }
    if (vDo["X0"]){
	model_params$pcmbase_model_box<-.set_pcmbase_model_box_X0(model_params$pcmbase_model_box,c(model_params$vY0[,1],model_params$vX0[,1]))
    }
    model_params$pcmbase_model_box
}

.par.transform<-function(params,EstimationParams,ModelParams=NULL){
## same parametrization for all models, since all are nested in the mvslouch
## if some parameters are not needed then one uses NA (NOT NULL)
    if (is.null(ModelParams)){
	ModelParams<-vector("list",13)
	names(ModelParams)<-c("A","B","mPsi","mPsi0","vY0","vX0","Syy","Syx","Sxy","Sxx","eigenSignsA","GivensQCsignsA","pcmbase_model_box")
	ModelParams$pcmbase_model_box<-EstimationParams$pcmbase_model_box
    }
    if (!is.null(EstimationParams$Fixed$A)){ModelParams$A<-EstimationParams$Fixed$A}
    else{
	if (EstimationParams$kY==1){
	    if (is.element("A",names(params))){ModelParams$A<-matrix(params["A"],ncol=1,nrow=1)}
	    else{ModelParams$A<-matrix(NA,ncol=1,nrow=1)}
	}
	else{
	    lAparams=switch(EstimationParams$Atype,
		SingleValueDiagonal={diag(params["A"],ncol=EstimationParams$kY,nrow=EstimationParams$kY)},
		Diagonal={diag(params[(which(names(params)=="Astart")):(which(names(params)=="Aend"))],ncol=EstimationParams$kY,nrow=EstimationParams$kY)},
		SymmetricPositiveDefinite={.sym.par(params[(which(names(params)=="Astart")):(which(names(params)=="Aend"))],EstimationParams$kY)},
		Symmetric={.par.transform.symmetric.matrix(params[(which(names(params)=="Astart")):(which(names(params)=="Aend"))],EstimationParams$kY)},
		TwoByTwo={.par.transform.twobytwo.matrix(params[(which(names(params)=="Astart")):(which(names(params)=="Aend"))])},
		UpperTri={.par.transform.uppertri.matrix(params[(which(names(params)=="Astart")):(which(names(params)=="Aend"))],EstimationParams$kY)},
		LowerTri={.par.transform.lowertri.matrix(params[(which(names(params)=="Astart")):(which(names(params)=="Aend"))],EstimationParams$kY)},
		DecomposablePositive={.par.transform.decomp.pos.real.matrix(params[(which(names(params)=="Astart")):(which(names(params)=="Aend"))],EstimationParams$kY,NULL,ifelse(is.null(EstimationParams$minAeigen),0.01,EstimationParams$minAeigen))},
		DecomposableNegative={.par.transform.decomp.neg.real.matrix(params[(which(names(params)=="Astart")):(which(names(params)=="Aend"))],EstimationParams$kY,NULL,ifelse(is.null(EstimationParams$maxAeigen),0.01,EstimationParams$maxAeigen))},
		DecomposableReal={
		    .par.transform.decomp.real.matrix(params[(which(names(params)=="Astart")):(which(names(params)=="Aend"))],EstimationParams$kY)
		}, 
		Invertible={.par.transform.invert.matrix(params[(which(names(params)=="Astart")):(which(names(params)=="Aend"))],EstimationParams$kY,"qr",ifelse(is.null(EstimationParams$minRdiag),0.01,EstimationParams$minRdiag))},
		Any={
		    if (!is.null(EstimationParams$signsA)){
			tmpsignsA<-EstimationParams$signsA
			vToNA<-c(which(tmpsignsA=="-"),which(tmpsignsA=="+"))
			Aval<-matrix(NA,ncol=ncol(tmpsignsA),nrow=nrow(tmpsignsA))
			if (length(vToNA)>0){tmpsignsA[vToNA]<-NA}			
			if (length(which(!is.na(tmpsignsA)))>0){
			    Aval[which(!is.na(tmpsignsA))]<-as.numeric(tmpsignsA[which(!is.na(tmpsignsA))])
			}
			if (length(which(is.na(tmpsignsA)))>0){
			    if (length(which(is.na(tmpsignsA)))==1){Aval[which(is.na(tmpsignsA))]<-params[which(names(params)=="A")]}
			    else{Aval[which(is.na(tmpsignsA))]<-params[(which(names(params)=="Astart")):(which(names(params)=="Aend"))]}
			}
		    }else{
			Aval<-matrix(params[(which(names(params)=="Astart")):(which(names(params)=="Aend"))],ncol=EstimationParams$kY,nrow=EstimationParams$kY,byrow=TRUE)
		    }
		    Aval
		}
	    )
	    if (is.element(EstimationParams$Atype,c("DecomposablePositive","DecomposableNegative","DecomposableReal","Invertible"))){
		ModelParams$A<-lAparams$A	    
		ModelParams$GivensQCsignsA<-lAparams$QcSigns
		if (is.element(EstimationParams$Atype,c("DecomposablePositive","DecomposableNegative","DecomposableReal"))){ 
		    ModelParams$eigenSignsA<-lAparams$eigenSigns
		}
	    }
	    else{ModelParams$A<-lAparams}
	}
	if ((!is.null(EstimationParams$diagA))&&(EstimationParams$Atype!="SymmetricPositiveDefinite")){
	    diag(ModelParams$A)=switch(EstimationParams$diagA,
		Positive={exp(diag(ModelParams$A))+ifelse(is.null(EstimationParams$minAdiag),0,EstimationParams$minAdiag)},
		Negative={(-1)*exp(diag(ModelParams$A))-ifelse(is.null(EstimationParams$maxAdiag),0,EstimationParams$maxAdiag)}
	    )
	}	  
	if (!is.null(EstimationParams$signsA)){ ## The user is allowed to specify signs in A but NO check is done whether now A will remain in the desired matrix class
	    ModelParams$A[which(EstimationParams$signsA=="0")]<-0
	    EstimationParams$signsA[which(EstimationParams$signsA=="0")]<-NA
	    ModelParams$A[which(EstimationParams$signsA==0)]<-0
	    EstimationParams$signsA[which(EstimationParams$signsA==0)]<-NA
	    ModelParams$A[which(EstimationParams$signsA=="-")]<- (-1)*exp(ModelParams$A[which(EstimationParams$signsA=="-")])
	    EstimationParams$signsA[which(EstimationParams$signsA=="-")]<-NA
	    ModelParams$A[which(EstimationParams$signsA=="+")]<- exp(ModelParams$A[which(EstimationParams$signsA=="+")])
	    EstimationParams$signsA[which(EstimationParams$signsA=="+")]<-NA
	    if (length(which(!is.na(EstimationParams$signsA)))>0){
	    ## we allow for setting any given pre-specified value here
		    ModelParams$A[which(!is.na(EstimationParams$signsA))]<-as.numeric(EstimationParams$signsA[which(!is.na(EstimationParams$signsA))])
	    }
	}
#	if (is.element("maxAabsval",names(EstimationParams))){
#	    mSignsA<-sign(ModelParams$A)
#	    ModelParams$A[which(abs(ModelParams$A)>EstimationParams$maxAabsval)]<-EstimationParams$maxAabsval
#	    ModelParams$A<-abs(ModelParams$A)*mSignsA
#	}
	ModelParams$A[which(abs(ModelParams$A)<1e-15)]<-0
    }

    if (!is.null(EstimationParams$Fixed$B)){ModelParams$B<-EstimationParams$Fixed$B}
    else{
	if ((EstimationParams$kX==1)&&(EstimationParams$kY==1)){
	    if (is.element("B",names(params))){ModelParams$B<-matrix(params["B"],ncol=1,nrow=1)}
	    else{ModelParams$B<-matrix(NA,ncol=1,nrow=1)}
	}
	else{
	    ModelParams$B=switch(EstimationParams$Btype,
		MinusA={(-1)*ModelParams$A},
		SingleValue={matrix(params["B"],ncol=EstimationParams$kX,nrow=EstimationParams$kY)},
		SingleValueDiagonal={diag(params["B"],ncol=EstimationParams$kX,nrow=EstimationParams$kY)},
	        Diagonal={diag(params[(which(names(params)=="Bstart")):(which(names(params)=="Bend"))],ncol=EstimationParams$kX,nrow=EstimationParams$kY)},
		Any={
		    if (!is.null(EstimationParams$signsB)){
			tmpsignsB<-EstimationParams$signsB
			vToNA<-c(which(tmpsignsB=="-"),which(tmpsignsB=="+"))
		        Bval<-matrix(NA,ncol=ncol(tmpsignsB),nrow=nrow(tmpsignsB))
			if (length(vToNA)>0){tmpsignsA[vToNA]<-NA}
			if (length(which(!is.na(tmpsignsB)))>0){
			    Bval[which(!is.na(tmpsignsB))]<-as.numeric(tmpsignsB[which(!is.na(tmpsignsB))])
			}
			if (length(which(is.na(tmpsignsB)))>0){
			    if (length(which(is.na(tmpsignsB)))==1){Bval[which(is.na(tmpsignsB))]<-params[which(names(params)=="B")]}
			    else{Bval[which(is.na(tmpsignsB))]<-params[(which(names(params)=="Bstart")):(which(names(params)=="Bend"))]}
			}
		    }else{
			Bval<- matrix(params[(which(names(params)=="Bstart")):(which(names(params)=="Bend"))],ncol=EstimationParams$kX,nrow=EstimationParams$kY,byrow=TRUE)
		    }
		    Bval
		}
	    )
	}
	if (!is.null(EstimationParams$signsB)){ ## The user is allowed to specify signs in B but NO check is done whether now B will remain in the desired matrix class
	    ModelParams$B[which(EstimationParams$signsB=="0")]<-0
	    EstimationParams$signsB[which(EstimationParams$signsB=="0")]<-NA
	    ModelParams$B[which(EstimationParams$signsB==0)]<-0
	    EstimationParams$signsB[which(EstimationParams$signsB==0)]<-NA
	    ModelParams$B[which(EstimationParams$signsB=="-")]<- (-1)*exp(ModelParams$B[which(EstimationParams$signsB=="-")])
    	    EstimationParams$signsB[which(EstimationParams$signsB=="-")]<-NA
	    ModelParams$B[which(EstimationParams$signsB=="+")]<- exp(ModelParams$B[which(EstimationParams$signsB=="+")])
    	    EstimationParams$signsB[which(EstimationParams$signsB=="+")]<-NA
	    if (length(which(!is.na(EstimationParams$signsB)))>0){
	    ## we allow for setting any given pre-specified value here
		    ModelParams$B[which(!is.na(EstimationParams$signsB))]<-as.numeric(EstimationParams$signsB[which(!is.na(EstimationParams$signsB))])
	    }
	}
#    	if (is.element("maxAabsval",names(EstimationParams))){
#	    mSignsB<-sign(ModelParams$B)
#	    ModelParams$B[which(abs(ModelParams$B)>EstimationParams$maxAabsval)]<-EstimationParams$maxAabsval
#	    ModelParams$B<-abs(ModelParams$B)*mSignsB
#	}
    	ModelParams$B[which(abs(ModelParams$B)<1e-15)]<-0
    }
    
    if (!is.null(EstimationParams$Fixed$mPsi)){ModelParams$mPsi<-EstimationParams$Fixed$mPsi}
    else{
	if ((EstimationParams$kY==1)&&((EstimationParams$mPsitype=="Global")||((EstimationParams$mPsitype=="Regimes")&&(length(EstimationParams$RegimeTypes)==1)))){
	    if (is.element("Psi",names(params))){ModelParams$mPsi<-matrix(params["Psi"],ncol=1,nrow=1)}
	    else{ModelParams$mPsi<-matrix(NA,ncol=1,nrow=1)}
	    if(EstimationParams$mPsitype=="Regimes"){colnames(ModelParams$mPsi)<-EstimationParams$RegimeTypes}
	}
	else{
	    ModelParams$mPsi=switch(EstimationParams$mPsitype,
		Global={matrix(params[(which(names(params)=="Psistart")):(which(names(params)=="Psiend"))],ncol=1)},
		Regimes={
		    mPsi<-matrix(params[(which(names(params)=="Psistart")):(which(names(params)=="Psiend"))],nrow=EstimationParams$kY,ncol=length(EstimationParams$RegimeTypes),byrow=FALSE)
		    names(mPsi)<-EstimationParams$RegimeTypes
		    mPsi
		}	    
	    )    
	}
    	if (!is.null(EstimationParams$signsmPsi)){ ## The user is allowed to specify signs in A but NO check is done whether now A will remain in the desired matrix class
	    ModelParams$mPsi[which(EstimationParams$signsmPsi=="0")]<-0
	    EstimationParams$signsmPsi[which(EstimationParams$signsmPsi=="0")]<-NA
	    ModelParams$mPsi[which(EstimationParams$signsmPsi==0)]<-0
	    EstimationParams$signsmPsi[which(EstimationParams$signsmPsi==0)]<-NA
	    ModelParams$mPsi[which(EstimationParams$signsmPsi=="-")]<- (-1)*exp(ModelParams$mPsi[which(EstimationParams$signsmPsi=="-")])
    	    EstimationParams$signsmPsi[which(EstimationParams$signsmPsi=="-")]<-NA
	    ModelParams$mPsi[which(EstimationParams$signsmPsi=="+")]<- exp(ModelParams$mPsi[which(EstimationParams$signsmPsi=="+")])
    	    EstimationParams$signsmPsi[which(EstimationParams$signsmPsi=="+")]<-NA
	    if (length(which(!is.na(EstimationParams$signsmPsi)))>0){
	    ## we allow for setting any given pre-specified value here
		    ModelParams$mPsi[which(!is.na(EstimationParams$signsmPsi))]<-as.numeric(EstimationParams$signsmPsi[which(!is.na(EstimationParams$signsmPsi))])
	    }
	}
    	ModelParams$mPsi[which(abs(ModelParams$mPsi)<1e-15)]<-0
    }
    if (!is.null(EstimationParams$Fixed$mPsi0)){ModelParams$mPsi0<-EstimationParams$Fixed$mPsi0}
    else{
	if (EstimationParams$kY==1){
	    if (is.element("Psi0",names(params))){ModelParams$mPsi0<-matrix(params["Psi0"],nrow=1,ncol=1)}
	    else{ModelParams$mPsi0<-matrix(NA,nrow=1,ncol=1)}
	}
	else{ModelParams$mPsi0<-params[(which(names(params)=="Psi0start")):(which(names(params)=="Psi0end"))]}
    	if (!is.null(EstimationParams$signsmPsi0)){ ## The user is allowed to specify signs in A but NO check is done whether now A will remain in the desired matrix class
	    ModelParams$mPsi0[which(EstimationParams$signsmPsi0=="0")]<-0
	    EstimationParams$signsmPsi0[which(EstimationParams$signsmPsi0=="0")]<-NA
	    ModelParams$mPsi0[which(EstimationParams$signsmPsi0==0)]<-0
	    EstimationParams$signsmPsi0[which(EstimationParams$signsmPsi0==0)]<-NA
	    ModelParams$mPsi0[which(EstimationParams$signsmPsi0=="-")]<- (-1)*exp(ModelParams$mPsi0[which(EstimationParams$signsmPsi0=="-")])
	    EstimationParams$signsmPsi0[which(EstimationParams$signsmPsi0=="-")]<-NA
	    ModelParams$mPsi0[which(EstimationParams$signsmPsi0=="+")]<- exp(ModelParams$mPsi0[which(EstimationParams$signsmPsi0=="+")])
	    EstimationParams$signsmPsi0[which(EstimationParams$signsmPsi0=="+")]<-NA
	    if (length(which(!is.na(EstimationParams$signsmPsi0)))>0){
	    ## we allow for setting any given pre-specified value here
		    ModelParams$mPsi0[which(!is.na(EstimationParams$signsmPsi0))]<-as.numeric(EstimationParams$signsmPsi0[which(!is.na(EstimationParams$signsmPsi0))])
	    }
	}
    	ModelParams$mPsi0[which(abs(ModelParams$mPsi0)<1e-15)]<-0
    }
    if (!is.null(EstimationParams$Fixed$vY0)){ModelParams$vY0<-EstimationParams$Fixed$vY0}
    else{
	if (EstimationParams$kY==1){
	    if (is.element("vY0",names(params))){ModelParams$vY0<-matrix(params["vY0"],ncol=1,nrow=1)}
	    else{ModelParams$vY0<-matrix(NA,ncol=1,nrow=1)}
	}
	else{ModelParams$vY0<-params[(which(names(params)=="vY0start")):(which(names(params)=="vY0end"))]}
    	if (!is.null(EstimationParams$signsvY0)){ ## The user is allowed to specify signs in A but NO check is done whether now A will remain in the desired matrix class
	    ModelParams$vY0[which(EstimationParams$signsvY0=="0")]<-0
	    EstimationParams$signsvY0[which(EstimationParams$signsvY0=="0")]<-NA
	    ModelParams$vY0[which(EstimationParams$signsvY0==0)]<-0
	    EstimationParams$signsvY0[which(EstimationParams$signsvY0==0)]<-NA
	    ModelParams$vY0[which(EstimationParams$signsvY0=="-")]<- (-1)*exp(ModelParams$vY0[which(EstimationParams$signsvY0=="-")])
	    EstimationParams$signsvY0[which(EstimationParams$signsvY0=="-")]<-NA
	    ModelParams$vY0[which(EstimationParams$signsvY0=="+")]<- exp(ModelParams$vY0[which(EstimationParams$signsvY0=="+")])
	    EstimationParams$signsvY0[which(EstimationParams$signsvY0=="+")]<-NA
	    if (length(which(!is.na(EstimationParams$signsvY0)))>0){
	    ## we allow for setting any given pre-specified value here
		    ModelParams$vY0[which(!is.na(EstimationParams$signsvY0))]<-as.numeric(EstimationParams$signsvY0[which(!is.na(EstimationParams$signsvY0))])
	    }
	}
    	ModelParams$vY0[which(abs(ModelParams$vY0)<1e-15)]<-0
    }
    if (!is.null(EstimationParams$Fixed$vX0)){ModelParams$vX0<-EstimationParams$Fixed$vX0}
    else{
	if (EstimationParams$kX==1){
	    if (is.element("vX0",names(params))){ModelParams$vX0<-matrix(params["vX0"],nrow=1,ncol=1)}
	    else{ModelParams$vX0<-matrix(NA,nrow=1,ncol=1)}
	}
	else{ModelParams$vX0<-params[(which(names(params)=="vX0start")):(which(names(params)=="vX0end"))]}
    	if (!is.null(EstimationParams$signsvX0)){ ## The user is allowed to specify signs in A but NO check is done whether now A will remain in the desired matrix class
	    ModelParams$vX0[which(EstimationParams$signsvX0=="0")]<-0
	    EstimationParams$signsvX0[which(EstimationParams$signsvX0=="0")]<-NA
	    ModelParams$vX0[which(EstimationParams$signsvX0==0)]<-0
	    EstimationParams$signsvX0[which(EstimationParams$signsvX0==0)]<-NA
	    ModelParams$vX0[which(EstimationParams$signsvX0=="-")]<- (-1)*exp(ModelParams$vX0[which(EstimationParams$signsvX0=="-")])
	    EstimationParams$signsvX0[which(EstimationParams$signsvX0=="-")]<-NA
	    ModelParams$vX0[which(EstimationParams$signsvX0=="+")]<- exp(ModelParams$vX0[which(EstimationParams$signsvX0=="+")])
	    EstimationParams$signsvX0[which(EstimationParams$signsvX0=="+")]<-NA
	    if (length(which(!is.na(EstimationParams$signsvX0)))>0){
	    ## we allow for setting any given pre-specified value here
		    ModelParams$vX0[which(!is.na(EstimationParams$signsvX0))]<-as.numeric(EstimationParams$signsvX0[which(!is.na(EstimationParams$signsvX0))])
	    }
	}
    	ModelParams$vX0[which(abs(ModelParams$vX0)<1e-15)]<-0
    }
    if (!is.null(EstimationParams$Fixed$Syy)){ModelParams$Syy<-EstimationParams$Fixed$Syy }
    else{ ## Generally symmetric structures for these matrices should suffice since we do Ct(C) anyway but just in case general 
	if (EstimationParams$kY==1){
	    if (is.element("Syy",names(params))){ModelParams$Syy<-matrix(params["Syy"],ncol=1,nrow=1)}
	    else{ModelParams$Syy<-matrix(NA,ncol=1,nrow=1)}
	}
	else{
	    ModelParams$Syy=switch(EstimationParams$Syytype,
		SingleValueDiagonal={diag(params["Syy"],ncol=EstimationParams$kY,nrow=EstimationParams$kY)},
		Diagonal={diag(params[(which(names(params)=="Syystart")):(which(names(params)=="Syyend"))],ncol=EstimationParams$kY,nrow=EstimationParams$kY)},
		Symmetric={.sym.par(params[(which(names(params)=="Syystart")):(which(names(params)=="Syyend"))],EstimationParams$kY)},
		UpperTri={.par.transform.uppertri.matrix(params[(which(names(params)=="Syystart")):(which(names(params)=="Syyend"))],EstimationParams$kY)},
		LowerTri={.par.transform.lowertri.matrix(params[(which(names(params)=="Syystart")):(which(names(params)=="Syyend"))],EstimationParams$kY)},
		Any={
		    if (!is.null(EstimationParams$signsSyy)){
			tmpsignsSyy<-EstimationParams$signsSyy
			vToNA<-c(which(tmpsignsSyy=="-"),which(tmpsignsSyy=="+"))
			Syyval<-matrix(NA,ncol=ncol(tmpsignsSyy),nrow=nrow(tmpsignsSyy))
			if (length(vToNA)>0){tmpsignsSyy[vToNA]<-NA}
			if (length(which(!is.na(tmpsignsSyy)))>0){
			    Syyval[which(!is.na(tmpsignsSyy))]<-as.numeric(tmpsignsSyy[which(!is.na(tmpsignsSyy))])
			}
			if (length(which(is.na(tmpsignsSyy)))>0){
			    if (length(which(is.na(tmpsignsSyy)))==1){Syyval[which(is.na(tmpsignsSyy))]<-params[which(names(params)=="Syy")]}
			    else{Syyval[which(is.na(tmpsignsSyy))]<-params[(which(names(params)=="Syystart")):(which(names(params)=="Syyend"))]}
			}
		    }else{Syyval<- matrix(params[(which(names(params)=="Syystart")):(which(names(params)=="Syyend"))],ncol=EstimationParams$kY,nrow=EstimationParams$kY,byrow=TRUE)}
		    Syyval
		    }
	    )	
	}
	if (!is.null(EstimationParams$diagSyy)){
		diag(ModelParams$Syy)=switch(EstimationParams$diagSyy,
		    Positive={exp(diag(ModelParams$Syy))},
		    Negative={(-1)*exp(diag(ModelParams$Syy))}
		)
	}
	if (!is.null(EstimationParams$signsSyy)){ ## The user is allowed to specify signs in A but NO check is done whether now A will remain in the desired matrix class
	    ModelParams$Syy[which(EstimationParams$signsSyy=="0")]<-0
	    EstimationParams$signsSyy[which(EstimationParams$signsSyy=="0")]<-NA
	    ModelParams$Syy[which(EstimationParams$signsSyy==0)]<-0
	    EstimationParams$signsSyy[which(EstimationParams$signsSyy==0)]<-NA
	    ModelParams$Syy[which(EstimationParams$signsSyy=="-")]<- (-1)*exp(ModelParams$Syy[which(EstimationParams$signsSyy=="-")])
	    EstimationParams$signsSyy[which(EstimationParams$signsSyy=="-")]<-NA
	    ModelParams$Syy[which(EstimationParams$signsSyy=="+")]<- exp(ModelParams$Syy[which(EstimationParams$signsSyy=="+")])
	    EstimationParams$signsSyy[which(EstimationParams$signsSyy=="+")]<-NA
	    if (length(which(!is.na(EstimationParams$signsSyy)))>0){
	    ## we allow for setting any given pre-specified value here
		    ModelParams$Syy[which(!is.na(EstimationParams$signsSyy))]<-as.numeric(EstimationParams$signsSyy[which(!is.na(EstimationParams$signsSyy))])
	    }
	}
    	if (is.element("maxSyyabsval",names(EstimationParams))){
	    mSignsSyy<-sign(ModelParams$Syy)
	    ModelParams$Syy[which(abs(ModelParams$Syy)>EstimationParams$maxSyyabsval)]<-EstimationParams$maxSyyabsval
	    ModelParams$Syy<-abs(ModelParams$Syy)*mSignsSyy
	}    	
    	ModelParams$Syy[which(abs(ModelParams$Syy)<1e-15)]<-0
    }
    if (!is.null(EstimationParams$Fixed$Syx)){ModelParams$Syx<-EstimationParams$Fixed$Syx} 
    else{
	if ((EstimationParams$kX==1)&&(EstimationParams$kY==1)){
	    if (is.element("Syx",names(params))){ModelParams$Syx<-matrix(params["Syx"],ncol=1,nrow=1)}
	    else{ModelParams$Syx<-matrix(NA,ncol=1,nrow=1)}
	}
	else{
	    if (!is.null(EstimationParams$signsSyx)){
			tmpsignsSyx<-EstimationParams$signsSyx
			vToNA<-c(which(tmpsignsSyx=="-"),which(tmpsignsSyx=="+"))
			if (length(vToNA)>0){tmpsignsSyx[vToNA]<-NA}
			if (length(which(!is.na(tmpsignsSyx)))>0){
			    ModelParams$Syx<-matrix(NA,ncol=ncol(tmpsignsSyx),nrow=nrow(tmpsignsSyx))
			    ModelParams$Syx[which(!is.na(tmpsignsSyx))]<-as.numeric(tmpsignsSyx[which(!is.na(tmpsignsSyx))])
			    if (length(which(is.na(tmpsignsSyx)))>0){
				if (length(which(is.na(tmpsignsSyx)))==1){ModelParams$Syx[which(is.na(tmpsignsSyx))]<-params[which(names(params)=="Syx")]}
				else{ModelParams$Syx[which(is.na(tmpsignsSyx))]<-params[(which(names(params)=="Syxstart")):(which(names(params)=="Syxend"))]}
			    }
			}
	    }else{ModelParams$Syx<-matrix(params[(which(names(params)=="Syxstart")):(which(names(params)=="Syxend"))],ncol=EstimationParams$kX,nrow=EstimationParams$kY,byrow=TRUE)}
	}
    	if (!is.null(EstimationParams$signsSyx)){ ## The user is allowed to specify signs in A but NO check is done whether now A will remain in the desired matrix class
	    ModelParams$Syx[which(EstimationParams$signsSyx=="0")]<-0
	    EstimationParams$signsSyx[which(EstimationParams$signsSyx=="0")]<-NA
	    ModelParams$Syx[which(EstimationParams$signsSyx==0)]<-0
	    EstimationParams$signsSyx[which(EstimationParams$signsSyx==0)]<-NA
	    ModelParams$Syx[which(EstimationParams$signsSyx=="-")]<- (-1)*exp(ModelParams$Syx[which(EstimationParams$signsSyx=="-")])
	    EstimationParams$signsSyx[which(EstimationParams$signsSyx=="-")]<-NA
	    ModelParams$Syx[which(EstimationParams$signsSyx=="+")]<- exp(ModelParams$Syx[which(EstimationParams$signsSyx=="+")])
	    EstimationParams$signsSyx[which(EstimationParams$signsSyx=="+")]<-NA
	    if (length(which(!is.na(EstimationParams$signsSyx)))>0){
	    ## we allow for setting any given pre-specified value here
		    ModelParams$Syx[which(!is.na(EstimationParams$signsSyx))]<-as.numeric(EstimationParams$signsSyx[which(!is.na(EstimationParams$signsSyx))])
	    }
	}
    	ModelParams$Syx[which(abs(ModelParams$Syx)<1e-15)]<-0

    }
    if (!(is.null(EstimationParams$Fixed$Sxy))){ModelParams$Sxy<-EstimationParams$Fixed$Sxy} 
    else{
    	if ((EstimationParams$kX==1)&&(EstimationParams$kY==1)){
    	    if (is.element("Sxy",names(params))){ModelParams$Sxy<-matrix(params["Sxy"],ncol=1,nrow=1)}
    	    else{ModelParams$Sxy<-matrix(NA,ncol=1,nrow=1)}
    	}
	else{
	    if (!is.null(EstimationParams$signsSxy)){
			tmpsignsSxy<-EstimationParams$signsSxy
			vToNA<-c(which(tmpsignsSxy=="-"),which(tmpsignsSxy=="+"))
			if (length(vToNA)>0){tmpsignsSxy[vToNA]<-NA}
			if (length(which(!is.na(tmpsignsSxy)))>0){
			    ModelParams$Sxy<-matrix(NA,ncol=ncol(tmpsignsSxy),nrow=nrow(tmpsignsSxy))
			    ModelParams$Sxy[which(!is.na(tmpsignsSxy))]<-as.numeric(tmpsignsSxy[which(!is.na(tmpsignsSxy))])
			    if (length(which(is.na(tmpsignsSxy)))>0){
				if (length(which(is.na(tmpsignsSxy)))==1){ModelParams$Sxy[which(is.na(tmpsignsSxy))]<-params[which(names(params)=="Sxy")]}
				else{ModelParams$Sxy[which(is.na(tmpsignsSxy))]<-params[(which(names(params)=="Sxystart")):(which(names(params)=="Sxyend"))]}
			    }
			}
	    }else{ModelParams$Sxy<-matrix(params[(which(names(params)=="Sxystart")):(which(names(params)=="Sxyend"))],ncol=EstimationParams$kY,nrow=EstimationParams$kX,byrow=TRUE)}
	}    
	if (!is.null(EstimationParams$signsSxy)){ ## The user is allowed to specify signs in A but NO check is done whether now A will remain in the desired matrix class
	    ModelParams$Sxy[which(EstimationParams$signsSxy=="0")]<-0
	    EstimationParams$signsSxy[which(EstimationParams$signsSxy=="0")]<-NA
	    ModelParams$Sxy[which(EstimationParams$signsSxy==0)]<-0
	    EstimationParams$signsSxy[which(EstimationParams$signsSxy==0)]<-NA
	    ModelParams$Sxy[which(EstimationParams$signsSxy=="-")]<- (-1)*exp(ModelParams$Sxy[which(EstimationParams$signsSxy=="-")])
	    EstimationParams$signsSxy[which(EstimationParams$signsSxy=="-")]<-NA
	    ModelParams$Sxy[which(EstimationParams$signsSxy=="+")]<- exp(ModelParams$Sxy[which(EstimationParams$signsSxy=="+")])
	    EstimationParams$signsSxy[which(EstimationParams$signsSxy=="+")]<-NA
	    if (length(which(!is.na(EstimationParams$signsSxy)))>0){
	    ## we allow for setting any given pre-specified value here
		    ModelParams$Sxy[which(!is.na(EstimationParams$signsSxy))]<-as.numeric(EstimationParams$signsSxy[which(!is.na(EstimationParams$signsSxy))])
	    }
	}
    	ModelParams$Sxy[which(abs(ModelParams$Sxy)<1e-15)]<-0

    }
    if (!is.null(EstimationParams$Fixed$Sxx)){ModelParams$Sxx<-EstimationParams$Fixed$Sxx}
    else{  ## Generally symmetric structures for these matrices should suffice since we do Ct(C) anyway but just in case general 
	if (EstimationParams$kX==1){
	    if (is.element("Sxx",names(params))){ModelParams$Sxx<-matrix(params["Sxx"],ncol=1,nrow=1)}
	    else{ModelParams$Sxx<-matrix(NA,ncol=1,nrow=1)}
	}
	else{
	    ModelParams$Sxx=switch(EstimationParams$Sxxtype,
		SingleValueDiagonal={diag(params["Sxx"],ncol=EstimationParams$kX,nrow=EstimationParams$kX)},
		Diagonal={diag(params[(which(names(params)=="Sxxstart")):(which(names(params)=="Sxxend"))],ncol=EstimationParams$kX,nrow=EstimationParams$kX)},
		Symmetric={.sym.par(params[(which(names(params)=="Sxxstart")):(which(names(params)=="Sxxend"))],EstimationParams$kX)},
		Any={
		    if (!is.null(EstimationParams$signsSxx)){
			tmpsignsSxx<-EstimationParams$signsSxx
			vToNA<-c(which(tmpsignsSxx=="-"),which(tmpsignsSxx=="+"))
			Sxxval<-matrix(NA,ncol=ncol(tmpsignsSxx),nrow=nrow(tmpsignsSxx))
			if (length(vToNA)>0){tmpsignsSxx[vToNA]<-NA}
			if (length(which(!is.na(tmpsignsSxx)))>0){
			    Sxxval[which(!is.na(tmpsignsSxx))]<-as.numeric(tmpsignsSxx[which(!is.na(tmpsignsSxx))])
			}
    			if (length(which(is.na(tmpsignsSxx)))>0){
			    if (length(which(is.na(tmpsignsSxx)))==1){Sxxval[which(is.na(tmpsignsSxx))]<-params[which(names(params)=="Sxx")]}
			    else{Sxxval[which(is.na(tmpsignsSxx))]<-params[(which(names(params)=="Sxxstart")):(which(names(params)=="Sxxend"))]}
			}
		    }else{Sxxval<-matrix(params[(which(names(params)=="Sxxstart")):(which(names(params)=="Sxxend"))],ncol=EstimationParams$kX,nrow=EstimationParams$kX,byrow=TRUE)}
		    Sxxval
		}
	    )
	}
	if (!is.null(EstimationParams$signsSxx)){ ## The user is allowed to specify signs in A but NO check is done whether now A will remain in the desired matrix class
	    ModelParams$Sxx[which(EstimationParams$signsSxx=="0")]<-0
	    EstimationParams$signsSxx[which(EstimationParams$signsSxx=="0")]<-NA
	    ModelParams$Sxx[which(EstimationParams$signsSxx==0)]<-0
	    EstimationParams$signsSxx[which(EstimationParams$signsSxx==0)]<-NA
	    ModelParams$Sxx[which(EstimationParams$signsSxx=="-")]<- (-1)*exp(ModelParams$Sxx[which(EstimationParams$signsSxx=="-")])
	    EstimationParams$signsSxx[which(EstimationParams$signsSxx=="-")]<-NA
	    ModelParams$Sxx[which(EstimationParams$signsSxx=="+")]<- exp(ModelParams$Sxx[which(EstimationParams$signsSxx=="+")])
	    EstimationParams$signsSxx[which(EstimationParams$signsSxx=="+")]<-NA
	    if (length(which(!is.na(EstimationParams$signsSxx)))>0){
	    ## we allow for setting any given pre-specified value here
		    ModelParams$Sxx[which(!is.na(EstimationParams$signsSxx))]<-as.numeric(EstimationParams$signsSxx[which(!is.na(EstimationParams$signsSxx))])
	    }
	}
	ModelParams$Sxx[which(abs(ModelParams$Sxx)<1e-15)]<-0
    }
    ModelParams    
}

.par.inv.transform<-function(ModelParams,EstimationParams){
## called in PhyloSDEestim.R
## same parametrization for all models, since all are nested in the mvslouch
## if some parameters are not needed then one uses NA (NOT NULL)
    params<-c() ## empty vector with which to start paramtrizing
    if ((is.element("A",names(ModelParams))) &&(is.null(EstimationParams$Fixed$A))){
	if (!is.null(EstimationParams$signsA)){ ## The user is allowed to specify signs in A but NO check is done whether now A will remain in the desired matrix class
		ModelParams$A[which(EstimationParams$signsA=="-")]<- log((-1)*ModelParams$A[which(EstimationParams$signsA=="-")])
		ModelParams$A[which(EstimationParams$signsA=="+")]<- log(ModelParams$A[which(EstimationParams$signsA=="+")])
	}
	if ((!is.null(EstimationParams$diagA))&&(EstimationParams$Atype!="SymmetricPositiveDefinite")){	
	    diag(ModelParams$A)=switch(EstimationParams$diagA,
		Positive={log(diag(ModelParams$A)-ifelse(is.null(EstimationParams$minAdiag),0,EstimationParams$minAdiag))},
		Negative={log((-1)*diag(ModelParams$A)-ifelse(is.null(EstimationParams$maxAdiag),0,EstimationParams$maxAdiag))}
	    )
	}
	if (EstimationParams$kY==1){
	    if ((is.null(EstimationParams$signsA))||(EstimationParams$signsA[1,1]=="-")||(EstimationParams$signsA[1,1]=="+")||is.na(EstimationParams$signsA[1,1])){Aparam<-c("A"=ModelParams$A[1,1])}
	    else{Aparam<-c()}
	}
	else{
	    lAparams=switch(EstimationParams$Atype,
		SingleValueDiagonal<-{ModelParams$A[1,1]},
		Diagonal={diag(ModelParams$A)},
		SymmetricPositiveDefinite={.sym.unpar(ModelParams$A)},
		Symmetric={.par.inv.transform.symmetric.matrix(ModelParams$A)},
		TwoByTwo={.par.inv.transform.twobytwo.matrix(ModelParams$A)},
		UpperTri={.par.inv.transform.uppertri.matrix(ModelParams$A)},
		LowerTri={.par.inv.transform.lowertri.matrix(ModelParams$A)},
		DecomposablePositive={.par.inv.transform.decomp.pos.real.matrix(ModelParams$A,ModelParams$eigenSignsA,ModelParams$GivensQCsignsA,1e-15,ifelse(is.null(EstimationParams$minAeigen),0.01,EstimationParams$minAeigen))},
		DecomposableNegative={.par.inv.transform.decomp.neg.real.matrix(ModelParams$A,ModelParams$eigenSignsA,ModelParams$GivensQCsignsA,1e-15,ifelse(is.null(EstimationParams$maxAeigen),0.01,EstimationParams$maxAeigen))},
		DecomposableReal={.par.inv.transform.decomp.real.matrix(ModelParams$A,ModelParams$eigenSignsA,ModelParams$GivensQCsignsA)}, 
		Invertible={.par.inv.transform.invert.matrix(ModelParams$A,1e-15,"qr",NULL,ifelse(is.null(EstimationParams$minRdiag),0.01,EstimationParams$minRdiag))},
		Any={
		    if (!is.null(EstimationParams$signsA)){
			tmpsignsA<-EstimationParams$signsA
			Atmp<-ModelParams$A
			vToNA<-c(which(tmpsignsA=="-"),which(tmpsignsA=="+"))
			if (length(vToNA)>0){tmpsignsA[vToNA]<-NA}
			if (length(which(!is.na(tmpsignsA)))>0){
			    Atmp[which(!is.na(tmpsignsA))]<-NA
			}
			Apars<-c(t(Atmp))
			if (length(which(is.na(Apars)))>0){Apars<-Apars[-which(is.na(Apars))]}
		    }else{Apars<-c(t(ModelParams$A))}
		    Apars
		}
	    )
	    if (is.element(EstimationParams$Atype,c("DecomposablePositive","DecomposableNegative","DecomposableReal","Invertible"))){
		Aparam<-lAparams$vParams
	    }
	    else{Aparam<-lAparams}		    
	    if (length(Aparam)>0){
		if((EstimationParams$Atype=="SingleValue")||(EstimationParams$Atype=="SingleValueDiagonal")||(length(Aparam)==1)){Aparam<-c("A"=Aparam)}
		else{
		    names(Aparam)<-sapply(1:length(Aparam),function(x){paste("A_",x,sep="")})
		    names(Aparam)[1]<-"Astart"
		    names(Aparam)[length(Aparam)]<-"Aend"	 
		}
	    }else{Aparam<-c()}
	}
	params<-c(params,Aparam)
    }
    
    if ((is.element("B",names(ModelParams))) &&(is.null(EstimationParams$Fixed$B))){
	if (!is.null(EstimationParams$signsB)){ ## The user is allowed to specify signs in B but NO check is done whether now A will remain in the desired matrix class
		ModelParams$B[which(EstimationParams$signsB=="-")]<- log((-1)*ModelParams$B[which(EstimationParams$signsB=="-")])
		ModelParams$B[which(EstimationParams$signsB=="+")]<- log(ModelParams$B[which(EstimationParams$signsB=="+")])
	}
	if ((EstimationParams$kX==1)&&(EstimationParams$kY==1)){
	    if ((is.null(EstimationParams$signsB))||(EstimationParams$signsB[1,1]=="-")||(EstimationParams$signsB[1,1]=="+")||is.na(EstimationParams$signsB[1,1])){params<-c(params,"B"=ModelParams$B[1,1])}
	}
	else{	
	    Bparam=switch(EstimationParams$Btype,
		MinusA={c(NA)},
		SingleValue={ModelParams$B[1,1]},
		SingleValueDiagonal<-{ModelParams$B[1,1]},
		Diagonal={diag(ModelParams$B)},
		Symmetric={.sym.unpar(ModelParams$B)},
		Any={
		    if (!is.null(EstimationParams$signsB)){
			tmpsignsB<-EstimationParams$signsB
			Btmp<-ModelParams$B
			vToNA<-c(which(tmpsignsB=="-"),which(tmpsignsB=="+"))
			if (length(vToNA)>0){tmpsignsB[vToNA]<-NA}
			if (length(which(!is.na(tmpsignsB)))>0){
			    Btmp[which(!is.na(tmpsignsB))]<-NA
			}
			Bpars<-c(t(Btmp))
			if (length(which(is.na(Bpars)))>0){Bpars<-Bpars[-which(is.na(Bpars))]}
		    }else{Bpars<-c(t(ModelParams$B))}
		    Bpars
		}
	    )
	    if ((length(Bparam)>0)&&((!is.na(Bparam[1])))){
		if (length(Bparam)>0){
		    if((EstimationParams$Btype=="SingleValue")||(EstimationParams$Btype=="SingleValueDiagonal")||(length(Bparam)==1)){Bparam<-c("B"=Bparam)}
		    else{
			names(Bparam)<-sapply(1:length(Bparam),function(x){paste("B_",x,sep="")})
			names(Bparam)[1]<-"Bstart"
			names(Bparam)[length(Bparam)]<-"Bend"	 
		    }
		}else{Bparam<-c()}
		params<-c(params,Bparam)
	    }
	}
    }
    
    if ((is.element("mPsi",names(ModelParams))) &&(is.null(EstimationParams$Fixed$mPsi))){
    	if (!is.null(EstimationParams$signsmPsi)){ ## The user is allowed to specify signs in A but NO check is done whether now A will remain in the desired matrix class
		ModelParams$mPsi[which(EstimationParams$signsmPsi=="-")]<- log((-1)*ModelParams$mPsi[which(EstimationParams$signsmPsi=="-")])
		ModelParams$mPsi[which(EstimationParams$signsmPsi=="+")]<- log(ModelParams$mPsi[which(EstimationParams$signsmPsi=="+")])
	}
    	if ((EstimationParams$kY==1)&&((EstimationParams$mPsitype=="Global")||((EstimationParams$mPsitype=="Regimes")&&(length(EstimationParams$RegimeTypes)==1)))){
	    if ((is.null(EstimationParams$signsmPsi))||(EstimationParams$signsmPsi[1,1]=="-")||(EstimationParams$signsmPsi[1,1]=="+")||is.na(EstimationParams$signsmPsi[1,1])){Psiparam<-c("Psi"=ModelParams$mPsi[,1])}
	    else{Psiparam<-c()}
	}else{
	    Psiparam=switch(EstimationParams$mPsitype,
		Global={ModelParams$mPsi[,1]},
		Regimes={c(ModelParams$mPsi)}
	    )    
	    names(Psiparam)<-sapply(1:length(Psiparam),function(x){paste("Psi_",x,sep="")})
	    names(Psiparam)[1]<-"Psistart"
	    names(Psiparam)[length(Psiparam)]<-"Psiend"	 
	}
	params<-c(params,Psiparam)
    }
    if ((is.element("mPsi0",names(ModelParams))) &&(is.null(EstimationParams$Fixed$mPsi0))){
    	if (!is.null(EstimationParams$signsmPsi0)){ ## The user is allowed to specify signs in A but NO check is done whether now A will remain in the desired matrix class
		ModelParams$mPsi0[which(EstimationParams$signsmPsi0=="-")]<- log((-1)*ModelParams$mPsi0[which(EstimationParams$signsmPsi0=="-")])
		ModelParams$mPsi0[which(EstimationParams$signsmPsi0=="+")]<- log(ModelParams$mPsi0[which(EstimationParams$signsmPsi0=="+")])
	}    	
    	if (EstimationParams$kY==1){
    	    if ((is.null(EstimationParams$signsmPsi0))||(EstimationParams$signsmPsi0[1,1]=="-")||(EstimationParams$signsmPsi0[1,1]=="+")||is.na(EstimationParams$signsmPsi0[1,1])){Psi0param<-c("Psi0"=ModelParams$mPsi0)}
    	    else{Psi0param<-c()}
    	}else{
	    Psi0param<-ModelParams$mPsi0
	    names(Psi0param)<-sapply(1:length(Psi0param),function(x){paste("Psi0_",x,sep="")})
	    names(Psi0param)[1]<-"Psi0start"
	    names(Psi0param)[length(Psi0param)]<-"Psi0end"	 
	}
	params<-c(params,Psi0param)
    }
    if ((is.element("vY0",names(ModelParams))) &&(is.null(EstimationParams$Fixed$vY0))){
	if (!is.null(EstimationParams$signsvY0)){ ## The user is allowed to specify signs in A but NO check is done whether now A will remain in the desired matrix class
		ModelParams$vY0[which(EstimationParams$signsvY0=="-")]<- log((-1)*ModelParams$vY0[which(EstimationParams$signsvY0=="-")])
		ModelParams$vY0[which(EstimationParams$signsvY0=="+")]<- log(ModelParams$vY0[which(EstimationParams$signsvY0=="+")])
	}
    	if (EstimationParams$kY==1){
    	    if ((is.null(EstimationParams$signsvY0))||(EstimationParams$signsvY0[1,1]=="-")||(EstimationParams$signsvY0[1,1]=="+")||is.na(EstimationParams$signsvY0[1,1])){vY0param<-c("vY0"=ModelParams$vY0)}
    	    else{vY0param<-c()}
    	}else{
    	    vY0param<-ModelParams$vY0
	    names(vY0param)<-sapply(1:length(vY0param),function(x){paste("vY0_",x,sep="")})
	    names(vY0param)[1]<-"vY0start"
	    names(vY0param)[length(vY0param)]<-"vY0end"	 
	}
	params<-c(params,vY0param)
    }
    if ((is.element("vX0",names(ModelParams))) &&(is.null(EstimationParams$Fixed$vX0))){
        if (!is.null(EstimationParams$signsvX0)){ ## The user is allowed to specify signs in A but NO check is done whether now A will remain in the desired matrix class
		ModelParams$vX0[which(EstimationParams$signsvX0=="-")]<- log((-1)*ModelParams$vX0[which(EstimationParams$signsvX0=="-")])
		ModelParams$vX0[which(EstimationParams$signsvX0=="+")]<- log(ModelParams$vX0[which(EstimationParams$signsvX0=="+")])
	}	
    	if (EstimationParams$kX==1){
    	    if ((is.null(EstimationParams$signsvX0))||(EstimationParams$signsvX0[1,1]=="-")||(EstimationParams$signsvX0[1,1]=="+")||is.na(EstimationParams$signsvX0[1,1])){vX0param<-c("vX0"=ModelParams$vX0)}
    	    else{vX0param<-c()}
    	}else{
	    vX0param<-ModelParams$vX0
	    names(vX0param)<-sapply(1:length(vX0param),function(x){paste("vX0_",x,sep="")})
	    names(vX0param)[1]<-"vX0start"
	    names(vX0param)[length(vX0param)]<-"vX0end"	 
	}
	params<-c(params,vX0param)
    }
    
    if ((is.element("Syy",names(ModelParams))) &&(is.null(EstimationParams$Fixed$Syy))){ 
    ## Generally symmetric structures for these matrices should suffice since we do Ct(C) anyway but just in case general 
	if (!is.null(EstimationParams$signsSyy)){ ## The user is allowed to specify signs in A but NO check is done whether now A will remain in the desired matrix class
		ModelParams$Syy[which(EstimationParams$signsSyy=="-")]<- log((-1)*ModelParams$Syy[which(EstimationParams$signsSyy=="-")])
		ModelParams$Syy[which(EstimationParams$signsSyy=="+")]<- log(ModelParams$Syy[which(EstimationParams$signsSyy=="+")])
	}    
        if (!is.null(EstimationParams$diagSyy)){
		diag(ModelParams$Syy)=switch(EstimationParams$diagSyy,
		    Positive={log(diag(ModelParams$Syy))},
		    Negative={log((-1)*diag(ModelParams$Syy))}
		)
	}
	if (EstimationParams$kY==1){
	    if ((is.null(EstimationParams$signsSyy))||(EstimationParams$signsSyy[1,1]=="-")||(EstimationParams$signsSyy[1,1]=="+")||is.na(EstimationParams$signsSyy[1,1])){Syyparam<-c("Syy"=ModelParams$Syy[1,1])}
	    else{Syyparam<-c()}
	}
	else{	
	    Syyparam=switch(EstimationParams$Syytype,
		SingleValueDiagonal={ModelParams$Syy[1,1]},
		Diagonal={diag(ModelParams$Syy)},
		Symmetric={.sym.unpar(ModelParams$Syy)},
		UpperTri={.par.inv.transform.uppertri.matrix(ModelParams$Syy)},
		LowerTri={.par.inv.transform.lowertri.matrix(ModelParams$Syy)},
		Any={
		    if (!is.null(EstimationParams$signsSyy)){
			tmpsignsSyy<-EstimationParams$signsSyy
	 	        Syytmp<-ModelParams$Syy			
			vToNA<-c(which(tmpsignsSyy=="-"),which(tmpsignsSyy=="+"))
			if (length(vToNA)>0){tmpsignsSyy[vToNA]<-NA}
			if (length(which(!is.na(tmpsignsSyy)))>0){
			    Syytmp[which(!is.na(tmpsignsSyy))]<-NA
			}
			Syypars<-c(t(Syytmp))
			if (length(which(is.na(Syypars)))>0){Syypars<-Syypars[-which(is.na(Syypars))]}
		    }else{Syypars<-c(t(ModelParams$Syy))}
		    Syypars
		}
	    )
	    if(length(Syyparam)>0){
		if ((EstimationParams$Syytype!="SingleValueDiagonal")&&(length(Syyparam)>1)){
		    names(Syyparam)<-sapply(1:length(Syyparam),function(x){paste("Syy_",x,sep="")})
		    names(Syyparam)[1]<-"Syystart"
		    names(Syyparam)[length(Syyparam)]<-"Syyend"	 
		}else{Syyparam<-c("Syy"=Syyparam)}
	    }else{Syyparam<-c()}
	}
	params<-c(params,Syyparam)	
    }
    
    if ((is.element("Syx",names(ModelParams))) &&(is.null(EstimationParams$Fixed$Syx))){
    	if ((EstimationParams$kX==1)&&(EstimationParams$kY==1)){ 
    	    if ((is.null(EstimationParams$signsSyx))||(EstimationParams$signsSyx[1,1]=="-")||(EstimationParams$signsSyx[1,1]=="+")||is.na(EstimationParams$signsSyx[1,1]))
    	    {Syxparam<-c("Syx"=ModelParams$Syx[1,1])}
	    else{Syxparam<-c()}
	}else{
	    if (!is.null(EstimationParams$signsSyx)){ ## The user is allowed to specify signs in A but NO check is done whether now A will remain in the desired matrix class
		ModelParams$Syx[which(EstimationParams$signsSyx=="-")]<- log((-1)*ModelParams$Syx[which(EstimationParams$signsSyx=="-")])
		ModelParams$Syx[which(EstimationParams$signsSyx=="+")]<- log(ModelParams$Syx[which(EstimationParams$signsSyx=="+")])
	    }
	    if ((!is.null(EstimationParams$signsSyx))&&(EstimationParams$Syxtype=="Any")){
			tmpsignsSyx<-EstimationParams$signsSyx
    		        Syxtmp<-ModelParams$Syx			
			vToNA<-c(which(tmpsignsSyx=="-"),which(tmpsignsSyx=="+"))
			if (length(vToNA)>0){tmpsignsSyx[vToNA]<-NA}
			if (length(which(!is.na(tmpsignsSyx)))>0){
			    Syxtmp[which(!is.na(tmpsignsSyx))]<-NA
			}
			Syxpars<-c(t(Syxtmp))
			if (length(which(is.na(Syxpars)))>0){Syxpars<-Syxpars[-which(is.na(Syxpars))]}
	    }else{Syxparam<-c(t(ModelParams$Syx))}
	    if (length(Syxparam)>0){
		if (length(Syxparam)==1){Syxparam<-c("Syx"=Syxparam)}
		else{
		    names(Syxparam)<-sapply(1:length(Syxparam),function(x){paste("Syx_",x,sep="")})
		    names(Syxparam)[1]<-"Syxstart"
    		    names(Syxparam)[length(Syxparam)]<-"Syxend"	 
    		}
    	    }else{Syxparam<-c()}
	}
	params<-c(params,Syxparam)	
    }

    if ((is.element("Sxy",names(ModelParams))) &&(is.null(EstimationParams$Fixed$Sxy))){
	if ((EstimationParams$kX==1)&&(EstimationParams$kY==1)){ 
	    if ((is.null(EstimationParams$signsSxy))||(EstimationParams$signsSxy[1,1]=="-")||(EstimationParams$signsSxy[1,1]=="+")||is.na(EstimationParams$signsSxy[1,1]))
	    {Sxyparam<-c("Sxy"=ModelParams$Sxy[1,1])}
	    else{Sxyparam<-c()}
	}else{
	    if (!is.null(EstimationParams$signsSxy)){ ## The user is allowed to specify signs in A but NO check is done whether now A will remain in the desired matrix class
		ModelParams$Sxy[which(EstimationParams$signsSxy=="-")]<- log((-1)*ModelParams$Sxy[which(EstimationParams$signsSxy=="-")])
		ModelParams$Sxy[which(EstimationParams$signsSxy=="+")]<- log(ModelParams$Sxy[which(EstimationParams$signsSxy=="+")])
	    }
	    if ((!is.null(EstimationParams$signsSxy))&&(EstimationParams$Sxytype=="Any")){
			tmpsignsSxy<-EstimationParams$signsSxy
			Sxytmp<-ModelParams$Sxy
			vToNA<-c(which(tmpsignsSxy=="-"),which(tmpsignsSxy=="+"))
			if (length(vToNA)>0){tmpsignsSxy[vToNA]<-NA}
			if (length(which(!is.na(tmpsignsSxy)))>0){
			    Sxytmp[which(!is.na(tmpsignsSxy))]<-NA
			}
			Sxypars<-c(t(Sxytmp))
			if (length(which(is.na(Sxypars)))>0){Sxypars<-Sxypars[-which(is.na(Sxypars))]}
	    }else{Sxyparam<-c(t(ModelParams$Sxy))}

	    if (length(Syxparam)>0){
	    	if (length(Sxyparam)==1){Sxyparam<-c("Sxy"=Sxyparam)}
		else{
		    names(Sxyparam)<-sapply(1:length(Sxyparam),function(x){paste("Sxy_",x,sep="")})
		    names(Sxyparam)[1]<-"Sxystart"
    		    names(Sxyparam)[length(Sxyparam)]<-"Sxyend"	 
    		}
    	    }else{Sxyparam<-c()}
	}
	params<-c(params,Sxyparam)	
    }

    if ((is.element("Sxx",names(ModelParams))) &&(is.null(EstimationParams$Fixed$Sxx))){ 
    ## Generally symmetric structures for these matrices should suffice since we do Ct(C) anyway but just in case general 
	if (EstimationParams$kX==1){
	    if ((is.null(EstimationParams$signsSxx))||(EstimationParams$signsSxx[1,1]=="-")||(EstimationParams$signsSxx[1,1]=="+")||is.na(EstimationParams$signsSxx[1,1])){Sxxparam<-c("Sxx"=ModelParams$Sxx[1,1])}
	    else{Sxxparam<-c()}
	}
	else{
	    if (!is.null(EstimationParams$signsSxx)){ ## The user is allowed to specify signs in A but NO check is done whether now A will remain in the desired matrix class
		ModelParams$Sxx[which(EstimationParams$signsSxx=="-")]<- log((-1)*ModelParams$Sxx[which(EstimationParams$signsSxx=="-")])
		ModelParams$Sxx[which(EstimationParams$signsSxx=="+")]<- log(ModelParams$Sxx[which(EstimationParams$signsSxx=="+")])
	    }
	    Sxxparam=switch(EstimationParams$Sxxtype,
		SingleValueDiagonal={ModelParams$Sxx[1,1]},
		Diagonal={diag(ModelParams$Sxx)},
		Symmetric={.sym.unpar(ModelParams$Sxx)},
		Any={
		    if (!is.null(EstimationParams$signsSxx)){
			tmpsignsSxx<-EstimationParams$signsSxx
			Sxxtmp<-ModelParams$Sxx
			vToNA<-c(which(tmpsignsSxx=="-"),which(tmpsignsSxx=="+"))
			if (length(vToNA)>0){tmpsignsSxx[vToNA]<-NA}
			if (length(which(!is.na(tmpsignsSxx)))>0){
			    Sxxtmp[which(!is.na(tmpsignsSxx))]<-NA
			}
			Sxxpars<-c(t(Sxxtmp))
			if (length(which(is.na(Sxxpars)))>0){Sxxpars<-Sxxpars[-which(is.na(Sxxpars))]}
		    }else{Sxxpars<-c(t(ModelParams$Sxx))}
		    Sxxpars
		}
	    )
	    if ((EstimationParams$Sxxtype!="SigleValueDiagonal")&&(length(Sxxparam)>1)){
		names(Sxxparam)<-sapply(1:length(Sxxparam),function(x){paste("Sxx_",x,sep="")})
		names(Sxxparam)[1]<-"Sxxstart"
		names(Sxxparam)[length(Sxxparam)]<-"Sxxend"	 
	    }else{Sxxparam<-c("Sxx"=Sxxparam)}
	}
	params<-c(params,Sxxparam)	
    }
    params    
}

.mvslouch_to_ouch_model<-function(model_params){
        newA<-rbind(cbind(model_params$A,model_params$B),matrix(0,nrow=ncol(model_params$B),ncol=ncol(model_params$A)+ncol(model_params$B)))
        colnames(newA)<-c(colnames(model_params$A),colnames(model_params$B))
        rownames(newA)<-c(colnames(model_params$A),colnames(model_params$B))
        newvY0<-matrix(c(model_params$vY0,model_params$vX0),ncol=1)
        rownames(newvY0)<-c(rownames(model_params$vY0),rownames(model_params$vX0))
        newSyy<-rbind(cbind(model_params$Syy,model_params$Syx),cbind(model_params$Sxy,model_params$Sxx))
        colnames(newSyy)<-c(colnames(model_params$Syy),colnames(model_params$Sxx))
        rownames(newSyy)<-c(colnames(model_params$Syy),colnames(model_params$Sxx))
        newmPsi<-rbind(model_params$mPsi,matrix(0,ncol=ncol(model_params$mPsi),nrow=nrow(model_params$vX0)))
        colnames(newmPsi)<-colnames(model_params$mPsi)
        rownames(newmPsi)<-c(rownames(model_params$mPsi),model_params$vX0)
        newmPsi0<-matrix(c(model_params$mPsi0,rep(0,nrow(model_params$vX0))),ncol=1)
        rownames(newmPsi0)<-c(rownames(model_params$newmPsi0),rownames(model_params$vX0))
        model_params$A<-newA;
        model_params$mPsi<-newmPsi
        model_params$mPsi0<-newmPsi0
        model_params$vY0<-newvY0
        model_params$Syy<-newSyy
        model_params$B<-NULL;model_params$Sxx<-NULL;model_params$Sxy<-NULL;model_params$Syx<-NULL;model_params$Sxx<-NULL;model_params$vX0<-NULL
        model_params
}
