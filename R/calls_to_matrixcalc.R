## This file is included in the mvSLOUCH R software package.

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this package or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .
## =========================================================================================================


.matrixcalc_is.diagonal.matrix <- function( x, tol=1e-8 )
{
## called in evolmodelest.R
    matrixcalc::is.diagonal.matrix( x, tol )
}

.matrixcalc_is.positive.semi.definite <- function( x, tol=1e-8 )
{
## called in PhyloSDE.R
    matrixcalc::is.positive.semi.definite( x, tol )
}

.matrixcalc_is.symmetric.matrix <- function( x )
{
## called in PhyloSDE.R, estimBM.R, evolmodelest.R, getESS.R, loglik.R, sdecovariancephyl.R, sdemoments.R
    matrixcalc::is.symmetric.matrix( x )
}

.matrixcalc_is.positive.definite <- function( x, tol=1e-8 )
{
## called in estimBM.R, evolmodelest.R, getESS.R, loglik.R, sdecovariancephyl.R, sdemoments.R
    matrixcalc::is.positive.definite( x, tol)
}


.matrixcalc_is.singular.matrix <- function( x, tol=1e-8 )
{
    matrixcalc::is.singular.matrix( x, tol)
}



