## This file is part of mvSLOUCH

## This software comes AS IS in the hope that it will be useful WITHOUT ANY WARRANTY, 
## NOT even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
## Please understand that there may still be bugs and errors. Use it at your own risk. 
## We take no responsibility for any errors or omissions in this package or for any misfortune 
## that may befall you or others as a result of its use. Please send comments and report 
## bugs to Krzysztof Bartoszek at krzbar@protonmail.ch .

## Function from the slouch package. Author : Jason Pienaar

`.make.states.mvsl` <-
function(pars.object){
x<-pars.object$Final.states
n<-length(x)
n.states<-length(pars.object$niche.code[,1])

#which nodes are ambiguous

count<-0

ambig<-NA
for(i in 1:n)
 {
  if(length(x[[i]])>=2)
   {
    count=count+1;
    ambig<-c(ambig,i)	
   }	
 }
ambig=ambig[-1]
n.ambig<-length(ambig)

# choose first character state

tmp<-matrix(data=NA, nrow=n, ncol=1) 
for(i in 1:n)
{
 tmp[i,1]<-x[[i]][1]
  for(j in 1:n.states)
  {
  if (tmp[i,1]==j) tmp[i,1]<-pars.object$niche.code[,1][[j]] 
  }
 }
 
# encode ambiguous states as character ambiguous

if(n.ambig!=0)
 {
  for(i in 1:n)
   {
    for(j in 1:n.ambig)
     {
      if(i==ambig[j]) tmp[i,1]="ambiguous"	
     }
   }
 }

return(as.factor(tmp))
}

