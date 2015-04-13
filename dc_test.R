#This function is used to compute the DC by formula (H-L)/(H+L)
#input is a matrix
#output is a list, DC is the DC
#phi is the skew-symmetric index(0 means completely symmetric, 0.5 means completely not symmetric)
#si is the symmetric index
dc_compute=function(Matrix){
  Matrix=as.matrix(Matrix)
  diag(Matrix)=0
  N=sum(Matrix)/2
  DC=sum(abs(Matrix-t(Matrix)))/2/sum(Matrix)
  S=(Matrix+t(Matrix))/2
  K=(Matrix-t(Matrix))/2
  phi=sum(diag(t(K)%*%K))/sum(diag(t(Matrix)%*%Matrix))
  #phi is the skew-symmetrical index
  si=1-phi
  #si is the symmetrical index
  result=list(DC=DC,S=S,K=K,phi=phi,si=si)
  return(result)
}

#This function is used to compute the p value of the symmetric test for a give matrix
#input is a relation matrix, p is the probability matrix used to do simulation
#N is the total results between each dyad
#ntimes is how many time we should simulate to get the empirical p.value
#output is the p value of the skew-symmetric index and DC in their
#empirical quantile
#PS  this is the one side test, that means is DC and the skew-symmetric index are too large
#the p value will be small and rejects the null hypothesis(This matrix is symmetric)
dc_test=function(matrix,p,N=20,ntimes=10000)
{
  matrix=as.matrix(matrix)
  n=nrow(matrix)
  phi=matrix(0,ntimes)
  DC=matrix(0,ntimes)
  result=rep(0,n^2*ntimes)
  dim(result)=c(n,n,ntimes)
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      number=rbinom(ntimes,N,p[i,j])
      for (k in 1:ntimes){
        result[i,j,k]=number[k]
        result[j,i,k]=N-number[k]
      }
    }
  }
for (i in 1:ntimes){
  a=dc_compute(result[,,k])
  phi[i]=a$phi;DC[i]=a$DC
}
a=dc_compute(matrix);phi_0=a$phi;DC_0<-a$DC
phi1=sort(c(phi,phi_0));DC1=sort(c(DC,DC_0))
phi.pvalue=min(1-which(phi1==phi_0)/(ntimes+1),1-which(phi1==phi_0)/(ntimes+1))
DC.pvalue=min(1-which(DC1==DC_0)/(ntimes+1),1-which(DC1==DC_0)/(ntimes+1))
if (phi.pvalue==0) phi.pvalue=1/ntimes
if (DC.pvalue==0) DC.pvalue=1/ntimes
mean_phi=mean(phi);mean_DC=mean(DC)
variance_phi=var(phi);variance_DC=var(DC)
p.value=list(DC.pvalue=DC.pvalue,phi.pvalue=phi.pvalue,mean_phi=mean_phi,mean_DC=mean_DC
            , variance_phi=variance_phi,variance_DC=variance_DC)
return(p.value)
}

