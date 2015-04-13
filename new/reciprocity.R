reciprocity<-function(M,N){
  Z_statistics=sum(M*N)
  total<-nrow(M)^2
  M1=M;N1=N
  diag(M1)=-1
  diag(N1)=-1
  m<-as.vector(M1);index=which(m==-1);m=m[-index]
  n<-as.vector(N1);n=n[-index]
  R_statistics=0
  for(i in 1:length(m)){
  R_statistics=R_statistics+(1+length(which(m<m[i]))+length(which(m<=m[i])))*(1+length(which(n<n[i]))+length(which(n<=n[i])))/4
  }
  n<-nrow(M)
  K_statistics=0
  for(i in 1:n){
    for(j in 1:(n-1)){
      for(k in (j+1):n){
        K_statistics<-K_statistics+sign((M[i,j]-M[i,k])*(N[i,j]-N[i,k]))
      }
    }
  }

  result=list(Z_statistics=Z_statistics,R_statistics=R_statistics,K_statistics=K_statistics)
return (result)
}

swap<-function(M,i,j){
  result=M
  k=result[i,];result[i,]=result[j,];result[j,]=k
  k=result[,i];result[,i]=result[,j];result[,j]=k
  return(result)
}

test<-function(M,N){
  n=nrow(M)
  emperical<-reciprocity(M,N)
  z=emperical$Z_statistics;r=emperical$R_statistics;k=emperical$K_statistics
  z_vec=r_vec=k_vec=rep(0,n*(n-1)/2)
  n1=1
  for (i in 1:(n-1)){
    for(j in (i+1):n){
      temp=swap(M,i,j)
      simulate=reciprocity(temp,N)
      z_vec[n1]=simulate$Z_statistics
      r_vec[n1]=simulate$R_statistics
      k_vec[n1]=simulate$K_statistics
      n1=n1+1
    }
  } 
n1=n1-1
z_pvalue=length(which(z_vec>=z))/n
r_pvalue=length(which(r_vec>=r))/n
k_pvalue=length(which(k_vec>=k))/n
result=list(z=z,z_pvalue=z_pvalue,r=r,r_pvalue=r_pvalue,k=k,k_pvalue=k_pvalue,k_vec)
return(result)
}