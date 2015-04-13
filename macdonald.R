#calculate the liearity by different methods  
calculate_linearity=function(M,method){
    M1=as.matrix(M)
    M=as.matrix(M)
    M1=M-t(M)
    n=nrow(M)
    matrix=M-t(M)
    matrix[which(matrix<=0)]=0
    M=sign(matrix)
    if (method=='Kendall'){
      d=n*(n-1)*(2*n-1)/12-0.5*sum(rowSums(M)^2)
      if (n %%2 ==0) dmax=(n^3-4*n)/24
      if (n%%2==1) dmax=(n^3-n)/24
      k=1-d/dmax
    }
else if(method=='Landau'){
  k=12/(n^3-n)*sum((rowSums(M)-(n-1)/2)^2)
}
else if(method=='modified-Landau'){
  result=rep(0,1000)
  u=M1[upper.tri(M1)]
  len=length(which(u==0))
  random=matrix(rbinom(len*1000,1,.5),1000,len)
  M2=M1
  M3=matrix(0,n,n)
  for (i in 1:1000){
    u=M1[upper.tri(M1)]
    u[u==0]=random[i,]
    M2[upper.tri(M2)]=u
    M2[lower.tri(M2)]=0
    M2[M2<0]=0
    M2=sign(M2)
    u[u<0]=0;u1=1-sign(u)
    M3[upper.tri(M3)]=u1
    M2=M2+t(M3)
    diag(M2)=0
    result[i]=12/(n^3-n)*sum((rowSums(M2)-(n-1)/2)^2)
  }
  k=mean(result)
}
return (k)
  }
# Calculate the transitive matrix proportion
  Ttri=function(M){
    n=nrow(M)
    M=as.matrix(M)
    matrix=M-t(M)
    matrix[matrix<0]=0
    M=sign(matrix)
    cycle=0;transitive=0;
    for (i in 1:(n-2)){
      for (j in (i+1):(n-1)){
        for (k in (j+1):n){
          s=M[i,j]+M[i,k]+M[j,k]+M[j,i]+M[k,i]+M[k,j]
          if (s==3){
            if ((M[i,j]+M[i,k]==1)&(M[j,k]+M[j,i]==1)&(M[k,i]+M[k,j])){
              cycle=cycle+1
            }else{
              transitive=transitive+1
            }
          }
        }
      }
    }
  Ttri=4*(transitive/(cycle+transitive)-.75)
  return(Ttri)
  }
  
 #calculate the pvalue of the matrix for simulation 
  tritest=function(M){
    M=as.matrix(M)
    tri_em=Ttri(M)
    n=nrow(M)
    M=M-t(M)
    total=length(which(M[upper.tri(M)]!=0))
    random=matrix(sample(c(1,-1),total*1000,replace=T),1000,total)
    n1=n*(n-1)/2
    result=rep(0,1000)
    for (i in 1:1000){
      simulation=matrix(0,n,n)
      u=simulation[upper.tri(simulation)]
    u[sort(sample(c(1:n1),total,replace=F))]=random[i,]
    simulation[upper.tri(simulation)]=u
    simulation=simulation-t(simulation)
    simulation[simulation<0]=0
    result[i]=Ttri(simulation)
    }
    pvalue=length(which(result>tri_em))/1000
    list1=list(pvalue=pvalue,tri=tri_em,emdistribution=result)
    return(list1)
  }