data1<-read.table('sheep.txt',header=T)
a=fun_2(data1,1000,1000)
# betasigma is the prior variance for the regression parameter beta
# variance_1 is the prior variance for the regression parameter sigma
# You should have rstan package before running this function
fun_2<-function(data1,betasigma=1000,variance_1=1000){
  library(rstan)
  regression<-"
  data{
  int n1;
  int n2;
  int y[n1,n1];
  matrix [n1,n2] cov;
  real betasigma;
  real variance_1;
}
  parameters{
  real <lower=-15,upper=15> d[n1];
  vector[n2] beta;
  real <lower=0> sigma_2;
  } 
  model{
  beta~normal(0,betasigma);
  sigma_2~normal(0,variance_1);
  d~normal(cov*beta,sigma_2);
  for (i in 1:(n1-1)){
  for(j in (i+1):n1){
  y[i,j]~binomial(y[i,j]+y[j,i],1/(1+exp(d[j]-d[i])));
  }
  }
  }
  "
  n1=nrow(data1);n2=ncol(data1)
  cov=scale(data1[,(n1+1):n2],center=T,scale=F)
  data2=data1[,1:n1]
  n2=n2-n1
  data=list(y=data2,cov=cov,n1=n1,n2=n2,betasigma=betasigma,variance_1=variance_1)
  fit <- stan(model_code = regression, model_name = "example11", 
              data = data, iter = 10000, chains = 2, verbose = FALSE) 
  
  combine<-function(string1){
    n=length(string1)
    save=''
    for (i in 1:n) save=paste0(save,string1[i])
    return (save)
  }
  
  
  b=matrix(0,n1,10000)
  newletter=c(letters,LETTERS)
  for (j in 1:2){
    a=fit@sim$samples[[j]]
    for (i in 1:n1){
      if (j==1) {
        b[i,1:5000]=a[[i]][5001:10000]}
      else{
        b[i,5001:10000]=a[[i]][5001:10000]
      }
    }
  }
  save1=NULL
  number=NULL
  for (i in 1:10000){
    index=combine(newletter[sort(b[,i],index.return=TRUE,decreasing=TRUE)$ix])
    if(any(save1==index)==FALSE){
      save1=c(save1,index)
      number=c(number,1)
    }else if (any(save1==index)==TRUE){
      number[which(save1==index)]=number[which(save1==index)]+1
    }
  }
  result=sort(number,decreasing=TRUE,index.return=TRUE)
  prob=result$x[1:8]/10000
  ranking=save1[result$ix[1:8]]
  
  kk=list(model=fit,ranking=ranking,prob=prob)
return(kk)
}

data<-read.table('sheep.txt')
