d 9 12 6 27
B 10 d 9 12 12
C 2 5 d 0 2
D 3 3 0 d 2
E 2 3 0 4 0


data1=matrix(c(0,10,2,3,2,9,0,5,3,3,12,9,0,0,0,6,12,0,0,4,27,12,2,2,0),5,5)
fun_1(data1,1000,3,1)

# para is the variance prior for the random position d
# choose is choosing which one as the focal to let it's d value be 0 as standard
# dist can be 1 or 2
# 1 is normal distribution and 2 is student t distribution
# If you use t-distribution,the degree of freedom will be n-2, n is the total amount of species 
# Install the rstan package from stan website before running this function
fun_1<-function(data1,para=1000,choose,dist=1,simulation=1000){
  library(rstan)
  bayesiani_si<-"
  data{
  int n;
  int y[n,n];
  real sigma1;
  int focal;
  int dist;
  }
  
  parameters{
  real <lower=-15,upper=15> d[n];
  } 
  
  transformed parameters{
  real <lower=-15,upper=15> d1[n];
  d1<-d;
  d1[focal]<-0;
  }
  model{
  if (dist==2){
  for  (i in 1:n)
{d[i]~student_t(n-2,0,sigma1);}
  }else{
  for  (i in 1:n)
{d[i]~normal(0,sigma1);}
  }
  for (i in 1:(n-1)){
  for(j in (i+1):n){
  y[i,j]~binomial(y[i,j]+y[j,i],1/(1+exp(d1[j]-d1[i])));
  }
  }
  }
  "
  n=nrow(data1)
  data=list(y=data1,n=n,focal=choose,dist=dist,sigma1=para)
  fit <- stan(model_code = bayesiani_si, model_name = "example11", 
              data = data, iter = simulation, chains = 2, verbose = FALSE) 
  
  combine<-function(string1){
    n=length(string1)
    save=''
    for (i in 1:n) save=paste0(save,string1[i])
    return (save)
  }
  getinverse<-function(a,b){
    n=length(a)
    a_1<-sort(a,decreasing=TRUE,index.return=TRUE)
    after<-b[a_1$ix]
    sum1=0
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        if (after[j]>after[i]) sum1=sum1+1
      }
    }
  return (sum1)
  }
  
  getarea<-function(a,b){
    n=length(a)
    a=a+rnorm(n,0,0.1)
    b=b+rnorm(n,0,0.1)
    result<-rep(0,100)
    area<-0
    for (i in 1:100){
      a1=quantile(b,0.01*i)
      area=area+length(which(a<=a1))/n
      result[i]=length(which(a<=a1))/n
    }
  answer=list(area=area/100,result=result)
  }
  
  
  
  b=matrix(0,n,simulation)
  newletter=c(letters,LETTERS)
  for (j in 1:2){
    a=fit@sim$samples[[j]]
    for (i in 1:n){
      if (j==1) {
        b[i,1:(simulation/2)]=a[[n+i]][(simulation/2+1):simulation]
      }else if(j==2){
        b[i,(simulation/2+1):simulation]=a[[n+i]][(simulation/2+1):simulation]
      }
    }
}

save1=NULL
number=NULL
times=simulation
empirical<-rep(0,times)
for (i in 1:times){
  for(j in 1:times){
    empirical[i]=empirical[i]+getinverse(b[,i],b[,j])
  }
}
newmatrix<-matrix(0,n,n)
mm<-rowMeans(b)
ntime<-data1+t(data1)
for(i in 1:n){
  for(j in 1:n){
    if(i!=j){
    if(mm[i]>mm[j]){
    newmatrix[i,j]=floor(1/(1+exp(mm[j]-mm[i]))*ntime[i,j])+1
    }else{
      newmatrix[i,j]=floor(1/(1+exp(mm[j]-mm[i]))*ntime[i,j])  
    }
    }
  }
}
data=list(y=newmatrix,n=n,focal=choose,dist=dist,sigma1=para)
fit1 <- stan(model_code = bayesiani_si, model_name = "example11", 
            data = data, iter = simulation, chains = 2, verbose = FALSE) 
b1=matrix(0,n,simulation)
for (j in 1:2){
  a=fit1@sim$samples[[j]]
  for (i in 1:n){
    if (j==1) {
      b1[i,1:(simulation/2)]=a[[n+i]][(simulation/2+1):simulation]
    }else if(j==2){
      b1[i,(simulation/2+1):simulation]=a[[n+i]][(simulation/2+1):simulation]
    }
  }
}

nullhypo<-rep(0,times)
for (i in 1:times){
  for(j in 1:times){
    nullhypo[i]=nullhypo[i]+getinverse(b1[,i],b1[,j])
  }
}

area1=getarea(nullhypo,empirical)
##do the simulation area
fit2 <- stan(model_code = bayesiani_si, model_name = "example11", 
             data = data, iter = simulation*2, chains = 2, verbose = FALSE) 
b2=matrix(0,n,simulation*2)
for (j in 1:2){
  a=fit2@sim$samples[[j]]
  for (i in 1:n){
    if (j==1) {
      b2[i,1:simulation]=a[[n+i]][(simulation+1):simulation]
    }else if(j==2){
      b2[i,(simulation+1):(2*simulation)]=a[[n+i]][(simulation+1):(2*simulation)]
    }
  }
}
em_area<-rep(0,5000)
reference<-rep(0,2*simulation)
for (i in 1:(2*times)){
  for(j in 1:(2*times)){
    reference[i]=reference[i]+getinverse(b2[,i],b2[,j])
  }
}

for (i in 1:5000){
  matrix1<-matrix(sample(2*times,2*times),2,times)
  ssss=getarea(reference[matrix1[1,]],reference[matrix1[2,]])
  em_area[i]=ssss$area
}

area_pvalue<-length(which(em_area>area1$area))/5000

for (i in 1:simulation){
  index=combine(newletter[sort(b[,i],index.return=TRUE,decreasing=TRUE)$ix])
  if(any(save1==index)==FALSE){
    save1=c(save1,index)
    number=c(number,1)
  }else if (any(save1==index)==TRUE){
    number[which(save1==index)]=number[which(save1==index)]+1
  }
}
result=sort(number,decreasing=TRUE,index.return=TRUE)
prob=result$x[1:8]/simulation
ranking=save1[result$ix[1:8]]
kk=list(model=fit,ranking=ranking,prob=prob,lineartest_pvalue=area_pvalue)
return (kk)
  }

