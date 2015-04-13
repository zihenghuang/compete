data1=matrix(c(0,10,2,3,2,9,0,5,3,3,12,9,0,0,0,6,12,0,0,4,27,12,2,2,0),5,5)
fun_1(data1,1000,3,1)

# para is the variance prior for the random position d
# choose is choosing which one as the focal to let it's d value be 0 as standard
# dist can be 1 or 2
# 1 is normal distribution and 2 is student t distribution
# If you use t-distribution,the degree of freedom will be n-2, n is the total amount of species 
# Install the rstan package from stan website before running this function
fun_1<-function(data1,para=1000,choose,dist=1){
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
            data = data, iter = 10000, chains = 2, verbose = FALSE) 

combine<-function(string1){
  n=length(string1)
  save=''
  for (i in 1:n) save=paste0(save,string1[i])
  return (save)
}


b=matrix(0,n,10000)
newletter=c(letters,LETTERS)
for (j in 1:2){
a=fit@sim$samples[[j]]
  for (i in 1:n){
if (j==1) {
  b[i,1:5000]=a[[n+i]][5001:10000]}
else{
  b[i,5001:10000]=a[[n+i]][5001:10000]
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
kk=list(model=fit,ranking=ranking,prob=prob,mean=(attr(fit@sim$samples[[1]],"mean_pars")+attr(fit@sim$samples[[2]],"mean_pars"))/2)
return (kk)
}
