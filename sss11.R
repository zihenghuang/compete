#' Transforms a frequency interaction sociomatrix (valued data) into a dichotomized 1/0 matrix
#'
#' @param m A matrix with individuals ordered identically in rows and columns.
#' @param type Determines the type of dichotomized matrix to be returned.
#'  \strong{\code{type}="wl"} is the default which returns a win-loss matrix
#'  with a '1' representing a consistent winner and a '0' representing a
#'  consistent loser for each dyad of the matrix. A consistent winner is
#'  defined as being the individual in each dyad that has absolutely more
#'  wins than defeats.   In the default condition if competitors have the
#'  same number of wins each, they both receive a 0.
#'  If \strong{\code{type}="wlties"} the default dichotomized win-loss
#'  matrix will be returned but it will also return 0.5 into cells for tied
#'  relationships.
#'  If \strong{\code{type}="wlties0"} the default dichotomized win-loss
#'  matrix will be returned but it will also return 0.5 into cells for tied
#'  relationships. Additionally, if two competitors never interacted both
#'  cells for that relationship will be returned with a 0.
#'  If \strong{\code{type}="wlbinom"} every relationship within the win-loss
#'  matrix is assessed for whether one competitor significantly wins more
#'  competitive interactions than the other competitor.  Significance is
#'  calculated using a binomial test with probability of p=0.05. A '1' is
#'  given to significant winners within a relationship and a '0' is given
#'  to significant losers or if neither individual is a winner.
#'  If \strong{\code{type}="wlbinomties"} The same procedure is done as for
#'  \strong{\code{type}="wlbinom"}, but if no signficiant winner/loser can
#'  be determined then a 0.5 is returned rather than a 0.
#'  If \strong{\code{type}="pa"} the inputted matrix will be turned into a
#'  dichotomized presence-absence matrix, with a '1' indicating that the
#'  competitor in a the row of the matrix beat the competitor in the column
#'  at least once. A '0' indicates that that competitor never beat the
#'  other competitor.
#'  If \strong{\code{type}="dom"} the inputted matrix will be turned into a
#'  dominance score matrix, with a '1' indicating that the
#'  competitor in a the row of the matrix dominates the competitor in the
#'  column. A '-1' indicates that that competitor in a row is subordinate
#'  to the competitor in the column. A '0.5' indicates a tie.  A '0'
#'  indicates an observational or structural zero.
#' @return A dichotomized win/loss or presence/absence matrix.
#' @examples
#' get_di_matrix(bonobos)
#' get_di_matrix(mouse)
#' @section References:
#' Appleby, Shizuka website
#' @section Further details:
#' fill this information in when have time.
#' @export
I_IS_20133<-function(Mk,p=c(.3,.5,.1,.1),a_max=50,nTries=30){
    original=Mk
    Mk=transform(Mk)
    attempt=1
    n<-ncol(Mk)
    index=nature(Mk)
     best=index
    M=matrix_change(Mk,index)
    initial=M
    index_initial=index
    l=fun(M)
    Imin=l[1];SImin=l[2]
    matrix1 = (M-t(M))/2
    while (attempt<a_max){
        t=1
        random_p = rmultinom(1,1,prob=p)
        random_p=which(random_p==1)
        stopIteration1=FALSE;stopIteration2=FALSE
        while (stopIteration1==FALSE){
            judge=0
            while (stopIteration2==FALSE){
                stopIteration2<-TRUE
                ####strategy 1
                if (attempt %% 2==1){
                    for (i in 1:(n-1)){
                        for (j in (i+1):n){
                            if (matrix1[i,j]<0){
                                sum1=sum(matrix1[i,(i+1):j])
                                if (sum1<0){
                                    matrix1=swap(matrix1,i,j)
                                    change=index[i];index[i]=index[j];index[j]=change
                                    stopIteration2=FALSE
                                }
                            }
                        }
                    }
                }

                #####strategy2
                else if (random_p==1){
                    for (i in 1:(n-1)){
                        for (j in (i+1):n){
                            sum1=sum(matrix1[i,(i+1):j])
                            if (sum1<0){
                                matrix1=swap(matrix1,i,j)
                                change=index[i];index[i]=index[j];index[j]=change
                                stopIteration2=FALSE
                            }
                        }
                    }
                }
                #####strategy4
               
           else if (random_p==2){
                   if (judge==0){
                        for (i in 1:n){
                        record=rep(0,n)
                        for (j in 1:n){
                            if (i==j){
                                record[j]=0
                            }
                            else if(j>i){
                                record[j]=sum(matrix1[i,(i+1):j])
                            }
                            else if(j<i){
                                record[j]=-sum(matrix1[i,j:(i-1)])
                            }
                        }
                        if (min(record)<0){
                           
                            min1<-sample(which(record==min(record)),1)
                            matrix1<-shift(matrix1,i,min1)
                            index<-shift_index(index,i,min1)
                            stopIteration2=FALSE
                        }
                    }
                    }
                    if (stopIteration2==TRUE){  
                    judge=1
                    record=matrix(0,2,n);lll=fun(matrix1)
                    for (i in 1:n){
                        for (j in 1:n){
                            if (i==j){
                                record[1,j]=1
                            }
                            else if(j>i){
                                record[1,j]=sum(matrix1[i,(i+1):j])
                            }
                            else if(j<i){
                                record[1,j]=-sum(matrix1[i,j:(i-1)])
                            }  
                            record[2,j]=delta_si(matrix1,i,j)[2]-lll[2]
                        }
                        a1=record[1,];a2=record[2,]
                        if (min(a1)==0){
                            a2[-which(a1==0)]=0
                            if (min(a2)<0){
                                final=sample(which(a2==min(a2)),1)
                                matrix1<-shift(matrix1,i,final)
                                index<-shift_index(index,i,final);lll=fun(matrix1)
                            stopIteration2=FALSE
                            }
                        }
                  }
                    }
                
            }
   #####strategy5
               
           else if (random_p==3){
                     if (judge==0){
                     count=0;pair=NULL
                    for (i in 1:n){
                        for (j in 1:n){
                            
                             if(j>i){
                                
                                if (sum(matrix1[i,(i+1):j])<0) {
                                      pair=rbind(pair,c(i,j));count=count+1
                          }    
                                     }
                             else if(j<i){
                              if (sum(matrix1[i,j:(i-1)])>0) {
                              pair=rbind(pair,c(i,j));count=count+1
                                                             }                            
                                          }
                                       }
                                   }
            if (count>0){                           
                            rand<-sample(count,1)                          
                            matrix1<-shift(matrix1,pair[rand,1],pair[rand,2])
                            index<-shift_index(index,pair[rand,1],pair[rand,2])
                            stopIteration2=FALSE
                        }
                        }
                    if (stopIteration2==TRUE){  
            judge=1      
            pair=NULL;count=0;lll=fun(matrix1)          
              for (i in 1:n){        
                        for (j in 1:n){
                            if (i==j){
                                delta_i=0
                            }
                            else if(j>i){
                                delta_i=sum(matrix1[i,(i+1):j])
                            }
                            else {
                                delta_i=-sum(matrix1[i,j:(i-1)])
                            }          
             if((delta_si(matrix1,i,j)[2]-lll[2]<0)&(delta_i==0)&(i!=j)){
                        pair=rbind(pair,c(i,j));count=count+1}
                        }
                    }
                        if (count>0){
                            rand=sample(count,1)
                                matrix1<-shift(matrix1,pair[rand,1],pair[rand,2])
                            index<-shift_index(index,pair[rand,1],pair[rand,2])
                            lll=fun(matrix1)
                            stopIteration2=FALSE
                            
                        }
                  
                    }
                }
   #####strategy3
               
           else if (random_p==4){
                   if (judge==0){
                    for (i in 1:n){
                        for (j in 1:n){
                            
                          if(j>i){
                                
                                if (sum(matrix1[i,(i+1):j])<0) {
matrix1<-shift(matrix1,i,j)
                            index<-shift_index(index,i,j)
                            stopIteration2=FALSE

}    
                            }
                            else if(j<i){
   if (sum(matrix1[i,j:(i-1)])>0) {
                             matrix1<-shift(matrix1,i,j)
                            index<-shift_index(index,i,j)
                            stopIteration2=FALSE

}                            
}
                        }
                         }
                           }
                    if (stopIteration2==TRUE){ 
                    judge=1
                    lll=fun(matrix1) 
                    for (i in 1:n){
                        for (j in 1:n){
                          delta_i=0
                            if(j>i){
                                delta_i=sum(matrix1[i,(i+1):j])
                            }
                            else if(j<i){
                                delta_i=-sum(matrix1[i,j:(i-1)])
                            }
                                  if((delta_si(matrix1,i,j)[2]-lll[2]<0)&(delta_i==0)&(i!=j)){
                         matrix1<-shift(matrix1,i,j)
                                index<-shift_index(index,i,j)
lll=fun(matrix1)                            
stopIteration2=FALSE}
}
}
                    }
                   }     
                    
                




}
            l<-fun(matrix1)
            I=l[1];SI=l[2]
            if ((I<Imin)|((I==Imin)&(SI<SImin))){
                best=list()
                best[[1]]=index
                Imin=I
                SImin=SI
                stopIteration1=FALSE
            }
            else if((I==Imin)&(SI==SImin)){
                t=t+1
                best[[length(best)+1]]=index
            }
            else{
                t=t+1
            }
            if ((SImin>0)&(t<nTries)){
                for (j in 2:n){
                kk=matrix1[j,1:(j-1)]
                kk=(kk+1)/2
                if (sum(kk)>0){
                    rand=sample((j-1),1)
                    matrix1=swap(matrix1,rand,j)
                    change=index[rand];index[rand]=index[j];index[j]=change
                    stopIteration1=FALSE
                               }
                               }
            }
            else{
                stopIteration1=TRUE
            }
        }
        attempt=attempt+1
    }
    best=unique(best)
    correlation=rep(0,length(best))
    for (k in 1:length(best)){
        correlation[k]=cor.test(c(1:n),best[[k]])$p.value
    }
        num<-which(correlation==min(correlation))
        result_1=list()
        result_2=list()
        for (i in 1:length(num)){
            result_1[[i]]=matrix_change(original,best[[num[i]]])
            result_2[[i]]=best[[num[i]]]
        }
    answer=list(best_matrix=result_1,best_order=result_2,I=Imin,SI=SImin)
    return(answer)
}


####
#THE FOLLOWING ARE ALL REQUIRED FUNCTION FOR 'fun_2013',please load them first
###


##
#This function is used to calculate the I and SI given the dominance matrix
##
fun<-function(M){
    n=ncol(M)
    result=rep(0,2)
    k=M-t(M)
    k[upper.tri(k)]=0
    result[1] = sum(k>0)
    a=length(which(k>0))
    y=which(k>0)%%n
    x=(which(k>0)-1)%/%n+1
    y[y==0]=y[y==0]+n
    if (a>0){
        result[1]=a
        result[2]=sum(y-x)
    }
return(result)
}

###
#This function is to modify the dominance matrix when ith individual and
#jth individual exchange their positions.
###
swap<-function(M,i,j){
    result=M
    k=result[i,];result[i,]=result[j,];result[j,]=k
    k=result[,i];result[,i]=result[,j];result[,j]=k
    return(result)
}

###
#'nature' is used to offer the initial sequence depends on the porpotion of dominance
###
nature<-function(M){
    n=ncol(M)
    index=c(1:n)
    D<-rep(0,n)
    S<-rep(0,n)
    for (i in 1:n){
        fast<-M-t(M)
        D[i]=D[i]+length(which(sign(fast[i,])==1))
        S[i]=S[i]+length(which(sign(fast[i,])==-1))
    }
    result=D/(D+S)
    save=sort(result,decreasing=TRUE,index.return=TRUE)
    sequence=save$ix
    result=save$x
    Dom_sub=D-S
    for (i in 1:(n-1)){
        if (result[i]==result[i+1]){
            if (Dom_sub[i]<Dom_sub[i+1]){
                mid<-sequence[i];sequence[i]=sequence[i+1];sequence[i+1]=mid
            }
        }
    }
    return(sequence)
}

####
#'matrix_change' is used to modify the dominance matrix according to the give sequence
###
matrix_change<-function(M,sequence){
    n=ncol(M)
    new_matrix<-NULL
    for(i in 1:n){
        temp<-M[sequence[i],]
        temp<-temp[sequence]
        new_matrix<-rbind(new_matrix,temp)
    }
    return(new_matrix)
}

####
#'transform' function is the function to transform the original matrix to the fictitious matrix
#M is the original matrix
####
transform<-function(M){
    n=ncol(M)
    for (i in 1:(n-1)){
        for (j in (i+1):n){
            if (M[j,i]>M[i,j]){
                M[j,i]=1;M[i,j]=-1
            }
            else if((M[j,i]>0)&(M[i,j]==M[j,i])){
                M[i,j]=M[j,i]=.5
            }
            else if(M[j,i]<M[i,j]){
            M[j,i]=-1;M[i,j]=1
            }
        }
    }
return(M)
}

####
#'shift' function is used to put the rank ith individual at the jth position
#M is the matrix
####
shift<-function(M,i,j){
    n=ncol(M)
    kk=M[,-i]
    if (j==1){
        left<-rep(0,n)
    }
    else{
        left<-cbind(rep(0,n),kk[,1:(j-1)])
    }
    if (j==n){
        right<-rep(0,n)}else{
        right<-cbind(kk[,j:(n-1)],rep(0,n))
    }
    matrix1=cbind(left,M[,i],right)
    matrix1=matrix1[,-c(1,n+2)]
    kk=matrix1[-i,]
    if (j==1){
        top<-rep(0,n)}else{
        top<-rbind(rep(0,n),kk[1:(j-1),])
    }
    if (j==n){
        down<-rep(0,n)
    }
    else{
        down<-rbind(kk[j:(n-1),],rep(0,n))
    }
    matrix2<-rbind(top,matrix1[i,],down)
    matrix2<-matrix2[-c(1,n+2),]
    return(matrix2)
}

####
#'shift_index' is a function change the rank
####
shift_index<-function(index,i,j){
    n<-length(index)
    memory<-index[-i]
    if (j==1){
        left=NULL
    }
    else{
        left<-memory[1:(j-1)]
    }
    if(j==n){
        right=NULL}
    else{
        right<-memory[j:(n-1)]
    }
    new<-c(left,index[i],right)
    return(new)
}
#### calculate the delta_si#####

delta_si<-function(m,i,j){
a=fun(shift(m,i,j))
return(a)}

