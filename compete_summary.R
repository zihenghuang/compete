cor_test <- function(M){
    n = nrow(M)
    d_score = ds(M)
    d_score_p = d_score$rank_p
    d_score_d = d_score$rank_d
    print("Dvide Score done!")
    isi = isi98(M,1000)
    isi_r = isi$best_order[[1]]
    print("I&SI 1998 done!")
    bay = bay_isi(M,1000,3,1)
    bay_ranking = bay$mean[(n+1):(2*n)]
    bay_r = order(bay_ranking)[n:1]
    print("Bayesian methods done!")
    total = matrix(c(d_score_p,d_score_d,isi_r,bay_r),nrow=4,ncol=n,byrow = TRUE)
    cor_matrix = matrix(1,4,4)
    cor.test(isi_r,bay_r,method = "spearman")
    for(i in 1:4){
        for(j in 1:4){
            cor_matrix[i,j] = cor.test(total[i,],total[j,],method = "spearman")$estimate
        }
    }
    name = c("DS_p","DS_d","I$SI","Bay")
    colnames(cor_matrix) = name
    rownames(cor_matrix) = name
    return(list(cor_matrix = cor_matrix, ds_p = d_score_p, ds_d = d_score_d, isi_r = isi_r, bay_r = bay_r))
}
