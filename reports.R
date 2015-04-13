# library(devtools)
# install_github('jalapic/engsoccerdata', username = "jalapic")
library(engsoccerdata)

data(package="engsoccerdata")
data = spainliga
attach(spainliga)
season_list = unique(Season)
n = length(season_list)
r = matrix(NA,6,n)
rownames(r) = c("year","I", "SI","t_r_hat","t_o_hat","dc_p_value")
r[1,] = season_list
print (n)
for (i in 84:n){
    print (i)
    sn = r[1,i]
    temp_data = data[data$Season == sn,]
    temp_data1 = temp_data[,c(3,4,7,8)]
    result = temp_data1[,3] - temp_data1[,4]
    result[result>0] = "W"
    result[result<0] = "L"
    result[result==0] = "T"
    temp_data1 = cbind(temp_data1,result)
    temp_data = temp_data1[,c(1,2,5)]
    temp_data2 = get_wl_dataframe(temp_data)
    temp_data2 = get_wl_matrix(temp_data2)
    isi = isi98(temp_data2,100)
    pp = matrix(0.5,nrow(temp_data2),nrow(temp_data2))
    dc = dc_test(temp_data2, pp)
    print("dc test done")
    r[6,i] = dc$DC.pvalue
    print("isi done")
    r[2,i] = isi$I
    r[3,i] = isi$SI
    tsi = tsi(temp_data2,3,datatype = "matrix")
    print("tsi done")
    r[4,i] = tsi$t_r_hat
    r[5,i] = tsi$t_o_hat


}

par(mfrow=c(2,2))
plot(season_list,r[3,],type = "l", xlab = "year", col = "blue", ylab = "I&SI")
points(season_list,r[2,],type = "l", xlab = "year", col = "red")
plot(season_list,r[4,],type = "l", xlab = "year", col = "blue", ylab = "TSI")
points(season_list,r[5,],type = "l", xlab = "year", col = "red")
plot(season_list,r[6,],type = "l", xlab = "year", col = "blue", ylab = "DC")
