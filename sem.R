#SEM
rm(list = ls(all.names = TRUE))
setwd("F:/R_code")
library(piecewiseSEM)
library(bestNormalize)
library(lme4)
library(mgcv) 


metanew <- read.csv("metamean.csv", row.names=1)
data=metanew



## 正态性检验
# 假设您的数据位于一个名为"data"的数据框中

# 对每一列进行Shapiro-Wilk正态性检验
normality_tests <- lapply(data[,7:50], function(x) {
  if(is.numeric(x)) {
    shapiro_test <- shapiro.test(x)
    return(shapiro_test)
  } else {
    return(NULL)
  }
})

# 打印检验结果
for (i in seq_along(normality_tests)) {
  if (!is.null(normality_tests[[i]])) {
    col_name <- names(data[,7:50])[i]
    print(paste("Shapiro-Wilk test for column", col_name))
    print(normality_tests[[i]])
  }
}
###p<0.05不符合正态性
##正态性变换，对每一列进行

BN_obj<- lapply(data[, 7:50], function(x) bestNormalize(x, allow_orderNorm = FALSE, out_of_sample = FALSE))
BN_obj####自动选择正态化方法根据pearson/df值，越小越好，选择相应的方法
# 提取每个变量的x.t值
x.t_list <- lapply(BN_obj, function(var) var$x.t)
#hist(x.t_list)
# 将提取的x.t值汇总成一个数据框
data <- as.data.frame(do.call(cbind, x.t_list))

data <-cbind(metanew[,1:6],data)


sem.model1 =psem(
  lm(C_stock~ Available_N+Water+density+NPP+Texture+MAT_wc2+MAP_V2+TSEA_V2+PSEA_V2+ Soil_pH+cycfac+soil_age,na.action = na.fail, data = data),
  lm(cycfac ~ Available_N+ Water+density+ NPP +Texture+ MAT_wc2+MAP_V2+TSEA_V2+PSEA_V2+ Soil_pH+ soil_age, data = data),
  lm(Available_N ~ Water+density+ NPP +Texture+ MAT_wc2+MAP_V2+TSEA_V2+PSEA_V2+ Soil_pH+ soil_age, data = data),
  lm(NPP~ Water+density+ Texture+ MAT_wc2+MAP_V2+TSEA_V2+PSEA_V2+ Soil_pH+ soil_age, data = data),
  lm(Water ~ soil_age +density+ Texture+ MAT_wc2+MAP_V2+TSEA_V2+PSEA_V2, data = data),
  lm(density~ soil_age + Texture+ MAT_wc2+MAP_V2+TSEA_V2+PSEA_V2, data = data),
  lm(Soil_pH ~ density +Water+ soil_age + Texture+ MAT_wc2+MAP_V2+TSEA_V2+PSEA_V2, data = data),
  lm(Texture~ MAT_wc2+MAP_V2+TSEA_V2+PSEA_V2+ soil_age, data = data),
  MAP_V2%~~%TSEA_V2,
  MAP_V2%~~%MAT_wc2,
  MAP_V2%~~%PSEA_V2,
  MAT_wc2%~~%PSEA_V2,
  MAT_wc2%~~%TSEA_V2
  
)


dSep(sem.model1)


#模型中没有不显著路径
summary(sem.model1, .progressBar = T)


fisherC(sem.model1)
A<-summary(sem.model1, .progressBar = T)
#coefs(sem.model1,standardize = "none")
coefs(sem.model1,standardize = "scale")#可显示标准化
AIC(sem.model1)

B<-A$coefficients
write.csv(B,file= "cycfacC.csv")
plot(sem.model1)
