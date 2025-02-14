rm(list = ls())
library(devtools)  
#install_github("ramellose/CoNetinR")  
#install_github("hallucigenia-sparsa/seqtime") 
#?getNetwork
#install.packages("gdata")


library(gdata) 
library(CoNetinR)
library(dplyr)
library(seqtime)

otu1<-read.delim('WA1hou.txt',sep='\t',row.names = 1)#格
otu1<-as.matrix(otu1)
##不需要进行随机抽样时,计算score
#get the correlation relationships between taxon-taxon by Spearman method using CoNetR package
result1 = getNetwork(mat = otu1, method="spearman", permutandboot= T,norm=F,rarefy=F,T.up=0.2, iters=999,T.down=-0.2,bh=T, min.occ=0, keep.filtered=F, plot=T, report.full=T, verbose=T,shuffle.samples = F, stand.rows = FALSE, pval.cor = F, permut = F, renorm = F
)
#因为我之前进行了数据过滤和标准化
#计算score无关permu等信息的影响
scores1 = result1$scores#score就是相关性


result1$pospercent
result1$negpercent
result1$posnumber
result1$negnumber
result1$total
result1$pvalues

#####计算P值
pmatrix1= getNetwork(mat = otu1, method="spearman",pval.cor = T, permut =F, renorm =F,permutandboot= F,norm=F,rarefy=F,T.up=0.2, iters=999,T.down=-0.2,bh=T, min.occ=0, keep.filtered=F, plot=T, report.full=T, verbose=T,shuffle.samples = F, stand.rows = FALSE
)###OTU数量必须超过样品数量
pmatrix1 = pmatrix1$pvalues

#adjmatrix = matrix(nrow = N, ncol = N)
adjmatrix1 = matrix(nrow =  nrow(otu1), ncol =  nrow(otu1))
adjmatrix1[lower.tri(adjmatrix1)] = scores1
adjmatrix1 = t(adjmatrix1)
adjmatrix1[lower.tri(adjmatrix1)] = scores1
for (i in 1:nrow(otu1)){
  for (j in 1:nrow(otu1)){
    if (is.na(adjmatrix1[i,j])){
      adjmatrix1[i,j] = 0
    }
    else if (abs(adjmatrix1[i,j] ) < 0.6){
      adjmatrix1[i,j] = 0
    }
  }
}

rownames(adjmatrix1)=rownames(otu1)
colnames(adjmatrix1)=rownames(otu1)

write.csv(adjmatrix1,"adjmatrix1.csv")


###利用bray
result2 = getNetwork(mat = otu1, method="bray", permutandboot= T,norm=F,rarefy=F,T.up=0.2, iters=999,T.down=-0.2,bh=T, min.occ=0, keep.filtered=F, plot=T, report.full=T, verbose=T,shuffle.samples = F, stand.rows = FALSE, pval.cor = F, permut = F, renorm = F)
scores2 = result2$scores
#shuffle.samples随机化处理： 对每个样本的数据进行随机打乱。例如，如果数据以矩阵形式给出，每列代表一个样本，每行代表一个特征（如OTU、基因表达量等），则在不改变每列（样本）内部数据顺序的情况下随机打乱每行的值
#####计算P值
pmatrix2= getNetwork(mat = otu1, method="bray",pval.cor = F, permut =T, permutandboot= T,renorm =T,norm=F,rarefy=F,T.up=0.2, iters=999,T.down=-0.2,bh=T, min.occ=0, keep.filtered=F, plot=T, report.full=T, verbose=T,shuffle.samples = F, stand.rows = FALSE
)
pmatrix2 = pmatrix2$pvalues

#adjmatrix = matrix(nrow = N, ncol = N)
adjmatrix2 = matrix(nrow =  nrow(otu1), ncol =  nrow(otu1))
adjmatrix2[lower.tri(adjmatrix2)] = scores2
adjmatrix2 = t(adjmatrix2)
adjmatrix2[lower.tri(adjmatrix2)] = scores2
for (i in 1:nrow(otu1)){
  for (j in 1:nrow(otu1)){
    if (is.na(adjmatrix2[i,j])){
      adjmatrix2[i,j] = 0
    }
    else if (pmatrix2[i,j] > 0.05){
      adjmatrix2[i,j] = 0.5
    }
  }
}
rownames(adjmatrix2)=rownames(otu1)
colnames(adjmatrix2)=rownames(otu1)
write.csv(adjmatrix2,"adjmatrix2.csv")
####两个矩阵取交集
#如果你的意思是找出两个矩阵中对应位置都不为零（或都满足某条件）的元素，你可以通过逻辑运算来实现。下面是一个例子，
#展示如何处理两个矩阵，只保留两个矩阵在相同位置上都大于0的元素，其余位置设为0（这可以视为一种“交集”的操作）。
##假设我们有两个矩阵mat1和mat2，我们想要创建一个新的矩阵，这个矩阵只保留mat1和mat2在相同位置上都大于0的元素，其余位置设为0。
# 示例矩阵
#set.seed(123) # 确保示例可重现
#mat1 <- matrix(sample(-2:2, 25, replace = TRUE), nrow = 5)
#mat2 <- matrix(sample(-2:2, 25, replace = TRUE), nrow = 5)

# 查看矩阵
adjmatrix1=read.csv("adjmatrix1.csv",row.names=1)
#mat2
adjmatrix2=read.csv("adjmatrix2.csv",row.names=1)
#mat2
# 创建一个新矩阵，只包含两个矩阵对应位置都不等于的元素.spearman关注r和P值，braycurtis只关注P值是否显著
mat_common <- ifelse(adjmatrix1 != 0 & adjmatrix2 != 0.5,adjmatrix1, 0)

# 查看结果
mat_common

write.csv(mat_common,"adjrs1.csv")
#这个矩阵输入igraph得出网络

#RMT
#following by an RMT-based approach that determines the correlation cut-off threshold in an automatic fashion. 
#Random matrix theory was initially proposed in the 1960s as a procedure to identify phase transitions associated with noise in physics and material science,
#and was later adopted for studying the behaviours of many other complex systems, including gene co-expression network construction for
#predicting gene functions and molecular ecological network construction.


#construct function, revised by microeco R package
nnsd = function(sp){
  nns <- NULL
  for(j in 2:length(sp)){
    nn <- abs(sp[j] - sp[j-1])
    nns <- c(nns, nn)
  }
  nns
}
rmt = function(cormat,lcor=0.6, hcor=0.8){
  s <- seq(0, 3, 0.1)
  pois <- exp(-s)
  geo <- 0.5 * pi * s * exp(-0.25 * pi * s^2)
  ps <- NULL  
  for(i in seq(lcor, hcor, 0.01)){
    cormat1 <- abs(cormat)
    cormat1[cormat1 < i] <- 0  
    eigen <- sort(eigen(cormat1)$value)
    ssp <- smooth.spline(eigen, control.spar = list(low = 0, high = 3)) 
    nnsd1 <- density(nnsd(ssp$y))
    nnsdpois <- density(nnsd(pois))
    chival1 <- sum((nnsd1$y - nnsdpois$y)^2/nnsdpois$y/512)
    ps <- rbind(ps, chival1)
    print(i*100)
  }
  ps <- cbind(ps, c(seq(lcor, hcor, 0.01)))
  tc <- ps[ps[,1] == min(ps[,1]), 2]
  tc
}


tc1 <- rmt(adj_mat1)
message("The optimized COR threshold: ", tc1, "...\n")

adj_mat1[abs(adj_mat1)< tc1] = 0

#save the adjacency matrix data
write.csv(adj_mat1,"adj_mat1.tc.csv")