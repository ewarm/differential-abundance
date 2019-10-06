library("edgeR")
library("DESeq")
library("DESeq2")
library("ggplot2")
library("grid")
library("scales")
library("metagenomeSeq")
library("phyloseq")
library("plyr")
library("reshape2")
library("ROCR")
library("MGLM")
library("dirmult")
library("compositions")
library("ancom.R")

s=Sys.time()

#set parameters
theta = c(0.05,0.1,0.15,0.2,0.25,0.3,0.35)
beta = 0.3
setseed = 12453210
#define a function to simulate data sets
simlist = function(beta, thetai) {
  set.seed(setseed)
  phi1 <- rnorm(40, mean = 1000, sd = 300)
  #标准化到单位值
  phi1 <- phi1/sum(phi1)
  #phi1 = phi2
  phi2 <- phi1
  
  if (length(which(phi1<0)) != 0) {
    print(phi1)
  }
  if (length(which(phi2<0)) != 0) {
    print(phi2)
  }
  
  s = c(order(phi1, decreasing = T)[10:12])
  t = c(order(phi2, decreasing = T)[13:15])
  effotu = c(s,t)
  
  #对概率密度施加效应量
  for (i in s) {
    phi1[i] = phi1[i] + 0.05 * beta
    phi2[i] = phi2[i] - 0.05 * beta
  }
  
  for (i in t) {
    phi1[i] = phi1[i] - 0.05 * beta
    phi2[i] = phi2[i] + 0.05 * beta
  }
  
  #随机生成测序深度
  N <- round(runif(100, min = 1000, max = 10000))
  
  #产生数据
  #raw/None
  mat <- matrix(NA, 100, 40)
  for (i in 1:100) {
    if (i < 51) {
      mat[i,] <- simPop(1, 40, N[i], phi1, thetai)$data
    } else {
      mat[i,] <- simPop(1, 40, N[i], phi2, thetai)$data
    }
  }
  mylist <- list(effotu = effotu, mat = mat, N = N)
  return(mylist)
}



#------------------------------------------------------------------#
#---------------------------plot1(t-test)--------------------------#
#------------------------------------------------------------------#
####Plot1(t-test)
#define a function to perform t-test
t_test = function(mat) {
  pvalue = apply(mat, 2, function(x) t.test(x[1:50], x[51:100], paired = TRUE)$p.value)
  pvalue = p.adjust(pvalue, method = "BH")
  TPs = sum(which(pvalue<0.05) %in% effotu)
  P = length(which(pvalue<0.05))
  recall = TPs / length(effotu)
  if (P == 0) {
    precision = 0
  } else {
    precision = TPs / P
  }
  res = c(recall, precision)
  return(res)
}


#初始化数据框
df = data.frame(Recall=c(), Precision=c(), beta=c(), theta = c(), type=c(),stringsAsFactors = F)
for (thetai in theta) {
  Recall_1 = Recall_2 = Recall_3 = Recall_4 = 0
  Precision1 = Precision2 = Precision3 = Precision4 = 0
  Recall = c()
  Precision = c()
  setseed <- 12453210  
  
  for (j in 1:100) {
    list = simlist(beta, thetai)
    setseed <- setseed + 1
    mat = list[["mat"]]
    N = list[["N"]]
    effotu = list[["effotu"]]
    
    mat_tss <- matrix(NA, 100, 40)
    for (k in 1:100) {
      mat_tss[k,] <- mat[k,] / N[k]
    }
    mat_tss <- data.frame(clr(mat_tss))
    
    mat_rarefying <- otu_table(mat, taxa_are_rows = FALSE)
    mat_rarres <- rarefy_even_depth(mat_rarefying, sample.size = 0.9 * min(rowSums(mat)))@.Data
    
    mat_eBay <- matrix(NA, 100, 40)
    fit <- MGLMfit(data = mat, dist = "DM")
    alpha <- fit@estimate
    tot = sum(alpha)
    for (m in 1:100) {
      for (n in 1:40) {
        mat_eBay[m, n] <- (mat[m, n] + alpha[n]) / (N[m] + tot)
      }
    }
    mat_eBay <- data.frame(clr(mat_eBay))
    
    #求解
	ttest_mat <- t_test(mat)
	ttest_mat_tss <- t_test(mat_tss)
	ttest_mat_rarres <- t_test(mat_rarres)
	ttest_mat_eBay <- t_test(mat_eBay)
	
	#召回率和准确率加和
    Recall_1 = Recall_1 + ttest_mat[1]
    Recall_2 = Recall_2 + ttest_mat_tss[1]
    Recall_3 = Recall_3 + ttest_mat_rarres[1]
    Recall_4 = Recall_4 + ttest_mat_eBay[1]
    
    Precision1 = Precision1 + ttest_mat[2]
    Precision2 = Precision2 + ttest_mat_tss[2]
    Precision3 = Precision3 + ttest_mat_rarres[2]
    Precision4 = Precision4 + ttest_mat_eBay[2]
    
    print(j)
    
  }
  Recall_1 = Recall_1 / 100
  Recall_2 = Recall_2 / 100
  Recall_3 = Recall_3 / 100
  Recall_4 = Recall_4 / 100
  
  Precision1 = Precision1 / 100
  Precision2 = Precision2 / 100
  Precision3 = Precision3 / 100
  Precision4 = Precision4 / 100
  
  df = rbind(df, cbind(Recall = c(Recall_1), Precision = c(Precision1), beta = c(beta), theta = c(thetai), type = c("none-t")))
  df = rbind(df, cbind(Recall = c(Recall_2), Precision = c(Precision2), beta = c(beta), theta = c(thetai), type = c("tss-t")))
  df = rbind(df, cbind(Recall = c(Recall_3), Precision = c(Precision3), beta = c(beta), theta = c(thetai), type = c("rarefying-t")))
  df = rbind(df, cbind(Recall = c(Recall_4), Precision = c(Precision4), beta = c(beta), theta = c(thetai), type = c("eBay-t")))
}
e=Sys.time()
print(e-s)

write.csv(df, file = "E://QQPCMgr(1)/Desktop/final/df_t_theta.csv")




#------------------------------------------------------------------#
#--------------------------plot2(wilcoxon)-------------------------#
#------------------------------------------------------------------#

setseed = 12453210

#define a function to do wilcoxon-test
wilcoxon_test = function(physeq) {
  pvalue = apply(physeq, 2, function(x) wilcox.test(x[1:50], x[51:100], paired = TRUE)$p.value)
  pvalue = p.adjust(pvalue, method = "BH")
  TPs = sum(which(pvalue<0.05) %in% effotu)
  P = length(which(pvalue<0.05))
  recall = TPs / length(effotu)
  if (P == 0) {
    precision = 0
  } else {
    precision = TPs / P
  }
  res = c(recall, precision)
  return(res)
}


#初始化数据框
df = data.frame(Recall=c(), Precision=c(), beta=c(), theta = c(), type=c(),stringsAsFactors = F)

#进行循环
for (thetai in theta) {
  Recall_1 = Recall_2 = Recall_3 = Recall_4 = 0
  Precision1 = Precision2 = Precision3 = Precision4 = 0
  Recall = c()
  Precision = c()
  setseed <- 12453210
  
  for (j in 1:100) {
    list = simlist(beta, thetai)
	setseed <- setseed + 1
    mat = list[["mat"]]
    N = list[["N"]]
    effotu = list[["effotu"]]
    
    mat_tss <- matrix(NA, 100, 40)
    for (k in 1:100) {
      mat_tss[k,] <- mat[k,] / N[k]
    }
    mat_tss <- data.frame(clr(mat_tss))
    
    mat_rarefying <- otu_table(mat, taxa_are_rows = FALSE)
    mat_rarres <- rarefy_even_depth(mat_rarefying, sample.size = 0.9 * min(rowSums(mat)))@.Data
    
    mat_eBay <- matrix(NA, 100, 40)
    fit <- MGLMfit(data = mat, dist = "DM")
    alpha <- fit@estimate
    tot = sum(alpha)
    for (m in 1:100) {
      for (n in 1:40) {
        mat_eBay[m, n] <- (mat[m, n] + alpha[n]) / (N[m] + tot)
      }
    }
    mat_eBay <- data.frame(clr(mat_eBay))
    
   #求解
	wilcoxon_mat <- wilcoxon_test(mat)
	wilcoxon_mat_tss <- wilcoxon_test(mat_tss)
	wilcoxon_mat_rarres <- wilcoxon_test(mat_rarres)
	wilcoxon_mat_eBay <- wilcoxon_test(mat_eBay)
	
	#召回率和准确率加和
    Recall_1 = Recall_1 + wilcoxon_mat[1]
    Recall_2 = Recall_2 + wilcoxon_mat_tss[1]
    Recall_3 = Recall_3 + wilcoxon_mat_rarres[1]
    Recall_4 = Recall_4 + wilcoxon_mat_eBay[1]
    
    Precision1 = Precision1 + wilcoxon_mat[2]
    Precision2 = Precision2 + wilcoxon_mat_tss[2]
    Precision3 = Precision3 + wilcoxon_mat_rarres[2]
    Precision4 = Precision4 + wilcoxon_mat_eBay[2]
	
    print(j)
  }
  
  Recall_1 = Recall_1 / 100
  Recall_2 = Recall_2 / 100
  Recall_3 = Recall_3 / 100
  Recall_4 = Recall_4 / 100
  
  Precision1 = Precision1 / 100
  Precision2 = Precision2 / 100
  Precision3 = Precision3 / 100
  Precision4 = Precision4 / 100
  
  df = rbind(df, cbind(Recall = c(Recall_1), Precision = c(Precision1), beta = c(beta), theta = c(thetai), type = c("none-wilcoxon")))
  df = rbind(df, cbind(Recall = c(Recall_2), Precision = c(Precision2), beta = c(beta), theta = c(thetai), type = c("tss-wilcoxon")))
  df = rbind(df, cbind(Recall = c(Recall_3), Precision = c(Precision3), beta = c(beta), theta = c(thetai), type = c("rarefying-wilcoxon")))
  df = rbind(df, cbind(Recall = c(Recall_4), Precision = c(Precision4), beta = c(beta), theta = c(thetai), type = c("eBay-wilcoxon")))
}

e=Sys.time()
print(e-s)
write.csv(df, file = "E://QQPCMgr(1)/Desktop/final/df_wilcoxon_theta.csv")






#------------------------------------------------------------------#
#-------------------------plot3&4(多种方法)------------------------#
#------------------------------------------------------------------#

setseed = 12453210

#define DESeq2 function
DESeq2_test = function(physeq) {  
  countData = t(physeq)
  countData = round(countData, digits = 0)
  countData = countData + 1L
  p = c(paste0('sample', 1:100))
  condition = factor(c(rep("type1", 50), rep("type2", 50)), levels = c("type1", "type2"))
  colData = data.frame(row.names = p, condition)
  dds = DESeqDataSetFromMatrix(countData, colData, design = ~ condition)
  dds <- DESeq(dds)
  res = results(dds, contrast = c("condition", "type1", "type2"))
  pvalue = res$pvalue
  pvalue = p.adjust(pvalue, method = "BH")
  TPs = sum(which(pvalue<0.05) %in% effotu)
  P = length(which(pvalue<0.05))
  recall = TPs / length(effotu)
  if (P == 0) {
    precision = 0
  } else {
    precision = TPs / P
  }
  res = c(recall, precision)
  return(res)
}

#define ANCOM function
ANCOM_test = function(physeq) {
  physeq <- physeq + 1L
  prac_otu <- data.frame(physeq)
  prac_otu[,41] <- c(rep("Control",50), rep("Treatment",50))
  colnames(prac_otu) <- c( paste0("OTU_", 1:40 ), "Group" )
  ancom.out <- ANCOM( OTUdat = prac_otu, sig = 0.20, multcorr = 2 )
  b <- ancom.out$detected
  b <- as.numeric(gsub("OTU_", "", b))
  TPs = sum(b %in% effotu)
  P = length(b)
  recall = TPs / length(effotu)
  if (P == 0) {
    precision = 0
  } else {
    precision = TPs / P
  }
  res = c(recall, precision)
  return(res)
}

#define metagenomeSeq function
metagenomeSeq_test = function(physeq) {
  OTU = as(otu_table(t(mat), taxa_are_rows = T), "matrix")
  OTU = OTU + 1
  rownames(OTU) <- c(1:40)
  colnames(OTU) <- c(1:100)
  type = c(rep("type1", 50), rep("type2", 50))
  samdataframe <- data.frame(type)
  ADF = AnnotatedDataFrame(samdataframe)
  TDF = AnnotatedDataFrame(data.frame(taxa = c(1:40), row.names = c(1:40)))
  MGS = newMRexperiment(counts = OTU, phenoData = ADF, featureData = TDF)

  p = cumNormStatFast(MGS)
  MGS = cumNorm(MGS, p = p)

  pd <- pData(MGS)[,1]
  mod = model.matrix(~pd)
  fit = fitZig(obj = MGS, mod=mod)
  y <- MRcoefs(fit, number = 40)
  
  TPs = sum(as.numeric(rownames(y)[which(y$adjPvalues<0.05)]) %in% effotu)
  P = length(which(y$adjPvalues<0.05))
  recall = TPs / length(effotu)
  if (P == 0) {
    precision = 0
  } else {
    precision = TPs / P
  }
  res = c(recall, precision)
  return(res)
}



df = data.frame(Recall=c(), Precision=c(), beta=c(), theta = c(), type=c(),stringsAsFactors = F)
for (thetai in theta) {
  Recall_1 = Recall_2 = Recall_3 = Recall_4 = 0
  Precision1 = Precision2 = Precision3 = Precision4 = 0
  Recall = c()
  Precision = c()
  setseed <- 12453210
  
  for (j in 1:100) {
    list = simlist(beta, thetai)
    setseed <- setseed + 1
    mat = list[["mat"]]
    N = list[["N"]]
    effotu = list[["effotu"]]
    
    DESeq2_mat <- DESeq2_test(mat)
	ANCOM_mat <- ANCOM_test(mat)
	metagenomeSeq_mat <- metagenomeSeq_test(mat)
	
	
    Recall_1 = Recall_1 + DESeq2_mat[1]
    Recall_2 = Recall_2 + ANCOM_mat[1]
    Recall_3 = Recall_3 + metagenomeSeq_mat[1]
    
    Precision1 = Precision1 + DESeq2_mat[2]
    Precision2 = Precision2 + ANCOM_mat[2]
    Precision3 = Precision3 + metagenomeSeq_mat[2]
	
    print(j)
  }
  
  Recall_1 = Recall_1 / 100
  Recall_2 = Recall_2 / 100
  Recall_3 = Recall_3 / 100

  Precision1 = Precision1 / 100
  Precision2 = Precision2 / 100
  Precision3 = Precision3 / 100
  
  df = rbind(df, cbind(Recall = c(Recall_1), Precision = c(Precision1), beta = c(beta), theta = c(thetai), type = c("DESeq2")))
  df = rbind(df, cbind(Recall = c(Recall_2), Precision = c(Precision2), beta = c(beta), theta = c(thetai), type = c("ANCOM")))
  df = rbind(df, cbind(Recall = c(Recall_3), Precision = c(Precision3), beta = c(beta), theta = c(thetai), type = c("metagenomeSeq")))

}

e=Sys.time()
print(e-s)
write.csv(df, file = "E://QQPCMgr(1)/Desktop/final/df_other_theta.csv")