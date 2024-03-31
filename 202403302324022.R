# Downloading and cleaning Depmap data https://depmap.org/portal/download/all/
# Assuming you have already downloaded and cleaned the Depmap data

## load libraries
library(dplyr)
library(tibble)
library(qs)
library(openxlsx)
library(purrr)

## Read data files
gene_effect <- data.table::fread("../public_databases/depmap/DepMap_Public_23Q4/CRISPRGeneEffect.csv", data.table = FALSE)
gene_dependency <- data.table::fread("../public_databases/depmap/DepMap_Public_23Q4/CRISPRGeneDependency.csv", data.table = FALSE)
expr_set <- data.table::fread("../public_databases/depmap/DepMap_Public_23Q4/CCLE_expression.csv", data.table = FALSE)
cell_info <- data.table::fread("../public_databases/depmap/DepMap_Public_23Q4/sample_info.csv", data.table = FALSE)
mut_data <- data.table::fread("../public_databases/depmap/DepMap_Public_23Q4/CCLE_mutations.csv", data.table = FALSE)


## Trimming cell line information
commonindex <- intersect(gene_effect$V1, expr_set$V1)

## gene_effect
gene_effect <- gene_effect %>%
  mutate(rowname = gene_effect[, 1]) %>%
  select(-1) %>%
  filter(rowname %in% commonindex) %>%
  rename_all(~gsub("\\s+\\(\\d+\\)", "", .)) %>%
  column_to_rownames(var = "rowname")
gene_effect[1:5, 1:5]

## gene_dependency
gene_dependency <- gene_dependency %>%
  mutate(rowname = gene_dependency[, 1]) %>%
  select(-1) %>%
  filter(rowname %in% commonindex) %>%
  rename_all(~gsub("\\s+\\(\\d+\\)", "", .)) %>%
  column_to_rownames(var = "rowname")
gene_dependency[1:5, 1:5]

## expr_set
expr_set <- expr_set %>%
  mutate(rowname = expr_set[, 1]) %>%
  select(-1) %>%
  filter(rowname %in% commonindex) %>%
  rename_all(~gsub("\\s+\\(\\d+\\)", "", .)) %>%
  column_to_rownames(var = "rowname")
expr_set[1:5, 1:5]

## cell_info
cell_info <- cell_info %>%
  mutate(rowname = cell_info[, 1]) %>%
  select(-1) %>%
  filter(rowname %in% commonindex) %>%
  column_to_rownames(var = "rowname")
table(cell_info$primary_disease)
nrow(cell_info)

## commonGenes
commonGenes <- intersect(colnames(gene_effect), colnames(expr_set))
gene_effect <- gene_effect[, commonGenes]
gene_dependency <- gene_dependency[, commonGenes]
expr_set <- expr_set[, commonGenes]

## save qs
qs::qsave(gene_effect, file = "data/raw/depmap_gene_effect.qs")
qs::qsave(gene_dependency, file = "data/raw/depmap_gene_dependency.qs")
qs::qsave(expr_set, file = "data/raw/ccle_expr_Set.qs")
qs::qsave(cell_info, file = "data/raw/ccle_cell_info.qs")
qs::qsave(mut_data, file = "data/raw/ccle_mut_data.qs")


# CCLE中基因和基因的相关性
rm(list = ls());gc()
exprSet <- qs::qread("data/raw/ccle_expr_Set.qs")
exprSet <- exprSet[, colSums(is.na(exprSet)) == 0]
# 批量计算KBTBD2与其它所有基因的相关性
corData <- colnames(exprSet) %>%
  map_df(function(gene2) {
    gene1 <- "KBTBD2"
    genedata <- exprSet[, gene1]
    dd <- cor.test(genedata, exprSet[, gene2])
    tibble(Gene1 = gene1, Gene2 = gene2, cor = dd$estimate, pvalue = dd$p.value)
  }) %>%
  arrange(desc(cor))
# Convert 'cor' column to numeric
corData$cor <- as.numeric(as.character(corData$cor))
# Save corData as a CSV file in the 'results/tables' folder
write.csv(corData, file = "results/tables/KBTBD2_CCLE_corgenes.csv", row.names = FALSE)
## 把所有疾病的细胞系提取出来
# Load cell info data
cell_info <- qs::qread("data/raw/ccle_cell_info.qs")
# Get unique primary diseases
unique_diseases <- unique(cell_info$primary_disease)
#去除少于7种细胞系的疾病
unique_diseases <- unique_diseases[table(cell_info$primary_disease) > 7]
# Loop through each primary disease
# Create an Excel file to store the results
output_file <- "results/tables/KBTBD2_CCLE_corgenes.xlsx"
wb <- createWorkbook()
# Use lapply with parallel processing
library(parallel)
cl <- makeCluster(3)  # Create a cluster with 4 cores
# Define a function to calculate correlation for each primary disease
calculate_correlation <- function(primary_disease) {
  # Filter exprSet for the current primary disease
  exprSet_filtered <- exprSet[rownames(cell_info[cell_info$primary_disease == primary_disease, ]), ]
  exprSet_filtered <- exprSet_filtered[, colSums(is.na(exprSet_filtered)) == 0]
  # Calculate correlation between KBTBD2 and all genes
  corData <- data.frame()
  gene1 <- "KBTBD2"
  genedata <- exprSet_filtered[, gene1]
  # Initialize corData as a data frame
  corData <- data.frame(Gene1 = character(), Gene2 = character(), cor = numeric(), pvalue = numeric())
  # lappy用起来
  corData <- do.call(rbind, lapply(seq_len(ncol(exprSet_filtered)), function(j) {
    ## 1. Indication
    print(j)
    ## 2. Calculation
    gene2 <- colnames(exprSet_filtered)[j]
    ## Check if there are enough finite observations
    if (sum(!is.na(genedata), !is.na(exprSet_filtered[, gene2]), na.rm = TRUE) > 2) {
      dd <- cor.test(genedata, exprSet_filtered[, gene2], use = "pairwise.complete.obs")
    ## 3. Storage
      if (!is.null(dd$estimate) && length(dd$estimate) == 1 && 
          !is.null(dd$p.value) && length(dd$p.value) == 1) {
        return(data.frame(Gene1 = gene1, Gene2 = gene2, cor = dd$estimate, pvalue = dd$p.value))
      }
    } else {
      return(data.frame(Gene1 = character(), Gene2 = character(), cor = numeric(), pvalue = numeric()))
    }
  }))
  # 添加列名
  colnames(corData) <- c("Gene1", "Gene2", "cor", "pvalue")
  ## Sort by correlation
  corData <- corData[order(corData$cor, decreasing = TRUE), ]
  
  return(list(primary_disease = primary_disease, corData = corData))
}
clusterExport(cl, c("exprSet", "cell_info"))
## Apply the function to each primary disease using parallel processing
results <- parLapply(cl, unique_diseases, calculate_correlation)
## Stop the cluster
stopCluster(cl)
## Write the results to the workbook
for (result in results) {
  primary_disease <- result$primary_disease
  corData <- result$corData
  
  addWorksheet(wb, primary_disease)
  writeData(wb, primary_disease, corData)
}
saveWorkbook(wb, output_file)


# Co-effect的探索
rm(list = ls());gc()
exprSet <- qs::qread("data/raw/depmap_gene_effect.qs")
exprSet <- exprSet[, colSums(is.na(exprSet)) == 0]
corData <- colnames(exprSet) %>%
  map_df(function(gene2) {
    gene1 <- "KBTBD2"
    genedata <- exprSet[, gene1]
    dd <- cor.test(genedata, exprSet[, gene2])
    tibble(Gene1 = gene1, Gene2 = gene2, cor = dd$estimate, pvalue = dd$p.value)
  }) %>%
  arrange(desc(cor))
# Convert 'cor' column to numeric
corData$cor <- as.numeric(as.character(corData$cor))
# Save corData as a CSV file in the 'results/tables' folder
write.csv(corData, file = "results/tables/KBTBD2_gene_effect_corgenes.csv", row.names = FALSE)
## 把所有疾病的细胞系提取出来
# Load cell info data
cell_info <- qs::qread("data/raw/ccle_cell_info.qs")
# Get unique primary diseases
unique_diseases <- unique(cell_info$primary_disease)
#去除少于7种细胞系的疾病
unique_diseases <- unique_diseases[table(cell_info$primary_disease) > 7]
# Loop through each primary disease
# Create an Excel file to store the results
output_file <- "results/tables/KBTBD2_gene_effect_corgenes.xlsx"
wb <- createWorkbook()
# Use lapply with parallel processing
library(parallel)
cl <- makeCluster(3)  # Create a cluster with 4 cores
# Define a function to calculate correlation for each primary disease
calculate_correlation <- function(primary_disease) {
  # Filter exprSet for the current primary disease
  exprSet_filtered <- exprSet[rownames(cell_info[cell_info$primary_disease == primary_disease, ]), ]
  exprSet_filtered <- exprSet_filtered[, colSums(is.na(exprSet_filtered)) == 0]
  # Calculate correlation between KBTBD2 and all genes
  corData <- data.frame()
  gene1 <- "KBTBD2"
  genedata <- exprSet_filtered[, gene1]
  # Initialize corData as a data frame
  corData <- data.frame(Gene1 = character(), Gene2 = character(), cor = numeric(), pvalue = numeric())
  # lappy用起来
  corData <- do.call(rbind, lapply(seq_len(ncol(exprSet_filtered)), function(j) {
    ## 1. Indication
    print(j)
    ## 2. Calculation
    gene2 <- colnames(exprSet_filtered)[j]
    ## Check if there are enough finite observations
    if (sum(!is.na(genedata), !is.na(exprSet_filtered[, gene2]), na.rm = TRUE) > 2) {
      dd <- cor.test(genedata, exprSet_filtered[, gene2], use = "pairwise.complete.obs")
    ## 3. Storage
      if (!is.null(dd$estimate) && length(dd$estimate) == 1 && 
          !is.null(dd$p.value) && length(dd$p.value) == 1) {
        return(data.frame(Gene1 = gene1, Gene2 = gene2, cor = dd$estimate, pvalue = dd$p.value))
      }
    } else {
      return(data.frame(Gene1 = character(), Gene2 = character(), cor = numeric(), pvalue = numeric()))
    }
  }))
  # 添加列名
  colnames(corData) <- c("Gene1", "Gene2", "cor", "pvalue")
  ## Sort by correlation
  corData <- corData[order(corData$cor, decreasing = TRUE), ]
  
  return(list(primary_disease = primary_disease, corData = corData))
}
clusterExport(cl, c("exprSet", "cell_info"))
## Apply the function to each primary disease using parallel processing
results <- parLapply(cl, unique_diseases, calculate_correlation)
## Stop the cluster
stopCluster(cl)
## Write the results to the workbook
for (result in results) {
  primary_disease <- result$primary_disease
  corData <- result$corData
  
  addWorksheet(wb, primary_disease)
  writeData(wb, primary_disease, corData)
}
saveWorkbook(wb, output_file)


# gene dependency的使用
rm(list = ls());gc()
exprSet <- qs::qread("data/raw/depmap_gene_dependency.qs")
exprSet <- exprSet[, colSums(is.na(exprSet)) == 0]
corData <- colnames(exprSet) %>%
  map_df(function(gene2) {
    gene1 <- "KBTBD2"
    genedata <- exprSet[, gene1]
    dd <- cor.test(genedata, exprSet[, gene2])
    tibble(Gene1 = gene1, Gene2 = gene2, cor = dd$estimate, pvalue = dd$p.value)
  }) %>%
  arrange(desc(cor))
# Convert 'cor' column to numeric
corData$cor <- as.numeric(as.character(corData$cor))
# Save corData as a CSV file in the 'results/tables' folder
write.csv(corData, file = "results/tables/KBTBD2_gene_dependency_corgenes.csv", row.names = FALSE)
## 把所有疾病的细胞系提取出来
# Load cell info data
cell_info <- qs::qread("data/raw/ccle_cell_info.qs")
# Get unique primary diseases
unique_diseases <- unique(cell_info$primary_disease)
#去除少于7种细胞系的疾病
unique_diseases <- unique_diseases[table(cell_info$primary_disease) > 7]
# Loop through each primary disease
# Create an Excel file to store the results
output_file <- "results/tables/KBTBD2_gene_dependency_corgenes.xlsx"
wb <- createWorkbook()
# Use lapply with parallel processing
library(parallel)
cl <- makeCluster(3)  # Create a cluster with 4 cores
# Define a function to calculate correlation for each primary disease
calculate_correlation <- function(primary_disease) {
  # Filter exprSet for the current primary disease
  exprSet_filtered <- exprSet[rownames(cell_info[cell_info$primary_disease == primary_disease, ]), ]
  exprSet_filtered <- exprSet_filtered[, colSums(is.na(exprSet_filtered)) == 0]
  # Calculate correlation between KBTBD2 and all genes
  corData <- data.frame()
  gene1 <- "KBTBD2"
  genedata <- exprSet_filtered[, gene1]
  # Initialize corData as a data frame
  corData <- data.frame(Gene1 = character(), Gene2 = character(), cor = numeric(), pvalue = numeric())
  # lappy用起来
  corData <- do.call(rbind, lapply(seq_len(ncol(exprSet_filtered)), function(j) {
    ## 1. Indication
    print(j)
    ## 2. Calculation
    gene2 <- colnames(exprSet_filtered)[j]
    ## Check if there are enough finite observations
    if (sum(!is.na(genedata), !is.na(exprSet_filtered[, gene2]), na.rm = TRUE) > 2) {
      dd <- cor.test(genedata, exprSet_filtered[, gene2], use = "pairwise.complete.obs")
    ## 3. Storage
      if (!is.null(dd$estimate) && length(dd$estimate) == 1 && 
          !is.null(dd$p.value) && length(dd$p.value) == 1) {
        return(data.frame(Gene1 = gene1, Gene2 = gene2, cor = dd$estimate, pvalue = dd$p.value))
      }
    } else {
      return(data.frame(Gene1 = character(), Gene2 = character(), cor = numeric(), pvalue = numeric()))
    }
  }))
  # 添加列名
  colnames(corData) <- c("Gene1", "Gene2", "cor", "pvalue")
  ## Sort by correlation
  corData <- corData[order(corData$cor, decreasing = TRUE), ]
  
  return(list(primary_disease = primary_disease, corData = corData))
}
clusterExport(cl, c("exprSet", "cell_info"))
## Apply the function to each primary disease using parallel processing
results <- parLapply(cl, unique_diseases, calculate_correlation)
## Stop the cluster
stopCluster(cl)
## Write the results to the workbook
for (result in results) {
  primary_disease <- result$primary_disease
  corData <- result$corData
  
  addWorksheet(wb, primary_disease)
  writeData(wb, primary_disease, corData)
}
saveWorkbook(wb, output_file)


# 基因A的Dependdency受到哪些基因调控
## 或者哪些基因的表达能够预测抑制剂疗效
rm(list = ls());gc()
geneEffect <- qs::qread("data/raw/depmap_gene_effect.qs")
exprSet <- qs::qread("data/raw/ccle_expr_Set.qs")
gene2 <- "KBTBD2"
genedata <- geneEffect[,"KBTBD2"]
corData <- do.call(rbind, lapply(colnames(exprSet), function(gene1) {
  dd <- cor.test(genedata, exprSet[, gene1])
  data.frame(Gene1 = gene1, Gene2 = gene2, cor = dd$estimate, pvalue = dd$p.value)
}))
corData <- corData[order(corData$cor, decreasing = TRUE), ]
#导出为预测抑制剂疗效基因列表
write.csv(corData, file = "results/tables/KBTBD2_inhibitor_effect_predictor.csv", row.names = FALSE)


### 基因A的表达量调控哪些基因的Dependdency
rm(list = ls())
exprSet <- qs::qread(file = "data/raw/ccle_expr_Set.qs")
geneEffect <- qs::qread(file = "data/raw/depmap_gene_effect.qs")
gene1 <- "KBTBD2"
genedata <- exprSet[,"KBTBD2"]
corData <- do.call(rbind, lapply(colnames(geneEffect), function(gene2) {
  dd <- cor.test(genedata, geneEffect[, gene2])
  data.frame(Gene1 = gene1, Gene2 = gene2, cor = dd$estimate, pvalue = dd$p.value)
}))
corData <- corData[order(corData$cor, decreasing = TRUE), ]
#导出为预测抑制剂疗效基因列表
write.csv(corData, file = "results/tables/KBTBD2_regulate_Genes_Dependdency.csv", row.names = FALSE)


### GSEA
library(clusterProfiler)
mygeneList <- -corData$cor
names(mygeneList) <- corData$Gene2
mygeneList <- sort(mygeneList,decreasing = T)
head(mygeneList)

geneSet <- read.gmt("resource/geneSets/h.all.v2023.2.Hs.symbols.gmt")

mygsea <- GSEA(geneList = mygeneList,TERM2GENE = geneSet)
data <- as.data.frame(mygsea)
library(ggplot2)
dotplot(mygsea,showCategory=30,
        split=".sign",
        font.size = 8,
        label_format = 60)+facet_grid(~.sign)

library(enrichplot)
gseaplot2(mygsea,"HALLMARK_MYC_TARGETS_V2",color = "red",pvalue_table = T)

########################################################
### IL6R
rm(list = ls())
exprSet <- readRDS(file = "output/ccle_exprSet_1005_17285.rds")
genedata <- exprSet[,"IL6R"]

geneEffect <- readRDS(file = "output/depmap_geneEffect_1005_17285.rds")

corData <- data.frame()
gene1 = "IL6R"
genedata <- exprSet[,"IL6R"]
for (i in seq_len(ncol(exprSet))) {
  
  ## 1指示
  print(i)
  
  ## 2计算
  gene2 = colnames(geneEffect)[i]
  dd = cor.test(genedata,geneEffect[,gene2])
  
  ## 3.储存
  corData[i,1] = gene1
  corData[i,2] = gene2
  corData[i,3] = dd$estimate
  corData[i,4] = dd$p.value
}

colnames(corData) <- c("exp","Dependency","cor","pvalue")

### GSEA
library(clusterProfiler)
mygeneList <- -corData$cor
names(mygeneList) <- corData$Dependency
mygeneList <- sort(mygeneList,decreasing = T)
head(mygeneList)

geneSet <- read.gmt("resource/geneSets/c2.cp.reactome.v2023.2.Hs.symbols.gmt")

mygsea <- GSEA(geneList = mygeneList,TERM2GENE = geneSet)
data <- as.data.frame(mygsea)
library(ggplot2)
dotplot(mygsea,showCategory=30,
        split=".sign",
        font.size = 8,
        label_format = 60)+facet_grid(~.sign)

library(enrichplot)
gseaplot2(mygsea,"REACTOME_CHROMOSOME_MAINTENANCE",color = "red",pvalue_table = T)

### 推荐阅读:
### cGAS–STING drives the IL-6-dependent survival of chromosomally instable cancers
### https://www.bioconductor.org/packages/release/bioc/vignettes/decoupleR/inst/doc/pw_bk.html
### https://saezlab.github.io/progeny/articles/progeny.html

### 拓展方向
### 机器学习
### 协同致死
### 药物敏感性
### 基因功能模块
################################################
################################################
### 作者：果子
### 更新时间：2023-11-22
### 微信公众号:果子学生信
### 私人微信：guotosky
### 个人邮箱：hello_guozi@126.com

############################################################
## 数据预处理
## 1.突变信息
rm(list = ls())
mutData_raw <- data.table::fread("data/DepMap_Public_20Q4v2/CCLE_mutations.csv",data.table = F)
library(dplyr)
library(tidyr)
library(tibble)
mutData <- mutData_raw %>% 
  dplyr::select(DepMap_ID,Hugo_Symbol,isDeleterious) %>% 
  mutate(isDeleterious= as.numeric(isDeleterious)) %>%
  group_by(DepMap_ID,Hugo_Symbol) %>% 
  summarise(isDeleterious=sum(isDeleterious)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "Hugo_Symbol",
              values_from = "isDeleterious",
              values_fill = 0) %>% 
  as.data.frame() %>% 
  column_to_rownames("DepMap_ID")

saveRDS(mutData,file = "data/DepMap_Public_20Q4v2/mutData.rds")

## 2.DepMap 基因Dependency 数据
rm(list = ls())
geneDependency <- data.table::fread("data/DepMap_Public_20Q4v2/Achilles_gene_dependency.csv",data.table = F)
test <- geneDependency[1:10,1:10]
## 修改列名
colnames(geneDependency) <- gsub("\\s+\\(\\d+\\)","",colnames(geneDependency))
## 第一列变行名
rownames(geneDependency) <- geneDependency[,1]
geneDependency <- geneDependency[,-1]
## 保存数据
saveRDS(geneDependency,file = "data/DepMap_Public_20Q4v2/geneDependency.rds")


## 3.细胞系信息
rm(list = ls())
cellinfor <- data.table::fread("data/DepMap_Public_20Q4v2/sample_info.csv",data.table = F)
rownames(cellinfor) <- cellinfor$DepMap_ID

saveRDS(cellinfor,file = "data/DepMap_Public_20Q4v2/cellinfor.rds")

###############################################################
###############################################################
### 图复现: ARID1A 和HMGCR 的例子
rm(list = ls())

### 加载R包和数据
library(dplyr)
library(tidyr)
library(tibble)
mutData <- readRDS(file = "data/DepMap_Public_20Q4v2/mutData.rds")
geneDependency <- readRDS(file = "data/DepMap_Public_20Q4v2/geneDependency.rds")
cellinfor <- readRDS(file = "data/DepMap_Public_20Q4v2/cellinfor.rds")

### 细胞系取交集，789
coID <- intersect(rownames(geneDependency),rownames(cellinfor)) %>% 
  intersect(rownames(mutData))

### 提取数据 ARID1A的突变信息，以及HMGCR的Dependency 信息
data <- data.frame(
  DepmapID = coID,
  Tissue = cellinfor[coID,"primary_disease"],
  Mutation = ifelse(mutData[coID,"ARID1A"]==0,"WT","Mut"),
  Dependency = geneDependency[coID,"HMGCR"]
)

## 数据筛选，Mut >=2,WT>=5(跟原文不同)
TissueStatus <- data %>% 
  select(Tissue,Mutation) %>%
  group_by(Tissue,Mutation) %>% 
  summarise(n =n()) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "Mutation",
              values_from = n,
              values_fill = 0) %>% 
  filter(Mut >=2,WT>=5) %>%
  mutate(Sum=Mut+WT) %>% 
  mutate(Mut_Number = sum(Mut),WT_Number=sum(WT)) 

  

mydata <- data %>% 
  filter(Tissue %in% TissueStatus$Tissue)

## 作图数据整理
Mean_Dep <- mydata %>% 
  filter(Mutation=="Mut") %>% 
  group_by(Tissue) %>% 
  summarise(Sensitivity = mean(Dependency))

Diff_Dep <- mydata %>%
  group_by(Tissue) %>%
  summarize(Difference = mean(Dependency[Mutation == "Mut"]) - mean(Dependency[Mutation == "WT"]))

## 数据合并
PlotData <- merge(Mean_Dep,Diff_Dep,by="Tissue") %>% 
  merge(TissueStatus,by="Tissue")


## 作图展示
library(ggplot2)
library(ggrepel)
library(stringr)
library(ggfun)

ggplot(PlotData, aes(x = Sensitivity, y = Difference)) +
  geom_point(aes(size = Sum), shape = 21, color = "black", fill = "green", alpha = 0.6, stroke = 1.5) + 
  geom_text_repel(aes(label = Tissue), vjust = 1.5, hjust = 1.5) + 
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf,
    fill = "yellow", alpha = 0.1, colour = "green", linetype = "dashed", linewidth = 1) +
  labs(x = "Sensitivity Score for Cell Lines with Mutations", 
  y = "Mutant vs. Wildtype Sensitivity", 
  title = paste0("Mutant (n=", unique(PlotData$Mut_Number),
          " cell lines) vs. Wildtype (n=",
          unique(PlotData$WT_Number), " cell lines)")) +
  theme_bw() + 
  guides(size = guide_legend(title = str_wrap("Number of cells", width = 8))) +
  theme(panel.grid.major = element_line(linetype = "dashed"),
    panel.grid.minor = element_line(linetype = "dashed", size = 0.5),
    panel.background = element_blank(),
    legend.position = c(.95, .05),
    legend.justification = c("right", "bottom"),
    legend.box.just = "right",
    legend.margin = margin(4, 4, 4, 4),
    legend.background = element_roundrect(color = "#808080", linetype = 1)
  )


##############################################################
### 换个基因 HMGCS1

data <- data.frame(
  DepmapID = coID,
  Tissue = cellinfor[coID,"primary_disease"],
  Mutation = ifelse(mutData[coID,"ARID1A"]==0,"WT","Mut"),
  Dependency = geneDependency[coID,"HMGCS1"]
)

## 数据筛选，Mut >=2,WT>=5(跟原文不同)
TissueStatus <- data %>% 
  select(Tissue,Mutation) %>%
  group_by(Tissue,Mutation) %>% 
  summarise(n =n()) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "Mutation",
              values_from = n,values_fill = 0) %>% 
  filter(Mut >=2,WT>=5) %>%
  mutate(Sum=Mut+WT) %>% 
  mutate(Mut_Number = sum(Mut),WT_Number=sum(WT)) 



mydata <- data %>% 
  filter(Tissue %in% TissueStatus$Tissue)

## 作图数据整理
Mean_Dep <- mydata %>% 
  filter(Mutation=="Mut") %>% 
  group_by(Tissue) %>% 
  summarise(Sensitivity = mean(Dependency))

Diff_Dep <- mydata %>%
  group_by(Tissue) %>%
  summarize(Difference = mean(Dependency[Mutation == "Mut"]) - mean(Dependency[Mutation == "WT"]))

## 数据合并
PlotData <- merge(Mean_Dep,Diff_Dep,by="Tissue")
PlotData <- merge(PlotData,TissueStatus,by="Tissue")

## 作图展示
library(ggplot2)
library(ggrepel)
library(stringr)

ggplot(PlotData, aes(x = Sensitivity, y = Difference)) +
  geom_point(aes(size = Sum), shape = 16, color = "blue", fill = "lightblue", alpha = 0.8, stroke = 1) + 
  geom_text_repel(aes(label = Tissue), vjust = -1, hjust = 0) + 
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf, 
       fill = "yellow", alpha = 0.1, colour = "green", linetype = "dashed", size = 1) +
  labs(x = "Sensitivity Score for Cell Lines with Mutations", 
   y = "Mutant vs. Wildtype Sensitivity", 
   title = paste0("Mutant (n=", unique(PlotData$Mut_Number),
          " cell lines) vs. Wildtype (n=",
          unique(PlotData$WT_Number), " cell lines)")) +
  theme_bw() + 
  guides(size = guide_legend(title = str_wrap("Number of cells", width = 8))) +
  theme(panel.grid.major = element_line(linetype = "dashed"),
    panel.grid.minor = element_line(linetype = "dashed", size = 0.5),
    panel.background = element_blank(),
    legend.position = c(.95, .05),
    legend.justification = c("right", "bottom"),
    legend.box.just = "right",
    legend.margin = margin(4, 4, 4, 4),
    legend.background = element_rect(color = "#808080", linetype = 1, fill = "white")
  )


###################################################
### 提取典型组织单独画图

library(ggplot2)
df = mydata[mydata$Tissue=="Endometrial/Uterine Cancer",]
df$Mutation <- factor(df$Mutation, levels = c("WT", "Mut"),
                      labels = c("ARID1A(wt) \n (n=11)", "ARID1A(mut) \n (n=15)"))
ggplot(df, aes(x = Mutation, y = Dependency, fill = Mutation)) +
  geom_violin(alpha = 0.5) +
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               dotsize = 0.5) +
  labs(y = "Dependency Score")+
  theme_bw() +
  theme(legend.position = "none")


################################################
################################################
### 作者：果子
### 更新时间：2023-11-22
### 微信公众号:果子学生信
### 私人微信：guotosky
### 个人邮箱：hello_guozi@126.com

############################################################
## 数据预处理
## 1.突变信息
rm(list = ls())
mutData_raw <- data.table::fread("data/DepMap_Public_22Q2/CCLE_mutations.csv",data.table = F)
library(dplyr)
library(tidyr)
library(tibble)
mutData <- mutData_raw %>% 
  dplyr::select(DepMap_ID,Hugo_Symbol,isDeleterious) %>% 
  mutate(isDeleterious= as.numeric(isDeleterious)) %>%
  group_by(DepMap_ID,Hugo_Symbol) %>% 
  summarise(isDeleterious=sum(isDeleterious)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "Hugo_Symbol",
              values_from = "isDeleterious",
              values_fill = 0) %>% 
  as.data.frame() %>% 
  column_to_rownames("DepMap_ID")

saveRDS(mutData,file = "data/DepMap_Public_22Q2/mutData.rds")

## 2.DepMap 基因Dependency 数据
rm(list = ls())
geneDependency <- data.table::fread("data/DepMap_Public_22Q2/CRISPR_gene_dependency.csv",data.table = F)
test <- geneDependency[1:10,1:10]
## 修改列名
colnames(geneDependency) <- gsub("\\s+\\(\\d+\\)","",colnames(geneDependency))
## 第一列变行名
rownames(geneDependency) <- geneDependency[,1]
geneDependency <- geneDependency[,-1]
## 保存数据
saveRDS(geneDependency,file = "data/DepMap_Public_22Q2/geneDependency.rds")


## 3.细胞系信息
rm(list = ls())
cellinfor <- data.table::fread("data/DepMap_Public_22Q2/sample_info.csv",data.table = F)
rownames(cellinfor) <- cellinfor$DepMap_ID

saveRDS(cellinfor,file = "data/DepMap_Public_22Q2/cellinfor.rds")

###############################################################
###############################################################
### 图复现: ARID1A 和HMGCR 的例子
rm(list = ls())

### 加载R包和数据
library(dplyr)
library(tidyr)
library(tibble)
mutData <- readRDS(file = "data/DepMap_Public_22Q2/mutData.rds")
geneDependency <- readRDS(file = "data/DepMap_Public_22Q2/geneDependency.rds")
cellinfor <- readRDS(file = "data/DepMap_Public_22Q2/cellinfor.rds")

### 细胞系取交集，789
coID <- intersect(rownames(geneDependency),rownames(cellinfor)) %>% 
  intersect(rownames(mutData))

### 提取数据 ARID1A的突变信息，以及HMGCR的Dependency 信息
data <- data.frame(
  DepmapID = coID,
  Tissue = cellinfor[coID,"primary_disease"],
  Mutation = ifelse(mutData[coID,"ARID1A"]==0,"WT","Mut"),
  Dependency = geneDependency[coID,"HMGCR"]
)

## 数据筛选，Mut >=2,WT>=5(跟原文不同)
TissueStatus <- data %>% 
  select(Tissue,Mutation) %>%
  group_by(Tissue,Mutation) %>% 
  summarise(n =n()) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "Mutation",
              values_from = n,
              values_fill = 0) %>% 
  filter(Mut >=2,WT>=5) %>%
  mutate(Sum=Mut+WT) %>% 
  mutate(Mut_Number = sum(Mut),WT_Number=sum(WT)) 



mydata <- data %>% 
  filter(Tissue %in% TissueStatus$Tissue)

## 作图数据整理
Mean_Dep <- mydata %>% 
  filter(Mutation=="Mut") %>% 
  group_by(Tissue) %>% 
  summarise(Sensitivity = mean(Dependency))

Diff_Dep <- mydata %>%
  group_by(Tissue) %>%
  summarize(Difference = mean(Dependency[Mutation == "Mut"]) - mean(Dependency[Mutation == "WT"]))

## 数据合并
PlotData <- merge(Mean_Dep,Diff_Dep,by="Tissue") %>% 
  merge(TissueStatus,by="Tissue")


## 作图展示
library(ggplot2)
library(ggrepel)
library(stringr)
library(ggfun)

ggplot(PlotData, aes(x = Sensitivity, y = Difference)) +
  geom_point(aes(size = Sum), shape = 21, color = "black", fill = "green", alpha = 0.6, stroke = 1.5) + 
  geom_text_repel(aes(label = Tissue), vjust = 1.5, hjust = 1.5) + 
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf,
           fill = "yellow", alpha = 0.1, colour = "green", linetype = "dashed", linewidth = 1) +
  labs(x = "Sensitivity Score for Cell Lines with Mutations", 
       y = "Mutant vs. Wildtype Sensitivity", 
       title = paste0("Mutant (n=", unique(PlotData$Mut_Number),
                      " cell lines) vs. Wildtype (n=",
                      unique(PlotData$WT_Number), " cell lines)")) +
  theme_bw() + 
  guides(size = guide_legend(title = str_wrap("Number of cells", width = 8))) +
  theme(panel.grid.major = element_line(linetype = "dashed"),
        panel.grid.minor = element_line(linetype = "dashed", size = 0.5),
        panel.background = element_blank(),
        legend.position = c(.95, .05),
        legend.justification = c("right", "bottom"),
        legend.box.just = "right",
        legend.margin = margin(4, 4, 4, 4),
        legend.background = element_roundrect(color = "#808080", linetype = 1)
  )


##############################################################
### 换个基因 HMGCS1

data <- data.frame(
  DepmapID = coID,
  Tissue = cellinfor[coID,"primary_disease"],
  Mutation = ifelse(mutData[coID,"ARID1A"]==0,"WT","Mut"),
  Dependency = geneDependency[coID,"HMGCS1"]
)

## 数据筛选，Mut >=2,WT>=5(跟原文不同)
TissueStatus <- data %>% 
  select(Tissue,Mutation) %>%
  group_by(Tissue,Mutation) %>% 
  summarise(n =n()) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "Mutation",
              values_from = n,values_fill = 0) %>% 
  filter(Mut >=2,WT>=5) %>%
  mutate(Sum=Mut+WT) %>% 
  mutate(Mut_Number = sum(Mut),WT_Number=sum(WT)) 



mydata <- data %>% 
  filter(Tissue %in% TissueStatus$Tissue)

## 作图数据整理
Mean_Dep <- mydata %>% 
  filter(Mutation=="Mut") %>% 
  group_by(Tissue) %>% 
  summarise(Sensitivity = mean(Dependency))

Diff_Dep <- mydata %>%
  group_by(Tissue) %>%
  summarize(Difference = mean(Dependency[Mutation == "Mut"]) - mean(Dependency[Mutation == "WT"]))

## 数据合并
PlotData <- merge(Mean_Dep,Diff_Dep,by="Tissue")
PlotData <- merge(PlotData,TissueStatus,by="Tissue")

## 作图展示
library(ggplot2)
library(ggrepel)
library(stringr)

ggplot(PlotData, aes(x = Sensitivity, y = Difference)) +
  geom_point(aes(size = Sum), shape = 21, color = "black", fill = "green", alpha = 0.6, stroke = 1.5) + 
  geom_text_repel(aes(label = Tissue), vjust = 1.5, hjust = 1.5) + 
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf,
           fill = "yellow", alpha = 0.1, colour = "green", linetype = "dashed", size = 1) +
  labs(x = "Sensitivity Score for Cell Lines with Mutations", 
       y = "Mutant vs. Wildtype Sensitivity", 
       title = paste0("Mutant (n=", unique(PlotData$Mut_Number),
                      " cell lines) vs. Wildtype (n=",
                      unique(PlotData$WT_Number), " cell lines)")) +
  theme_bw() + 
  guides(size = guide_legend(title = str_wrap("Number of cells", width = 8))) +
  theme(panel.grid.major = element_line(linetype = "dashed"),
        panel.grid.minor = element_line(linetype = "dashed", size = 0.5),
        panel.background = element_blank(),
        legend.position = c(.95, .05),
        legend.justification = c("right", "bottom"),
        legend.box.just = "right",
        legend.margin = margin(4, 4, 4, 4),
        legend.background = element_roundrect(color = "#808080", linetype = 1)
  )


###################################################
### 提取典型组织单独画图

library(ggplot2)
df = mydata[mydata$Tissue=="Endometrial/Uterine Cancer",]
df$Mutation <- factor(df$Mutation, levels = c("WT", "Mut"),
                      labels = c("ARID1A(wt) \n (n=11)", "ARID1A(mut) \n (n=15)"))
ggplot(df, aes(x = Mutation, y = Dependency, fill = Mutation)) +
  geom_violin(alpha = 0.5) +
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               dotsize = 0.5) +
  labs(y = "Dependency Score")+
  theme_bw() +
  theme(legend.position = "none")

###################################################
### 写成函数

mutPlot <- function(mutGene,targerGene){
  data <- data.frame(
    DepmapID = coID,
    Tissue = cellinfor[coID,"primary_disease"],
    Mutation = ifelse(mutData[coID,mutGene]==0,"WT","Mut"),
    Dependency = geneDependency[coID,targerGene]
  )
  
  ## 数据筛选，Mut >=2,WT>=5(跟原文不同)
  TissueStatus <- data %>% 
    select(Tissue,Mutation) %>%
    group_by(Tissue,Mutation) %>% 
    summarise(n =n()) %>% 
    ungroup() %>% 
    pivot_wider(names_from = "Mutation",
                values_from = n,values_fill = 0) %>% 
    filter(Mut >=2,WT>=5) %>%
    mutate(Sum=Mut+WT) %>% 
    mutate(Mut_Number = sum(Mut),WT_Number=sum(WT)) 
  
  
  
  mydata <- data %>% 
    filter(Tissue %in% TissueStatus$Tissue)
  
  ## 作图数据整理
  Mean_Dep <- mydata %>% 
    filter(Mutation=="Mut") %>% 
    group_by(Tissue) %>% 
    summarise(Sensitivity = mean(Dependency))
  
  Diff_Dep <- mydata %>%
    group_by(Tissue) %>%
    summarize(Difference = mean(Dependency[Mutation == "Mut"]) - mean(Dependency[Mutation == "WT"]))
  
  ## 数据合并
  PlotData <- merge(Mean_Dep,Diff_Dep,by="Tissue")
  PlotData <- merge(PlotData,TissueStatus,by="Tissue")
  
  ## 作图展示
  library(ggplot2)
  library(ggrepel)
  library(stringr)
  
  ggplot(PlotData, aes(x = Sensitivity, y = Difference)) +
    geom_point(aes(size = Sum), shape = 21, color = "black", fill = "green", alpha = 0.6, stroke = 1.5) + 
    geom_text_repel(aes(label = Tissue), vjust = 1.5, hjust = 1.5) + 
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf,
             fill = "yellow", alpha = 0.1, colour = "green", linetype = "dashed", size = 1) +
    labs(x = "Sensitivity Score for Cell Lines with Mutations", 
         y = "Mutant vs. Wildtype Sensitivity", 
         title = paste0("Mutant (n=", unique(PlotData$Mut_Number),
                        " cell lines) vs. Wildtype (n=",
                        unique(PlotData$WT_Number), " cell lines)")) +
    theme_bw() + 
    guides(size = guide_legend(title = str_wrap("Number of cells", width = 8))) +
    theme(panel.grid.major = element_line(linetype = "dashed"),
          panel.grid.minor = element_line(linetype = "dashed", size = 0.5),
          panel.background = element_blank(),
          legend.position = c(.95, .05),
          legend.justification = c("right", "bottom"),
          legend.box.just = "right",
          legend.margin = margin(4, 4, 4, 4),
          legend.background = element_roundrect(color = "#808080", linetype = 1)
    )
  
} 


mutPlot(mutGene = "ARID1A",targerGene = "HMGCR")
mutPlot(mutGene = "TP53",targerGene = "KBTBD2")

### 能不能批量提取？
### GZ07
### 1. 已知突变基因，找靶基因
### 2. 已知感兴趣的基因，找哪个基因突变后对他影响最大
### 批量的指标是什么？

################################################
################################################
### 作者：果子
### 更新时间：2023-12-06
### 微信公众号:果子学生信
### 私人微信：guotosky
### 个人邮箱：hello_guozi@126.com
### 有个问题, BRAF的V600E数据没有

## 数据预处理
## CCLE 究竟如何界定 突变的？
library(dplyr)
library(tidyr)
library(tibble)

#######################################################
### 22Q4V2版本
rm(list = ls())
mutData_raw <- data.table::fread("data/DepMap_Public_20Q4v2/CCLE_mutations.csv",data.table = F)
test <- mutData_raw[mutData_raw$Hugo_Symbol=="BRAF",]
table(mutData_raw$Variant_Classification,mutData_raw$isDeleterious)
table(mutData_raw$Variant_Classification,mutData_raw$Variant_annotation)
table(mutData_raw$Variant_annotation,mutData_raw$isDeleterious)

index <- c("De_novo_Start_OutOfFrame", "Frame_Shift_Del", 
           "Frame_Shift_Ins", "IGR", "Intron", 
           "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", 
           "Start_Codon_Del", "Start_Codon_Ins",
           "Stop_Codon_Del", "Stop_Codon_Ins")

### 1749,19543
mutData <- mutData_raw %>% 
  mutate(isDeleterious = Variant_Classification %in% index) %>% 
  dplyr::select(DepMap_ID,Hugo_Symbol,isDeleterious) %>% 
  mutate(isDeleterious= as.numeric(isDeleterious)) %>%
  group_by(DepMap_ID,Hugo_Symbol) %>% 
  summarise(isDeleterious=sum(isDeleterious)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "Hugo_Symbol",
              values_from = "isDeleterious",
              values_fill = 0) %>% 
  as.data.frame() %>% 
  column_to_rownames("DepMap_ID")

saveRDS(mutData,file = "data/DepMap_Public_20Q4v2/mutData_with_Missense.rds")

### 测试
### 修改前,85,704
### 修改后,118,671
### 原文,133,656

mutData <- readRDS(file = "data/DepMap_Public_20Q4v2/mutData.rds")
#mutData <- readRDS(file = "data/DepMap_Public_20Q4v2/mutData_with_Missense.rds")
geneDependency <- readRDS(file = "data/DepMap_Public_20Q4v2/geneDependency.rds")
cellinfor <- readRDS(file = "data/DepMap_Public_20Q4v2/cellinfor.rds")

### 细胞系取交集，789
coID <- intersect(rownames(geneDependency),rownames(cellinfor)) %>% 
  intersect(rownames(mutData))

### 提取数据 ARID1A的突变信息，以及HMGCR的Dependency 信息
data <- data.frame(
  DepmapID = coID,
  Tissue = cellinfor[coID,"primary_disease"],
  Mutation = ifelse(mutData[coID,"ARID1A"]==0,"WT","Mut"),
  Dependency = geneDependency[coID,"HMGCR"]
)

## 数据筛选，Mut >=2,WT>=5(跟原文不同)
TissueStatus <- data %>% 
  select(Tissue,Mutation) %>%
  group_by(Tissue,Mutation) %>% 
  summarise(n =n()) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "Mutation",
              values_from = n,
              values_fill = 0) %>% 
  filter(Mut >=0,WT>=0) %>%
  mutate(Sum=Mut+WT) %>% 
  mutate(Mut_Number = sum(Mut),WT_Number=sum(WT)) 

#######################################################
### 22Q2版本
rm(list = ls())
mutData_raw <- data.table::fread("data/DepMap_Public_22Q2/CCLE_mutations.csv",data.table = F)
test <- mutData_raw[mutData_raw$Hugo_Symbol=="BRAF",]
table(mutData_raw$Variant_Classification,mutData_raw$isDeleterious)

index <- c("De_novo_Start_OutOfFrame", "Frame_Shift_Del", 
           "Frame_Shift_Ins", "IGR", "Intron", 
           "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", 
           "Start_Codon_Del", "Start_Codon_Ins",
           "Stop_Codon_Del", "Stop_Codon_Ins")

mutData <- mutData_raw %>% 
  mutate(isDeleterious = Variant_Classification %in% index) %>% 
  dplyr::select(DepMap_ID,Hugo_Symbol,isDeleterious) %>% 
  mutate(isDeleterious= as.numeric(isDeleterious)) %>%
  group_by(DepMap_ID,Hugo_Symbol) %>% 
  summarise(isDeleterious=sum(isDeleterious)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "Hugo_Symbol",
              values_from = "isDeleterious",
              values_fill = 0) %>% 
  as.data.frame() %>% 
  column_to_rownames("DepMap_ID")

saveRDS(mutData,file = "data/DepMap_Public_22Q2/mutData_with_Missense.rds")

#mutData <- readRDS(file = "data/DepMap_Public_20Q4v2/mutData.rds")
#mutData <- readRDS(file = "data/DepMap_Public_20Q4v2/mutData_with_Missense.rds")
mutData <- readRDS(file = "data/DepMap_Public_22Q2/mutData_with_Missense.rds")
### 22Q4V2
# geneDependency <- readRDS(file = "data/DepMap_Public_20Q4v2/geneDependency.rds")
# cellinfor <- readRDS(file = "data/DepMap_Public_20Q4v2/cellinfor.rds")
### 22Q2
geneDependency <- readRDS(file = "data/DepMap_Public_22Q2/geneDependency.rds")
cellinfor <- readRDS(file = "data/DepMap_Public_22Q2/cellinfor.rds")

### 细胞系取交集，789
coID <- intersect(rownames(geneDependency),rownames(cellinfor)) %>% 
  intersect(rownames(mutData))

### 提取数据 ARID1A的突变信息，以及HMGCR的Dependency 信息
data <- data.frame(
  DepmapID = coID,
  Tissue = cellinfor[coID,"primary_disease"],
  Mutation = ifelse(mutData[coID,"ARID1A"]==0,"WT","Mut"),
  Dependency = geneDependency[coID,"HMGCR"]
)

## 数据筛选，Mut >=2,WT>=5(跟原文不同)
TissueStatus <- data %>% 
  select(Tissue,Mutation) %>%
  group_by(Tissue,Mutation) %>% 
  summarise(n =n()) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "Mutation",
              values_from = n,
              values_fill = 0) %>% 
  filter(Mut >=0,WT>=0) %>%
  mutate(Sum=Mut+WT) %>% 
  mutate(Mut_Number = sum(Mut),WT_Number=sum(WT)) 

##################################################################
## 新的版本 23Q2 
rm(list = ls())
library(dplyr)
library(tidyr)
library(tibble)
mutData_raw <- data.table::fread("data/DepMap_Public_23Q2/OmicsSomaticMutations.csv",data.table = F)

table(mutData_raw$VariantInfo,mutData_raw$CCLEDeleterious)

test <- mutData_raw[mutData_raw$HugoSymbol=="BRAF",]

index <- c("FRAME_SHIFT_INS", 
           "IN_FRAME_DEL", 
           "MISSENSE", 
           "NONSENSE",
           "NONSTOP", 
           "START_CODON_INS")

## ModelID 就是以前的DepMap_ID
mutData <- mutData_raw %>% 
  mutate(isDeleterious = VariantInfo %in% index) %>% 
  dplyr::select(ModelID,HugoSymbol,isDeleterious) %>% 
  rename(DepMap_ID=ModelID,Hugo_Symbol=HugoSymbol) %>% 
  mutate(isDeleterious= as.numeric(isDeleterious)) %>%
  group_by(DepMap_ID,Hugo_Symbol) %>% 
  summarise(isDeleterious=sum(isDeleterious)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "Hugo_Symbol",
              values_from = "isDeleterious",
              values_fill = 0) %>% 
  as.data.frame() %>% 
  column_to_rownames("DepMap_ID")

saveRDS(mutData,file = "data/DepMap_Public_23Q2/mutData.rds")

### 23Q2 有自己的细胞系信息
cellinfor <- data.table::fread(file = "data/DepMap_Public_23Q2/Model.csv",data.table = F)
colnames(cellinfor)[1] = "DepMap_ID"

cellinfor22Q2 <- readRDS(file = "data/DepMap_Public_22Q2/cellinfor.rds")
cellinfor22Q2 <- cellinfor22Q2[,c("DepMap_ID","lineage","primary_disease")]
cellinfor <- merge(cellinfor,cellinfor22Q2,by = "DepMap_ID")
rownames(cellinfor) <- cellinfor$DepMap_ID
saveRDS(cellinfor,file = "data/DepMap_Public_23Q2/cellinfor.rds")

### 也有自己的geneDependency
## 2.DepMap 基因Dependency 数据
rm(list = ls())
### 1095 vs 1086
geneDependency <- data.table::fread("data/DepMap_Public_23Q2/CRISPRGeneDependency.csv",data.table = F)
test <- geneDependency[1:10,1:10]

## 修改列名
colnames(geneDependency) <- gsub("\\s+\\(\\d+\\)","",colnames(geneDependency))
## 第一列变行名
rownames(geneDependency) <- geneDependency[,1]
geneDependency <- geneDependency[,-1]
## 保存数据
saveRDS(geneDependency,file = "data/DepMap_Public_23Q2/geneDependency.rds")

###########################################################################
rm(list = ls())
mutData <- readRDS(file = "data/DepMap_Public_23Q2/mutData.rds")
cellinfor <- readRDS(file = "data/DepMap_Public_23Q2/cellinfor.rds")
geneDependency <- readRDS(file = "data/DepMap_Public_23Q2/geneDependency.rds")

### 细胞系取交集，1072
coID <- intersect(rownames(geneDependency),rownames(cellinfor)) %>% 
  intersect(rownames(mutData))

### 提取数据 ARID1A的突变信息，以及HMGCR的Dependency 信息
### "lineage","primary_disease","OncotreeCode","OncotreeSubtype","OncotreePrimaryDisease","OncotreeLineage"
Tissue = "lineage"
mutGene = "ARID1A"
targetGene = "HMGCR"
data <- data.frame(
  DepmapID = coID,
  Tissue = cellinfor[coID,Tissue],
  Mutation = ifelse(mutData[coID,mutGene]==0,"WT","Mut"),
  Dependency = geneDependency[coID,targetGene]
)

## 数据筛选，Mut >=3,WT>=5(跟原文不同)
TissueStatus <- data %>% 
  select(Tissue,Mutation) %>%
  group_by(Tissue,Mutation) %>% 
  summarise(n =n()) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "Mutation",
              values_from = n,
              values_fill = 0) %>% 
  filter(Mut >=3,WT>=5) %>%
  mutate(Sum=Mut+WT) %>% 
  mutate(Mut_Number = sum(Mut),WT_Number=sum(WT)) 

mydata <- data %>% 
  filter(Tissue %in% TissueStatus$Tissue)

## 作图数据整理
Mean_Dep <- mydata %>% 
  filter(Mutation=="Mut") %>% 
  group_by(Tissue) %>% 
  summarise(Sensitivity = mean(Dependency))

Diff_Dep <- mydata %>%
  group_by(Tissue) %>%
  summarize(Difference = mean(Dependency[Mutation == "Mut"]) - mean(Dependency[Mutation == "WT"]))

## 数据合并
PlotData <- merge(Mean_Dep,Diff_Dep,by="Tissue") %>% 
  merge(TissueStatus,by="Tissue")


## 作图展示
library(ggplot2)
library(ggrepel)
library(stringr)
library(ggfun)

ggplot(PlotData, aes(x = Sensitivity, y = Difference)) +
  geom_point(aes(size = Sum), shape = 21, color = "black", fill = "green", alpha = 0.6, stroke = 1.5) + 
  geom_text_repel(aes(label = Tissue), vjust = 1.5, hjust = 1.5) + 
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf,
           fill = "yellow", alpha = 0.1, colour = "green", linetype = "dashed", linewidth = 1) +
  labs(x = "Sensitivity Score for Cell Lines with Mutations", 
       y = "Mutant vs. Wildtype Sensitivity", 
       title = paste0("Mutant (n=", unique(PlotData$Mut_Number),
                      " cell lines) vs. Wildtype (n=",
                      unique(PlotData$WT_Number), " cell lines)")) +
  theme_bw() + 
  guides(size = guide_legend(title = str_wrap("Number of cells", width = 8))) +
  theme(panel.grid.major = element_line(linetype = "dashed"),
        panel.grid.minor = element_line(linetype = "dashed", linewidth = 0.5),
        panel.background = element_blank(),
        legend.position = c(.95, .05),
        legend.justification = c("right", "bottom"),
        legend.box.just = "right",
        legend.margin = margin(4, 4, 4, 4),
        legend.background = element_roundrect(color = "#808080", linetype = 1)
  )

#############################################################
#############################################################
### 数据又更新了！2024年11月23日
### 等待确认
rm(list = ls())
library(dplyr)
library(tidyr)
library(tibble)
mutData_raw <- data.table::fread("data/DepMap_Public_23Q4/OmicsSomaticMutations.csv",data.table = F)
names(table(mutData_raw$VariantInfo))
names(table(mutData_raw$VariantType))

### 推荐阅读
### CCLE 如何标记突变
### https://forum.depmap.org/t/how-is-the-isdeleterious-column-in-the-ccle-mutations-csv-file-determined/129


################################################
################################################
### 作者：果子
### 更新时间：2023-12-06
### 微信公众号:果子学生信
### 私人微信：guotosky
### 个人邮箱：hello_guozi@126.com

### 批量找到ARID1A突变后哪个基因受到的影响最大

rm(list = ls())
### 加载数据分析三剑客
library(dplyr)
library(tidyr)
library(tibble)

mutData <- readRDS(file = "data/DepMap_Public_23Q2/mutData.rds")
cellinfor <- readRDS(file = "data/DepMap_Public_23Q2/cellinfor.rds")
geneDependency <- readRDS(file = "data/DepMap_Public_23Q2/geneDependency.rds")
### 细胞系取交集，1072
coID <- intersect(rownames(geneDependency),rownames(cellinfor)) %>% 
  intersect(rownames(mutData))

Tissue = "lineage"
mutGene = "ARID1A"
targetGene = "HMGCR"

# targetGene= "ARID1B"
# targetGene= "SCAP"
# targetGene= "YPEL5"

data <- data.frame(
  DepmapID = coID,
  Tissue = cellinfor[coID,Tissue],
  Mutation = ifelse(mutData[coID,mutGene]==0,"WT","Mut"),
  Dependency = geneDependency[coID,targetGene]
)

## 数据筛选，Mut >=2,WT>=5(跟原文不同)
TissueStatus <- data %>% 
  select(Tissue,Mutation) %>%
  group_by(Tissue,Mutation) %>% 
  summarise(n =n(),.groups = "drop") %>% 
  ungroup() %>% 
  pivot_wider(names_from = "Mutation",
              values_from = n,
              values_fill = 0) %>% 
  filter(Mut >=3,WT>=5) %>%
  mutate(Sum=Mut+WT) %>% 
  mutate(Mut_Number = sum(Mut),WT_Number=sum(WT)) 

mydata <- data %>% 
  filter(Tissue %in% TissueStatus$Tissue)

## 作图数据整理
Mean_Dep <- mydata %>% 
  filter(Mutation=="Mut") %>% 
  group_by(Tissue) %>% 
  summarise(Sensitivity = mean(Dependency))

Diff_Dep <- mydata %>%
  group_by(Tissue) %>%
  summarize(Difference = mean(Dependency[Mutation == "Mut"]) - mean(Dependency[Mutation == "WT"]))

## 数据合并
PlotData <- merge(Mean_Dep,Diff_Dep,by="Tissue") %>% 
  merge(TissueStatus,by="Tissue")

## 作图展示
library(ggplot2)
library(ggrepel)
library(stringr)
library(ggfun)

ggplot(PlotData, aes(x = Sensitivity, y = Difference)) +
  geom_point(aes(size = Sum), shape = 21, color = "black", fill = "green", alpha = 0.6, stroke = 1.5) + 
  geom_text_repel(aes(label = Tissue), vjust = 1.5, hjust = 1.5) + 
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf,
           fill = "yellow", alpha = 0.1, colour = "green", linetype = "dashed", linewidth = 1) +
  labs(x = "Sensitivity Score for Cell Lines with Mutations", 
       y = "Mutant vs. Wildtype Sensitivity", 
       title = paste0("Mutant (n=", unique(PlotData$Mut_Number),
                      " cell lines) vs. Wildtype (n=",
                      unique(PlotData$WT_Number), " cell lines)")) +
  theme_bw() + 
  guides(size = guide_legend(title = str_wrap("Number of cells", width = 8))) +
  theme(panel.grid.major = element_line(linetype = "dashed"),
        panel.grid.minor = element_line(linetype = "dashed", linewidth = 0.5),
        panel.background = element_blank(),
        legend.position = c(.95, .05),
        legend.justification = c("right", "bottom"),
        legend.box.just = "right",
        legend.margin = margin(4, 4, 4, 4),
        legend.background = element_roundrect(color = "#808080", linetype = 1)
  )

## 为了批量,定义一个分数
score = mean(PlotData$Sensitivity)*mean(PlotData$Difference)*100

# score = as.numeric(PlotData$Sensitivity %*% PlotData$Sum)*
#   as.numeric(PlotData$Difference %*% PlotData$Sum)/1000
#######################################################
### 写循环
genelist = colnames(geneDependency)
Tissue = "lineage"
mutGene = "ARID1A"

### 增加一个判断
### mutGene 如果 在所有细胞中不突变
### TissueStatus 没有行数，跳过
results = data.frame()
for(i in 1:length(genelist) ){
  print(i)
  targetGene = genelist[i]
  data <- data.frame(
    DepmapID = coID,
    Tissue = cellinfor[coID,Tissue],
    Mutation = ifelse(mutData[coID,mutGene]==0,"WT","Mut"),
    Dependency = geneDependency[coID,targetGene]
  )
  
  ## 数据筛选，Mut >=2,WT>=5(跟原文不同)
  TissueStatus <- data %>% 
    select(Tissue,Mutation) %>%
    group_by(Tissue,Mutation) %>% 
    summarise(n =n(),.groups = "drop") %>% 
    ungroup() %>% 
    pivot_wider(names_from = "Mutation",
                values_from = n,
                values_fill = 0) %>% 
    filter(Mut >=3,WT>=5) %>%
    mutate(Sum=Mut+WT) %>% 
    mutate(Mut_Number = sum(Mut),WT_Number=sum(WT)) 
  
  mydata <- data %>% 
    filter(Tissue %in% TissueStatus$Tissue)
  
  ## 作图数据整理
  Mean_Dep <- mydata %>% 
    filter(Mutation=="Mut") %>% 
    group_by(Tissue) %>% 
    summarise(Sensitivity = mean(Dependency))
  
  Diff_Dep <- mydata %>%
    group_by(Tissue) %>%
    summarize(Difference = mean(Dependency[Mutation == "Mut"]) - mean(Dependency[Mutation == "WT"]))
  
  ## 数据合并
  PlotData <- merge(Mean_Dep,Diff_Dep,by="Tissue") %>% 
    merge(TissueStatus,by="Tissue")
  
  score = mean(PlotData$Sensitivity)*mean(PlotData$Difference)*100
  # score = as.numeric(PlotData$Sensitivity %*% PlotData$Sum)*
  #   as.numeric(PlotData$Difference %*% PlotData$Sum)/1000
  
  results[i,1] = mutGene
  results[i,2] = targetGene
  results[i,3] = mean(PlotData$Sensitivity)
  results[i,4] = mean(PlotData$Difference)
  results[i,5] = score
}

colnames(results) = c("mutGene","targetGene","mena_Sensitivity","mean_Diff","Score")
saveRDS(results,file = "output/mut2target_results.rds")
#saveRDS(results,file = "output/mut2target_results_dot.rds")

results <- readRDS(file = "output/mut2target_results.rds")

### https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3954704/
### GZ07 批量课程


################################################
################################################
### 作者：果子
### 更新时间：2023-12-13
### 微信公众号:果子学生信
### 私人微信：guotosky
### 个人邮箱：hello_guozi@126.com

### 批量找到ARID1A突变后哪个基因受到的影响最大
### 加上细胞类型

rm(list = ls())
### 加载数据分析三剑客
library(dplyr)
library(tidyr)
library(tibble)

mutData <- readRDS(file = "data/DepMap_Public_23Q2/mutData.rds")
cellinfor <- readRDS(file = "data/DepMap_Public_23Q2/cellinfor.rds")
geneDependency <- readRDS(file = "data/DepMap_Public_23Q2/geneDependency.rds")
### 细胞系取交集，1072
coID <- intersect(rownames(geneDependency),rownames(cellinfor)) %>% 
  intersect(rownames(mutData))

Tissue = "lineage"
mutGene = "ARID1A"


### 优化流程
## 数据筛选，Mut >=2,WT>=5(跟原文不同)
TissueStatus <- data.frame(
  Tissue = cellinfor[coID,Tissue],
  Mutation = ifelse(mutData[coID,mutGene]==0,"WT","Mut")
) %>% 
  select(Tissue,Mutation) %>%
  group_by(Tissue,Mutation) %>% 
  summarise(n =n(),.groups = "drop") %>% 
  ungroup() %>% 
  pivot_wider(names_from = "Mutation",
              values_from = n,
              values_fill = 0) %>% 
  filter(Mut >=3,WT>=5) %>%
  mutate(Sum=Mut+WT) %>% 
  mutate(Mut_Number = sum(Mut),WT_Number=sum(WT)) 


targetGene = "HMGCR"
# targetGene= "ARID1B"
# targetGene= "PIK3CA"
targetGene= "GATA3"

data <- data.frame(
  DepmapID = coID,
  Tissue = cellinfor[coID,Tissue],
  Mutation = ifelse(mutData[coID,mutGene]==0,"WT","Mut"),
  Dependency = geneDependency[coID,targetGene]
)

mydata <- data %>% 
  filter(Tissue %in% TissueStatus$Tissue)

## 作图数据整理
Mean_Dep <- mydata %>% 
  filter(Mutation=="Mut") %>% 
  group_by(Tissue) %>% 
  summarise(Sensitivity = mean(Dependency))

Diff_Dep <- mydata %>%
  group_by(Tissue) %>%
  summarize(Difference = mean(Dependency[Mutation == "Mut"]) - mean(Dependency[Mutation == "WT"]))

## 数据合并
PlotData <- merge(Mean_Dep,Diff_Dep,by="Tissue") %>% 
  merge(TissueStatus,by="Tissue")

## 作图展示
library(ggplot2)
library(ggrepel)
library(stringr)
library(ggfun)

ggplot(PlotData, aes(x = Sensitivity, y = Difference)) +
  geom_point(aes(size = Sum), shape = 21, color = "black", fill = "green", alpha = 0.6, stroke = 1.5) + 
  geom_text_repel(aes(label = Tissue), vjust = 1.5, hjust = 1.5) + 
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf,
           fill = "yellow", alpha = 0.1, colour = "green", linetype = "dashed", linewidth = 1) +
  labs(x = "Sensitivity Score for Cell Lines with Mutations", 
       y = "Mutant vs. Wildtype Sensitivity", 
       title = paste0("Mutant (n=", unique(PlotData$Mut_Number),
                      " cell lines) vs. Wildtype (n=",
                      unique(PlotData$WT_Number), " cell lines)")) +
  theme_bw() + 
  guides(size = guide_legend(title = str_wrap("Number of cells", width = 8))) +
  theme(panel.grid.major = element_line(linetype = "dashed"),
        panel.grid.minor = element_line(linetype = "dashed", linewidth = 0.5),
        panel.background = element_blank(),
        legend.position = c(.95, .05),
        legend.justification = c("right", "bottom"),
        legend.box.just = "right",
        legend.margin = margin(4, 4, 4, 4),
        legend.background = element_roundrect(color = "#808080", linetype = 1)
  )

## 为了批量,定义总分数
## 以及每个类型的分数
score_all = mean(PlotData$Sensitivity)*mean(PlotData$Difference)*100
score_celltype = PlotData$Sensitivity*PlotData$Difference*100

#######################################################
### 写循环
genelist = colnames(geneDependency)
Tissue = "lineage"
mutGene = "ARID1A"

TissueStatus <- data.frame(
  Tissue = cellinfor[coID,Tissue],
  Mutation = ifelse(mutData[coID,mutGene]==0,"WT","Mut")
) %>% 
  select(Tissue,Mutation) %>%
  group_by(Tissue,Mutation) %>% 
  summarise(n =n(),.groups = "drop") %>% 
  ungroup() %>% 
  pivot_wider(names_from = "Mutation",
              values_from = n,
              values_fill = 0) %>% 
  filter(Mut >=3,WT>=5) %>%
  mutate(Sum=Mut+WT) %>% 
  mutate(Mut_Number = sum(Mut),WT_Number=sum(WT)) 


results = data.frame()
for(i in 1:length(genelist) ){
  print(i)
  targetGene = genelist[i]
  data <- data.frame(
    DepmapID = coID,
    Tissue = cellinfor[coID,Tissue],
    Mutation = ifelse(mutData[coID,mutGene]==0,"WT","Mut"),
    Dependency = geneDependency[coID,targetGene]
  )

  mydata <- data %>% 
    filter(Tissue %in% TissueStatus$Tissue)
  
  ## 作图数据整理
  Mean_Dep <- mydata %>% 
    filter(Mutation=="Mut") %>% 
    group_by(Tissue) %>% 
    summarise(Sensitivity = mean(Dependency))
  
  Diff_Dep <- mydata %>%
    group_by(Tissue) %>%
    summarize(Difference = mean(Dependency[Mutation == "Mut"]) - mean(Dependency[Mutation == "WT"]))
  
  ## 数据合并
  PlotData <- merge(Mean_Dep,Diff_Dep,by="Tissue") %>% 
    merge(TissueStatus,by="Tissue")
  
  score_all = mean(PlotData$Sensitivity)*mean(PlotData$Difference)*100
  score_celltype = PlotData$Sensitivity*PlotData$Difference*100
  
  results[i,1] = mutGene
  results[i,2] = targetGene
  results[i,3] = mean(PlotData$Sensitivity)
  results[i,4] = mean(PlotData$Difference)
  results[i,5:(4+nrow(PlotData))] = score_celltype
  results[i,(5+nrow(PlotData))] = score_all
}

colnames(results) = c("mutGene","targetGene","mena_Sensitivity","mean_Diff",TissueStatus$Tissue,"Score")
saveRDS(results,file = "output/mut2target_results_celltype.rds")

results <- readRDS(file = "output/mut2target_results_celltype.rds")

### 未来迭代版本
### 考虑突变基因是否适配流程
### 采用并行化提速, 50

### 新的话题
### mut target_gene：
### target_gene mut 

### https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3954704/
