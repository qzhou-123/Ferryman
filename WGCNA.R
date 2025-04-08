# 0. 清空环境并设置工作目录
rm(list = ls())
setwd("/Users/zhouquan/Desktop/WGCNA")

# 1. 加载必要的库
library(WGCNA)
library(data.table)
library(dplyr)
library(ggplot2)
library(GEOquery)
library(limma)
library(TCGAbiolinks)
library(clusterProfiler)
library(org.Hs.eg.db)

########################### 获取数据 ###########################
# 2. 设置GEO数据集编号
gse_number = "GSE248830"

# 3. 下载并加载GEO数据集
eSet <- getGEO(gse_number, destdir = '.', getGPL = FALSE)
eSet = eSet[[1]]

# 4. 提取表达数据并进行标准化
exp <- exprs(eSet)
if (TRUE) { exp = limma::normalizeBetweenArrays(exp) }

# 5. 查看标准化后的表达数据
head(exp)
nrow(exp)
boxplot(exp)

# 6. 创建包含组别信息的数据框并选择子集表达数据
coldata <- data.frame(group = factor(rep(c('treat', 'control'), each = 11)))
df <- exp[, c(34:44, 23:33)] %>% as.data.frame()
write.table(df, "gene_matrix.csv", quote = FALSE, row.names = TRUE, sep = "\t")
head(df)
dim(df)

########################### 差异分析 ###########################
# 7. 创建组别信息的因子变量并生成设计矩阵
f <- factor(c(rep("treat", 11), rep("control", 11)))
desigN <- model.matrix(~ 0 + f)
colnames(desigN) <- levels(f)
print(desigN)

# 8. 拟合线性模型并进行对比分析
fit = lmFit(df, desigN)
contrast.matrix <- makeContrasts(treat - control, levels = desigN)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# 9. 获取差异表达基因的顶表并筛选差异表达基因
tempOutput <- topTable(fit2, coef = 1, n = Inf, adjust = "BH")
up <- subset(tempOutput, P.Value < 0.05 & logFC > 1)
down <- subset(tempOutput, P.Value < 0.05 & logFC < -1)
diff <- subset(tempOutput, P.Value < 0.05 & abs(logFC) > 1)
write.table(diff, "diff.table", quote = FALSE, row.names = TRUE, sep = "\t")

########################### WGCNA分析 ###########################
# 10. 转置数据并选择高方差基因
data <- as.data.frame(t(df))
m.vars <- apply(data, 2, var)
expro.upper <- data[, which(m.vars > quantile(m.vars, probs = seq(0, 1, 0.25))[4])]
data <- as.data.frame(t(expro.upper))
dim(data)

# 11. 检查样本和基因质量
gsg <- goodSamplesGenes(data, verbose = 3)
gsg$allOK
dim(data)

# 12. 选择软阈值幂并生成图表
powers <- seq(1, 7.9, by = 0.5)
dat <- as.data.frame(t(data))
sft <- pickSoftThreshold(dat, powerVector = powers, verbose = 5)

par(mfrow = c(1, 2))
cex1 <- 0.95
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit,signed R^2", type = "n", main = "Scale independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels = powers, cex = cex1, col = "red")
abline(h = 0.95, col = "red")

plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n", main = "Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = cex1, col = "red")

# 13. 绘制样本聚类树
sample_tree <- hclust(dist(t(data)), method = "average")
par(cex = 0.7)
par(mar = c(1, 4, 2, 0))
plot(sample_tree, main = "Sample clustering", cex.axis = 1.5, cex.main = 2)

# 14. 进行模块检测
cor <- WGCNA::cor
net <- blockwiseModules(dat, power = 6, maxBlockSize = 1000, TOMType = 'unsigned', minModuleSize = 10, reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = TRUE, saveTOMFileBase = "drought", verbose = 3)
cor <- stats::cor
table(net$colors)

# 15. 提取模块标签和颜色
moduleLabels <- net$colors
moduleColors <- labels2colors(moduleLabels)

# 16. 绘制模块聚类树及模块颜色
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]], "Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)

# 17. 提取基因树
geneTree <- net$dendrograms[[1]]
moduleLabelsAutomatic = net$colors
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
moduleColorsWW = moduleColorsAutomatic

# 18. 计算模块特征基因
nGenes <- ncol(dat)
nSamples <- nrow(dat)
MEsO = moduleEigengenes(dat, moduleColorsWW)$eigengenes
MEsWW = orderMEs(MEsO)

# 19. 计算模块-性状相关性及其p值
modTraitCor = cor(MEsWW, desigN, use = "p")
modTRaitp = corPvalueStudent(modTraitCor, nSamples)

# 20. 生成相关性和p值矩阵的文本表示并绘制热图
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTRaitp, 1), ")", sep = "")
dim(textMatrix) <- dim(modTraitCor)

png("heatmap.png", width = 1500, height = 1000)
par(mar = c(8, 10, 4, 10))
customColors <- colorRampPalette(c("blue", "white", "red"))(50)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(desigN), yLabels = names(MEsWW), ySymbols = names(MEsWW), colorLabels = FALSE, colors = customColors, textMatrix = textMatrix, setStdMargins = FALSE, cex.lab = 1.8, cex.text = 1.5, zlim = c(-1, 1), legendSpace = 0.3, xColorWidth = 0.01, yColorWidth = 0.01, main = "Module-trait relationships")
dev.off()

# 21. 计算模块成员(MM)值和基因显著性(GS)值
modNames <- substring(names(MEsWW), 3)
geneModuleMembership <- as.data.frame(cor(dat, MEsWW, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) <- paste("MM", modNames, sep = "")
names(MMPvalue) <- paste("p.MM", modNames, sep = "")
geneTraitSignificance <- as.data.frame(cor(dat, desigN, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) <- paste("GS.", colnames(desigN), sep = "")
names(GSPvalue) <- paste("p.GS.", colnames(desigN), sep = "")

# 22. 选择"blue"模块并计算MM值和GS值
module <- "blue"
column <- match(module, modNames)
blue_moduleGenes <- names(net$colors)[which(moduleColors == module)]
MM <- abs(geneModuleMembership[blue_moduleGenes, column])
GS <- abs(geneTraitSignificance[blue_moduleGenes, 1])

# 23. 绘制模块成员和基因显著性关系图
png("Module_blue_membership_gene_significance.png", width = 800, height = 600)
par(mfrow = c(1, 1))
verboseScatterplot(MM, GS, xlab = paste("Module Membership in", module, "module"), ylab = "Gene significance for Basal", main = paste("Module membership vs. gene significance\n"), cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

# 24. 选择蓝色模块中具有高MM值和GS值的关键基因并保存到CSV文件
blue_hub <- blue_moduleGenes[(GS > 0.4 & MM > 0.7)]
write.csv(blue_hub, 'blue_hub_gene.csv', quote = TRUE)

########################### 富集分析 ###########################
# 25. 读取蓝色模块的关键基因并去重
brown_hub_gene <- read.csv('blue_hub_gene.csv', header = FALSE)
brown_hub_gene <- brown_hub_gene[-1, ]
brown_hub_gene <- brown_hub_gene[, -1]
gene <- unique(brown_hub_gene)

# 26. 基因符号转换为ENTREZ ID
sig_DP_entrezId <- mapIds(x = org.Hs.eg.db, keys = gene, keytype = "SYMBOL", column = "ENTREZID")
sig_DP_entrezId <- na.omit(sig_DP_entrezId)

# 27. 进行GO生物过程(BP)富集分析
go_bp <- enrichGO(gene = sig_DP_entrezId, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", pvalueCutoff = 0.5, qvalueCutoff = 0.5, readable = TRUE)

# 28. 绘制GO富集分析的柱状图
barplot(go_bp)

