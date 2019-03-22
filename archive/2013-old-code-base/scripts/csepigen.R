library(LiblineaR)
rr = read.table('hg19_H1hesc_chr1-2-3_200_Rep1_gene_ext3kb_Genc3c_strict_51.bed', header=T, comment.char='', sep='\t')
names(rr)
properties <- c('X.chrom', 'start', 'end', 'id', 'strand', 'length', 'RPKM', 'ensembl_gene_id', 'class', 'label', 'assembly', 'element', 'gene_name', 'gene_type', 'havana_gene_id', 'source', 'seq')
features <- names(rr) !%in% properties
features <- names(rr) %in% properties
features
features <- names(rr)[!(names(rr) %in% properties)]
features
summary(rr[,features])
?LiblineaR
?LiblineaR
svm <- LiblineaR(data=rr[,features], labels=rr$label, type=1, cost=1)
svm
?LiblineaR
svm <- LiblineaR(data=rr[,features], labels=rr$label, type=1, cost=1, cross=10)
svm
rr = read.table('hg19_H1hesc_chr1-2-3_200_Rep1_gene_ext3kb_Genc3c_normal_50.bed', header=T, comment.char='', sep='\t')
svm <- LiblineaR(data=rr[,features], labels=rr$label, type=1, cost=1, cross=10)
svm
library(randomForest)
rf <- randomForest(x=rr[,features], y=as.factor(rr$label), ntree=500, importance=T)
rf
rr = read.table('hg19_H1hesc_chr1-2-3_200_Rep1_gene_ext3kb_Genc3c_strict_51.bed', header=T, comment.char='', sep='\t')
rf <- randomForest(x=rr[,features], y=as.factor(rr$label), ntree=500, importance=T)
rf
varImpPlot(rf, sort=True, n.var=50)
varImpPlot(rf, sort=TRUE, n.var=50)
varImpPlot(rf, sort=TRUE, n.var=50, type=2, class=1)
varImpPlot(rf, sort=TRUE, n.var=50, type=2, class=as.factor(1))
varImpPlot(rf, sort=TRUE, n.var=50, class=as.factor(1))
varImpPlot(rf, sort=TRUE, n.var=50, class=as.factor(-1))
varImpPlot(rf, sort=TRUE, n.var=50, class=as.factor(2))
rf
varImpPlot(rf, sort=TRUE, n.var=50, class=1)
varImpPlot(rf, sort=TRUE, n.var=50, class=-1)
plot(margin(rf))
partialPlot(rf)
partialPlot(rf, 'M00803')
partialPlot(rf, x.var='M00803')
partialPlot(rf, rr[,features], x.var='M00803')
partialPlot(rf, rr[,features], x.var='M00803', plot=TRUE)
partialPlot(rf, rr[,features], x.var='num_CGIovl', plot=TRUE)
partialPlot(rf, rr[,features], x.var='num_CGIovl', plot=TRUE, which.class=-1)
partialPlot(rf, rr[,features], x.var='num_CGIovl', plot=TRUE, which.class=1)
varUsed(rf, FALSE, TRUE)
vc <- varUsed(rf, FALSE, TRUE)
vc
names(vc) <- features
vc
max(vc)
which(max(vc))
which(vc[vc > 800])
dim(vc)
vc[vc>800]
vc[vc>500]
vc[vc>400]
rf.51 <- rf
rr = read.table('hg19_H1hesc_chr1-2_Rep1_200_TSS_ext5kb_Genc3c_normal_70.bed', comment.char='', header=T, sep='\t')
names(rf)
names(rr)
properties.gene <- c('X.chrom', 'start', 'end', 'id', 'strand', 'length', 'RPKM', 'ensembl_gene_id', 'class', 'label', 'assembly', 'element', 'gene_name', 'gene_type', 'havana_gene_id', 'source', 'seq')
properties.tss <- c('X.chrom', 'start', 'end', 'id', 'strand', 'length', 'RPKM', 'ensembl_gene_id', 'ensembl_transcript_id', 'havana_transcript_id', 'transcript_type', 'class', 'label', 'assembly', 'element', 'gene_name', 'gene_type', 'havana_gene_id', 'source', 'seq')
features.tss <- names(rr)[!(names(rr) %in% properties.tss)]
features.tss
rf.tss <- randomForest(x=rr[,features.tss], y=as.factor(rr$label), ntree=1000, importance=T)
rf.tss
row.names(rr)
summary(rr$HIGH)
summary(rr$class)
train.high <- sample(row.names(rr[class=='HIGH',]), 4000, replace=F)
train.high <- sample(row.names(rr[rr$class=='HIGH',]), 4000, replace=F)
train.low <- sample(row.names(rr[rr$class=='LOW',]), 4000, replace=F)
train.idx <- c(train.high, train.low)
rf.tss.bal <- randomForest(x=rr[train.idx,features.tss], y=as.factor(rr[train.idx,]$label), ntree=500, importance=T)
rf.tss.bal
varImpPlot(rf, sort=T, n.var=20, type=2)
varImpPlot(rf, sort=T, n.var=20, type=2, class=-1)
varImpPlot(rf, sort=T, n.var=20, class=-1)
varImpPlot(rf, sort=T, n.var=20, class=1)
varImpPlot(rf, sort=T, n.var=20, type=2)
?LiblineaR
?LiblineaR
svm.tss.bal <- Liblinear(data=rr[train.idx,features], labels=rr[train.idx,]$label, bias=T, type=1, cost=1, cross=10)
svm.tss.bal <- LiblineaR(data=rr[train.idx,features], labels=rr[train.idx,]$label, bias=T, type=1, cost=1, cross=10)
svm.tss.bal
svm.tss.bal <- LiblineaR(data=rr[train.idx,features], labels=rr[train.idx,]$label, bias=T, type=1, cost=1)
library(caret)
svm.tss.bal
max(svm.tss.bal$W)
min(svm.tss.bal$W)
which(svm.tss.bal$w[max(svm.tss.bal$W)])
which(svm.tss.bal$w == max(svm.tss.bal$W)])
which(svm.tss.bal$w == max(svm.tss.bal$W))
max(svm.tss.bal$W)
names(max(svm.tss.bal$W))
str(max(svm.tss.bal$W))
?is.ma
?is.max
which(svm.tss.bal$W > 0.638
)
svm.tss.bal$W[728]
svm.tss.bal$W[1,728]
rr.feat.scaled <- scale(rr[,features],T,T)
which(is.na(rr.feat.scaled, arr.ind=T)
)
which(is.na(rr.feat.scaled), arr.ind=T)
?LiblineaR
?var
vars <- apply(rr[,features],2,var)
vars
max(vars)
min(vars)
?cor
cors <- cor(rr[,features])
cors
dim(cors)
which(cors > 0.3, arr.ind=T)
high.cor <- which(cors > 0.6, arr.ind=T)
dim(high.cor)
head(high.cor)
cors[1,1]
uns <- unique(row.names(cors))
length(uns)
uns <- unique(row.names(high.cor))
length(uns)
hist(rr[,'CTA'])
hist(rr[,'CTA'], plot=F)
var(rr[,'CTA'])
hist(rr[,'CTA'], plot=T)
hist(rr[,'X12comb_1'], plot=T)
hist(rr[,'X12comb_1'], plot=F)
var(rr[,'X12comb_1'])
dim(rr[,'X12comb_1'])
length(rr[,'X12comb_1'])
which(vars > 0 && var < 0.0001)
which(vars > 0 && var < 0.001)
which(vars > 0 && var < 0.01)
which(vars > 0 && var < 0.1)
which(vars > 0 & var < 0.1)
which(vars > 0 && vars < 0.1)
which(vars > 0 && vars < 0.01)
which(vars > 0 && vars < 0.001)
which(vars > 0 &1 vars < 0.1)
which(vars > 0 & vars < 0.1)
var(rr[,'num_CGIovl'])
hist(rr[,'num_CGIovl'])
table(rr[,'num_CGIovl'])
max(table(rr[,'num_CGIovl']))
?var
foo <- var(rr[,features])
dim(foo)
var.test
vars
str(vars)
which(vars == 0)
names(which(vars == 0))
90 / 100
table(rr[,'num_CGIovl'])
names(which(max(table(rr[,'num_CGIovl']))))
names(which(max(table(rr[,'num_CGIovl'])))
)
names(max(table(rr[,'num_CGIovl'])))
max.count <- apply(rr[,325:330], 2, FUN(X) {max(table(X))})
max.count <- apply(rr[,325:330], 2, FUN(X) max(table(X)))
max.count <- apply(rr[,325:330], 2, max(table(X)))
max.count <- apply(rr[,325:330], 2, FUN(X) max(table(X)))
max.count <- apply(rr[,325:330], 2, FUN(X) (max(table(X))))
?apply
max.count <- apply(rr[,325:330], 2, function(X) (max(table(X))))
max.count
names(max.count)
rem <- names(which(max.count > 100, arr.ind=F))
rem
removeNearZeroVarCorPred <- function(dataframe, thresholdfreq)
{
print(dim(dataframe))
# remove all features with zero variance
var.col <- apply(dataframe, 2, var)
zeroVarPred <- names(which(var.col == 0))
reduced.data <- dataframe[,-zeroVarPred]
# remove all features that are close to zero variance
# frequency of most common value >= threshold freq
threshold.abs <- nrow(dataframe) * (thresholdfreq / 100)
max.count <- apply(reduced.data, 2, FUN(X) {max(table(X))} )
nearZeroPred <- names(which(max.count >= threshold.abs, arr.ind=F))
reduced.data <- reduced.data[, -nearZeroPred]
print(dim(reduced.data))
return(reduced.data)
removeNearZeroVarCorPred <- function(dataframe, thresholdfreq)
{
print(dim(dataframe))
# remove all features with zero variance
var.col <- apply(dataframe, 2, var)
zeroVarPred <- names(which(var.col == 0))
reduced.data <- dataframe[,-zeroVarPred]
# remove all features that are close to zero variance
# frequency of most common value >= threshold freq
threshold.abs <- nrow(dataframe) * (thresholdfreq / 100)
max.count <- apply(reduced.data, 2, FUN(X) {max(table(X))} )
nearZeroPred <- names(which(max.count >= threshold.abs, arr.ind=F))
reduced.data <- reduced.data[, -nearZeroPred]
print(dim(reduced.data))
return(reduced.data)
removeNearZeroVarCorPred <- function(dataframe, thresholdfreq)
{
print(dim(dataframe))
# remove all features with zero variance
var.col <- apply(dataframe, 2, var)
zeroVarPred <- names(which(var.col == 0))
reduced.data <- dataframe[,-zeroVarPred]
# remove all features that are close to zero variance
# frequency of most common value >= threshold freq
threshold.abs <- nrow(dataframe) * (thresholdfreq / 100)
max.count <- apply(reduced.data, 2, FUN(X) {max(table(X))} )
nearZeroPred <- names(which(max.count >= threshold.abs, arr.ind=F))
reduced.data <- reduced.data[, -nearZeroPred]
print(dim(reduced.data))
return(reduced.data)
source('~/work/code/sandbox/csepigen/scripts/csepigen_ML.R')
source('~/work/code/sandbox/csepigen/scripts/csepigen_ML.R')
source('~/work/code/sandbox/csepigen/scripts/csepigen_ML.R')
redrr <- removeNearZeroVarCorPred(rr[,features], 90)
source('~/work/code/sandbox/csepigen/scripts/csepigen_ML.R')
redrr <- removeNearZeroVarCorPred(rr[,features], 90)
redrr.scaled <- scale(redrr, T, T)
is.na(redrr.scaled)
sum(is.na(redrr.scaled))
svm.tss.bal.scaled <- LiblineaR(data=redrr[train.idx,], labels=rr[train.idx,]$label, bias=T, type=1, cost=1, cross=10)
svm.tss.bal.scaled
rf.tss.bal.scaled <- randomForest(x=redrr[train.idx,], y=as.factor(redrr[train.idx,]$label), ntree=500, importance=T)
rf.tss.bal.scaled <- randomForest(x=redrr[train.idx,], y=as.factor(rr[train.idx,]$label), ntree=500, importance=T)
rf.tss.bal.scaled
library(e1071)
svm.tss.bal.scaled.rbf <- svm(x=redrr[train.idx,], y=as.factor(rr[train.idx,]$label), scale=F, cachesize=200, cross=10, kernel=radial, probability=T)
svm.tss.bal.scaled.rbf <- svm(x=redrr[train.idx,], y=as.factor(rr[train.idx,]$label), scale=F, cachesize=200, cross=10, kernel='radial basis', probability=T)
svm.tss.bal.scaled.rbf <- svm(x=redrr[train.idx,], y=as.factor(rr[train.idx,]$label), scale=F, cachesize=200, cross=10, kernel=rbf, probability=T, cost=1)
svm.tss.bal.scaled.rbf <- svm(x=redrr[train.idx,], y=as.factor(rr[train.idx,]$label), scale=F, cachesize=200, cross=10, kernel="radial", probability=T, cost=1)
svm.tss.bal.scaled.rbf
str(svm.tss.bal.scaled.rbf)
svm.tss.bal.scaled.rbf$tot_accuracy
svm.tss.bal.scaled.rbf$tot.accuracy
plot(svm.tss.bal.scaled.rbf)
plot(svm.tss.bal.scaled.rbf, redrr, )
dc <- svm.tss.bal.scaled.rbf$decision.values
dc
max(dc)
min(dc)
plot(dc)
which(dc == max(dc))
dc[685]
names(dc[685])
attr(dc[685])
str(dc[685])
dc[685]
str(dc)
attr(dc[685],"dimnames")
dc[685]$

dc[685][[1]]
dc[685][[2]]
dc[685][1]
dc[685][2]
dc[685][3]
dimnames(dc[685])
dimnames(dc)
dimnames(dc)[[1]]
dimnames(dc)[[1]][685]
redrr.scaled[685,]
svm.tss.bal.scaled.rbf
str(svm.tss.bal.scaled.rbf)
rbf.feat.w <- t(svm.tss.bal.scaled.rbf$coefs) %*% svm.tss.bal.scaled.rbf$SV
rbf.feat.w
rbf.feat.dc <- rbf.feat.w %*% t(as.matrix(redrr.scaled[train.idx,])) - svm.tss.bal.scaled.rbf$rho
rbf.feat.dc <- t(rbf.feat.w %*% t(as.matrix(redrr.scaled[train.idx,]))) - svm.tss.bal.scaled.rbf$rho
dim(redrr.scaled)
length(train.idx)

dim(rbf.feat.w)
redrr.scaled[train.idx,[
redrr.scaled[train.idx,]
head(redrr)
head(redrr.scaled)
dim(redrr.scaled)
head(train.idx)
redrr.scaled[train.idx,]
redrr.scaled[as.numeric(train.idx),]
savehistory('csepigen.R')
