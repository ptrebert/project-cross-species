#!/TL/opt/bin/Rscript

#library(e1071)
library(caret)
#library(randomForest)
library(matrixStats)
library(doMC)

ncpus <- 64
server_name <- as.character(Sys.info()['nodename'])
registerDoMC(cores=ncpus)

source('/home/pebert/work/code/sandbox/csepigen/trunk/scripts/generate_report.R')

get_regexps <- function(type, species)
{
	re.keep_all <- 'keep_all'
	re.no_peaks <- '(^X[0-9]+_peak$)'
	re.no_tfbs <- '(^(M0.+|X[12]{1,2}comb_[0-9]+)$)'
	re.no_signal <- '(^X[0-9]+_p[24689]{1}[05]{1}$)'
	re.no_ext_seq <- '(^(perc.+|num.+|obsexp_CpG)$)'
	re.baseline_seq <- paste(re.no_peaks, re.no_tfbs, re.no_signal, re.no_ext_seq, sep='|')
	re.baseline_ext <- paste(re.no_peaks, re.no_tfbs, re.no_signal, sep='|')
	if (species == 'human')
	{
		hg19 <- '(^(X3|X4|X5|X6|X7|X8)_p[24689]{1}[05]{1}$)'
		mm9_to_hg19 <- '(^(X20|X21|X22|X23|X24|X72)_p[24689]{1}[05]{1}$)'
		hg19_mm9_hg19 <- '(^(X33|X34|X35|X36|X37|X38)_p[24689]{1}[05]{1}$)'
		mm9_hg19_mm9_hg19 <- '(^(X45|X46|X48|X49|X50|X74)_p[24689]{1}[05]{1}$)'
		if (type == 'expression')
		{
			re.noTFBS_noPeak_HG19 <- paste(re.no_peaks, re.no_tfbs, mm9_to_hg19, hg19_mm9_hg19, mm9_hg19_mm9_hg19, sep='|')
			re.noTFBS_noPeak_1xMM9 <- paste(re.no_peaks, re.no_tfbs, hg19, hg19_mm9_hg19, mm9_hg19_mm9_hg19, sep='|')
			re.noTFBS_noPeak_2xHG19 <- paste(re.no_peaks, re.no_tfbs, hg19, mm9_to_hg19, mm9_hg19_mm9_hg19, sep='|')
			re.noTFBS_noPeak_3xMM9 <- paste(re.no_peaks, re.no_tfbs, hg19, mm9_to_hg19, hg19_mm9_hg19, sep='|')
			re.noTFBS_noSignal <- paste(re.no_tfbs, re.no_signal, sep='|')
			expr_features_human <- list('all'=re.keep_all,
							'no_epigen'=paste(re.no_peaks, re.no_signal, sep='|'),
							'no_tfbs_epigen'=re.baseline_ext,
							'only_seq'=re.baseline_seq,
							'hg19_noTFBS_noPeak'=re.noTFBS_noPeak_HG19,
							'mm91x_noTFBS_noPeak'=re.noTFBS_noPeak_1xMM9,
							'hg192x_noTFBS_noPeak'=re.noTFBS_noPeak_2xHG19,
							'mm93x_noTFBS_noPeak'=re.noTFBS_noPeak_3xMM9,
							'noTFBS_noSignal'=re.noTFBS_noSignal,
							'noSignal'=re.no_signal
						)

			#expr_features_human <- list('only_seq_reptest'=re.baseline_seq)

			return(expr_features_human)
		}
		else
		{
			stopifnot(type == 'histone')
			histone_features_human <- list('all'=re.keep_all,
								'no_epigen'=paste(re.no_peaks, re.no_signal, sep='|'),
								'only_seq'=re.baseline_seq,
								'no_tfbs_epigen'=re.baseline_ext,
								'hg19_noTFBS_noPeak'=re.noTFBS_noPeak_HG19,
								'mm91x_noTFBS_noPeak'=re.noTFBS_noPeak_1xMM9,
								'hg192x_noTFBS_noPeak'=re.noTFBS_noPeak_2xHG19
								)
			return(histone_features_human)
		}
	}
	else
	{
		stopifnot(species == 'mouse')
		mm9 <- '(^(X1|X13|X14|X15|X16|X17)_p[24689]{1}[05]{1}$)'
		hg19_to_mm9 <- '(^(X26|X27|X28|X29|X30|X31)_p[24689]{1}[05]{1}$)'
		mm9_hg19_mm9 <- '(^(X39|X40|X41|X42|X43|X73)_p[24689]{1}[05]{1}$)'
		if (type == 'histone')
		{
			re.noTFBS_MM9 <- paste(re.no_tfbs, hg19_to_mm9, mm9_hg19_mm9, sep='|')
			re.noTFBS_1xHG19 <- paste(re.no_tfbs, mm9, mm9_hg19_mm9, sep='|')
			re.noTFBS_2xMM9 <- paste(re.no_tfbs, mm9, hg19_to_mm9, sep='|')

			histone_features_mouse <- list('all'=re.keep_all,
								'no_epigen'=paste(re.no_peaks, re.no_signal, sep='|'),
								'only_seq'=re.baseline_seq,
								'no_tfbs_epigen'=re.baseline_ext,
								'mm9_no_tfbs'=re.noTFBS_MM9,
								'hg191x_no_tfbs'=re.noTFBS_1xHG19,
								'mm92x_no_tfbs'=re.noTFBS_2xMM9
							)
			histone_features_mouse <- list('mm9_no_tfbs'=re.noTFBS_MM9,
								'hg191x_no_tfbs'=re.noTFBS_1xHG19,
								'mm92x_no_tfbs'=re.noTFBS_2xMM9
							)
			return(histone_features_mouse)
		}
		else
		{
			stopifnot(type == 'expression')
		}

	}
	return(NULL)
}


init_reporting <- function(folder, name)
{
	template.folder <- '/home/pebert/work/code/sandbox/csepigen/trunk/reports'
	template.style <- paste(template.folder, 'styling_tmpl.css', sep='/')
	template.nav <- paste(template.folder, 'nav_tmpl.html', sep='/')
	template.index <- paste(template.folder, 'index_tmpl.html', sep='/')
	template.content <- paste(template.folder, 'content_tmpl.html', sep='/')
	report.main <- paste(paste(folder, name, sep='/'))
	report.sub <- paste(paste(folder, name, 'content_data', sep='/'))
	if (!file.exists(paste(folder, name, 'content_data', sep='/'))) {
		dir.create(path=paste(folder, name, 'content_data', sep='/'), showWarnings=TRUE, recursive=TRUE)
	}
	file.copy(template.style, paste(report.sub, 'style.css', sep='/'), overwrite=TRUE)
	file.copy(template.nav, paste(report.sub, 'nav.html', sep='/'), overwrite=TRUE)
	file.copy(template.index, paste(report.main, 'index.html', sep='/'), overwrite=TRUE)
	file.copy(template.content, paste(report.sub, 'content.html', sep='/'), overwrite=TRUE)
	reporting.files <- list('report_main_folder'=report.main,
							'report_sub_folder'=report.sub,
							'report_css'=paste(report.sub, 'style.css', sep='/'),
							'report_nav'=paste(report.sub, 'nav.html', sep='/'),
							'report_index'=paste(report.main, 'index.html', sep='/'),
							'report_content'=paste(report.sub, 'content.html', sep='/')
							)
	return(reporting.files)
}

descriptors <- c('chrom', 'X.chrom', 'start', 'end', 'strand', 'length',
				'id', 'class', 'label', 'assembly', 'matching_id', 'relax',
				'gene_name', 'gene_type', 'ensembl_gene_id', 'RPKM',
				'ensembl_transcript_id', 'transcript_type', 'source',
				'transcript_name', 'havana_gene_name', 'havana_gene_id',
				'havana_transcript_id', 'havana_transcript_name', 'seq', 'element',
				'promoter_class', 'RefSeq_mRNA')

# rev: order from dark blue to light blue
colors.blue <- rev(c('141d90', '121a81', '101773', '0e1464', '0c1156', '0a0e48'))

match.predictors <- function(var.names, regexp, boolean=FALSE)
{
	matched.items <- unlist(sapply(var.names, function(x) grepl(regexp, x)))
	if (boolean) { return(as.vector(matched.items)) }
	return(names(which(matched.items)))
}

outwrite <- function(outstring)
{
    cat(paste(Sys.time(), outstring, sep=':  '), sep='\n', file=stdout(), fill=F)
    return(NULL)
}

removeNearZeroVarPred <- function(dataframe, thresholdfreq, descr)
{
	outwrite(paste('Dimension of dataframe:', dim(dataframe), sep=' '))
	reduced.data <- dataframe[, !(names(dataframe) %in% descr)]
	# remove all features with zero variance
	var.col <- colVars(reduced.data)
	zeroVarPred <- names(which(var.col == 0))
	reduced.data <- reduced.data[, !(names(reduced.data) %in% zeroVarPred)]
	# remove all features that are close to zero variance
	# frequency of most common value >= threshold freq
	threshold.abs <- nrow(dataframe) * thresholdfreq
	max.count <- apply(reduced.data, 2, function(X) max(table(X)) )
	nearZeroPred <- names(which(max.count >= threshold.abs, arr.ind=F))
	reduced.data <- reduced.data[, !(names(reduced.data) %in% nearZeroPred)]
	outwrite(paste('Dimension of dataframe w/o n0var predictors:', dim(reduced.data), sep=' '))
	reduced.data <- cbind(reduced.data, dataframe[, names(dataframe) %in% descr])
	return(reduced.data)
}

createScaledDataList <- function(dataframe, dataindices, descr)
{
	traindata <- scale(dataframe[dataindices$training, features], center=TRUE, scale=TRUE)
	scalecenter <- attr(traindata, 'scaled:center')
	scalescale <- attr(traindata, 'scaled:scale')
	traindata <- data.frame(traindata)
	trainlabels <- as.factor(dataframe[dataindices$training, ]$label)
	testdata <- scale(dataframe[dataindices$testing, features], center=scalecenter, scale=scalescale)
	testdata <- data.frame(testdata)
	testlabels <- as.factor(dataframe[dataindices$testing, ]$label)
	dataset <- list(train.data=traindata, train.labels=trainlabels, test.data=testdata, test.labels=testlabels)
	stopifnot(!is.na(traindata), !is.na(testdata))
	return(dataset)
}

trainTuneRF <- function(tuninglist, tcontrol, use.feat)
{
	predictors <- floor(sqrt(ncol(tuninglist$train.data[,use.feat])))
	low.lim <- ceiling(predictors*0.5)
	high.lim <- floor(predictors*1.5)
	mtry <- seq(from=low.lim, to=high.lim, by=2)
	rfreport <- create_content_paragraph("RandomForest setup")
	rfreport <- paste(rfreport, create_list(c("Run: Classification",
						"Selection metric: Accuracy",
						"Number of trees: 500",
						paste("Default value for model paramter mtry: ", predictors, sep=""),
						paste("Values of mtry: ", paste(mtry, collapse=" | "), sep=""))),
				sep="\n")
	rf.grid <- expand.grid(.mtry=mtry)
	tuned.rf <- train(x=tuninglist$train.data[,use.feat], y=as.factor(tuninglist$train.labels),
					method='rf', metric='Accuracy', maximize=TRUE,
					trControl=tcontrol, tuneGrid=rf.grid,
					ntree=500, importance=TRUE)
	bestparams <- tuned.rf$bestTune
	paratable <- create_table("Parameters of selected model", names(bestparams), matrix(bestparams, byrow=TRUE, ncol=length(names(bestparams))) )
	rfreport <- paste(rfreport, paratable, sep="\n\n")
	tmp <- list(model=tuned.rf, reporting=rfreport)
	return(tmp)
}

trainTuneLinSVM <- function(tuninglist, tcontrol, use.feat)
{
	cvals <- 2^(seq(from=-6, to=6, by=0.5))
	lin.grid <- expand.grid(.C=cvals)
	svmreport <- create_content_paragraph("Linear SVM setup")
	svmreport <- paste(svmreport, create_list(c("Run: Classification",
						"Selection metric: Accuracy",
						"Default value of model parameter C: 1",
						paste("Tested values of C (exact): ", paste("2^", seq(-6,6,0.5), sep="", collapse=" | "), sep=""),
						paste("Tested values of C (rounded): ", paste(signif(cvals, 4), collapse=" | "), sep="")
						)),
				sep="\n")
	tuned.linsvm <- train(x=tuninglist$train.data[,use.feat], y=as.factor(tuninglist$train.labels),
						method='svmLinear', metric='Accuracy', maximize=TRUE,
						trControl=tcontrol, tuneGrid=lin.grid)
	bestparams <- tuned.linsvm$bestTune
	paratable <- create_table("Parameters of selected model", names(bestparams), matrix(bestparams, byrow=TRUE, ncol=length(names(bestparams))) )
	svmreport <- paste(svmreport, paratable, sep="\n\n")
	tmp <- list(model=tuned.linsvm, reporting=svmreport)
	return(tmp)
}

trainTuneRBFSVM <- function(tuninglist, tcontrol, use.feat)
{
	cvals <- 2^(seq(from=-2, to=6, by=0.5))
	predictors <- ncol(tuninglist$train.data)
	# default is 1/dim
	lim.low <- 1/(predictors * 1.4)
	lim.high <- 1/(predictors * 0.7)	
	sigmas <- seq(from=lim.low, to=lim.high, length.out=length(cvals))
	rbf.grid <- expand.grid(.C=cvals, .sigma=sigmas)
	tuned.rbfsvm <- train(x=tuninglist$train.data[,use.feat], y=as.factor(tuninglist$train.labels),
						method='svmRadial', metric='Accuracy', maximize=TRUE,
						trControl=tcontrol, tuneGrid=rbf.grid)
	return(tuned.rbfsvm)
}

reduceFeatureSet <- function(all.names, re.match)
{
	if (re.match == 'keep_all')
	{
		return(all.names)
	}
	rem.feat <- grepl(re.match, all.names, fixed=FALSE)
	return(all.names[!rem.feat])
}

################## BEGIN Histone

createHistTuningList <- function(dataset, low.limit)
{
	high.limit <- 1. - low.limit
	limits <- as.numeric(quantile(dataset$perc_CpG, c(low.limit, high.limit)))
	lowCpG <- limits[1]
	highCpG <- limits[2]
	num.train <- 500
	num.test <- 100
	all.low <- row.names(dataset[dataset$perc_CpG <= lowCpG, ])
	all.high <- row.names(dataset[dataset$perc_CpG >= highCpG, ])
	indexlist <- list(mixed=sampleHistonePeaks(dataset, '1', '-1', num.train, num.test),
					low_CpG=sampleHistonePeaks(dataset[all.low, ], '1', '-1', num.train, num.test),
					high_CpG=sampleHistonePeaks(dataset[all.high, ], '1', '-1', num.train, num.test))
	mixed <- buildDataset.Histone(dataset, indexlist$mixed)
	low_CpG <- buildDataset.Histone(dataset, indexlist$low_CpG)
	high_CpG <- buildDataset.Histone(dataset, indexlist$high_CpG)
	all.tuninglists <- list(mixed=mixed, low_CpG=low_CpG, high_CpG=high_CpG)
	result <- list(alltuninglists=all.tuninglists, indexlist=indexlist)
	return(result)
}

sampleHistonePeaks <- function(dataset, pos.label, neg.label, numtrain, numtest)
{
	outwrite('=== Sampling histone peaks')
	init_train <- numtrain
	init_test <- numtest
	idxlist <- list()
	outwrite('Class distribution')
	write.table(table(dataset$class))
	num_class <- table(dataset$label)
	outwrite('CpG composition of classes')
	write.table(aggregate(dataset$perc_CpG, by=list(dataset$class), FUN=summary), file=stdout(), quote=TRUE)
	if (num_class[pos.label] > numtrain + numtest & num_class[neg.label] > numtrain + numtest)
	{
		outwrite('Easy case, enough data')
		train.pos.idx <- sample(row.names(dataset[dataset$label == pos.label, ]), numtrain, replace=FALSE)
		train.neg.idx <- row.names(dataset[dataset$matching_id %in% (dataset[train.pos.idx,]$id),])
		if (length(train.neg.idx) < numtrain)
		{
			outwrite(paste(numtrain - length(train.neg.idx), 'matched regions (training) not included in dataset', sep=' '))
			tmp <- sample(row.names(dataset[dataset$label == neg.label & !(row.names(dataset) %in% train.neg.idx), ]), numtrain - length(train.neg.idx), replace=FALSE)
			train.neg.idx <- c(train.neg.idx, tmp)
		}
		idxlist$train.pos <- train.pos.idx
		idxlist$train.neg <- train.neg.idx
		dataset <- dataset[!(row.names(dataset) %in% c(idxlist$train.pos, idxlist$train.neg)), ]
		test.pos.idx <- sample(row.names(dataset[dataset$label == pos.label, ]), numtest, replace=FALSE)
		test.neg.idx <- row.names(dataset[dataset$matching_id %in% (dataset[test.pos.idx,]$id), ])
		if (length(test.neg.idx) < numtest)
		{
			outwrite(paste(numtest - length(test.neg.idx), 'matched regions (testing) not included in dataset', sep=' '))
			tmp <- sample(row.names(dataset[dataset$label == neg.label & !(row.names(dataset) %in% test.neg.idx), ]), numtest - length(test.neg.idx), replace=FALSE)
			test.neg.idx <- c(test.neg.idx, tmp)
		}
		idxlist$test.pos <- test.pos.idx
		idxlist$test.neg <- test.neg.idx
	}
	else
	{
		outwrite('Need to bootstrap')
		relative.amount <- numtrain %/% numtest
		if (num_class[pos.label] > numtrain + numtest)
		{
			train.pos.idx <- sample(row.names(dataset[dataset$label == pos.label, ]), numtrain, replace=FALSE)
			test.pos.idx <- sample(row.names(dataset[dataset$label == pos.label & !(row.names(dataset) %in% train.pos.idx), ]), numtest, replace=FALSE)
			idxlist$train.pos <- train.pos.idx
			idxlist$test.pos <- test.pos.idx
		}
		else
		{
			outwrite('Need to bootstrap for pos label')
			num.pos.test <- num_class[pos.label] %/% (relative.amount + 1)
			test.pos.idx <- sample(row.names(dataset[dataset$label == pos.label, ]), numtest, replace=ifelse(num.pos.test >= numtest, FALSE, TRUE))
			train.pos.idx <- sample(row.names(dataset[dataset$label == pos.label & !(row.names(dataset) %in% test.pos.idx), ]), numtrain, replace=ifelse(num_class[pos.label] - num.pos.test >= numtrain, FALSE, TRUE))
			idxlist$train.pos <- train.pos.idx
			idxlist$test.pos <- test.pos.idx
		}
		if (num_class[neg.label] > numtrain + numtest)
		{
			train.neg.idx <- sample(row.names(dataset[dataset$label == neg.label, ]), numtrain, replace=FALSE)
			test.neg.idx <- sample(row.names(dataset[dataset$label == neg.label & !(row.names(dataset) %in% train.neg.idx), ]), numtest, replace=FALSE)
			idxlist$train.neg <- train.neg.idx
			idxlist$test.neg <- test.neg.idx
		}
		else
		{
			outwrite('Need to bootstrap for neg label')
			num.neg.test <- num_class[neg.label] %/% (relative.amount + 1)
			test.neg.idx <- sample(row.names(dataset[dataset$label == neg.label, ]), numtest, replace=ifelse(num.neg.test >= numtest, FALSE, TRUE))
			train.neg.idx <- sample(row.names(dataset[dataset$label == neg.label & !(row.names(dataset) %in% test.neg.idx), ]), numtrain, replace=ifelse(num_class[neg.label] - num.neg.test >= numtrain, FALSE, TRUE))
			idxlist$train.neg <- train.neg.idx
			idxlist$test.neg <- test.neg.idx
		}
	}	
	stopifnot(length(c(idxlist$train.pos, idxlist$train.neg)) == 2*init_train, length(c(idxlist$test.pos, idxlist$test.neg)) == 2*init_test)
	return(idxlist)
}

buildDataset.Histone <- function(dataset, indices)
{
	trainset.pos <- dataset[indices$train.pos, ]
	trainset.neg <- dataset[indices$train.neg, ]
	trainset.pos$label <- 'peak'
	trainset.neg$label <- 'nopeak'
	trainset <- data.frame(rbind(trainset.pos, trainset.neg))
	train.data <- trainset[, !(names(trainset) %in% descriptors)]
	train.labels <- trainset$label
	train.data <- scale(train.data, center=TRUE, scale=TRUE)
	scalecenter <- attr(train.data, 'scaled:center')
	scalescale <- attr(train.data, 'scaled:scale')
	train.data <- data.frame(train.data)
	testset.pos <- dataset[indices$test.pos, ]
	testset.neg <- dataset[indices$test.neg, ]
	testset.pos$label <- 'peak'
	testset.neg$label <- 'nopeak'
	testset <- data.frame(rbind(testset.pos, testset.neg))
	test.data <- testset[, !(names(testset) %in% descriptors)]
	test.labels <- testset$label
	test.data <- scale(test.data, center=scalecenter, scale=scalescale)
	test.data <- data.frame(test.data)
	tuninglist <- list(train.raw=trainset, train.data=train.data, train.labels=train.labels,
					test.raw=testset, test.data=test.data, test.labels=test.labels)
	return(tuninglist)
}

tuneModels.Histone <- function(inputlists, red.feat)
{
	# tuninglists consists of 4 sublists with the four tuning scenarios
	save_list <- list()
	scenarios <- names(inputlists)
	for (i in 1:length(scenarios))
	{
		outwrite(paste('Tuning scenario:', scenarios[i], sep=' '))
		tuninglist <- inputlists[[scenarios[i]]]
		outwrite('CpG composition of regions')
		trainset <- tuninglist$train.raw
		write.table(aggregate(trainset$perc_CpG, by=list(trainset$label), FUN=summary), file=stdout(), quote=TRUE)
		all.feat <- names(tuninglist$train.data)
		use.feat <- reduceFeatureSet(all.feat, red.feat)
		outwrite('Using features...')
		outwrite(paste(use.feat, collapse=' '))
		sef <- 'oneSE' # or 'best'
		outwrite('Start tuning')
		tctl <- trainControl(method='cv', number=10, repeats=2,
						savePredictions=TRUE, classProbs=TRUE,
						summaryFunction=defaultSummary,
						selectionFunction=sef, allowParallel=TRUE)
		rf <- trainTuneRF(tuninglist, tctl, use.feat)
		outwrite('RandomForest done')
		svm.lin <- trainTuneLinSVM(tuninglist, tctl, use.feat)
		outwrite('Linear SVM done')
		svm.rbf <- trainTuneRBFSVM(tuninglist, tctl, use.feat)
		outwrite('RBF SVM done')
		models <- list(rf, svm.lin, svm.rbf)
		names(models) <- c('rf', 'svmLin', 'svmRBF')
		save_list[[length(save_list) + 1]] <- list(models=models, features=use.feat)
		outwrite(paste(rep('-', 80), collapse=''))
	}
	names(save_list) <- scenarios
	return(save_list)
}

runTestData.Histone <- function(lsts.model.feat, tuninglists)
{
	save_list <- list()
	scenarios <- names(tuninglists)
	for (i in 1:length(scenarios))
	{
		outwrite(paste('Running test data for scenario', scenarios[i], sep=' '))
		model.feat <- lsts.model.feat[[scenarios[i]]]
		tuninglist <- tuninglists[[scenarios[i]]]
		testset <- tuninglist$test.raw
		outwrite('CpG composition of regions')
		write.table(aggregate(testset$perc_CpG, by=list(testset$label), FUN=summary), file=stdout(), quote=TRUE)
		models <- model.feat$models
		feat <- model.feat$features
		preds <- extractPrediction(models, testX=tuninglist$test.data[,feat], testY=as.factor(tuninglist$test.labels))
		probs <- extractProb(models, testX=tuninglist$test.data[,feat], testY=as.factor(tuninglist$test.labels))
		rbf.pred <- subset(preds, model=='svmRadial' & dataType=='Test', select=c(obs,pred))
		lin.pred <- subset(preds, model=='svmLinear' & dataType=='Test', select=c(obs,pred))
		rf.pred <- subset(preds, model=='rf' & dataType=='Test', select=c(obs,pred))
		cmlist <- list(rf.test=confusionMatrix(rf.pred$pred, rf.pred$obs),
					 lin.test=confusionMatrix(lin.pred$pred, lin.pred$obs),
					 rbf.test=confusionMatrix(rbf.pred$pred, rbf.pred$obs))
		cm.names <- names(cmlist$rf.test$overall)
		outwrite('RandomForest confusion matrix for test data')
		outwrite(paste(cm.names, as.vector(cmlist$rf.test$overall), sep=': ', collapse='\n'))
		outwrite('Linear SVM confusion matrix for test data')
		outwrite(paste(cm.names, as.vector(cmlist$lin.test$overall), sep=': ', collapse='\n'))
		outwrite('RBF SVM confusion matrix for test data')
		outwrite(paste(cm.names, as.vector(cmlist$rbf.test$overall), sep=': ', collapse='\n'))
		save_list[[length(save_list) + 1]] <- list(testPred=preds, testProbs=probs, testConfMatrix=cmlist)
	}
	names(save_list) <- scenarios
	return(save_list)
}

################## END Histone

################## BEGIN Expression

runTestData.Expression <- function(lsts.model.feat, tuninglists)
{
	save_list <- list()
	scenarios <- names(tuninglists)
	reporting <- ""
	for (i in 1:length(scenarios))
	{
		report_section <- create_content_paragraph("Description of test data")
		outwrite(paste('Running test data for scenario', scenarios[i], sep=' '))
		model.feat <- lsts.model.feat[[scenarios[i]]]
		tuninglist <- tuninglists[[scenarios[i]]]
		outwrite('Composition of test data')
		testset <- tuninglist$test.raw
		testtable <- table(testset$promoter_class, testset$label)
		write.table(testtable, file=stdout(), quote=FALSE)
		#reporting
		names.row <- t(t(attr(testtable, "dimnames")[[1]]))
		names.col <- c("Promotertype", attr(testtable, "dimnames")[[2]])
		testtable <- data.frame(cbind(names.row, matrix(testtable, ncol=length(names.col)-1)) )
		#reporting
		report_section <- paste(report_section, create_table("Composition of test data", names.col, testtable), sep="\n")
		outwrite('CpG composition of promoters')
		cpgtable <- aggregate(testset$perc_CpG, by=list(testset$promoter_class), FUN=summary)
		write.table(cpgtable, file=stdout(), quote=FALSE)
		#reporting
		names.col <- c("Promotertype", attr(cpgtable$x, "dimnames")[[2]])
		names.row <- t(t(as.vector(cpgtable$Group.1)))
		cpgtable <- data.frame(cbind(names.row, as.matrix(cpgtable$x)) )
		#reporting
		cpgtable <- create_table("CpG composition of promoters", names.col, cpgtable)
		report_section <- paste(report_section, cpgtable, sep="\n")
		models <- model.feat$models
		feat <- model.feat$features
		preds <- extractPrediction(models, testX=tuninglist$test.data[,feat], testY=as.factor(tuninglist$test.labels))
		probs <- extractProb(models, testX=tuninglist$test.data[,feat], testY=as.factor(tuninglist$test.labels))
		#rbf.pred <- subset(preds, model=='svmRadial' & dataType=='Test', select=c(obs,pred))
		lin.pred <- subset(preds, model=='svmLinear' & dataType=='Test', select=c(obs,pred))
		rf.pred <- subset(preds, model=='rf' & dataType=='Test', select=c(obs,pred))
		#cmlist <- list(rf.test=confusionMatrix(rf.pred$pred, rf.pred$obs),
		#			 lin.test=confusionMatrix(lin.pred$pred, lin.pred$obs),
		#			 rbf.test=confusionMatrix(rbf.pred$pred, rbf.pred$obs))
		cmlist <- list(rf.test=confusionMatrix(rf.pred$pred, rf.pred$obs),
					lin.test=confusionMatrix(lin.pred$pred, lin.pred$obs))
		#cm.names <- names(cmlist$rf.test$overall[[c('Accuracy', 'AccuracyLower', 'AccuracyUpper', 'AccuracyPValue', 'AccuracyNull')]])
		report_section <- paste(report_section, create_content_paragraph("Model performance on test dataset"), sep="\n")
		perf.names <- c('Accuracy', 'AccuracyLower', 'AccuracyUpper', 'AccuracyPValue', 'AccuracyNull')
		cm.names <- c("Accuracy", "95% CI lower", "95% CI upper", "Acc. p-val", "NIR")
		rfacc <- matrix(signif(cmlist$rf.test$overall[perf.names], 4), byrow=TRUE, ncol=length(cm.names))
		rftable <- create_table("RandomForest performance", cm.names, rfacc)
		outwrite('RandomForest confusion matrix for test data')
		outwrite(paste(cm.names, as.vector(cmlist$rf.test$overall[perf.names]), sep=': ', collapse='\n'))
		outwrite('Linear SVM confusion matrix for test data')
		outwrite(paste(cm.names, as.vector(cmlist$lin.test$overall[perf.names]), sep=': ', collapse='\n'))
		svmacc <- matrix(signif(cmlist$lin.test$overall[perf.names], 4), byrow=TRUE, ncol=length(cm.names))
		svmtable <- create_table("Linear SVM performance", cm.names, svmacc)
		report_section <- paste(report_section, rftable, svmtable, sep="\n\n")
		report_section <- wrap_section(report_section, paste("Testing scenario: ", scenarios[i], sep=""), "h2")
		reporting <- paste(reporting, report_section, sep="\n\n")
		#outwrite('RBF SVM confusion matrix for test data')
		#outwrite(paste(cm.names, as.vector(cmlist$rbf.test$overall), sep=': ', collapse='\n'))
		save_list[[length(save_list) + 1]] <- list(testPred=preds, testProbs=probs, testConfMatrix=cmlist)
	}
	reporting <- wrap_section(reporting, "Model Performance on Test Data", "h1")
	reporting <- wrap_visibility_toggle(reporting, "div_model_testing", "Model performance report")
	names(save_list) <- scenarios
	tmp <- list(output=save_list, reporting=reporting)
	return(tmp)
}

tuneModels.Expression <- function(inputlists, red.feat)
{
	# tuninglists consists of 4 sublists with the four tuning scenarios
	save_list <- list()
	scenarios <- names(inputlists)
	reporting <- ""
	for (i in 1:length(scenarios))
	{
		report_section <- create_content_paragraph("Description of training data")
		outwrite(paste('Tuning scenario:', scenarios[i], sep=' '))
		tuninglist <- inputlists[[scenarios[i]]]
		outwrite('Composition of training data')
		trainset <- tuninglist$train.raw
		traintable <- table(trainset$promoter_class, trainset$label)
		write.table(traintable, file=stdout(), quote=FALSE)
		#reporting
		names.row <- t(t(attr(traintable, "dimnames")[[1]]))
		names.col <- c("Promotertype", attr(traintable, "dimnames")[[2]])
		traintable <- data.frame(cbind(names.row, matrix(traintable, ncol=length(names.col)-1)) )
		#reporting		
		traintable <- create_table("Composition of training data", names.col, traintable)
		report_section <- paste(report_section, traintable, sep='\n')
		outwrite('CpG composition of promoters')
		cpgtable <- aggregate(trainset$perc_CpG, by=list(trainset$promoter_class), FUN=summary)
		write.table(cpgtable, file=stdout(), quote=FALSE)
		#reporting
		names.col <- c("Promotertype", attr(cpgtable$x, "dimnames")[[2]])
		names.row <- t(t(as.vector(cpgtable$Group.1)))
		cpgtable <- data.frame(cbind(names.row, as.matrix(cpgtable$x)) )
		#reporting
		cpgtable <- create_table("CpG composition of promoters", names.col, cpgtable)
		report_section <- paste(report_section, cpgtable, sep='\n')
		all.feat <- names(tuninglist$train.data)
		use.feat <- reduceFeatureSet(all.feat, red.feat)
		if (i==1)
		{
			report_tuning <- create_content_paragraph("General parameters of the tuning process")
			tuning_params <- c("Method: CV", "Folds: 10", "CV repeats: 2", "Summary: default", "Selection rule: oneSE")
			report_tuning <- paste(report_tuning, create_list(tuning_params))
			report_tuning <- wrap_section(report_tuning, "Setup of tuning process", "h2")
			report_feat <- create_content_paragraph("The following features were used in all models")
			remainder <- 20 - (length(use.feat) %% 20)
			feat.table <- create_table("Used features", paste("FeatCol-", c(1:20), sep=""), matrix(c(use.feat, rep("-", remainder)), ncol=20, byrow=TRUE))
			report_feat <- paste(report_feat, feat.table, sep='\n')
			report_feat <- wrap_section(report_feat, "Model Features", "h2")
			report_feat <- wrap_visibility_toggle(report_feat, "div_rep_feat", "feature report")
			reporting <- paste(report_tuning, report_feat, reporting, sep='\n')
		}
		outwrite('Using features...')
		outwrite(paste(use.feat, collapse=' '))
		sef <- 'oneSE' # or 'best'
		outwrite('Start tuning')
		tctl <- trainControl(method='cv', number=10, repeats=2,
						savePredictions=TRUE, classProbs=TRUE,
						summaryFunction=defaultSummary,
						selectionFunction=sef, allowParallel=TRUE)
		tmp <- trainTuneRF(tuninglist, tctl, use.feat)
		rf <- tmp[['model']]
		report_section <- paste(report_section, tmp[['reporting']])
		outwrite('RandomForest done')
		tmp <- trainTuneLinSVM(tuninglist, tctl, use.feat)
		svm.lin <- tmp[['model']]
		report_section <- paste(report_section, tmp[['reporting']])
		outwrite('Linear SVM done')
		#svm.rbf <- trainTuneRBFSVM(tuninglist, tctl, use.feat)
		#outwrite('RBF SVM done')
		#models <- list(rf, svm.lin, svm.rbf)

		# reporting
		report_section <- wrap_section(report_section, paste('Tuning scenario: ', scenarios[i], sep=''), "h2")
		reporting <- paste(reporting, report_section, sep='\n\n')
		# reporting
		models <- list(rf, svm.lin)
		save_list[[length(save_list) + 1]] <- list(models=models, features=use.feat)
		outwrite(paste(rep('-', 80), collapse=''))
	}
	reporting <- wrap_section(reporting, "Model Tuning", "h1")
	reporting <- wrap_visibility_toggle(reporting, "div_model_tuning", "Model tuning report")
	names(save_list) <- scenarios
	tmp <- list(output=save_list, reporting=reporting)
	return(tmp)
}

samplePromoterClasses <- function(dataset, gold_label, data_label, numtrain, numtest)
{
	idxlist <- list()
	init_train <- numtrain
	init_test <- numtest
	outwrite(paste('Gold label', gold_label, sep=' '))
	outwrite(paste('Data label', data_label, sep=' '))
	num_promclasses <- table(dataset$promoter_class)
	write.table(num_promclasses, file=stdout())
	if (num_promclasses[gold_label] >= numtrain + numtest)
	{	# easy, enough gold labels
		outwrite('Perfect case, enough gold standard labels')
		train.idx <- sample(row.names(dataset[dataset$promoter_class == gold_label, ]), numtrain, replace=FALSE)
		test.idx <- sample(row.names(dataset[!(row.names(dataset) %in% train.idx) & dataset$promoter_class == gold_label, ]), numtest, replace=FALSE)
		idxlist$training <- train.idx
		idxlist$testing <- test.idx
	}
	else
	{
		outwrite('Sampling including data-based labels / UNDEF state promoters')
		relative.amount <- numtrain %/% numtest
		num.test.gold <- num_promclasses[gold_label] %/% (relative.amount + 1)
		num.test.data <- max(0, min(num_promclasses[data_label] %/% (relative.amount + 1), numtest - num.test.gold))
		numtest <- numtest - num.test.gold - num.test.data
		test.idx.gold <- sample(row.names(dataset[dataset$promoter_class == gold_label, ]), num.test.gold, replace=FALSE)
		test.idx.data <- sample(row.names(dataset[dataset$promoter_class == data_label, ]), num.test.data, replace=FALSE)
		test.idx.undef <- c()
		num.test.undef <- num_promclasses['UNDEF'] %/% (relative.amount + 1)
		sep.test.undef <- sample(row.names(dataset[dataset$promoter_class == 'UNDEF', ]), num.test.undef, replace=FALSE)
		if (numtest > 0)
		{
			if (num.test.undef < numtest)
			{
				outwrite('Bootstrapping for test data')
				test.idx.undef <- sample(sep.test.undef, numtest, replace=TRUE)
			}
			else
			{
				test.idx.undef <- sample(sep.test.undef, numtest, replace=FALSE)
			}
		}
		idxlist$testing <- c(test.idx.gold, test.idx.data, test.idx.undef)
		dataset <- dataset[!(row.names(dataset) %in% idxlist$testing), ]
		num_promclasses <- table(dataset$promoter_class)
		train.idx.gold <- row.names(dataset[dataset$promoter_class == gold_label, ])
		num.train.data <- max(0, min(numtrain - length(train.idx.gold), num_promclasses[data_label]))
		train.idx.data <- sample(row.names(dataset[dataset$promoter_class == data_label, ]), num.train.data, replace=FALSE)
		train.idx.undef <- c()
		numtrain <- numtrain - length(train.idx.gold) - length(train.idx.data)
		if (numtrain > 0)
		{
			if (num_promclasses['UNDEF'] < numtrain)
			{
				outwrite('Bootstrapping for train data')
				train.idx.undef <- sample(row.names(dataset[dataset$promoter_class == 'UNDEF', ]), numtrain, replace=TRUE)
			}
			else
			{
				train.idx.undef <- sample(row.names(dataset[dataset$promoter_class == 'UNDEF', ]), numtrain, replace=FALSE)
			}
		}
		idxlist$training <- c(train.idx.gold, train.idx.data, train.idx.undef)
	}
	outwrite(paste('No. of training samples:', length(idxlist$training), sep=' '))
	outwrite(paste('No. of test samples:', length(idxlist$testing), sep=' '))
	stopifnot(length(idxlist$training) == init_train, length(idxlist$testing) == init_test)
	outwrite('Returning index list')
	return(idxlist)
}

buildDataset.Expression <- function(dataset, low.idx, high.idx)
{
	trainset.low <- dataset[low.idx$training, ]
	trainset.high <- dataset[high.idx$training, ]
	trainset.low$label <- 'low'
	trainset.high$label <- 'high'
	trainset <- data.frame(rbind(trainset.low, trainset.high))
	train.data <- trainset[, !(names(trainset) %in% descriptors)]
	train.labels <- trainset$label
	train.data <- scale(train.data, center=TRUE, scale=TRUE)
	scalecenter <- attr(train.data, 'scaled:center')
	scalescale <- attr(train.data, 'scaled:scale')
	train.data <- data.frame(train.data)
	testset.low <- dataset[low.idx$testing, ]
	testset.high <- dataset[high.idx$testing, ]
	testset.low$label <- 'low'
	testset.high$label <- 'high'
	testset <- data.frame(rbind(testset.low, testset.high))
	test.data <- testset[, !(names(testset) %in% descriptors)]
	test.labels <- testset$label
	test.data <- scale(test.data, center=scalecenter, scale=scalescale)
	test.data <- data.frame(test.data)
	tuninglist <- list(train.raw=trainset, train.data=train.data, train.labels=train.labels,
					test.raw=testset, test.data=test.data, test.labels=test.labels)
	return(tuninglist)
}

createExprTuningList <- function(dataset, low.limit)
{
	high.limit <- 1. - low.limit
	limits <- as.numeric(quantile(dataset$RPKM, c(low.limit, high.limit)))
	lowRPKM <- limits[1]
	highRPKM <- limits[2]
	num.train <- 500
	num.test <- 100
	all.low <- row.names(dataset[dataset$RPKM <= lowRPKM, ])
	all.high <- row.names(dataset[dataset$RPKM >= highRPKM, ])
	indexlist <- list(low_LCG=samplePromoterClasses(dataset[all.low, ], 'gLCG', 'dLCG', num.train, num.test),
					low_HCG=samplePromoterClasses(dataset[all.low, ], 'gHCG', 'dHCG', num.train, num.test),
					high_LCG=samplePromoterClasses(dataset[all.high, ], 'gLCG', 'dLCG', num.train, num.test),
					high_HCG=samplePromoterClasses(dataset[all.high, ], 'gHCG', 'dHCG', num.train, num.test))
	lowLCG_highLCG <- buildDataset.Expression(dataset, indexlist$low_LCG, indexlist$high_LCG)
	lowLCG_highHCG <- buildDataset.Expression(dataset, indexlist$low_LCG, indexlist$high_HCG)
	lowHCG_highHCG <- buildDataset.Expression(dataset, indexlist$low_HCG, indexlist$high_HCG)
	lowHCG_highLCG <- buildDataset.Expression(dataset, indexlist$low_HCG, indexlist$high_LCG)
	all.tuninglists <- list(lowLCG_highLCG=lowLCG_highLCG, lowLCG_highHCG=lowLCG_highHCG,
					lowHCG_highHCG=lowHCG_highHCG, lowHCG_highLCG=lowHCG_highLCG)
	result <- list(alltuninglists=all.tuninglists, indexlist=indexlist)
	return(result)
}

############### END Expression

createTuningList <- function(datafile_loc, datatype)
{
	dataset <- read.table(datafile_loc, header=TRUE, sep='\t', comment.char='')
	outwrite(paste('Read file:', datafile_loc, sep=' '))
	dataset <- removeNearZeroVarPred(dataset, 0.95, descriptors)
	if (datatype == 'idx_expr')
	{
		qlims <- c(.05, .15, .25, .45)
		strq <- c('perc05', 'perc15', 'perc25', 'perc45')
		tunlists.allperc <- list()
		for (i in 1:length(qlims))
		{
			outwrite(paste('Working on percentile:', strq[i], sep=' '))
			comb_lists <- createExprTuningList(dataset, qlims[i])
			outwrite('Tuning lists created')
			tunlists.allperc[[i]] <- comb_lists
			outwrite(paste(rep('-', 20), collapse=''))
		}
		names(tunlists.allperc) <- strq
		rdat.file <- gsub('.bed', '.idxexpr.Rdat', basename(datafile_loc))
		rdat.loc <- paste(dirname(datafile_loc), rdat.file, sep='/')
		save(tunlists.allperc, file=rdat.loc)
	}
	else
	{
		rdat.file <- gsub('.bed', '.idxhist.Rdat', basename(datafile_loc))
		rdat.loc <- paste(dirname(datafile_loc), rdat.file, sep='/')
		comb_lists <- createHistTuningList(dataset, .3)
		outwrite('Tuning lists created')
		save(comb_lists, file=rdat.loc)
	}
	return(NULL)
}

run_histone_pipeline <- function(feature_sets)
{
	outwrite('Histone run started...')
	cargs <- commandArgs(TRUE)
	stopifnot(cargs[1] == 'tune_hist' | cargs[1] == 'idx_hist')
	if (cargs[1] == 'idx_hist')
	{
		datafile <- cargs[2]
		stopifnot(!is.na(datafile))
		createTuningList(datafile, cargs[1])
	}
	else
	{
		tune.rdat <- tail(cargs, -1)
    	stopifnot(!is.na(tune.rdat))
    	data.files <- unlist(strsplit(tune.rdat, ' ', fixed=TRUE))
    	for (i in 1:length(data.files))
    	{
    		outwrite(paste('Processing file:', data.files[i], sep=' '))
    		load(data.files[i])
    		for (k in 1:length(feature_sets))
    		{
    			fset <- names(feature_sets)[k]
    			outwrite(paste('Processing feature set:', fset, sep=' '))
    			save.list <- list()
    			lsts.model.feat <- tuneModels.Histone(comb_lists[["alltuninglists"]], feature_sets[[k]])
    			perf.lists <- runTestData.Histone(lsts.model.feat, comb_lists[["alltuninglists"]])
    			save.list <- list(modelsFeatures=lsts.model.feat, performances=perf.lists)
    			outwrite(paste(rep('#', 80), collapse=''))
    			save.file <- paste((unlist(strsplit(data.files[i], '.', fixed=TRUE)))[1], '.', fset, '.models.Rdat', sep='')
    			save(save.list, file=save.file)
    			outwrite('File saved')
    		}
    	} 		
	}
	outwrite('Run finished')
	return(0)
}

run_expression_pipeline <- function(feature_sets)
{
	outwrite('Expression run started...')
    cargs <- commandArgs(TRUE)
    stopifnot(cargs[1] == 'tune_expr' | cargs[1] == 'idx_expr')
    if (cargs[1] == 'idx_expr')
    {
    	datafile <- cargs[2]
    	stopifnot(!is.na(datafile))
    	createTuningList(datafile, cargs[1])
    }
    else
    {
    	tune.rdat <- tail(cargs, -1)
    	stopifnot(!is.na(tune.rdat))
    	data.files <- unlist(strsplit(tune.rdat, ' ', fixed=TRUE))
    	# data.files is vector file1 file2 fileN report1 report2 reportN
    	stopifnot(length(data.files) %% 2 == 0)
    	num.files <- length(data.files) %/% 2
    	for (i in 1:num.files)
    	{
    		# reporting
    		report.files <- init_reporting('/TL/epigenetics2/work/pebert/projects/validation/results/2013_Dec/reports', data.files[i+num.files])
    		descr_page <- add_new_content_page(report.files[['report_content']], report.files[['report_sub_folder']], 'description.html')
    		report_nav <- raw_read_template(report.files[['report_nav']])
    		report_css <- raw_read_template(report.files[['report_css']])
    		report_nav <- replace_toplevel_menu_nosub('description', 'description.html', 'Run Description', report_nav)
    		report_css <- replace_css_menu_top_entry('description', colors.blue[1], report_css)
    		descr_content <- raw_read_template(descr_page)
    		descr_formatted <- create_content_paragraph('This is a description of the current settings')
    		descr_formatted <- paste(descr_formatted, create_list(c(paste('Processing file: ', data.files[i], sep=""),
    														'Training samples: 1000',
    														'Testing samples: 200',
    														'Bootstrapping: if necessary')),
    							sep='\n')
    		descr_formatted <- paste(descr_formatted, create_content_paragraph("Testing scenario: binary prediction 
    									of transcript expression (high/low) where top 5% (15%, 25%, 45%) vs bottom 5% (15%, 25%, 45%) 
    									of data have been used to group low RPKM and high RPKM genes (promoters; quantified around TSS -3kb/+2kb). 
    									Since the CpG content of promoters can be very indicative of transcription levels (high CpG content promoters 
    									are often associated with ubiquitously expressed housekeeping genes), another level of separating 
    									groups of genes is necessary. Based on a published study (gold standard), promoters were classified 
    									as HCG and LCG (high/low CpG content). For promoters were no label was availabel in the study, a 
    									label was inferred from CpG content (promoters in the 30% group with highest/lowest CpG content; data based label) 
    									or marked as UNDEF otherwise."))
    		descr_formatted <- wrap_section(descr_formatted, 'Settings for this run', 'h1')
    		start <- Sys.time()
    		# reporting
    		outwrite(paste('Processing file:', data.files[i], sep=' '))
    		load(data.files[i])
    		for (k in 1:length(feature_sets))
    		{
    			fset <- names(feature_sets)[k]
    			# reporting
    			report_css <- replace_css_menu_top_entry(fset, sample(colors.blue, 1), report_css)
    			tmp <- replace_toplevel_menu_sub(fset, fset, fset, report_nav)
    			report_nav <- tmp[['navstr']]
    			sub_id <- tmp[['subid']]
    			# reporting
    			outwrite(paste('Processing feature set:', fset, sep=' '))
    			save.list <- list()
    			list.names <- c()
    			strq <- c('perc05', 'perc15', 'perc25', 'perc45')
    			for (j in 1:length(strq))
    			{
    				# reporting
    				extend <- TRUE
    				if (j == length(strq)) { extend <- FALSE }
    				new_id <- paste(fset, '_', strq[j], sep='')
    				new_page <- paste(new_id, '.html', sep='')
    				report_page <- add_new_content_page(report.files[['report_content']], report.files[['report_sub_folder']], new_page)
    				page_content <- raw_read_template(report_page)
    				report_nav <- replace_submenu(new_id, new_page, strq[j], sub_id, extend, report_nav)
    				# reporting
    				outwrite(paste('Processing entry:', strq[j], sep=' '))
    				comb_lists <- tunlists.allperc[[j]]
    				tmp_lst <- tuneModels.Expression(comb_lists[["alltuninglists"]], feature_sets[[k]])
    				lsts.model.feat <- tmp_lst[['output']]
    				tuning_reports <- tmp_lst[['reporting']]
    				tmp_lst <- runTestData.Expression(lsts.model.feat, comb_lists[["alltuninglists"]])
    				perf.lists <- tmp_lst[['output']]
    				testing_reports <- tmp_lst[['reporting']]
    				entry.name <- paste(basename(data.files[i]), strq[j], fset, sep='_')
    				list.names <- c(list.names, entry.name)
    				save.list[[j]] <- list(modelsFeatures=lsts.model.feat, performances=perf.lists)
    				outwrite(paste(rep('#', 80), collapse=''))
    				# reporting
    				page_formatted <- paste(tuning_reports, testing_reports, sep="\n\n")
    				page_content <- add_formatted_content(page_formatted, page_content)
    				raw_write_page(report_page, page_content)
    				# reporting
    			}
    			names(save.list) <- list.names
    			save.file <- paste((unlist(strsplit(data.files[i], '.', fixed=TRUE)))[1], '.', fset, '.models.Rdat', sep='')
    			save(save.list, file=save.file)
    			outwrite('File saved')
    		}
    		# reporting
    		end <- Sys.time()
    		run.time <- round(as.numeric(difftime(end, start, units="mins")), 2)
    		run.stats <- create_table('Execution info', c('Server', '# processes', 'Runtime (mins)'), matrix(c(server_name, ncpus, run.time), nrow=1, ncol=3))
    		run.stats <- wrap_section(run.stats, 'Run statistics', 'h3')
    		descr_formatted <- paste(descr_formatted, run.stats, sep='\n\n')
    		descr_content <- add_formatted_content(descr_formatted, descr_content)
    		raw_write_page(descr_page, descr_content)
    		raw_write_page(report.files[['report_nav']], report_nav)
    		raw_write_page(report.files[['report_css']], report_css)
    		# reporting
    	} 		
    }
    outwrite('Run finished')
    return(0)
}

run <- function()
{
	cargs <- commandArgs(TRUE)
	stopifnot(!is.na(cargs[1]))
	if (grepl('hist', cargs[1], fixed=TRUE))
	{
		retcode <- run_histone_pipeline(get_regexps('histone', 'mouse'))
		return(retcode)
	}
	else
	{
		if (grepl('expr', cargs[1], fixed=TRUE))
		{
			retcode <- run_expression_pipeline(get_regexps('expression', 'human'))
			return(retcode)
		}
		else
		{
			outwrite('Do not know what to do')
			return(2)
		}
	}
	return(0)
}

retcode <- run()
quit(save='no', status=retcode)
