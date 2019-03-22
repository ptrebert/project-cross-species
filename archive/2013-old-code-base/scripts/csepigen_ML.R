#!/TL/opt/bin/Rscript

#library(e1071)
library(caret)
#library(randomForest)
library(matrixStats)
library(doMC)

registerDoMC(cores=16)



descriptors <- c('chrom', 'X.chrom', 'start', 'end', 'strand', 'length',
				'id', 'class', 'label', 'assembly', 'matching_id', 'relax',
				'gene_name', 'gene_type', 'ensembl_gene_id', 'RPKM',
				'ensembl_transcript_id', 'transcript_type', 'source',
				'transcript_name', 'havana_gene_name', 'havana_gene_id',
				'havana_transcript_id', 'havana_transcript_name', 'seq', 'element')

signal.features <- c('p20', 'p40', 'p60', 'p80', 'p90', 'p95')

re.signal <- '^X1_p[24689]{1}[05]{1}$'
re.tfbs <- '^M0'

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
	features <- names(dataframe)[!(names(dataframe) %in% descr)]
	#
	# create ctrl-me1
	#signal.feat <- match.predictors(features, re.signal)
	#features <- features[!(features %in% signal.feat)]
	#
	traindata <- scale(dataframe[dataindices$training, features], center=TRUE, scale=TRUE)
	scalecenter <- attr(traindata, 'scaled:center')
	scalescale <- attr(traindata, 'scaled:scale')
	traindata <- data.frame(traindata)
	trainlabels <- as.factor(dataframe[dataindices$training, ]$label)
	testdata <- scale(dataframe[dataindices$testing, features], center=scalecenter, scale=scalescale)
	testdata <- data.frame(testdata)
	testlabels <- as.factor(dataframe[dataindices$testing, ]$label)
	dataset <- list(train.data=traindata, train.labels=trainlabels, test.data=testdata, test.labels=testlabels)
	return(dataset)
}

trainTuneRF <- function(dataset, tcontrol)
{
	predictors <- floor(sqrt(ncol(dataset$train.data)))
	low.lim <- ceiling(predictors*0.8)
	high.lim <- floor(predictors*1.3)
	tuned.rf <- tune.randomForest(tunecontrol=tcontrol, x=dataset$train.data, y=dataset$train.labels,
								validation.x=dataset$test.data, validation.y=dataset$test.labels,
								nodesize=1, mtry=seq(low.lim, high.lim, 2), ntree=seq(50,600,50), importance=TRUE)
	return(tuned.rf)
}

trainTuneLinSVM <- function(dataset, tcontrol)
{
	tuned.linsvm <- tune.svm(tunecontrol=tcontrol, x=dataset$train.data, y=dataset$train.labels,
							validation.x=dataset$test.data, validation.y=dataset$test.labels,
							type='C-classification', kernel='linear', cost=seq(1,11,2), scale=FALSE)
	return(tuned.linsvm)
}

trainTuneRBFSVM <- function(dataset, tcontrol)
{
	predictors <- ncol(dataset$train.data)
	# default is 1/dim
	lim.low <- 1/(predictors * 1.3)
	lim.high <- 1/(predictors * 0.8)
	tuned.rbfsvm <- tune.svm(tunecontrol=tcontrol, x=dataset$train.data, y=dataset$train.labels,
							validation.x=dataset$test.data, validation.y=dataset$test.labels,
							type='C-classification', kernel='radial', cost=seq(1,11,2),
							gamma=seq(from=lim.low, to=lim.high, length.out=10), scale=FALSE)
	return(tuned.rbfsvm)
}

tuneML <- function(tuninglist)
{
	outwrite('Start tuning')
	tctl <- trainControl(method='cv', number=10, repeats=1,
					savePredictions=TRUE, classProbs=TRUE)

	return(NULL)
}

createExprTuningList <- function(dataset, low.limit)
{
	high.limit <- 1. - low.limit
	limits <- as.numeric(quantile(dataset$RPKM, c(low.limit, high.limit)))
	lowRPKM <- limits[1]
	highRPKM <- limits[2]
	train.pos <- sample(row.names(dataset[dataset$RPKM >= highRPKM, ]), 500, replace=FALSE) 
	train.neg <- sample(row.names(dataset[dataset$RPKM <= lowRPKM, ]), 500, replace=FALSE)
	test.pos <- sample(row.names(dataset[!(row.names(dataset) %in% train.pos) & dataset$RPKM >= highRPKM, ]), 100, replace=FALSE)
	test.neg <- sample(row.names(dataset[!(row.names(dataset) %in% train.neg) & dataset$RPKM <= lowRPKM, ]), 100, replace=FALSE)
	dataset[c(train.pos, test.pos), ]$label <- 1
	dataset[c(train.neg, test.neg), ]$label <- -1
	traindata <- scale(dataset[c(train.pos, train.neg), !(names(dataset) %in% descriptors)], center=TRUE, scale=TRUE)
	scalecenter <- attr(traindata, 'scaled:center')
	scalescale <- attr(traindata, 'scaled:scale')
	traindata <- data.frame(traindata)
	trainlabels <- as.factor(dataset[c(train.pos, train.neg), ]$label)
	testdata <- scale(dataset[c(test.pos, test.neg), !(names(dataset) %in% descriptors)], center=scalecenter, scale=scalescale)
	testdata <- data.frame(testdata)
	testlabels <- as.factor(dataset[c(test.pos, test.neg), ]$label)
	tuninglist <- list(train.data=traindata, train.labels=trainlabels, test.data=testdata, test.labels=testlabels)
	return(tuninglist)
}

createTuningList <- function(datafile_loc, datatype)
{
	dataset <- read.table(datafile_loc, header=TRUE, sep='\t', comment.char='')
	outwrite(paste('Read file:', datafile_loc, sep=' '))
	dataset <- removeNearZeroVarPred(dataset, 0.9, descriptors)
	if (datatype == 'idxexpr')
	{
		qlims <- c(.05, .15, .25, .45)
		strq <- c('.q05.', '.q15.', '.q25.', '.q45.')
		for (i in 1:length(qlims))
		{
			rdat.file <- gsub('.bed', paste(strq[i], 'idxexpr.Rdat', sep=''), basename(datafile_loc))
			rdat.loc <- paste(dirname(datafile_loc), rdat.file, sep='/')
			tuninglist <- createExprTuningList(dataset, qlims[i])
			save(tuninglist, file=rdat.loc)
		}
	}
	else
	{
		rdat.file <- gsub('.bed', '.idxhist.Rdat', basename(datafile_loc))
		rdat.loc <- paste(dirname(datafile_loc), rdat.file, sep='/')
		tuninglist <- createHistTuningList(dataset)
		save(tuninglist, file=rdat.loc)
	}
	return(NULL)
}

run <- function()
{
	outwrite('Run started...')
    cargs <- commandArgs(TRUE)
    if (cargs[1] == 'idxhist' | cargs[1] == 'idxexpr')
    {
    	datafile <- cargs[2]
    	stopifnot(!is.na(datafile))
    	createTuningList(datafile, cargs[1])
    }
    else
    	if (cargs[1] == 'tune')
    	{
    		tune.rdat <- cargs[2]
    		stopifnot(!is.na(tune.rdat))
    		load(tune.rdat)
    		tuneML(tuninglist)
    	}
    	else { return(2) }
    outwrite('Run finished')
    return(0)
	"
    datafile <- cargs[1]
    file.dataindices <- cargs[2]
    outwrite('=========== Processing new file')
    outwrite(paste('File name:', datafile, sep=' '))
    load(file.dataindices)
    dataset <- read.table(datafile, header=TRUE, sep='\t', comment.char='')
    clean.dataset <- removeNearZeroVarPred(dataset, 0.9, descriptors)
    datalist <- createScaledDataList(clean.dataset, dataindices, descriptors)
    tcontrol <- tune.control(sampling='cross', cross=10, best.model=TRUE, performances=TRUE)
    outwrite('Tuning random forest')
    tune.obj.rf <- trainTuneRF(datalist, tcontrol)
    outwrite(paste('Best RF performance:', tune.obj.rf$best.performance, sep=' '))
    outwrite('Tuning linear SVM')
    tune.obj.linsvm <- trainTuneLinSVM(datalist, tcontrol)
    outwrite(paste('Best lin. SVM performance:', tune.obj.linsvm$best.performance, sep=' '))
    outwrite('Tuning radial basis SVM')
    tune.obj.rbfsvm <- trainTuneRBFSVM(datalist, tcontrol)
    outwrite(paste('Best RBF SVM performance:', tune.obj.rbfsvm$best.performance, sep=' '))
    outname <- gsub('.bed', '_ctrl-me1-9ac-me3.Rdat', datafile)
    save(tune.obj.rf, tune.obj.linsvm, tune.obj.rbfsvm, datafile, file=outname)
	outwrite('=========== Processing finished')
	"
}

retcode <- run()
quit(save='no', status=retcode)
