load('/TL/epigenetics2/work/pebert/projects/validation/featuredatasets/hg19_H1hesc_chr1-2_200_Rep1_TSS_ext5kb_Genc3c_normal_72.q05.idxexpr.Rdat')
source('csepigen_ML.R')
test <- tuneML(tuninglist)
source('csepigen_ML.R')
source('csepigen_ML.R')
source('csepigen_ML.R')
test <- tuneML(tuninglist)
source('csepigen_ML.R')
test <- tuneML(tuninglist)
source('csepigen_ML.R')
test <- tuneML(tuninglist)
source('csepigen_ML.R')
test.rf <- tuneML(tuninglist)
str(test.rf)
test.rf$finalModel$err.rate
test.rf$finalModel$Importance
test.rf$finalModel$importance
varImpPlot(test.rf$finalModel$importance)
plot(test.rf, top=20)
plotClassProb(test.rf)
plotClassProbs(test.rf)
extractProb(test.rf)
pred.rf <- predict(test.rf, testX=tuninglist$test.data, testY=as.factor(tuninglist$test.labels))
str(pred.rf)
pred.rf <- extractPrediction(test.rf, testX=tuninglist$test.data, testY=as.factor(tuninglist$test.labels))
models <- list()
models[[1]] <- test.rf
pred.rf <- extractPrediction(models, testX=tuninglist$test.data, testY=as.factor(tuninglist$test.labels))
str(pred.rf)
plotObsVsPred(pred.rf)
plotClassProbs(pred.rf)
pred.rf <- extractPrediction(models, testX=tuninglist$test.data, testY=as.factor(tuninglist$test.labels))
prob.rf <- extractProb(models, testX=tuninglist$test.data, testY=as.factor(tuninglist$test.labels))
plotClassProbs(prob.rf)
print(models[[1]])
confusionMatrix(test.rf)
plot.varImp.train(test.rf)
plot(test.rf)
plot(test.rf$importance)
str(test.rf)
plot(test.rf$finalModel)
plot(test.rf$finalModel$Importance)
plot(test.rf$finalModel$importance)
library(randomForest)
varImpPlot(test.rf$finalModel$importance)
varImpPlot(test.rf$finalModel)
confusionMatrix(pred.rf)
confusionMatrix(pred.rf$pred, pred.rf$obs)
str(pred.rf$pred)
str(pred.rf)
pred.rf$dataType == 'Test'
confusionMatrix(pred.rf$pred[pred.rf$dataType == 'Test'], pred.rf$obs[pred.rf$dataType == 'Training'])
confusionMatrix(pred.rf$pred[pred.rf$dataType == 'Test'], pred.rf$obs[pred.rf$dataType == 'Test'])
varImp(test.rf)
varImp(test.rf, scale=F)
plot(varImp(test.rf, scale=F))
plot(varImp(test.rf, scale=F), top=20)
plot(varImp(test.rf, scale=T), top=20)
plot(varImp(test.rf, scale=T, useModel=F), top=20)
?subset
subset(pred.rf, dataType=='Training')
subset(pred.rf, dataType=='Test')
2^2
10^0.25
2^0.25
2^0.1
2^-0.1
2^-0.5
2^-1
2^-1.5
2^-2
2^6
2^c(1,2,3)
2^(seq(-2,6,0.5))
source('csepigen_ML.R')
foo <- list(1,2,3)
foo
source('csepigen_ML.R')
models <- tuneML(tuninglist)
source('csepigen_ML.R')
models <- tuneML(tuninglist)
?extractPrediction
preds <- extractPrediction(models, testX=tuninglist$$test.data, testY=as.factor(tuninglist$test.labels))
preds <- extractPrediction(models, testX=tuninglist$test.data, testY=as.factor(tuninglist$test.labels))
str(preds_
)
str(preds)
rf.pred <- subset(preds, model=rf, dataType=Test, select=c(obs,pred))
confusionMatrix(rf.pred)
confusionMatrix(rf.pred$pred, rf.pred$obs)
str(rf.pred)
rf.pred <- subset(preds, model==rf & dataType==Test, select=c(obs,pred))
rf.pred <- subset(preds, model==rf && dataType==Test, select=c(obs,pred))
rf.pred <- subset(preds, model==rf ,dataType==Test, select=c(obs,pred))
rf.pred <- subset(preds, model=='' & dataType=='Test', select=c(obs,pred))
confusionMatrix(rf.pred$pred, rf.pred$obs)
confusionMatrix(as.factor(c(1,1,1,2,2,2)), as.factor(c(1,1,2,2,2,2))
)
confusionMatrix(as.factor(c(1,1,2,2,2,2)), as.factor(c(1,1,1,2,2,2))
)
preds
str(preds)
lin.pred <- subset(preds, model==svmLinear ,dataType==Test, select=c(obs,pred))
preds$model
lin.pred <- subset(preds, model==svmLinear ,dataType==Test, select=c(obs,pred))
lin.pred <- subset(preds, model=='svmLinear' ,dataType=='Test', select=c(obs,pred))
lin.pred <- subset(preds, model=='svmLinear' & dataType=='Test', select=c(obs,pred))
confusionMatrix(lin.pred$preds, lin.pred$obs)
confusionMatrix(lin.pred$pred, lin.pred$obs)
confusionMatrix(rf.pred$pred, rf.pred$obs)
lin.pred <- subset(preds, model=='svmRadial' & dataType=='Test', select=c(obs,pred))
lin.pred <- subset(preds, model=='svmLinear' & dataType=='Test', select=c(obs,pred))
rbf.pred <- subset(preds, model=='svmRadial' & dataType=='Test', select=c(obs,pred))
confusionMatrix(rbf.pred$pred, rbf.pred$obs)
?extractProb
probs <- extractProb(models, testX=tuninglist$test.data, testY=as.factor(tuninglist$test.labels))
plotClassProbs(probs)
?split
?strsplit
strsplit('file1,file2,file3', ',', fixed=T)
unlist(strsplit('file1,file2,file3', ',', fixed=T))
l <- list()
l$'foo' <- 1
l$'bar' <- 2
l
?save
?save.history
?history
savehistory(file='mlsession.r')
