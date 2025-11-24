library(Hmisc)
library(randomForest)
library(ggplot2)
library(splines)
library(reshape2)
library(pROC)
library(dplyr)
library(tidyverse)
library(hash)
library(rfPermute)
library(caret)
library(ggprism)

GetOptimRF = function(data_train,y_colname,m_try = seq(1,30),n_tree=seq(100,3000,100))
{
  ctrl = trainControl(method = "cv",number = 10)
  grid = expand.grid(mtry=m_try)
  x_train = data_train[, -which(names(data_train) == y_colname)]
  y_train = data_train[,y_colname]
  rf_model = caret::train(x = x_train,
                          y = y_train,
                          method = 'rf',
                          trControl = ctrl,
                          tuneGrid = grid)
  
  best_mtry = rf_model$bestTune$mtry

  grid<- expand.grid(mtry=c(best_mtry))
  modellist<- list()

  for(ntree in n_tree)
  {
    fit<-train(x = x_train,
               y = y_train,
               method="rf",
               metric="Accuracy",
               tuneGrid=grid,
               trControl=ctrl,
               ntree=ntree)
    key <- toString(ntree)
    modellist[[key]]<- fit
  }
  results<-resamples(modellist)
  
  best_ntree = summary(results)$models[ which.max(summary(results)$statistics$Accuracy[,'Mean']) ]
  
  print(list(mtry=best_mtry,
             ntree=best_ntree))
  
  return(list(
    mtry=best_mtry,
    ntree=best_ntree
  ))
}

otu = read.csv("OTU_ID.csv",row.names = 1,check.names=FALSE)
group = read.csv("OTUGroup.csv",row.names = 1,check.names=FALSE)
taxid = read.csv("OTU_TaxID.csv")

otu = data.frame(t(otu),check.names = FALSE)
otu = cbind(otu, group)
otu$Group = as.factor(otu$Group)
set.seed(114514)

HashTable = hash()
.set(HashTable,keys=make.keys(taxid$X.ID),values=make.keys(taxid$X.NAME))

trainIndex = createDataPartition(group$Group,p=0.7,list=F)
otu_train <- otu[trainIndex, ]
otu_test <- otu[-trainIndex, ]

###

best_params <- GetOptimRF(otu_train,"Group")

rf_train = randomForest(Group ~ .,
                        data=otu_train,
                        mtry = best_params$mtry,
                        ntree = best_params$ntree,
                        importance=TRUE,
                        proximity=TRUE,
                        na.action=na.omit)
rf_train

datamelt = melt(rf_train$err.rate,id.vars="OOB")
train_conf = rf_train[["confusion"]]
oob = sum(train_conf[,"class.error"] * rowSums(train_conf[,1:nrow(train_conf)])) / sum(train_conf[,1:nrow(train_conf)])
oob_title = paste("OOB estimate of error rate:",oob,seq="")
ggplot(datamelt, aes(x=Var1, y=value)) +
  geom_line(aes(color=Var2)) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(color = 'black', fill = 'transparent'),
        legend.position = c(0.99, 0.99), legend.justification = c(1, 1),
        legend.background = element_rect(color = "black")) +
  labs(title = oob_title,x = 'Number of Trees', y = 'Error')+
  scale_color_discrete(name="Group")
# plot(randomForest::margin(rf_train, otu_train$Group), main = '')
tiff('Top30.tiff',,res=300,width = 2000,height = 2000,units='px')
varImpPlot(rf_train,scale = FALSE,main = '')
dev.off()

train_predict = predict(rf_train, otu_train,type = "class")
compare_train = table(train_predict, otu_train$Group)
compare_train
sum(diag(compare_train)/sum(compare_train))

test_predict = predict(rf_train, otu_test,type = "class")
compare_test = table(otu_test$Group, test_predict, dnn = c('Actual', 'Predicted'))
compare_test
sum(diag(compare_test)/sum(compare_test))

imp_otu = as.data.frame(rf_train$importance,check.names = FALSE)
imp_otu = imp_otu[order(imp_otu[,'MeanDecreaseAccuracy'], decreasing = TRUE), ][1:30, ]

tiff('Top30 MDA.tiff',,res=300,width = 3000,height = 2000,units='px')
ggplot(data=imp_otu,aes(x=MeanDecreaseAccuracy,y=reorder(rownames(imp_otu),MeanDecreaseAccuracy)))+
  geom_bar(position = position_dodge(),
           width = 0.5,
           stat = "identity",
           fill="steelblue")+
  theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5),panel.background = element_rect(color = 'black', fill = 'transparent')) +
  labs(title = 'Top 30',x = 'MeanDecreaseAccuracy', y = '')
dev.off()

###############################ROC####################################

pred_p = predict(rf_train, otu_test,type = "prob")
case = 'T'
control = 'C'
roc_obj<-roc(otu_test$Group,pred_p[,case],levels = c(case,control),direction = '>')
auc(roc_obj)
ci(roc_obj)
plot.roc(roc_obj,
         print.auc=TRUE,
         main= paste(case,'(Case) vs ',control,'(Control)',sep=''),
         auc.polygon=TRUE,
         grid=c(0.1, 0.2),
         grid.col=c("black", "black"),
         max.auc.polygon=TRUE,
         auc.polygon.col="orange",
         max.auc.polygon.col="grey",
         print.thres=F
)

############################Multi ROC#################################

pred_p = predict(rf_train, otu_test,type = "prob")

roc_obj<-pROC::multiclass.roc(otu_test$Group,pred_p)
auc(roc_obj)
for(i in 1:length(roc_obj$rocs))
{
  for(j in 1:2)
  {
    roc = roc_obj$rocs[[i]][[j]]
    title = paste(roc$levels[1],'(Case) vs ',roc$levels[2],'(Control)',sep='')
    auc(roc)
    ci(roc)
    tiff(paste(title,'.tiff',sep=''),res=300,width = 1800,height = 1200,units='px')
    plot.roc(roc,
             print.auc=TRUE,
             main= title,
             auc.polygon=TRUE,
             grid=c(0.1, 0.2),
             grid.col=c("black", "black"),
             max.auc.polygon=TRUE,
             auc.polygon.col="orange",
             max.auc.polygon.col="grey",
             print.thres=F
    )
    dev.off()
  }
}

############################ConfMat#################################

pred_r = predict(rf_train, otu_test,type = "response")
confmat = caret::confusionMatrix(data = pred_r,reference = otu_test$Group)
confmat_df = as.data.frame.matrix( confmat$table )
confmat_normalized = round(confmat_df/rowSums(confmat_df),2)
confmat_normalized = replace(confmat_normalized,is.na(confmat_normalized),0)
confmat_normalized$Actual = rownames(confmat_normalized)
confmat_draw = melt(confmat_normalized)
ggplot(data=confmat_draw,aes(Actual,variable,fill=value))+
  geom_tile()+
  geom_text(aes(label=scales::percent(value)))+
  scale_fill_gradient(low="#F0F0F0",high="#3575b5")+
  labs(x="Actual",y="Predicted",title="Confusion Matrix")+
  theme_prism(border = T)+
  theme(legend.position = "none")
write.csv(confmat_df,"confusion_matrix_test.csv")

# conf_error = rf_train$confusion[,'class.error']
# confmat_df = data.frame(rf_train$confusion[,1:nrow(rf_train$confusion)])
# confmat_normalized = round(confmat_df/rowSums(confmat_df),2)
# confmat_normalized = replace(confmat_normalized,is.na(confmat_normalized),0)
# confmat_normalized$Actual = rownames(confmat_normalized)
# confmat_draw = melt(confmat_normalized)
# ggplot(data=confmat_draw,aes(Actual,variable,fill=value))+
#   geom_tile()+
#   geom_text(aes(label=scales::percent(value)))+
#   scale_fill_gradient(low="#F0F0F0",high="#3575b5")+
#   labs(x="Actual",y="Predicted",title="Confusion Matrix")+
#   theme_prism(border = T)+
#   theme(legend.position = "none")
# confmat_df[,'class.error'] = conf_error
# write.csv(confmat_df,"confusion_matrix_train.csv")

##################################Output###############################

name_select = rownames(imp_otu)
name_select = append(name_select,'Group')

otu_out = otu[,name_select]
otu_out = t(otu_out)
otu_out = data.frame(cbind(otu_out,row.names(otu_out)),check.names = FALSE)
colnames(otu_out)[ncol(otu_out)] = '#NAME'
for(i in 1:(nrow(otu_out)-1))
{
  if(!has.key(rownames(otu_out)[i],HashTable))
  {
    print(rownames(otu_out)[i])
  }
  otu_out[i,'#NAME'] = HashTable[[ rownames(otu_out)[i] ]]
}
write.csv(otu_out,"top30.csv")

### Feature Stability

set.seed(114514)
B <- 1000
n_top <- 30

features <- colnames(otu_train)[colnames(otu_train)!="Group"]

MDA_list <- list()
Gini_list <- list()
top_list <- list()

for (b in 1:B) {

  idx <- sample(1:nrow(otu_train), replace = TRUE)
  data_boot <- otu_train[idx, ]

  rf_boot <- randomForest(
    Group ~ .,
    data = data_boot,
    ntree = best_params$ntree,
    mtry = best_params$mtry,
    importance = TRUE
  )
  
  imp_df <- as.data.frame(rf_boot$importance)

  MDA_list[[b]]  <- imp_df$MeanDecreaseAccuracy
  Gini_list[[b]] <- imp_df$MeanDecreaseGini
  
  names(MDA_list[[b]]) <- features
  names(Gini_list[[b]]) <- features

  top_features <- names(sort(MDA_list[[b]], decreasing = TRUE))[1:n_top]
  top_list[[b]] <- top_features
  
  cat("Bootstrap:", b, "/", B, "\n")
}

MDA_mat  <- do.call(rbind, MDA_list)
Gini_mat <- do.call(rbind, Gini_list)

importance_mean_MDA <- colMeans(MDA_mat)
importance_sd_MDA   <- apply(MDA_mat, 2, sd)

importance_mean_Gini <- colMeans(Gini_mat)
importance_sd_Gini   <- apply(Gini_mat, 2, sd)

importance_df <- data.frame(
  Feature = features,
  Mean_MDA = importance_mean_MDA,
  SD_MDA = importance_sd_MDA,
  Mean_Gini = importance_mean_Gini,
  SD_Gini = importance_sd_Gini
)

stability <- sapply(features, function(f){
  sum(sapply(top_list, function(x) f %in% x))
}) / B

stability_df <- data.frame(
  Feature = features,
  Frequency = stability
)

stability_df <- stability_df[order(stability_df$Frequency, decreasing = TRUE), ]

write.csv(importance_df, "Bootstrap_importance_mean_sd.csv", row.names = FALSE)
write.csv(stability_df, "Bootstrap_feature_freq_top30.csv", row.names = FALSE)
