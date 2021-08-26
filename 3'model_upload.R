all2<-read.csv("hunglun_gtag_only_0_1.5_2-fold_candidates_maxent_snp_exct_alt_U2_novelGC_folds_motif5s_len_open_count_dG_pbp_cond_dsvm.csv", head=T)

filter<-read.csv("filter_wt.mt_unspliced.noncanonical_5_percent.csv", head=T)
all7<-all2[which((all2$wt_count.1 == ">100") & (all2$mt_count.1 == ">100")),] #4154
all<-all2[-which((all2$wt_count == "<100") | (all2$mt_count == "<100") | (all2$overall == "0x")| (all2$overall == "1.5x")),] #4350 ## for AG-model

all$AG<-as.character(all$AG)
all$AG[all$AG %in% NA] <- "ref"
all$AG[all$AG %in% "reduce/add_AG" |all$AG %in% "same_AG(indel)"  ] <- "ref"
all$AG<-factor(all$AG, levels=c("ref", "add_AG", "reduce_AG"))
all$database[all$database %in% NA]<-"manual"
all$database<-factor(all$database, levels=c("cosmic", "dbSNP", "ClinVar", "hgmd", "manual"))
all$svm_bp_ct <- as.numeric(as.character(all$svm_bp_ct))
all$X1stbp_position <- factor(all$X1stbp_position, levels = c("-3", "-2", "0"))
all$bp_position.svm <- factor(all$bp_position.svm, levels = c("no", "-5", "-4", "-1", "1", "-3", "-2", "0")) 

all$wt.dG <- as.numeric(all$wt.dG)
all$mt.dG <- as.numeric(all$mt.dG)
all$novel.maxent[all$novel.maxent %in% NA]<- mean(all$novel.maxent, na.rm=T)
all$novel.svm_scr[all$novel.svm_scr %in% NA]<- mean(all$novel.svm_scr, na.rm=T)
all$mt.ex.in._delta_gc.percent.[all$mt.ex.in._delta_gc.percent. %in% NA] <- mean(all$mt.ex.in._delta_gc.percent., na.rm=T)
all$mt_splicing_efficiency[all$mt_splicing_efficiency %in% NA] <- mean(all$mt_splicing_efficiency, na.rm=T)
all$maxent.m[all$maxent.m %in% NA]<- mean(all$maxent.m, na.rm=T)
all$mm.m[all$mm.m %in% NA]<- mean(all$mm.m, na.rm=T)
all$wmm.m[all$wmm.m %in% NA]<- mean(all$wmm.m, na.rm=T)
all$ss_dist.m[all$ss_dist.m %in% NA]<- mean(all$ss_dist.m, na.rm=T)
all$bp_scr.m[all$bp_scr.m %in% NA]<- mean(all$bp_scr.m, na.rm=T)
all$ppt_scr.m[all$ppt_scr.m %in% NA]<- mean(all$ppt_scr.m, na.rm=T)
all$svm_scr.m[all$svm_scr.m %in% NA]<- mean(all$svm_scr.m, na.rm=T)
all$dsvm[all$dsvm %in% NA]<- mean(all$dsvm, na.rm=T)
all$svm_scr.m[all$svm_scr.m %in% NA]<- mean(all$svm_scr.m, na.rm=T)
all$novel_mt.maxent<- as.numeric(all$novel_mt.maxent)
all$novel_mt.mm<- as.numeric(all$novel_mt.mm)
all$novel_mt.wmm<- as.numeric(all$novel_mt.wmm)
all$delta_svm<-as.numeric(all$delta_svm)
all$ex.count[all$ex.count %in% NA]<- mean(all$ex.count, na.rm=T)
all$PhyloP30[all$PhyloP30 %in% NA]<- mean(all$PhyloP30, na.rm=T)
all$PhyloP20[all$PhyloP20 %in% NA]<- mean(all$PhyloP20, na.rm=T)
all$PhyloP7[all$PhyloP7 %in% NA]<- mean(all$PhyloP7, na.rm=T)

all$variant_type[all$variant_type %in% "deletion"] <-"indel"
all$variant_type[all$variant_type %in% "insertion"] <-"indel"
all$variant_type[all$variant_type %in% NA] <-"indel"
all$variant_type<-factor(all$variant_type, levels = c("SNP", "indel"))
all$same_bp[all$same_bp %in% NA] <- "FALSE"
all$simp.variant_class[all$simp.variant_class %in% NA] <- "undefined"
all$simp.variant_class <- factor(all$simp.variant_class, levels= c("undefined", "Uncertain significance", "Benign/Likely benign", 
"Functional evidence", "Pathogenic/Likely pathogenic"))
all$Anumber.10[all$Anumber.10 %in% NA|all$Anumber.10 %in% "0"] <- "no"
all$wt.mt_3ss_distance.right. <- -all$wt.mt_3ss_distance.right.
all$bp.right<- -all$bp.right
all$overall[which(all$mt %in% filter$doubt_5per)] <- "FALSE"

AG1<-all[which(all$AG == "add_AG"),] #317 #5percent:315 #-100-multi 283
AG2<-AG1[,c(1:105, 315:321, 510:549)]
AG3<-AG1[,c(106:314, 322:509)] 
AG4<-AG3[, colSums(AG3, na.rm=T)>2]
AG5<-AG2[,-c(1:22, 29, 31:33, 35:50, 52, 63, 68, 89, 95:97, 99, 108:109, 113:114, 118, 122, 126, 130, 134, 147:152)] #with exp bp(34) and U2(90:91)
AG<-cbind(AG5, AG4)
AG$overall <- factor(AG$overall, levels=c("FALSE", "2x")) #267, 50 #269, 48
AG6<-AG[apply(is.na(AG), 1, sum)==0,]
dim(AG6)
dim(AG) #283 395

###################standardized AG model###############################################################
library(caret)
library(rsample)
library(glmnet)
library(dplyr)
library(ggplot2)

AG_model<-model.matrix(object = overall ~ ., data =  AG )[, -1]
sds<-apply(AG_model, 2, sd)
sds[sds==0]<-1
AG_model_std<-t(t(AG_model)/sds)

sampsize<-dim(AG_model)[1]
set.seed(37)
index<-sample(1:sampsize, replace=FALSE)
sel.train<-index[1:round(sampsize*2/3)]
sel.test<-index[(round(sampsize*2/3)+1):sampsize]
AG_train_x<-AG_model_std[sel.train,]
AG_test_x<-AG_model_std[sel.test,]
AG_train_y<-as.numeric(AG$overall[sel.train])
AG_test_y<-as.numeric(AG$overall[sel.test])

# create a lasso model
library(pROC)
set.seed(555)
cv_lasso <- cv.glmnet(x = AG_train_x, y = AG_train_y, alpha = 1, family = "binomial")
pred <- predict(cv_lasso, newx = AG_test_x, s = cv_lasso$lambda.min, type="response")
test_roc<-roc(AG_test_y ~ pred, plot=TRUE, print.auc=TRUE)  
pdf(file="AG_ROC881_save.pdf", width = 4, height = 4)

coef(cv_lasso, s = "lambda.min") %>% # 308 x 1 sparse Matrix of class "dgCMatrix"
  as.matrix() %>% 
  as.data.frame() %>% 
  add_rownames(var = "var") %>% 
  `colnames<-`(c("var","coef")) %>%
  filter(var != "(Intercept)") %>%  
  top_n(10, wt = abs(coef)) %>% 
  ggplot(aes(coef, reorder(var, coef))) +
  geom_point() +
  ggtitle("Top 10 influential variables") +
  xlab("Coefficient") +
  ylab(NULL)

#obtain coefficient and intercept
coef<-coef(cv_lasso, s = "lambda.min")
AG_model<-model.matrix(object = overall ~ ., data =  AG7 )[, -1] #only columns in AG120

###########non-standard AG model#################################
AG_model<-model.matrix(object = overall ~ ., data =  AG )[, -1]
sampsize<-dim(AG_model)[1]
set.seed(37)
index<-sample(1:sampsize, replace=FALSE)
sel.train<-index[1:round(sampsize*2/3)]
sel.test<-index[(round(sampsize*2/3)+1):sampsize]
AG_train_x<-AG_model[sel.train,]
AG_test_x<-AG_model[sel.test,]
AG_train_y<-as.numeric(AG7$overall[sel.train])
AG_test_y<-as.numeric(AG7$overall[sel.test])


# create a lasso model
library(pROC)
set.seed(555)
cv_lasso <- cv.glmnet(x = AG_train_x, y = AG_train_y, alpha = 1, family = "binomial")
pred <- predict(cv_lasso, newx = AG_test_x, s = cv_lasso$lambda.min, type="response")
test_roc<-roc(AG_test_y ~ pred, plot=TRUE, print.auc=TRUE)  
pdf(file="AG_ROC881_save.pdf", width = 4, height = 4)

#find the best cutoff
ag.index<-roc(AG_test_y ~ pred, plot=TRUE, print.auc=TRUE) %>%
coords(transpose = FALSE) %>%
filter(sensitivity > 0.6,
	 specificity > 0.6)
ag.index[(ag.index[,2]+ag.index[,3])== max(ag.index[,2]+ag.index[,3]),]

#obtain coefficient and intercept
coef<-coef(cv_lasso, s = "lambda.min")

################train nonAG model#################################
all<-all7[-which((all7$X1stbp_position %in% NA)| (all7$overall == "0x")| (all7$overall == "1.5x")),] #remove multiple bp #3970 551

all$AG<-as.character(all$AG)
all$AG[all$AG %in% NA] <- "ref"
all$AG[all$AG %in% "reduce/add_AG" |all$AG %in% "same_AG(indel)"  ] <- "ref"
all$AG<-factor(all$AG, levels=c("ref", "add_AG", "reduce_AG"))
all$database[all$database %in% NA]<-"manual"
all$database<-factor(all$database, levels=c("cosmic", "dbSNP", "ClinVar", "hgmd", "manual"))
all$svm_bp_ct <- as.numeric(as.character(all$svm_bp_ct))
all$X1stbp_position <- factor(all$X1stbp_position, levels = c("-3", "-2", "0"))
all$bp_position.svm <- factor(all$bp_position.svm, levels = c("no", "-5", "-4", "-1", "1", "-3", "-2", "0")) 

all$wt.dG <- as.numeric(all$wt.dG)
all$mt.dG <- as.numeric(all$mt.dG)
all$novel.maxent[all$novel.maxent %in% NA]<- mean(all$novel.maxent, na.rm=T)
all$novel.svm_scr[all$novel.svm_scr %in% NA]<- mean(all$novel.svm_scr, na.rm=T)
all$mt.ex.in._delta_gc.percent.[all$mt.ex.in._delta_gc.percent. %in% NA] <- mean(all$mt.ex.in._delta_gc.percent., na.rm=T)
all$mt_splicing_efficiency[all$mt_splicing_efficiency %in% NA] <- mean(all$mt_splicing_efficiency, na.rm=T)
all$maxent.m[all$maxent.m %in% NA]<- mean(all$maxent.m, na.rm=T)
all$mm.m[all$mm.m %in% NA]<- mean(all$mm.m, na.rm=T)
all$wmm.m[all$wmm.m %in% NA]<- mean(all$wmm.m, na.rm=T)
all$ss_dist.m[all$ss_dist.m %in% NA]<- mean(all$ss_dist.m, na.rm=T)
all$bp_scr.m[all$bp_scr.m %in% NA]<- mean(all$bp_scr.m, na.rm=T)
all$ppt_scr.m[all$ppt_scr.m %in% NA]<- mean(all$ppt_scr.m, na.rm=T)
all$svm_scr.m[all$svm_scr.m %in% NA]<- mean(all$svm_scr.m, na.rm=T)
all$dsvm[all$dsvm %in% NA]<- mean(all$dsvm, na.rm=T)
all$svm_scr.m[all$svm_scr.m %in% NA]<- mean(all$svm_scr.m, na.rm=T)
all$novel_mt.maxent<- as.numeric(all$novel_mt.maxent)
all$novel_mt.mm<- as.numeric(all$novel_mt.mm)
all$novel_mt.wmm<- as.numeric(all$novel_mt.wmm)
all$delta_svm<-as.numeric(all$delta_svm)
all$ex.count[all$ex.count %in% NA]<- mean(all$ex.count, na.rm=T)
all$PhyloP30[all$PhyloP30 %in% NA]<- mean(all$PhyloP30, na.rm=T)
all$PhyloP20[all$PhyloP20 %in% NA]<- mean(all$PhyloP20, na.rm=T)
all$PhyloP7[all$PhyloP7 %in% NA]<- mean(all$PhyloP7, na.rm=T)

all$variant_type[all$variant_type %in% "deletion"] <-"indel"
all$variant_type[all$variant_type %in% "insertion"] <-"indel"
all$variant_type[all$variant_type %in% NA] <-"indel"
all$variant_type<-factor(all$variant_type, levels = c("SNP", "indel"))
all$same_bp[all$same_bp %in% NA] <- "FALSE"
all$simp.variant_class[all$simp.variant_class %in% NA] <- "undefined"
all$simp.variant_class <- factor(all$simp.variant_class, levels= c("undefined", "Uncertain significance", "Benign/Likely benign", 
"Functional evidence", "Pathogenic/Likely pathogenic"))
all$Anumber.10[all$Anumber.10 %in% NA|all$Anumber.10 %in% "0"] <- "no"
all$wt.mt_3ss_distance.right. <- -all$wt.mt_3ss_distance.right.
all$bp.right<- -all$bp.right
all$overall[which(all$mt %in% filter$doubt_5per)] <- "FALSE"
AG1<-all[which(all$AG == "add_AG"),] 
nonAG<-setdiff(all, AG1)
all3<-nonAG
all3$overall <- factor(all3$overall, levels=c("FALSE", "2x"))
all4<-all3[,-c(1:23, 29, 31, 33, 35:50, 52, 63, 68, 73, 79:82, 84, 86, 89, 92:318, 510:511, 515, 519, 523, 524:540, 544:555)]
all5<-all4[apply(is.na(all4), 1, sum)==0,] # remove manually selected multiple bp candidates (14=7x2 introns)

nonAG_model<-model.matrix(object = overall ~ ., data =  all5 )[, -1]
sds<-apply(nonAG_model, 2, sd)
sds[sds==0]<-1
nonAG_model_std<-t(t(nonAG_model)/sds)

sampsize<-dim(nonAG_model)[1]
set.seed(1010)
index<-sample(1:sampsize, replace=FALSE)
sel.train<-index[1:round(sampsize*2/3)]
sel.test<-index[(round(sampsize*2/3)+1):sampsize]
nonAG_train_x<-nonAG_model_std[sel.train,]
nonAG_test_x<-nonAG_model_std[sel.test,]
nonAG_train_y<-as.numeric(all5$overall[sel.train])
nonAG_test_y<-as.numeric(all5$overall[sel.test])

library(pROC)
set.seed(27) #w/o wt_splicing_eff

cv_lasso <- cv.glmnet(x = nonAG_train_x, y = nonAG_train_y, alpha = 1, family = "binomial")
pred <- predict(cv_lasso, newx = nonAG_test_x, s = cv_lasso$lambda.min, type="response")
test_roc<-roc(nonAG_test_y ~ pred, plot=TRUE, print.auc=TRUE)

coef(cv_lasso, s = "lambda.min") %>% # 308 x 1 sparse Matrix of class "dgCMatrix"
  as.matrix() %>% 
  as.data.frame() %>% 
  add_rownames(var = "var") %>% 
  `colnames<-`(c("var","coef")) %>%
  filter(var != "(Intercept)") %>% 
  top_n(10, wt = abs(coef)) %>% 
  ggplot(aes(coef, reorder(var, coef))) +
  geom_point() +
  ggtitle("Top 10 influential variables") +
  xlab("Coefficient") +
  ylab(NULL)

###########non-standard nonAG model##################
all4<-all3[,-c(1:23, 28:29, 31:33, 35:50, 52, 63, 68, 73, 79:82, 84, 86,  89, 92:316, 322:511, 515, 519, 523, 524:540, 544:555)] # noRBP
all5<-all4[apply(is.na(all4), 1, sum)==0,]
nonAG_model<-model.matrix(object = overall ~ ., data =  all5)[, -1]

sampsize<-dim(nonAG_model)[1]
set.seed(1010)
index<-sample(1:sampsize, replace=FALSE)
sel.train<-index[1:round(sampsize*2/3)]
sel.test<-index[(round(sampsize*2/3)+1):sampsize]
nonAG_train_x<-nonAG_model[sel.train,]
nonAG_test_x<-nonAG_model[sel.test,]
nonAG_train_y<-as.numeric(all5$overall[sel.train])
nonAG_test_y<-as.numeric(all5$overall[sel.test])

library(pROC)
set.seed(27) #w/o wt_splicing_eff

cv_lasso <- cv.glmnet(x = nonAG_train_x, y = nonAG_train_y, alpha = 1, family = "binomial")
pred <- predict(cv_lasso, newx = nonAG_test_x, s = cv_lasso$lambda.min, type="response")
test_roc<-roc(nonAG_test_y ~ pred, plot=TRUE, print.auc=TRUE)

# Returns thresholds, sensitivities and specificities:
index<-roc(nonAG_test_y ~ pred, plot=TRUE, print.auc=TRUE) %>%
coords(transpose = FALSE) %>%
filter(sensitivity > 0.6,
	 specificity > 0.6)
index[(index[,2]+index[,3])== max(index[,2]+index[,3]),]
max(index[,2]+index[,3])

