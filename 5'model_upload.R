all1<-read.csv("snv_5ss_snp_gt_maxent_gc_patho_dG_alt_len_rbp_con_3ss.csv", head=T)
all1$in1.len[all1$in1.len %in% NA] <- round(mean(all1$in1.len, na.rm=T))
all<-na.omit(all1)

all$delta.gt<-factor(all$delta.gt, levels=c("Ref", "Add_GT", "Reduce_GT"))
all$patho.s <- factor(all$patho.s, levels= c("no", "Uncertain significance", "Benign/Likely benign", "Pathogenic/Likely pathogenic"))
all2 <- all[,92:1211]
all3 <- all2[, abs(colSums(all2, na.rm=T))>5]
all4 <- cbind(all[,c(1:91, 1212:1214)], all3)

##################################################################################
library(caret)
library(rsample)
library(glmnet)
library(dplyr)
library(ggplot2)
all4$is_sdv <- factor(all4$is_sdv, levels=c("FALSE", "TRUE"))

all5<-all4[,-c(1:17, 18:21, 23, 25:27, 33:39, 40, 53:65, 67:69, 92:915)]  
all6<-all5[apply(is.na(all5), 1, sum)==0,] # remove manually selected multiple bp candidates (14=7x2 introns)

ss5_model<-model.matrix(object = is_sdv ~ ., data =  all5 )[, -1]
sds<-apply(ss5_model, 2, sd)
sds[sds==0]<-1
ss5_model_std<-t(t(ss5_model)/sds)

sampsize<-dim(ss5_model)[1]
set.seed(555)

index<-sample(1:sampsize, replace=FALSE)
sel.train<-index[1:round(sampsize*2/3)]
sel.test<-index[(round(sampsize*2/3)+1):sampsize]
ss5_train_x<-ss5_model_std[sel.train,]
ss5_test_x<-ss5_model_std[sel.test,]
ss5_train_y<-as.numeric(all5$is_sdv[sel.train])
ss5_test_y<-as.numeric(all5$is_sdv[sel.test])
dim(ss5_test_x)
dim(ss5_train_x)
length(ss5_test_y)
length(ss5_train_y)

# create a lasso model
library(pROC)
set.seed(1010) #w/o wt_splicing_eff

cv_lasso <- cv.glmnet(x = ss5_train_x, y = ss5_train_y, alpha = 1, family = "binomial")
pred <- predict(cv_lasso, newx = ss5_test_x, s = cv_lasso$lambda.min, type="response")
test_roc<-roc(ss5_test_y ~ pred, plot=TRUE, print.auc=TRUE)

coef(cv_lasso, s = "lambda.min") %>% # 308 x 1 sparse Matrix of class "dgCMatrix"
  as.matrix() %>% 
  as.data.frame() %>% 
  add_rownames(var = "var") %>% 
  `colnames<-`(c("var","coef")) %>%
  filter(var != "(Intercept)") %>%  #­ç°£ºI¶Z¶µ
  top_n(10, wt = abs(coef)) %>% 
  ggplot(aes(coef, reorder(var, coef))) +
  geom_point() +
  ggtitle("Top 10 influential variables") +
  xlab("Coefficient") +
  ylab(NULL)

coef.non<-coef(cv_lasso, s = "lambda.min")

######non-STD model for validation#############################
all5<-all4[,-c(1:17, 18:23, 25:27, 33:39, 40, 53:65, 67:69, 95:915)] #22: nat_index == wt_splicing_eff ##remove RBP 0.946
non_model<-model.matrix(object = is_sdv ~ ., data =  all5 )[, -1]

sampsize<-dim(non_model)[1]
set.seed(555)

index<-sample(1:sampsize, replace=FALSE)
sel.train<-index[1:round(sampsize*2/3)]
sel.test<-index[(round(sampsize*2/3)+1):sampsize]
non_train_x<-non_model[sel.train,]
non_test_x<-non_model[sel.test,]
non_train_y<-as.numeric(all6$is_sdv[sel.train])
non_test_y<-as.numeric(all6$is_sdv[sel.test])
dim(non_test_x)
dim(non_train_x)
length(non_test_y)
length(non_train_y)

# create a lasso model
library(pROC)
set.seed(1010) #w/o wt_splicing_eff

cv_lasso <- cv.glmnet(x = non_train_x, y = non_train_y, alpha = 1, family = "binomial")
pred <- predict(cv_lasso, newx = non_test_x, s = cv_lasso$lambda.min, type="response")
test_roc<-roc(non_test_y ~ pred, plot=TRUE, print.auc=TRUE)

# Returns thresholds, sensitivities and specificities:
index<-roc(non_test_y ~ pred, plot=TRUE, print.auc=TRUE) %>%
coords(transpose = FALSE) %>%
filter(sensitivity > 0.6,
	 specificity > 0.6)
index[(index[,2]+index[,3])== max(index[,2]+index[,3]),]
max(index[,2]+index[,3])

coef.non<-coef(cv_lasso, s = "lambda.min")
