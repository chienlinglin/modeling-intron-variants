library(tidyr)
library(stringr)
library(caret)
library(rsample)  # data splitting 
library(glmnet)   # implementing regularized regression approaches
library(dplyr)    # basic data manipulation procedures
library(ggplot2)  # plotting

setwd("/Volumes/TOSHIBA EXT/yiting/Hung Lun/home2_cllin_20191223/aligncomrev_0707/fisher/summarize/gtag_only/")

####### new svm file
all <- read.csv(file = "/Volumes/TOSHIBA EXT/for_server/new_svm_0719/hunglun_gtag_only_0_1.5_2-fold_candidates_maxent_snp_exct_cons7_alt_U2_novelGC_folds_motif5s_len_open_count_dG_pbp_dsvm.csv",row.names = 1,stringsAsFactors = F)

## from dr.lin added attributes
#all <- read.csv(file = "by_cllin/20201026_factors_for_glm/hunglun_gtag_only_0_1.5_2-fold_candidates_maxent_svm_snp_ct_cons_alt_U2_novelgc_fold_corrected.csv",row.names = 1,stringsAsFactors = F)
all <- all[-which(is.na(all$wt.mt_3ss_distance.right.)),]

### remove objects with read counts less than 100
over100 <- read.csv(file = "Recheck_100count/hunglun_wt.mt.count_recheck_100_count.csv",row.names = 1)
over100$wt.mt <- NA
over100$wt.mt[which(over100$wt_count==">100" & over100$mt_count==">100")] <- TRUE  # 4215

all <- all[-which(all$mt %in% over100$mt[which(is.na(over100$wt.mt))]),]


mttype <- all[,c(1,2)]
mttype$l <- mttype$mt
mttype <- separate(mttype,"l", into=c("l","r"),sep="_-")
mttype <- separate(mttype,"r", into=c("ref","alt"),sep="to")
ref_nchar <- nchar(str_extract(mttype$ref,"[A-Z]+"), type = "chars", allowNA = FALSE, keepNA = NA)
alt_nchar <- nchar(str_extract(mttype$alt,"[A-Z]+"), type = "chars", allowNA = FALSE, keepNA = NA)
mttype$ref_nchar <- as.numeric(ref_nchar)
mttype$alt_nchar <- as.numeric(alt_nchar)
mttype$ref_nchar[is.na(mttype$ref_nchar)] <- 0
mttype$alt_nchar[is.na(mttype$alt_nchar)] <- 0
mttype$ref_nchar <- as.numeric(mttype$ref_nchar)
mttype$alt_nchar <- as.numeric(mttype$alt_nchar)
mttype$type <- NA
mttype$type[which(mttype$ref_nchar<mttype$alt_nchar)] <- "insertion"
mttype$type[which(mttype$ref_nchar>mttype$alt_nchar)] <- "deletion"
mttype$type[which(mttype$ref_nchar==mttype$alt_nchar)] <- "SNP"
mttype$ref <- str_replace_all(mttype$ref,pattern = "[0-9]",replacement = "")
mttype$wt.snp <- mttype$ref
mttype$mt.snp <- mttype$alt
mttype$wt.snp[which(mttype$ref_nchar<mttype$alt_nchar)] <- "insertion"
mttype$mt.snp[which(mttype$ref_nchar<mttype$alt_nchar)] <- "insertion"
mttype$wt.snp[which(mttype$ref_nchar>mttype$alt_nchar)] <- "deletion"
mttype$mt.snp[which(mttype$ref_nchar>mttype$alt_nchar)] <- "deletion"
all$wt.snp <- mttype$wt.snp
all$mt.snp <- mttype$mt.snp

all$bp.snp <- NA
all$bp.snp[which(all$bp_position=="0")] <- str_sub(all$zero[which(all$bp_position=="0")],start = 6,end = 6)   # 1747
all$bp.snp[which(all$bp_position=="-2")] <- str_sub(all$minus2[which(all$bp_position=="-2")],start = 6,end = 6) # 1137
all$bp.snp[which(all$bp_position=="-3")] <- str_sub(all$minus3[which(all$bp_position=="-3")],start = 6,end = 6) #1378
all$bp.snp[is.na(all$bp.snp)] <- "multiple" #776


## original file
all$AG<-as.character(all$AG)
all$AG[all$AG %in% NA] <- "ref"
all$AG[all$AG %in% "reduce/add_AG" | all$AG %in% "same_AG(indel)"  ] <- "ref"
all$AG<-factor(all$AG, levels=c("ref", "add_AG", "reduce_AG"))
all$database[all$database %in% NA]<-"manual"
all$database<-factor(all$database, levels=c("cosmic", "dbSNP", "ClinVar", "hgmd"))
all$svm_bp_ct <- as.numeric(as.character(all$svm_bp_ct))
all$X1stbp_position <- factor(all$X1stbp_position, levels = c("-3", "-2", "0"))

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
all$mt.Anumber[all$mt.Anumber %in% NA]<- mean(all$mt.Anumber, na.rm=T)
all$U2.mfe.cmt[all$U2.mfe.cmt %in% NA]<- mean(all$U2.mfe.cmt, na.rm=T)


all$wt.snp <- factor(all$wt.snp, levels=c("C", "insertion","deletion", "T", "G", "A"))
all$mt.snp <- factor(all$mt.snp, levels=c("C", "insertion","deletion","T", "G", "A"))
all$bp.snp <- factor(all$bp.snp, levels=c("C", "T", "G", "A","multiple"))
all$same_bp[all$same_bp %in% NA] <- "FALSE"
all$simp.variant_class[all$simp.variant_class %in% NA] <- "undefined"
all$simp.variant_class <- factor(all$simp.variant_class, levels= c("undefined", "Uncertain significance", "Benign/Likely benign", 
                                                                   "Functional evidence", "Pathogenic/Likely pathogenic"))
all$Anumber.10[all$Anumber.10 %in% NA] <- "no"
fil5per <- read.csv(file = "Check_count_percentage/filter_wt.mt_unspliced.noncanonical_5_percent.csv",row.names = 1)
all[which(all$mt %in% fil5per[,1]),48] <- NA


all1 <- all

all1$wt_splicing_efficiency2 <- all1$wt_splicing_efficiency/100
all1$mt_splicing_efficiency2 <- all1$mt_splicing_efficiency/100
all1$delta_splicing_efficency2 <- all1$mt_splicing_efficiency2 - all1$wt_splicing_efficiency2

### log transformation for normal distributed data for modeling
all1$wt_splicing_efficiency2[which(all1$wt_splicing_efficiency2==1)] <- 0.999
all1$wt_splicing_efficiency2[which(all1$wt_splicing_efficiency2==0)] <- 0.001
all1$mt_splicing_efficiency2[which(all1$mt_splicing_efficiency2==1)] <- 0.999
all1$mt_splicing_efficiency2[which(all1$mt_splicing_efficiency2==0)] <- 0.001
all1$mt_splicing_efficiency3 <- log(all1$mt_splicing_efficiency2/(1-all1$mt_splicing_efficiency2))


which(all$mt_splicing_efficiency== 0) #62
all1 <- all1[-which(all1$mt_splicing_efficiency2== 0.001 ),] # 4507
hist(all1$mt_splicing_efficiency2)
hist(all1$mt_splicing_efficiency3)


#################################################################################
### standardized non-AG model with wt efficiency
remove <- which(colSums(all1[,c(105:313,321:508)])<10)

# factors that have only one value
values_count <- sapply(lapply(all1, unique), length)  # Identify variables with 1 value
remove2 <- which(values_count==1)

# non ag have no novel 3'ss related factors
remove3 <- which(str_detect(colnames(all1),pattern = "novel"))

# 94:96 novel_mt/maxent/mm/wmm
allf <- all1[-which(all1$AG=="add_AG"),-unique(c(1:22,24,28:50,51,62,67,72,79,94:96,98,105:314,316,317,509,510,540,543:548,549:553,555:557,remove,remove2,remove3))]
which(is.na(allf), arr.ind=TRUE) 
allf <- na.omit(allf)
# lacking some phylop score
#row col
#1288  941 222
#1288  941 223
#3075 2148 223
which(is.na(allf), arr.ind=TRUE) 

# Create training (70%) and test (30%) sets 
# Use set.seed for reproducibility
set.seed(123)

eff_model <- model.matrix(object = mt_splicing_efficiency3 ~ ., data = allf )[, -1]
sds<-apply(eff_model, 2, sd)
sds[sds==0]<-1
eff_model_std<-t(t(eff_model)/sds)

sampsize<-dim(eff_model)[1]
set.seed(37)
index<-sample(1:sampsize, replace=FALSE)
sel.train<-index[1:round(sampsize*2/3)]
sel.test<-index[(round(sampsize*2/3)+1):sampsize]
eff_train_x<-eff_model_std[sel.train,]
eff_test_x<-eff_model_std[sel.test,]

eff_train_y<-as.numeric(allf$mt_splicing_efficiency3[sel.train])
eff_test_y<-as.numeric(allf$mt_splicing_efficiency3[sel.test])
dim(eff_test_x)  # 1280 264
dim(eff_train_x) # 2559 264

# create a lasso model
cv_lasso <- cv.glmnet(x = eff_train_x, y = eff_train_y, alpha = 1, family = "gaussian")

# Compute R^2 from true and predicted values
eval_results <- function(true, predicted, df) {
  SSE <- sum((predicted - true)^2)
  SST <- sum((true - mean(true))^2)
  R_square <- 1 - SSE / SST
  RMSE = sqrt(SSE/nrow(df))
  
  # Model performance metrics
  data.frame(
    RMSE = RMSE,
    Rsquare = R_square
  )
  
}

# Prediction and evaluation on train data
predictions_train <- predict(cv_lasso, s = cv_lasso$lambda.min, newx = eff_train_x)
eval_results(eff_train_y, predictions_train, eff_train_x)

# Prediction and evaluation on test data
predictions_test <- predict(cv_lasso, s = cv_lasso$lambda.min, newx = eff_test_x)
eval_results(eff_test_y, predictions_test, eff_test_x)
#RMSE  Rsquare
# 0.5384975 0.591142

coef(cv_lasso, s="lambda.min")
write.csv(as.matrix(coef(cv_lasso, s="lambda.min")),"non_ag_wt_eff_model_coef_0519.csv")



coef(cv_lasso, s = "lambda.min") %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "var") %>% 
  `colnames<-`(c("var","coef")) %>%
  filter(var != "(Intercept)") %>%  
  top_n(15, wt = abs(coef)) %>% 
  ggplot(aes(coef, reorder(var, coef))) +
  geom_point() +
  ggtitle("Top 15 influential variables in non-ag splicing efficiency(wt_eff_add)") +
  xlab("Coefficient") +
  ylab(NULL)
dev.off()


final <- data.frame(cbind(eff_test_y, predictions_test))
colnames(final)[2] <- "predicted"

pear_cor <- round(cor(final$predicted, final$eff_test_y), 3) #0.77
spear_cor <- round(cor(final$predicted, final$eff_test_y,method = "spearman"), 3) #0.75
rmselabel <- paste0("R square= 0.60", "\n RMSE = 0.54")
p <- ggplot(data=final, aes(x=predicted, y=eff_test_y)) + geom_point(colour = "gray") + ggtitle("non_ag_eff(+wt eff): real vs prediction") + geom_smooth(method=lm, se=T,colour="peachpuff4")+
  theme_light()+annotate(geom = "text",label=rmselabel,x=0,y=3.3,colour="lightsalmon",fontface = "bold") + annotate(geom = "text",label="Pearson's p = 0.77",x=1.9,y=-1.7,colour="peachpuff4",fontface = "bold")

ggsave("non_AG_no_wt_eff_standardized_efficiency_performance.pdf", width = 4, height = 4)

#################################################################################
### standardized non-AG model without wt efficiency
remove <- which(colSums(all1[,c(105:313,321:508)])<10)

# factors that have only one value
values_count <- sapply(lapply(all1, unique), length)  # Identify variables with 1 value
remove2 <- which(values_count==1)

# non ag have no novel 3'ss related factors
remove3 <- which(str_detect(colnames(all1),pattern = "novel"))

# 94:96 novel_mt/maxent/mm/wmm
allf <- all1[-which(all1$AG=="add_AG"),-unique(c(1:22,24,27:50,51,62,67,72,79,94:96,98,105:314,316,317,509,510,540,543:548,549:553,555:557,remove,remove2,remove3))]
which(is.na(allf), arr.ind=TRUE) 
allf <- na.omit(allf)
# lacking some phylop score
#row col
#1288  941 235
#1288  941 236
#3075 2148 236
which(is.na(allf), arr.ind=TRUE) 

### remove RBPs
allf <- allf[,-(42:191)]

# Create training (70%) and test (30%) sets 
# Use set.seed for reproducibility
set.seed(123)

eff_model <- model.matrix(object = mt_splicing_efficiency3 ~ ., data = allf )[, -1]
sds<-apply(eff_model, 2, sd)
sds[sds==0]<-1
eff_model_std<-t(t(eff_model)/sds)

sampsize<-dim(eff_model)[1]
set.seed(37)
index<-sample(1:sampsize, replace=FALSE)
sel.train<-index[1:round(sampsize*2/3)]
sel.test<-index[(round(sampsize*2/3)+1):sampsize]
eff_train_x<-eff_model_std[sel.train,]
eff_test_x<-eff_model_std[sel.test,]

eff_train_y<-as.numeric(allf$mt_splicing_efficiency3[sel.train])
eff_test_y<-as.numeric(allf$mt_splicing_efficiency3[sel.test])
dim(eff_test_x)  # 1279 244
dim(eff_train_x) # 2558 244


# create a lasso model
cv_lasso <- cv.glmnet(x = eff_train_x, y = eff_train_y, alpha = 1, family = "gaussian")


# Compute R^2 from true and predicted values
eval_results <- function(true, predicted, df) {
  SSE <- sum((predicted - true)^2)
  SST <- sum((true - mean(true))^2)
  R_square <- 1 - SSE / SST
  RMSE = sqrt(SSE/nrow(df))
  
  # Model performance metrics
  data.frame(
    RMSE = RMSE,
    Rsquare = R_square
  )
  
}

# Prediction and evaluation on train data
predictions_train <- predict(cv_lasso, s = cv_lasso$lambda.min, newx = eff_train_x)
eval_results(eff_train_y, predictions_train, eff_train_x)

# Prediction and evaluation on test data
predictions_test <- predict(cv_lasso, s = cv_lasso$lambda.min, newx = eff_test_x)
eval_results(eff_test_y, predictions_test, eff_test_x)
#RMSE   Rsquare
#0.7442244 0.2190684

coef(cv_lasso, s="lambda.min")
write.csv(as.matrix(coef(cv_lasso, s="lambda.min")),"non_ag_wt_eff_model_coef_0519.csv")


coef(cv_lasso, s = "lambda.min") %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "var") %>% 
  `colnames<-`(c("var","coef")) %>%
  filter(var != "(Intercept)") %>%  
  top_n(15, wt = abs(coef)) %>% 
  ggplot(aes(coef, reorder(var, coef))) +
  geom_point() +
  ggtitle("Top 15 influential variables in non-ag splicing efficiency(wt_eff_add)") +
  xlab("Coefficient") +
  ylab(NULL)
dev.off()


final <- data.frame(cbind(eff_test_y, predictions_test))
colnames(final)[2] <- "predicted"

pear_cor <- round(cor(final$predicted, final$eff_test_y), 3) #0.47
spear_cor <- round(cor(final$predicted, final$eff_test_y,method = "spearman"), 3) #0.46
rmselabel <- paste0("R square= 0.22", "\n RMSE = 0.74")
p <- ggplot(data=final, aes(x=predicted, y=eff_test_y)) + geom_point(colour = "gray") + ggtitle("non_ag_eff(+wt eff): real vs prediction") + geom_smooth(method=lm, se=T,colour="peachpuff4")+
  theme_light()+annotate(geom = "text",label=rmselabel,x=0,y=3.3,colour="lightsalmon",fontface = "bold") + annotate(geom = "text",label="Pearson's p = 0.47",x=1.9,y=-1.7,colour="peachpuff4",fontface = "bold")
ggsave("non_AG_no_wt_eff_standardized_efficiency_performance.pdf", width = 4, height = 4)

#################################################################################
### non-standardized non-AG model without wt efficiency for server 

### Remove RBPs motifs appears <10 times across all non-AG sequences
remove <- which(colSums(all1[,c(105:313,321:508)])<10)

# factors that have only one value
values_count <- sapply(lapply(all1, unique), length)  # Identify variables with 1 value
remove2 <- which(values_count==1)

# non ag have no novel 3'ss related factors
remove3 <- which(str_detect(colnames(all1),pattern = "novel"))

# 94:96 novel_mt/maxent/mm/wmm
allf <- all1[-which(all1$AG=="add_AG"),-unique(c(1:22,24,27:50,51,62,67,72,79,94:96,98,105:314,316,317,509,510,540,543:548,549:553,555:557,remove,remove2,remove3))]
which(is.na(allf), arr.ind=TRUE) 
allf <- na.omit(allf)
# lacking some phylop score
#row col
#1288  941 235
#1288  941 236
#3075 2148 236
which(is.na(allf), arr.ind=TRUE) 

### remove RBPs
allf <- allf[,-(42:191)]


# Create training (70%) and test (30%) sets 
# Use set.seed for reproducibility
set.seed(123)

eff_model <- model.matrix(object = mt_splicing_efficiency3 ~ ., data = allf )[, -1]

sampsize<-dim(eff_model)[1]
set.seed(37)
index<-sample(1:sampsize, replace=FALSE)
sel.train<-index[1:round(sampsize*2/3)]
sel.test<-index[(round(sampsize*2/3)+1):sampsize]
eff_train_x<-eff_model[sel.train,]
eff_test_x<-eff_model[sel.test,]

eff_train_y<-as.numeric(allf$mt_splicing_efficiency3[sel.train])
eff_test_y<-as.numeric(allf$mt_splicing_efficiency3[sel.test])
dim(eff_test_x) # 1279  94
dim(eff_train_x)# 2558  94


# create a lasso model
cv_lasso <- cv.glmnet(x = eff_train_x, y = eff_train_y, alpha = 1, family = "gaussian")

# Compute R^2 from true and predicted values
eval_results <- function(true, predicted, df) {
  SSE <- sum((predicted - true)^2)
  SST <- sum((true - mean(true))^2)
  R_square <- 1 - SSE / SST
  RMSE = sqrt(SSE/nrow(df))
  
  # Model performance metrics
  data.frame(
    RMSE = RMSE,
    Rsquare = R_square
  )
  
}

# Prediction and evaluation on train data
predictions_train <- predict(cv_lasso, s = cv_lasso$lambda.min, newx = eff_train_x)
eval_results(eff_train_y, predictions_train, eff_train_x)

# Prediction and evaluation on test data
predictions_test <- predict(cv_lasso, s = cv_lasso$lambda.min, newx = eff_test_x)
eval_results(eff_test_y, predictions_test, eff_test_x)
##
#RMSE   Rsquare
# 0.7442244 0.2190684

### obtain coefficient
coef(cv_lasso, s="lambda.min")
write.csv(as.matrix(coef(cv_lasso, s="lambda.min")),"non_ag_no_wt_eff_model_coef_0812_nostd2.csv")

