pacman::p_load(data.table,dplyr,tidyr,purrr,lightgbm,xgboost,ranger,fastshap,caret,kernlab,pROC,PRROC,writexl,ggplot2,patchwork,catboost)
setwd(dpath <- "D:/analysis/stroke/")
infile <- file.path(dpath,"dat_ml.csv")
Y <- "cvd_stroke_i"; target <- paste0(Y,".Yt2e")

norm01 <- \(x){s <- sum(x,na.rm=T); if(!is.finite(s)||s<=0) rep(0,length(x)) else x/s}
rank_desc <- \(x) rank(-x,ties.method="average")
logloss_vec <- \(y,p,eps=1e-15){p <- pmin(pmax(p,eps),1-eps); -mean(y*log(p)+(1-y)*log(1-p))}
get_select_k <- function(loss_mean, cutoff = .002) { # 把阈值从 0.005 收紧到 0.002
	n <- length(loss_mean)
	if (n <= 2) return(n)
	d <- c(NA, head(loss_mean, -1) - tail(loss_mean, -1))
	for (i in 2:(n - 1)) if (isTRUE(d[i] < cutoff && d[i + 1] < cutoff)) return(i - 1)
	n
}
make_folds <- function(y,k=5,seed=12345){set.seed(seed); i1 <- sample(which(y==1)); i0 <- sample(which(y==0)); fs <- vector("list",k); for(j in 1:k) fs[[j]] <- integer(0); for(j in seq_along(i1)) fs[[((j-1)%%k)+1]] <- c(fs[[((j-1)%%k)+1]],i1[j]); for(j in seq_along(i0)) fs[[((j-1)%%k)+1]] <- c(fs[[((j-1)%%k)+1]],i0[j]); fs}
set.seed(12345)

mydf <- fread(infile)
feat_lst <- grep("_res$", names(mydf), value=TRUE)
mydf <- mydf[, c("eid","region_code",target,feat_lst), with=FALSE][complete.cases(mydf[, c("eid","region_code",target,feat_lst), with=FALSE])]
fold_id_lst <- sort(unique(mydf$region_code)); stopifnot(length(fold_id_lst)==10)
mkxy <- \(idx, vars) list(X=as.matrix(mydf[idx, vars, with=FALSE]), y=mydf[[target]][idx])
threads <- max(1L, parallel::detectCores()-1L)

fit_lgb <- \(X,y) lightgbm::lgb.train(
	params=list(objective="binary",metric="auc",is_unbalance=TRUE,verbosity=-1L,seed=12345L,learning_rate=.01,num_leaves=10L,max_depth=10L,bagging_fraction=.7,feature_fraction=.7,bagging_freq=1L,num_threads=threads),
	data=lightgbm::lgb.Dataset(as.matrix(X),label=y), nrounds=300L)
pred_lgb <- \(m,X) predict(m,as.matrix(X))
fit_xgb <- \(X,y) xgb.train(
	params=list(objective="binary:logistic",eval_metric="logloss",eta=.01,max_depth=10L,subsample=.7,colsample_bytree=.7,seed=12345L,nthread=threads),
	data=xgb.DMatrix(as.matrix(X),label=y), nrounds=300L, verbose=0)
pred_xgb <- \(m,X) predict(m,as.matrix(X))
fit_rf  <- \(X,y){d <- as.data.frame(X); d$.y <- factor(ifelse(y==1,"case","control"), levels=c("control","case")); ranger::ranger(.y~., d, probability=TRUE, num.trees=300, importance="impurity", num.threads=threads, seed=12345)}
pred_rf <- \(m,X) predict(m,data=as.data.frame(X))$predictions[,"case"]
fit_glm  <- \(X,y) glm(y~., data=data.frame(y=y,as.data.frame(X)), family=binomial())
pred_glm <- \(m,X) as.numeric(predict(m,newdata=as.data.frame(X),type="response"))
fit_svm  <- \(X,y) kernlab::ksvm(x=as.matrix(X), y=factor(ifelse(y==1,"case","control"), levels=c("control","case")), kernel="rbfdot", prob.model=TRUE, scaled=TRUE)
pred_svm <- \(m,X) as.numeric(predict(m,as.matrix(X),type="probabilities")[,"case"])
fit_cat  <- \(X,y) catboost.train(catboost.load_pool(as.data.frame(X),label=y), NULL, params=list(loss_function="Logloss",eval_metric="AUC",iterations=300,depth=10,learning_rate=.01,subsample=.7,random_seed=12345,thread_count=threads,logging_level="Silent"))
pred_cat <- \(m,X) as.numeric(catboost.predict(m,catboost.load_pool(as.data.frame(X)),prediction_type="Probability"))

#-----------------------------#
# 3) feature selection
#-----------------------------#
## 3.1 LightGBM importance (Python原样主线：Gain + Split + SHAP)
lgb_gain <- lgb_cover <- lgb_shap <- setNames(rep(0,length(feat_lst)), feat_lst)
for(fd in fold_id_lst){
	tr <- mydf$region_code!=fd; te <- !tr; a1 <- mkxy(tr,feat_lst); a2 <- mkxy(te,feat_lst); m <- fit_lgb(a1$X,a1$y)
	imp <- lightgbm::lgb.importance(model=m); g0 <- c0 <- setNames(rep(0,length(feat_lst)), feat_lst)
	if(nrow(imp)){g0[imp$Feature] <- imp$Gain; c0[imp$Feature] <- imp$Frequency}
	sh <- predict(m,a2$X,type="contrib"); if(is.vector(sh)) sh <- matrix(sh,ncol=length(feat_lst)+1); s0 <- setNames(colMeans(abs(sh[,seq_along(feat_lst),drop=FALSE])), feat_lst)
	lgb_gain <- lgb_gain + norm01(g0); lgb_cover <- lgb_cover + norm01(c0); lgb_shap <- lgb_shap + norm01(s0)
}
FoodImportance <- data.frame(Food=feat_lst, ShapValues_cv=lgb_shap/10, TotalGain_cv=lgb_gain/10, TotalCover_cv=lgb_cover/10)
FoodImportance$Ensemble <- rowMeans(FoodImportance[,c("ShapValues_cv","TotalGain_cv","TotalCover_cv")]); FoodImportance <- FoodImportance[order(-FoodImportance$TotalGain_cv), ]
fwrite(FoodImportance, file.path(dpath,"FoodImportance.csv"))
saveRDS(list(mydf=mydf, feat_lst=feat_lst, fold_id_lst=fold_id_lst, target=target, threads=threads), file.path(dpath, "ml_state_after_31.rds"))


st <- readRDS(file.path(dpath, "ml_state_after_31.rds"))
mydf <- st$mydf; feat_lst <- st$feat_lst; fold_id_lst <- st$fold_id_lst; target <- st$target; threads <- st$threads
FoodImportance <- fread(file.path(dpath, "FoodImportance.csv")) %>% as.data.frame()



## 3.2 LGB多方法SFS：自动定主方法和k
lgb_m <- data.frame(Method=c("LGB_SHAP","LGB_Gain","LGB_Cover","LGB_Ensemble"), imp_col=c("ShapValues_cv","TotalGain_cv","TotalCover_cv","Ensemble"), file=c("AccLOSS_ShapValues.csv","AccLOSS_TotalGain.csv","AccLOSS_TotalCover.csv","AccLOSS_Ensemble.csv"))
run_lgb_sfs <- \(imp_col,meth,file){ord <- FoodImportance$Food[order(-FoodImportance[[imp_col]])]; yfull <- unlist(lapply(fold_id_lst,\(fd) mydf[[target]][mydf$region_code==fd])); tmp <- c(); out <- list(); for(f in ord){tmp <- c(tmp,f); lc <- c(); pfull <- c(); for(fd in fold_id_lst){tr <- mydf$region_code!=fd; a1 <- mkxy(tr,tmp); a2 <- mkxy(!tr,tmp); m <- fit_lgb(a1$X,a1$y); p2 <- pred_lgb(m,a2$X); lc <- c(lc,logloss_vec(a2$y,p2)); pfull <- c(pfull,p2)}; out[[length(out)+1]] <- data.frame(Food=f,loss_cv_mean=mean(lc),loss_cv_sd=sd(lc),loss_all=logloss_vec(yfull,pfull),check.names=FALSE)}; d <- bind_rows(out) %>% mutate(Method=meth,Importance=FoodImportance[[imp_col]][match(Food,FoodImportance$Food)],Delta_loss=lag(loss_cv_mean)-loss_cv_mean,Select=+(row_number()<=get_select_k(loss_cv_mean,cutoff=.002))); fwrite(d %>% mutate(across(c(loss_cv_mean,loss_cv_sd,loss_all,Delta_loss,Importance), ~round(.x,5))), file.path(dpath,file)); d}
tab_lgb_all <- bind_rows(lapply(seq_len(nrow(lgb_m)), \(i) run_lgb_sfs(lgb_m$imp_col[i],lgb_m$Method[i],lgb_m$file[i])))
lgb_choice <- tab_lgb_all %>% group_by(Method) %>% summarise(k=max(c(1,which(Select==1))), loss_at_k=loss_cv_mean[k], best_loss=min(loss_cv_mean), .groups="drop") %>% arrange(loss_at_k,k,best_loss)
lgb_main_method <- lgb_choice$Method[1]; k_main <- lgb_choice$k[1]
tab_lgb_main <- tab_lgb_all %>% filter(Method==lgb_main_method) %>% arrange(desc(Importance))
main_foods <- tab_lgb_main$Food[1:k_main]
write.csv(tab_lgb_all,file.path(dpath,"LGB_Methods_Summary.csv"),row.names=FALSE)
write.csv(lgb_choice,file.path(dpath,"LGB_MethodChoice.csv"),row.names=FALSE)
write.csv(tab_lgb_main,file.path(dpath,"LGB_MainMethod_Ranking.csv"),row.names=FALSE)
write.csv(data.frame(Food=main_foods,LGB_main_method=lgb_main_method,k_main=k_main),file.path(dpath,"LGB_MainFeatures.csv"),row.names=FALSE)

## 3.3 XGB model-based + SHAP 排名
xgb_gain <- xgb_cover <- xgb_shap <- setNames(rep(0,length(feat_lst)), feat_lst)
for(fd in fold_id_lst){tr <- mydf$region_code!=fd; a1 <- mkxy(tr,feat_lst); a2 <- mkxy(!tr,feat_lst); m <- fit_xgb(a1$X,a1$y); imp <- xgb.importance(model=m,feature_names=feat_lst); g0 <- c0 <- setNames(rep(0,length(feat_lst)),feat_lst); if(nrow(imp)){g0[imp$Feature] <- imp$Gain; c0[imp$Feature] <- imp$Cover}; sh <- predict(m,a2$X,predcontrib=TRUE); if(is.vector(sh)) sh <- matrix(sh,ncol=length(feat_lst)+1); xgb_gain <- xgb_gain + norm01(g0); xgb_cover <- xgb_cover + norm01(c0); xgb_shap <- xgb_shap + norm01(setNames(colMeans(abs(sh[,seq_along(feat_lst),drop=FALSE])),feat_lst))}
XGB_Importance <- data.frame(Food=feat_lst,XGB_Gain=xgb_gain/10,XGB_Cover=xgb_cover/10,XGB_SHAP=xgb_shap/10,XGB_Ensemble=rowMeans(cbind(xgb_gain/10,xgb_cover/10,xgb_shap/10))) %>% arrange(desc(XGB_Gain))
fwrite(XGB_Importance,file.path(dpath,"XGB_Importance.csv"))

## 3.4 RF model-based + SHAP 排名
rf_imp <- rf_shap <- setNames(rep(0,length(feat_lst)), feat_lst)
for(fd in fold_id_lst){tr <- mydf$region_code!=fd; d1 <- as.data.frame(mydf[tr,c(target,feat_lst),with=FALSE]); d1[[target]] <- factor(ifelse(d1[[target]]==1,"case","control"),levels=c("control","case")); d2 <- as.data.frame(mydf[!tr,feat_lst,with=FALSE]); m <- fit_rf(d1[,feat_lst,drop=FALSE], +(d1[[target]]=="case")); rf_imp <- rf_imp + norm01(setNames(m$variable.importance[feat_lst],feat_lst)); sh <- fastshap::explain(m,X=d1[,feat_lst,drop=FALSE],newdata=d2,pred_wrapper=\(obj,newdata) predict(obj,data=as.data.frame(newdata))$predictions[,"case"],nsim=20,adjust=TRUE); rf_shap <- rf_shap + norm01(setNames(colMeans(abs(sh)),colnames(sh)))}
RF_Importance <- data.frame(Food=feat_lst,RF_Impurity=rf_imp/10,RF_SHAP=rf_shap/10) %>% arrange(desc(RF_Impurity))
fwrite(RF_Importance,file.path(dpath,"RF_Importance.csv"))

## 3.4b CatBoost model-based + SHAP 排名
library(catboost)
cat_imp <- cat_shap <- setNames(rep(0,length(feat_lst)), feat_lst)
for(fd in fold_id_lst){tr <- mydf$region_code!=fd; a1 <- mkxy(tr,feat_lst); a2 <- mkxy(!tr,feat_lst); m <- fit_cat(a1$X,a1$y); p_tr <- catboost.load_pool(as.data.frame(a1$X),label=a1$y); p_te <- catboost.load_pool(as.data.frame(a2$X),label=a2$y); i0 <- setNames(as.numeric(catboost.get_feature_importance(model=m,pool=p_tr,type="FeatureImportance",thread_count=threads)),feat_lst); sh <- catboost.get_feature_importance(model=m,pool=p_te,type="ShapValues",thread_count=threads); cat_imp <- cat_imp + norm01(i0); cat_shap <- cat_shap + norm01(setNames(colMeans(abs(sh[,seq_along(feat_lst),drop=FALSE])),feat_lst))}
CatBoost_Importance <- data.frame(Food=feat_lst,CatBoost_Imp=cat_imp/10,CatBoost_SHAP=cat_shap/10) %>% arrange(desc(CatBoost_Imp))
fwrite(CatBoost_Importance,file.path(dpath,"CatBoost_Importance.csv"))

## 3.5 RFE-based
y01 <- \(z) if(is.factor(z) || is.character(z)) as.integer(z == "case") else as.integer(z == 1)
rfe_rank <- \(type, vars, dat_use, ycol) {
	dat_use <- as.data.frame(dat_use); v <- vars; drop <- character(0)
	tr_ctrl <- caret::trainControl(method = "cv", number = 3, classProbs = TRUE, summaryFunction = twoClassSummary, allowParallel = TRUE)
	while (length(v) > 1) {
		if (type == "lgb") {m <- fit_lgb(as.matrix(dat_use[, v, drop = FALSE]), y01(dat_use[[ycol]])); imp <- lightgbm::lgb.importance(model = m); s <- setNames(rep(0, length(v)), v); if (nrow(imp)) s[imp$Feature] <- imp$Gain}
		else if (type == "xgb") {m <- fit_xgb(as.matrix(dat_use[, v, drop = FALSE]), y01(dat_use[[ycol]])); imp <- xgb.importance(model = m, feature_names = v); s <- setNames(rep(0, length(v)), v); if (nrow(imp)) s[imp$Feature] <- imp$Gain}
		else if (type == "rf") {m <- fit_rf(dat_use[, v, drop = FALSE], y01(dat_use[[ycol]])); s <- setNames(m$variable.importance[v], v)}
		else if (type == "glm") {yy <- y01(dat_use[[ycol]]); m <- glm(yy ~ ., data = data.frame(yy = yy, dat_use[, v, drop = FALSE]), family = binomial()); s <- setNames(rep(0, length(v)), v); cc <- abs(coef(m)[-1]); s[names(cc)] <- cc}
		else {m <- caret::train(x = dat_use[, v, drop = FALSE], y = dat_use[[ycol]], method = "svmLinear", metric = "ROC", trControl = tr_ctrl, tuneLength = 1); s <- caret::varImp(m)$importance; s <- if (ncol(s) > 1) rowMeans(s) else s[, 1]}
		s[is.na(s)] <- 0; d <- names(sort(s))[1]; drop <- c(drop, d); v <- setdiff(v, d)
	}
	data.frame(Food = c(v, rev(drop)), rank = seq_along(vars))
}
dat_cls <- as.data.frame(mydf[, c(target, feat_lst), with = FALSE]); dat_cls[[target]] <- factor(ifelse(dat_cls[[target]] == 1, "case", "control"), levels = c("control", "case"))
i1 <- sample(which(dat_cls[[target]] == "case"), min(3000, sum(dat_cls[[target]] == "case"))); i0 <- sample(which(dat_cls[[target]] == "control"), min(3000, sum(dat_cls[[target]] == "control"))); dat_svm <- dat_cls[c(i1, i0), ]
RFE_Ranks <- Reduce(\(a,b) left_join(a,b,by="Food"), list(
	rfe_rank("lgb", feat_lst, mydf, target) %>% rename(RFE_LGB_rank = rank),
	rfe_rank("xgb", feat_lst, mydf, target) %>% rename(RFE_XGB_rank = rank),
	rfe_rank("rf",  feat_lst, mydf, target) %>% rename(RFE_RF_rank  = rank),
	rfe_rank("glm", feat_lst, dat_cls, target) %>% rename(RFE_GLM_rank = rank),
	rfe_rank("svm", feat_lst, dat_svm, target) %>% rename(RFE_SVM_rank = rank)
))
write.csv(RFE_Ranks, file.path(dpath, "RFE_Ranks.csv"), row.names = FALSE)


## 3.6 汇总：只看LGB主结果在其他方法里的排名是否也靠前
FoodImportance <- fread(file.path(dpath,"FoodImportance.csv")) %>% as.data.frame()
XGB_Importance <- fread(file.path(dpath,"XGB_Importance.csv")) %>% as.data.frame()
RF_Importance <- fread(file.path(dpath,"RF_Importance.csv")) %>% as.data.frame()
CatBoost_Importance <- fread(file.path(dpath,"CatBoost_Importance.csv")) %>% as.data.frame()
RFE_Ranks <- fread(file.path(dpath,"RFE_Ranks.csv")) %>% as.data.frame()
lgb_choice <- fread(file.path(dpath,"LGB_MethodChoice.csv")) %>% as.data.frame()
tab_lgb_main <- fread(file.path(dpath,"LGB_MainMethod_Ranking.csv")) %>% as.data.frame()

rk <- \(x) rank(-x, ties.method = "min")
lgb_main_method <- lgb_choice$Method[1]; k_main <- lgb_choice$k[1]
main_foods <- if(file.exists(file.path(dpath,"LGB_MainFeatures.csv"))) fread(file.path(dpath,"LGB_MainFeatures.csv"))$Food else head(tab_lgb_main$Food, k_main)
Rank_All <- FoodImportance %>% 
	transmute(Food, LGB_SHAP_rank=rk(ShapValues_cv), LGB_Gain_rank=rk(TotalGain_cv), LGB_Cover_rank=rk(TotalCover_cv), LGB_Ensemble_rank=rk(Ensemble)) %>%
	left_join(XGB_Importance %>% transmute(Food, XGB_Gain_rank=rk(XGB_Gain), XGB_Cover_rank=rk(XGB_Cover), XGB_SHAP_rank=rk(XGB_SHAP), XGB_Ensemble_rank=rk(XGB_Ensemble)), by="Food") %>%
	left_join(RF_Importance %>% transmute(Food, RF_Impurity_rank=rk(RF_Impurity), RF_SHAP_rank=rk(RF_SHAP)), by="Food") %>%
	left_join(CatBoost_Importance %>% transmute(Food, CatBoost_Imp_rank=rk(CatBoost_Imp), CatBoost_SHAP_rank=rk(CatBoost_SHAP)), by="Food") %>%
	left_join(RFE_Ranks, by="Food")
other_rank_cols <- grep("_rank$", names(Rank_All), value=TRUE); other_rank_cols <- setdiff(other_rank_cols, grep("^LGB_", other_rank_cols, value=TRUE))
Main_Stability <- Rank_All %>% filter(Food %in% main_foods) %>% mutate(LGB_main_method=lgb_main_method, k_main=k_main, LGB_main_rank=match(Food, main_foods))
Main_Stability[paste0(other_rank_cols,"_topk")] <- lapply(Main_Stability[other_rank_cols], \(z) +(z <= k_main))
Main_Stability$n_topk_other <- rowSums(Main_Stability[paste0(other_rank_cols,"_topk")], na.rm=TRUE)
Main_Stability <- Main_Stability %>% arrange(LGB_main_rank)

write.csv(Rank_All, file.path(dpath,"ModelRanks_AllFoods.csv"), row.names=FALSE)
write.csv(Main_Stability, file.path(dpath,"MainFeatures_Stability.csv"), row.names=FALSE)
write_xlsx(list(
	LGB_MethodChoice=lgb_choice,
	LGB_MainMethod_Ranking=tab_lgb_main,
	MainFeatures_LGB=data.frame(Food=main_foods),
	XGB_Importance=XGB_Importance,
	RF_Importance=RF_Importance,
	CatBoost_Importance=CatBoost_Importance,
	RFE_Ranks=RFE_Ranks,
	ModelRanks_AllFoods=Rank_All,
	MainFeatures_Stability=Main_Stability
), file.path(dpath,"ML_Stability_Summary.xlsx"))

# TableS6
lgb_mod  <- tab_lgb_main %>% arrange(desc(Importance))
xgb_mod  <- XGB_Importance %>% arrange(desc(XGB_Gain))
rf_mod   <- RF_Importance %>% arrange(desc(RF_Impurity))
cat_mod  <- CatBoost_Importance %>% arrange(desc(CatBoost_Imp))
lgb_shap <- FoodImportance %>% arrange(desc(ShapValues_cv))
xgb_shap <- XGB_Importance %>% arrange(desc(XGB_SHAP))
rf_shap  <- RF_Importance %>% arrange(desc(RF_SHAP))
cat_shap <- CatBoost_Importance %>% arrange(desc(CatBoost_SHAP))
rfe_lgb  <- RFE_Ranks %>% arrange(RFE_LGB_rank) %>% pull(Food)
rfe_xgb  <- RFE_Ranks %>% arrange(RFE_XGB_rank) %>% pull(Food)
rfe_rf   <- RFE_Ranks %>% arrange(RFE_RF_rank) %>% pull(Food)
rfe_glm  <- RFE_Ranks %>% arrange(RFE_GLM_rank) %>% pull(Food)
rfe_svm  <- RFE_Ranks %>% arrange(RFE_SVM_rank) %>% pull(Food)
nmax <- max(nrow(lgb_mod), nrow(xgb_mod), nrow(rf_mod), nrow(cat_mod), nrow(lgb_shap), nrow(xgb_shap), nrow(rf_shap), nrow(cat_shap), length(rfe_lgb), length(rfe_xgb), length(rfe_rf), length(rfe_glm), length(rfe_svm))
pad <- \(x,n) c(x, rep(NA, n - length(x)))
TableS6 <- data.frame(
	Ranking = 1:nmax,
	LGB_Food = pad(lgb_mod$Food, nmax), LGB_Importance = pad(round(lgb_mod$Importance, 8), nmax), LGB_Loss = pad(round(lgb_mod$loss_cv_mean, 5), nmax), LGB_Loss_SD = pad(round(lgb_mod$loss_cv_sd, 5), nmax), LGB_Delta_loss = pad(round(lgb_mod$Delta_loss, 5), nmax), LGB_Select = pad(lgb_mod$Select, nmax),
	XGB_ModelBased = pad(xgb_mod$Food, nmax), RF_ModelBased = pad(rf_mod$Food, nmax), CatBoost_ModelBased = pad(cat_mod$Food, nmax),
	LGB_SHAP = pad(lgb_shap$Food, nmax), XGB_SHAP = pad(xgb_shap$Food, nmax), RF_SHAP = pad(rf_shap$Food, nmax), CatBoost_SHAP = pad(cat_shap$Food, nmax),
	RFE_LGB = pad(rfe_lgb, nmax), RFE_XGB = pad(rfe_xgb, nmax), RFE_RF = pad(rfe_rf, nmax), RFE_GLM = pad(rfe_glm, nmax), RFE_SVM = pad(rfe_svm, nmax)
)
write.csv(TableS6, file.path(dpath,"TableS6.csv"), row.names = FALSE)
write_xlsx(list(
	LGB_MethodChoice = lgb_choice,
	LGB_MainMethod_Ranking = tab_lgb_main,
	FoodImportance = FoodImportance,
	XGB_Importance = XGB_Importance,
	RF_Importance = RF_Importance,
	CatBoost_Importance = CatBoost_Importance,
	RFE_Ranks = RFE_Ranks,
	TableS6 = TableS6
), file.path(dpath,"TableS6.xlsx"))


#-----------------------------#
# 4) prediction analysis
#-----------------------------#
rd <- \(f) fread(file.path(dpath, f)) %>% as.data.frame()
topk <- \(x,k) unique(na.omit(x))[1:min(k,length(unique(na.omit(x))))]
as_num <- \(x) as.numeric(unlist(x, use.names = FALSE))
fmt_ci <- \(x) sprintf("%.3f [%.3f - %.3f]", mean(x), quantile(x,.025), quantile(x,.975))
append_csv <- \(x,f) fwrite(x, f, append=file.exists(f), col.names=!file.exists(f))
calc_metrics <- \(y,p){
	y <- as_num(y); p <- as_num(p); stopifnot(length(y)==length(p))
	ro <- pROC::roc(y, p, quiet=TRUE, direction="<")
	cd <- pROC::coords(ro, x="best", best.method="youden",
		ret=c("threshold","sensitivity","specificity"), transpose=FALSE)
	th <- as_num(cd[["threshold"]])[1]; se <- as_num(cd[["sensitivity"]])[1]; sp <- as_num(cd[["specificity"]])[1]
	yh <- +(p >= th)
	c(AUC=as.numeric(pROC::auc(ro)), Sensitivity=se, Specificity=sp,
	  Accuracy=mean(yh==y), Youden_Index=se+sp-1,
	  APR=PRROC::pr.curve(scores.class0=p[y==1], scores.class1=p[y==0], curve=FALSE)$auc.integral)
}
cv_iso_pred <- \(X1,y1,X2,fit_fun,pred_fun,seed=12345){
	y1 <- as_num(y1); fs <- make_folds(y1,5,seed); oof <- rep(NA_real_, length(y1))
	for(i in seq_along(fs)){va <- fs[[i]]; tr <- setdiff(seq_along(y1), va); oof[va] <- as_num(pred_fun(fit_fun(X1[tr,,drop=FALSE], y1[tr]), X1[va,,drop=FALSE]))}
	iso <- isoreg(oof, y1); f <- approxfun(iso$x[iso$ord], iso$yf, ties=mean, rule=2)
	p <- pmin(pmax(as_num(f(as_num(pred_fun(fit_fun(X1,y1), X2)))), 0), 1); stopifnot(length(p)==nrow(X2)); p
}

tab_lgb_all <- rd("LGB_Methods_Summary.csv")
lgb_choice <- rd("LGB_MethodChoice.csv")
xgb <- rd("XGB_Importance.csv")
rf <- rd("RF_Importance.csv")
rfe <- rd("RFE_Ranks.csv")
catb <- rd("CatBoost_Importance.csv")
k_sel <- lgb_choice$k[1]; lgb_main_method <- lgb_choice$Method[1]
get_lgb <- \(m) tab_lgb_all %>% filter(Method==m) %>% arrange(desc(Importance)) %>% pull(Food) %>% topk(k_sel)
feature_sets <- list(# svm太费时间了
	All15=feat_lst,
	LGB_ModelBased=get_lgb(lgb_main_method), XGB_ModelBased=topk(xgb$Food[order(-xgb$XGB_Gain)], k_sel),
	RF_ModelBased=topk(rf$Food[order(-rf$RF_Impurity)], k_sel), CatBoost_ModelBased=topk(catb$Food[order(-catb$CatBoost_Imp)], k_sel),
	LGB_SHAP=get_lgb("LGB_SHAP"), XGB_SHAP=topk(xgb$Food[order(-xgb$XGB_SHAP)], k_sel),
	RF_SHAP=topk(rf$Food[order(-rf$RF_SHAP)], k_sel), CatBoost_SHAP=topk(catb$Food[order(-catb$CatBoost_SHAP)], k_sel),
	RFE_LGB=topk(rfe$Food[order(rfe$RFE_LGB_rank)], k_sel), RFE_XGB=topk(rfe$Food[order(rfe$RFE_XGB_rank)], k_sel),
	RFE_RF=topk(rfe$Food[order(rfe$RFE_RF_rank)], k_sel), RFE_GLM=topk(rfe$Food[order(rfe$RFE_GLM_rank)], k_sel),
	LGB_Plus_WholeGrains=unique(c(get_lgb(lgb_main_method),"Whole_grains_res")),
	LGB_Plus_WG_GLV=unique(c(get_lgb(lgb_main_method),"Whole_grains_res","Green_leafy_vegetables_res")),
	LGB_Plus_WG_GLV_PM=unique(c(get_lgb(lgb_main_method),"Whole_grains_res","Green_leafy_vegetables_res","Processed_meats_res"))
)
feature_sets <- lapply(feature_sets, \(x) intersect(unique(x), feat_lst))
jobs <- bind_rows(
	data.frame(Block=paste("All",length(feat_lst),"food groups"), Method=c("LightGBM","XGBoost","Random forest","CatBoost"), Set="All15", Model=c("lgb","xgb","rf","cat")),
	data.frame(Block=paste0("Model-based top-",k_sel," selected food groups"),
		Method=c("LightGBM","XGBoost","Random forest","CatBoost","LightGBM + Whole grains","LightGBM + Whole grains + Green leafy vegetables","LightGBM + Whole grains + Green leafy vegetables + Processed meats"),
		Set=c("LGB_ModelBased","XGB_ModelBased","RF_ModelBased","CatBoost_ModelBased","LGB_Plus_WholeGrains","LGB_Plus_WG_GLV","LGB_Plus_WG_GLV_PM"),
		Model=c("lgb","xgb","rf","cat","lgb","lgb","lgb")),
	data.frame(Block=paste0("SHAP value-based top-",k_sel," selected food groups"),
		Method=c("LightGBM","XGBoost","Random forest","CatBoost"),
		Set=c("LGB_SHAP","XGB_SHAP","RF_SHAP","CatBoost_SHAP"),
		Model=c("lgb","xgb","rf","cat")),
	data.frame(Block=paste0("RFE-based top-",k_sel," selected food groups"),
		Method=c("LightGBM","XGBoost","Random forest","Logistic regression"),
		Set=c("RFE_LGB","RFE_XGB","RFE_RF","RFE_GLM"),
		Model=c("lgb","xgb","rf","glm"))
)
fit_map <- list(lgb=fit_lgb,xgb=fit_xgb,rf=fit_rf,glm=fit_glm,cat=fit_cat)
pred_map <- list(lgb=pred_lgb,xgb=pred_xgb,rf=pred_rf,glm=pred_glm,cat=pred_cat)
write.csv(bind_rows(lapply(names(feature_sets), \(nm) data.frame(Set=nm, Food=feature_sets[[nm]]))), file.path(dpath, paste0("Prediction_FeatureSets_top",k_sel,".csv")), row.names=FALSE)
pred_file <- file.path(dpath, paste0("PredProbs_TableS7_top",k_sel,"_extended.csv"))
metric_file <- file.path(dpath, paste0("MetricsByFold_TableS7_top",k_sel,"_extended.csv"))
if(file.exists(pred_file)) file.remove(pred_file); if(file.exists(metric_file)) file.remove(metric_file)
for(i in seq_len(nrow(jobs))){
	v <- feature_sets[[jobs$Set[i]]]; ff <- fit_map[[jobs$Model[i]]]; pf <- pred_map[[jobs$Model[i]]]
	for(fd in fold_id_lst){
		tr <- mydf$region_code!=fd; te <- !tr
		X1 <- as.data.frame(mydf[tr,..v]); X2 <- as.data.frame(mydf[te,..v]); y1 <- mydf[[target]][tr]; y2 <- mydf[[target]][te]
		p2 <- tryCatch(cv_iso_pred(X1,y1,X2,ff,pf,12345), error=\(e){cat("\nFAILED:",jobs$Method[i],"|",jobs$Set[i],"| fold",fd,"\n"); stop(e)})
		append_csv(data.frame(eid=mydf$eid[te], region_code=mydf$region_code[te], y_true=y2, Block=jobs$Block[i], Method=jobs$Method[i], Set=jobs$Set[i], y_pred=as_num(p2)), pred_file)
		append_csv(data.frame(Block=jobs$Block[i], Method=jobs$Method[i], Fold=fd, t(calc_metrics(y2,p2)), check.names=FALSE), metric_file)
		cat("done:", jobs$Method[i], "|", jobs$Set[i], "| fold", fd, "\n")
	}
}
PredProbs <- rd(basename(pred_file)); metrics_long <- rd(basename(metric_file))
TableS7_tidy <- metrics_long %>% group_by(Block,Method) %>% summarise(AUC=fmt_ci(AUC),Sensitivity=fmt_ci(Sensitivity),Specificity=fmt_ci(Specificity),Accuracy=fmt_ci(Accuracy),Youden_Index=fmt_ci(Youden_Index),APR=fmt_ci(APR),.groups="drop")
write.csv(TableS7_tidy, file.path(dpath,"TableS7_tidy.csv"), row.names=FALSE)

metric_order <- c("AUC","Sensitivity","Specificity","Accuracy","Youden_Index","APR")
tabs7 <- lapply(split(TableS7_tidy, TableS7_tidy$Block), \(d){out <- data.frame(Metric=metric_order); for(m in d$Method) out[[m]] <- unlist(d[d$Method==m, metric_order], use.names=FALSE); out})
mk_block <- \(nm){d <- tabs7[[nm]]; h <- as.data.frame(setNames(as.list(rep(NA_character_, ncol(d))), names(d))); h$Metric <- nm; b <- h; b$Metric <- NA_character_; dplyr::bind_rows(h, d, b)}
TableS7_csv <- dplyr::bind_rows(lapply(names(tabs7), mk_block))
write.csv(TableS7_csv, file.path(dpath,"TableS7.csv"), row.names=FALSE)


#-----------------------------#
# 4b) incremental prediction: age + sex only
#-----------------------------#
covs.base2 <- c("age","sex.b")
grp <- as.data.table(readRDS(file.path(dpath,"dat.grp.rds"))); grp[, eid := as.character(eid)]
mydf2 <- copy(mydf); mydf2[, eid := as.character(eid)]
idx <- match(mydf2$eid, grp$eid); mydf2[, (covs.base2) := grp[idx, ..covs.base2]]
stopifnot(nrow(mydf2) == nrow(mydf), !anyNA(mydf2[, ..covs.base2]), identical(sort(unique(mydf2$region_code)), fold_id_lst))

diet7 <- unique(fread(file.path(dpath,"LGB_MainFeatures.csv"))$Food)
diet8 <- unique(c(diet7, "Whole_grains_res")) 
mkx <- \(v){x <- model.matrix(~ . - 1, data = as.data.frame(mydf2[, ..v])); storage.mode(x) <- "double"; x}
Xlst <- list(
	Diet7 = mkx(diet7),
	Diet8 = mkx(diet8),
	Baseline = mkx(covs.base2),
	Baseline_Diet7 = mkx(c(covs.base2, diet7)),
	Baseline_Diet8 = mkx(c(covs.base2, diet8))
)
y2 <- mydf2[[target]]
ev <- \(X, lab){
	z <- lapply(fold_id_lst, \(fd){
		tr <- mydf2$region_code != fd; te <- !tr
		p <- cv_iso_pred(X[tr,,drop=FALSE], y2[tr], X[te,,drop=FALSE], fit_lgb, pred_lgb, 12345)
		list(
			pred = data.frame(eid = mydf2$eid[te], region_code = mydf2$region_code[te], Fold = fd, Block = lab, y_true = y2[te], y_pred = p),
			met  = data.frame(Block = lab, Fold = fd, t(calc_metrics(y2[te], p)), check.names = FALSE)
		)
	})
	list(pred = bind_rows(lapply(z, `[[`, "pred")), met = bind_rows(lapply(z, `[[`, "met")))
}
res_inc <- lapply(setNames(names(Xlst), names(Xlst)), \(nm){cat("running:", nm, "\n"); ev(Xlst[[nm]], nm)})
pred_inc <- bind_rows(lapply(res_inc, `[[`, "pred"))
met_inc  <- bind_rows(lapply(res_inc, `[[`, "met"))
sum_inc  <- met_inc %>% group_by(Block) %>% summarise(
	AUC = fmt_ci(AUC), Sensitivity = fmt_ci(Sensitivity), Specificity = fmt_ci(Specificity),
	Accuracy = fmt_ci(Accuracy), Youden_Index = fmt_ci(Youden_Index), APR = fmt_ci(APR), .groups = "drop"
)
mk_delta <- \(b){
	x <- met_inc %>% filter(Block %in% c("Baseline", b)) %>%
		select(Block, Fold, AUC, APR, Youden_Index) %>%
		tidyr::pivot_wider(names_from = Block, values_from = c(AUC, APR, Youden_Index))
	data.frame(
		compare = paste0(b, " - Baseline"),
		n_fold = nrow(x),
		mean_dAUC = mean(x[[paste0("AUC_", b)]] - x$AUC_Baseline),
		median_dAUC = median(x[[paste0("AUC_", b)]] - x$AUC_Baseline),
		n_pos = sum((x[[paste0("AUC_", b)]] - x$AUC_Baseline) > 0),
		p_wilcox = wilcox.test(x[[paste0("AUC_", b)]] - x$AUC_Baseline, mu = 0, alternative = "greater", exact = TRUE)$p.value,
		p_sign = binom.test(sum((x[[paste0("AUC_", b)]] - x$AUC_Baseline) > 0), nrow(x), p = 0.5, alternative = "greater")$p.value,
		mean_dAPR = mean(x[[paste0("APR_", b)]] - x$APR_Baseline),
		median_dAPR = median(x[[paste0("APR_", b)]] - x$APR_Baseline),
		mean_dYouden = mean(x[[paste0("Youden_Index_", b)]] - x$Youden_Index_Baseline),
		median_dYouden = median(x[[paste0("Youden_Index_", b)]] - x$Youden_Index_Baseline)
	)
}
delta_inc <- bind_rows(mk_delta("Baseline_Diet7"), mk_delta("Baseline_Diet8"))
sum_inc; delta_inc

write.csv(pred_inc,  file.path(dpath,"Prediction_incremental_age_sex_PredProbs.csv"), row.names = FALSE)
write.csv(met_inc,   file.path(dpath,"Prediction_incremental_age_sex_ByFold.csv"), row.names = FALSE)
write.csv(sum_inc,   file.path(dpath,"Prediction_incremental_age_sex_Summary.csv"), row.names = FALSE)
write.csv(delta_inc, file.path(dpath,"Prediction_incremental_age_sex_delta.csv"), row.names = FALSE)
write_xlsx(list(Pred = pred_inc, ByFold = met_inc, Summary = sum_inc, Delta = delta_inc), file.path(dpath,"Prediction_incremental_age_sex.xlsx"))


#==============================#
# 5) Main Plots Fig3ab
#==============================#
pacman::p_load(data.table, dplyr, tidyr, ggplot2, patchwork, pROC, forcats, scales)
rd <- \(f) as.data.frame(fread(file.path(dpath, f)))
d <- sapply(c("LGB_Methods_Summary.csv","LGB_MethodChoice.csv","LGB_MainMethod_Ranking.csv",
              "Prediction_incremental_age_sex_ByFold.csv","Prediction_incremental_age_sex_Summary.csv",
              "Prediction_incremental_age_sex_PredProbs.csv","Prediction_incremental_age_sex_delta.csv"), rd, simplify=F)
names(d) <- c("tab_all","choice","tab_main","met_plot","sum_plot","pred_plot","delta_sum")
lgb <- d$choice$Method[1]; k <- d$choice$k[1]
m_blks <- c("Baseline", "Baseline_Diet7") # "Diet7"
m_labs <- c(Diet7="DOSE", Baseline="Baseline", Baseline_Diet7="Baseline + DOSE")
thm <- theme_bw(base_size=10) + theme(axis.title=element_text(size=11, face="bold"), axis.text=element_text(size=9), legend.title=element_blank(), panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(), plot.margin=margin(10,10,10,10))

# Fig3a
xgb_imp <- rd("XGB_Importance.csv")
rf_imp <- rd("RF_Importance.csv")
rfe_rk <- rd("RFE_Ranks.csv")
tk <- \(x, k) {x <- unique(na.omit(x)); x[1:min(k, length(x))]}; food_lab <- \(x) { x <- gsub("Apples Pears", "Apples/Pears", trimws(gsub("_res$|_", " ", x))); gsub("High fat diary", "High-fat dairy", x) }
fs <- list(LGB_MB = tk((d$tab_all %>% filter(Method == lgb) %>% arrange(desc(Importance)) %>% pull(Food)), k), XGB_MB = tk(xgb_imp$Food[order(-xgb_imp$XGB_Gain)], k), RF_MB = tk(rf_imp$Food[order(-rf_imp$RF_Impurity)], k), LGB_SHAP = tk((d$tab_all %>% filter(Method == "LGB_SHAP") %>% arrange(desc(Importance)) %>% pull(Food)), k), XGB_SHAP = tk(xgb_imp$Food[order(-xgb_imp$XGB_SHAP)], k), RF_SHAP = tk(rf_imp$Food[order(-rf_imp$RF_SHAP)], k), RFE_LGB = tk(rfe_rk$Food[order(rfe_rk$RFE_LGB_rank)], k), RFE_XGB = tk(rfe_rk$Food[order(rfe_rk$RFE_XGB_rank)], k), RFE_RF = tk(rfe_rk$Food[order(rfe_rk$RFE_RF_rank)], k), RFE_GLM = tk(rfe_rk$Food[order(rfe_rk$RFE_GLM_rank)], k), RFE_SVM = tk(rfe_rk$Food[order(rfe_rk$RFE_SVM_rank)], k))
info <- data.frame(Method = c("LGB_MB","XGB_MB","RF_MB", "LGB_SHAP","XGB_SHAP","RF_SHAP", "RFE_LGB","RFE_XGB","RFE_RF","RFE_GLM","RFE_SVM"), lab = c("LightGBM","XGBoost","Random forest", "LightGBM","XGBoost","Random forest", "LightGBM","XGBoost","Random forest","Logistic regression","Support vector machine"), grp = c(rep("Model-based approach", 3), rep("SHAP value-based approach", 3), rep("RFE-based approach", 5)), y = c(12,11,10, 7.5,6.5,5.5, 3,2,1,0,-1))
grp_header_df <- info %>% group_by(grp) %>% summarise(y = max(y) + 0.8, .groups = "drop")
dfAB <- d$tab_main %>% slice(1:15) %>% mutate(F_lab = factor(food_lab(Food), levels = food_lab(Food)), food_idx = row_number(), topk = food_idx <= k, loss_lower = loss_cv_mean - loss_cv_sd, loss_upper = loss_cv_mean + loss_cv_sd)
track_df <- tidyr::expand_grid(Method = info$Method, Food = dfAB$Food) %>% left_join(bind_rows(lapply(names(fs), \(nm) data.frame(Method = nm, Food = fs[[nm]], sel = 1L))), by = c("Method", "Food")) %>% left_join(info, by = "Method") %>% mutate(sel = tidyr::replace_na(sel, 0L), F_lab = factor(food_lab(Food), levels = levels(dfAB$F_lab)), topk = match(Food, dfAB$Food) <= k)
C <- list(hi="#2F5D8A", lo="#D9E1EA", line="#3B668F", rib="#B8C9D9", cut="#7B8FA6", ora="#D55E00", dot="#4A7EBB", bg=scales::alpha("#EAF1F7", 0.35)) # 降低 alpha
col_frame <- "grey35"; lw_frame <- 0.55; n_food <- nlevels(dfAB$F_lab); imp_max <- max(dfAB$Importance, na.rm = TRUE) * 1.05; a <- imp_max / (max(dfAB$loss_upper) - min(dfAB$loss_lower)); b <- -a * min(dfAB$loss_lower); ptAB <- dfAB[dfAB$food_idx == k, ]; left_margin_size <- 130 
pAB_top <- ggplot(dfAB, aes(F_lab)) + geom_rect(aes(xmin = 0.5, xmax = k + 0.5, ymin = -Inf, ymax = Inf), fill = C$bg, inherit.aes = FALSE) + geom_col(aes(y = Importance, fill = topk), width = .64) + geom_ribbon(aes(x = food_idx, ymin = a * loss_lower + b, ymax = a * loss_upper + b, group = 1), fill = C$rib, alpha = .18, inherit.aes = FALSE) + geom_line(aes(x = food_idx, y = a * loss_cv_mean + b, group = 1), colour = C$line, linewidth = 1, inherit.aes = FALSE) + geom_point(aes(x = food_idx, y = a * loss_cv_mean + b, color = topk), size = 2.2, inherit.aes = FALSE) + geom_vline(xintercept = k + .5, linetype = "22", colour = C$cut, linewidth = .7) + geom_vline(xintercept = c(0.5, n_food + 0.5), colour = col_frame, linewidth = lw_frame) + annotate("text", x = k + .58, y = a * ptAB$loss_cv_mean + b + imp_max * .02, label = paste0("n = ", k), colour = C$ora, fontface = "bold", size = 3.8, hjust = 0) + scale_fill_manual(values = c(`TRUE` = C$hi, `FALSE` = C$lo), guide = "none") + scale_color_manual(values = c(`TRUE` = C$ora, `FALSE` = C$line), guide = "none") + scale_y_continuous(name = "Food importance", expand = expansion(c(0, .05)), sec.axis = sec_axis(~(. - b) / a, name = "CV log loss")) + thm + theme(panel.border = element_blank(), axis.line.x.top = element_line(colour = col_frame, linewidth = lw_frame), axis.line.x.bottom = element_blank(), axis.line.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_text(size = 12, face = "bold"), plot.margin = margin(6, 10, -1, left_margin_size))
count_df <- track_df %>% group_by(Food, F_lab) %>% summarise(SelectionCount = sum(sel), .groups = "drop") %>% mutate(y = -2.5)
line_df <- track_df %>% filter(sel == 1) %>% group_by(F_lab) %>% summarise(y_min = min(y), y_max = max(y), count = n()) %>% filter(count > 1)
pAB_track <- ggplot(track_df, aes(F_lab, y)) + geom_rect(aes(xmin = 0.5, xmax = k + 0.5, ymin = -Inf, ymax = Inf), fill = C$bg, inherit.aes = FALSE) + geom_point(colour = "grey92", size = 1.5) + 
	geom_segment(data = line_df, aes(x = F_lab, xend = F_lab, y = y_min, yend = y_max), colour = C$dot, linewidth = 0.7) +
	geom_point(data = filter(track_df, sel == 1), shape = 21, size = 3, fill = C$dot, colour = "white", stroke = 0.4) +
	geom_text(data = count_df, aes(x = F_lab, y = -2.5, label = SelectionCount), size = 3.5, fontface = "bold", colour = C$hi) +
	annotate("text", x = 0.45, y = -2.5, label = "Selection count", hjust = 1, fontface = "bold", size = 3.2, colour = "grey35") +
	geom_hline(yintercept = c(8.75, 4.25), colour = "grey88", linewidth = .4) +
	geom_vline(xintercept = c(0.5, n_food + 0.5), colour = col_frame, linewidth = lw_frame) +
	geom_text(data = grp_header_df, aes(x = 0.45, y = y, label = grp), hjust = 1, fontface = "bold", size = 3.2, colour = "grey35", inherit.aes = FALSE) +
	scale_y_continuous(breaks = c(-2.5, info$y), labels = c("", info$lab), expand = expansion(add = c(1.0, 1.2))) +
	coord_cartesian(clip = "off") +
	theme_minimal(base_size = 8.8) +
	theme(panel.grid = element_blank(), axis.line.x.bottom = element_line(colour = col_frame, linewidth = lw_frame), axis.text.x = element_text(angle = 45, hjust = 1, colour = "black", size = 10, face = "bold"), axis.text.y = element_text(size = 8.8, colour = "black"), plot.margin = margin(-1, 10, 8, left_margin_size))
pAB <- pAB_top / pAB_track + plot_layout(heights = c(2.2, 1.8)) & labs(x = NULL, y = NULL); pAB
ggsave("Fig3a.pdf", plot = pAB, width = 14, height = 8, device = cairo_pdf)

# Fig3b Inset Incremental AUC
roc_df <- d$pred_plot %>% filter(Block %in% m_blks) %>% split(.$Block) %>% lapply(\(x) data.frame(FPR = 1 - (ro <- roc(x$y_true, x$y_pred, quiet = TRUE))$specificities, TPR = ro$sensitivities, Block = x$Block[1])) %>% bind_rows()
sum2 <- d$sum_plot %>% filter(Block %in% m_blks) %>% mutate(Block = factor(Block, m_blks)) %>% arrange(Block)
lvec <- setNames(paste0(m_labs[as.character(sum2$Block)], " (AUC = ", sub(" .*", "", sum2$AUC), ")"), sum2$Block)
dd <- d$delta_sum %>% filter(compare == "Baseline_Diet7 - Baseline")
pv <- dd$p_wilcox[1]; pv_sci <- format(pv, scientific = TRUE, digits = 3); pv_split <- strsplit(pv_sci, "e")[[1]]
pv_m <- as.numeric(pv_split[1]); pv_e <- as.integer(pv_split[2])
lab_p <- sprintf("'Paired Wilcoxon ' * italic(P) * ' = %.2f' ~ '\\u00D7' ~ 10^{%d}", pv_m, pv_e)
pC <- ggplot(roc_df, aes(FPR, TPR, colour = factor(Block, m_blks))) +
    geom_path(linewidth = 1.2) + geom_abline(slope = 1, linetype = 2, colour = "grey80", linewidth = 0.7) + coord_equal() +
    scale_colour_manual(values = c(Baseline = C$line, Baseline_Diet7 = C$ora), labels = lvec) +
    annotate("text", x = .40, y = .16, label = lab_p, hjust = 0, vjust = 0, size = 4, colour = "grey20", fontface = "bold", parse = TRUE) +
    labs(x = "False positive rate", y = "True positive rate") + thm +
    theme(panel.grid = element_blank(), panel.grid.minor = element_blank(), legend.position = c(.34, .84), legend.background = element_blank(), legend.key = element_blank(), legend.text = element_text(size = 8.5, face = "bold"), plot.margin = margin(10, 20, 10, 10))
ggsave("Fig3b.pdf", pC, width = 5, height = 5, device = "pdf", useDingbats = FALSE)