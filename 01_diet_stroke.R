pacman::p_load(readxl, dplyr, data.table, tidyverse, ggrepel, lubridate, purrr, survival, ggplot2, forcats)

dir0 <- ifelse(Sys.info()[["sysname"]] == "Windows", "D:/", "/work/sph-huangj")
indir = paste0(dir0, "/data/ukb/phe")
invisible(lapply(c("phe.f.R", "assoc.f.R", "plot.f.R"), function(f) source(file.path(dir0, "scripts", "f", f))))

setwd("D:/analysis/stroke")

## ------------------------------------------------------------------
## 🚩 饮食 QC
## ------------------------------------------------------------------
std_10_90 <- function(x){  # 统一标准化到同一尺度（10th→90th 的增量） Per 10th-to-90th percentile increment
	x <- as.numeric(x); q <- quantile(x, c(0.10, 0.90), na.rm = TRUE, names = FALSE); if(!all(is.finite(q)) || q[2] == q[1]) rep(NA_real_, length(x)) else (x - q[1]) / (q[2] - q[1])
}
score_q5 <- function(x){ # 真正五分位 -> 0/25/50/75/100
	x <- as.numeric(x); out <- rep(NA_real_, length(x)); ok <- is.finite(x); out[ok] <- c(0, 25, 50, 75, 100)[dplyr::ntile(x[ok], 5)]; out
}
score_q4 <- function(x){ # modern的分组：先算四分位切点，再按真正的分位点阈值划区间（不拆相同分值，但可能因为分数太离散，最后组数变少）
	x <- as.numeric(x); q <- unique(quantile(x, c(0,.25,.5,.75,1), na.rm=TRUE, names=FALSE))
	if(length(q) < 3) factor(rep(NA_character_, length(x))) else cut(x, breaks=q, include.lowest=TRUE, right=TRUE, ordered_result=TRUE)
}
score_0_100 <- function(x){
	x <- as.numeric(x); ok <- is.finite(x); if(sum(ok) < 2) rep(NA_real_, length(x)) else {mn <- min(x[ok]); mx <- max(x[ok]); if(mx == mn) rep(NA_real_, length(x)) else {out <- rep(NA_real_, length(x)); out[ok] <- (x[ok] - mn) / (mx - mn) * 100; out}}
}
score_0_100_q5 <- function(x){ # 先缩放到 0-100，再按固定宽度切成 5 档
	x <- as.numeric(x); ok <- is.finite(x); if(sum(ok) < 2) rep(NA_real_, length(x)) else {mn <- min(x[ok]); mx <- max(x[ok]); if(mx == mn) rep(NA_real_, length(x)) else {out <- rep(NA_real_, length(x)); out[ok] <- (x[ok] - mn) / (mx - mn) * 100; as.numeric(as.character(cut(out, breaks = c(-Inf, 20, 40, 60, 80, Inf), labels = c(0, 25, 50, 75, 100), right = FALSE)))}}
}

dat0 <- readRDS(paste0(indir, "/Rdata/all.rds"))
# dat0 <- dat0 %>% left_join(readRDS("dose_score.rds"), by = "eid")
sum_cols <- names(dat0)[grepl("^diet\\..+\\.sum$", names(dat0))]
dat <- dat0 %>% mutate(
	across(all_of(sum_cols), std_10_90, .names = "{sub('\\\\.sum$', '', .col)}.std"),
    across(all_of(sum_cols), score_0_100, .names = "{sub('\\\\.sum$', '', .col)}.s100"),
	across(all_of(sum_cols), score_0_100_q5, .names = "{sub('\\\\.sum$', '', .col)}.pts"),
	across(all_of(sum_cols), score_q5, .names = "{sub('\\\\.sum$', '', .col)}.q5"),
	across(all_of(sum_cols), f3c, .names = "{sub('\\\\.sum$', '', .col)}.3c")
)
lapply(dat %>% select(matches("\\.sum$")), summary)
lapply(dat %>% select(matches("\\.(pts|q5)$")), table, useNA = "always")


## ------------------------------------------------------------------
## 🚩 生成变量
## ------------------------------------------------------------------
# p6154	drug.aspirin1	all
# p22036 pa_met_recom; p20107 fhist_father; p20110 fhist_mother; p20111 fhist_siblings
# Stroke，Ischemic stroke，Intracerebral hemorrhage，Subarachnoid hemorrhage
# early-onset ischemic stroke ([EOS]; onset 18–59 years) and late-onset ischemic stroke ([LOS]; onset ≥60 years)
dat <- dat %>% filter(enroll >= 1) %>% mutate( # %>% filter(gen_ethnicity==1) ; & as.vector(ethnic.c) == "White"
	ethnic.c3 = factor(ifelse(ethnic.c=="White", "White", "Others"), levels=c("White", "Others")),
	edu.c = factor(case_when(grepl("(^|\\|)1(\\||$)", as.character(edu)) ~ "high", grepl("(^|\\|)(2|5|6)(\\||$)", as.character(edu)) ~ "middle", grepl("(^|\\|)(3|4|7)(\\||$)", as.character(edu)) ~ "low", TRUE ~ NA_character_), levels = c("low", "middle", "high")),
	job.b = factor(ifelse(is.na(job), NA, ifelse(grepl("[1267]", as.character(job)), 1, 0)), levels = c(0, 1), labels = c("unemployed", "employed")),
	living_alone.b = factor(ifelse(householdCt == 1, 1, ifelse(householdCt > 1, 0, NA)), levels = c(0, 1), labels = c("no", "yes")),
	center = as.character(center),
	region = case_when(center %in% c("c11012", "c11020", "c11018") ~ "London", center %in% c("c11001", "c11016", "c10003", "c11024", "c11025", "c11008") ~ "North West", center %in% c("c11017", "c11009", "c11027") ~ "North East", center %in% c("c11010", "c11014") ~ "Yorkshire", center %in% c("c11021", "c11006") ~ "West Midlands", center %in% c("c11013") ~ "East Midlands", center %in% c("c11002", "c11007", "c11026") ~ "South East", center %in% c("c11011", "c11028") ~ "South West", center %in% c("c11005", "c11004") ~ "Scotland", center %in% c("c11003", "c11022", "c11023") ~ "Wales", TRUE ~ NA_character_),
	region = factor(region, levels = c("London", "North West", "North East", "Yorkshire", "West Midlands", "East Midlands", "South East", "South West", "Scotland", "Wales")), 
	region_code = as.integer(region) - 1,
	smoke.c = factor(smoke_status, levels = 0:2, labels = c("never", "previous", "current")),
	alcohol.c = factor(alcohol_status, levels = 0:2, labels = c("never", "previous", "current")),
	days_pa_mod.c = f3c(days_pa_mod),
	sleep_duration.b = factor(ifelse(!is.na(sleep_duration) & sleep_duration >= 7 & sleep_duration <= 9, 1, 0), levels = c(0, 1), labels = c("unhealthy", "healthy")),#7–9 小时为健康
	age.b = factor(ifelse(age < 60, "<60", "≥60"), levels = c("<60", "≥60")), # 用于亚组分析
	bmi.b = factor(ifelse(bmi < 25, "<25", ">=25"), levels = c("<25", ">=25")), # 用于亚组分析
	drug.aspirin = factor(as.integer(rowSums(across(starts_with("drug.aspirin1"), ~ grepl("1", as.character(.x)))) > 0)),	
	fhist_stroke = factor(rowSums(across(c(fhist_father, fhist_mother, fhist_siblings), ~ !is.na(.x) & grepl("(^|\\|)2(\\||$)", as.character(.x)))) > 0),
)
	drugs <- c("drug.lipid", "drug.htn", "drug.dm", "drug.aspirin")
	dxs <- c("primary_ahtn", "t2dm", "cvd") #  "cancer"
	date.cvd <- paste0("fod_icd10_", c("cvd_cad", "cvd_hfail", "cvd_afib", "cvd_pad"))
	dat$fod_icd10_cvd <- do.call(pmin, c(dat[, date.cvd], list(na.rm = TRUE)))
	for(n in dxs) { dat[[n]] <- !is.na(dat[[paste0("fod_icd10_", n)]]) & !is.na(dat$start_date) & (dat[[paste0("fod_icd10_", n)]] <= dat$start_date) }
dat <- dat %>% mutate(
    prior_qi_stroke = !is.na(fod_icd10_qi_stroke) & !is.na(start_date) & (fod_icd10_qi_stroke <= start_date),
	across(all_of(dxs), ~ factor(ifelse(., "Yes", "No"), levels = c("No", "Yes"))),
	across(all_of(drugs), ~factor(., levels = c(0, 1), labels = c("No", "Yes"))),
) %>% filter(!prior_qi_stroke) 
follow_end <- as.Date(date_follow_end)
dat <- dat %>% mutate(end0 = pmin(date_lost, date_death, follow_end, na.rm = TRUE)) %>% filter(is.na(end0) | end0 > start_date)
# I6[0-4]	qi_stroke;	I63	cvd_stroke_i	Ischemic stroke;	I61	cvd_stroke_ih	Intracerebral hemorrhage;	I60	cvd_stroke_sh	Subarachnoid hemorrhage	
Ys <- c("qi_stroke", "cvd_stroke_i", "cvd_stroke_ih", "cvd_stroke_sh") 
out_map <- c("qi_stroke" = "Overall stroke", "cvd_stroke_i" = "Ischemic stroke", "cvd_stroke_ih" = "Intracerebral hemorrhage", "cvd_stroke_sh" = "Subarachnoid hemorrhage")
for (Y in Ys) {
	dat[grep(paste0("^", Y, "\\.Y[?]t2e"), names(dat))] <- NULL
	dat <- t2e(dat, NA, paste0("fod_icd10_", Y), "birth_date", "start_date", "date_lost", "date_death", date_follow_end, Y, "year")
}
sapply(paste0(Ys, ".Yt2e"), \(x) table(dat[[x]], useNA = "ifany"))
dat <- dat %>% mutate(across(ends_with("Yt2e"), as.integer))
nrow(dat) # 202331


## ------------------------------------------------------------------
## covs
## ------------------------------------------------------------------
pacman::p_load(tidyr, readr)
covs.base <- c("age", "sex.b", paste0("PC", 1:10))
covs.add <- c("ethnic.c3", "edu.c", "tdi", "job.b", "living_alone.b", "enroll", "energy.kJ", "smoke.c", "days_pa_mod.c", "sleep_duration.b", "bmi") # strata(region)
diet_with_alcohol <- c("medi24", "mind", "ahei", "redii", "fds", "upf"); diet_without_alcohol <- c("dash", "hpdi", "phdi", "digm", "hlcd", "hlfd", "modern") #🏮
diet_lst <- c(diet_with_alcohol, diet_without_alcohol)
	diet_std <- unlist(lapply(diet_lst, \(d) grep(paste0(d, ".*\\.std$"), names(dat), value = TRUE)))
	diet_sum <- unlist(lapply(diet_lst, \(d) grep(paste0(d, ".*\\.sum$"), names(dat), value = TRUE)))
	diet_q5  <- unlist(lapply(diet_lst, \(d) grep(paste0(d, ".*\\.q5$"), names(dat), value = TRUE)))
get_model_list <- function(diet) { 
	alc_in_score <- diet %in% diet_with_alcohol
	m1 <- c(covs.base, covs.add, if (!alc_in_score) "alcohol.c", "fhist_stroke")
	m2 <- c(m1, dxs)  
	m3 <- c(m2, drugs) 
	list(m1 = unique(m1), m2 = unique(m2), m3 = unique(m3))
}
get_covs <- function(diet, model = c("m1", "m2", "m3")) {model <- match.arg(model); get_model_list(diet)[[model]]}
covs.all <- unique(c(covs.base, covs.add, dxs, drugs, "fhist_stroke", "alcohol.c")) 

dat1 <- dat %>% filter(if_all(all_of(paste0("PC", 1:10)), ~ !is.na(.)))
covs_miss <- tibble(variable = covs.all) %>% mutate(
	class = sapply(variable, \(v) class(dat1[[v]])[1]), n_total = nrow(dat1),
	n_missing = sapply(variable, \(v) sum(is.na(dat1[[v]]))), pct_missing = round(n_missing / n_total * 100, 2), n_nonmissing = n_total - n_missing,
	n_level = sapply(variable, \(v) {x <- dat1[[v]]; if(is.factor(x) || is.character(x)) dplyr::n_distinct(x, na.rm = TRUE) else NA_integer_})
) %>% arrange(desc(pct_missing), variable); covs_miss
nrow(dat1) # 198456
saveRDS(dat1, "dat1.rds") 
write.csv(covs_miss, "covs_miss.csv", row.names = FALSE)
table_s <- covs_miss %>% filter(n_missing > 0) %>% select(variable, n_missing) %>% rename(`Covariates` = variable, `No. of sample with missing covariates` = n_missing)
writexl::write_xlsx(table_s, "Table_S_Missing_Covariates.xlsx")

tab_event <- lapply(Ys, function(Y) {
	y <- paste0(Y, ".Yt2e"); t <- paste0(Y, ".t2e")
	case_n <- sum(dat1[[y]] == 1, na.rm = TRUE); total_n <- sum(!is.na(dat1[[y]]))
	py <- sum(dat1[[t]], na.rm = TRUE)
	ir <- case_n / py * 1e5
	tibble(outcome = Y, 
		total = total_n, case = case_n, case_pct = round(case_n / total_n * 100, 1), 
		person_years = round(py, 1), ir_1e5 = round(ir, 1), ir_label = paste0(round(ir, 1), " per 100,000 person-years"),
		mean_ageonset = round(mean(dat1$age[dat1[[y]] == 1] + dat1[[t]][dat1[[y]] == 1], na.rm = TRUE), 1), 
		mean_fu = round(mean(dat1[[t]], na.rm = TRUE), 1), sd_fu = round(sd(dat1[[t]], na.rm = TRUE), 1), median_fu = round(median(dat1[[t]], na.rm = TRUE), 1)
	)
}) %>% bind_rows()
tab_event <- tab_event %>% mutate(outcome = recode(outcome, !!!out_map)); tab_event 
write.csv(tab_event, "tab_event.csv", row.names = FALSE)

# 数值变量 → 中位数填补；分类变量 → Missing 单独一类
cont_vars <- c("bmi", "tdi", "enroll", "energy.kJ")
cat_vars  <- c("days_pa_mod.c", "job.b", "edu.c", "living_alone.b", "ethnic.c3", "smoke.c", "alcohol.c")
dat.imp <- dat1 %>% mutate(
	across(all_of(cont_vars), ~ replace(.x, is.na(.x), median(.x, na.rm = TRUE))),
	across(all_of(cat_vars), ~ forcats::fct_na_value_to_level(.x, level = "Missing"))
)
sapply(dat.imp[, c(cont_vars, cat_vars)], \(x) sum(is.na(x)))
saveRDS(dat.imp, "dat.imp.rds") 


## ------------------------------------------------------------------
## 分析
## ------------------------------------------------------------------
dat1 <- readRDS("dat1.rds")
dat.imp <- readRDS("dat.imp.rds")
dat <- dat1 %>% select(-all_of(covs.all)) %>% left_join(dat.imp %>% select(eid, all_of(covs.all)), by = "eid") # 换imp的协变量

run_cox <- function(dat, Y, X, covs, strata_var = NULL) {
	y <- paste0(Y, ".Yt2e"); t <- paste0(Y, ".t2e")
	stopifnot(all(c(y, t, X, covs) %in% names(dat)))
	rhs <- paste(c(X, covs), collapse = " + ")
	if (!is.null(strata_var)) {rhs <- paste(rhs, paste0("strata(", strata_var, ")"), sep = " + ")}
	fml <- as.formula(paste0("Surv(", t, ", ", y, ") ~ ", rhs))
	fit <- survival::coxph(fml, data = dat)
	est <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>% dplyr::filter(term == X)
	tibble(outcome = Y, exposure = X, total = fit$n, case = fit$nevent, HR = est$estimate, LCI = est$conf.low, HCI = est$conf.high, P = est$p.value)
}
run_cox_cat <- function(dat, Y, X, covs, strata_var = NULL, ref_level = NULL, ord_levels = NULL) {
	y <- paste0(Y, ".Yt2e"); t <- paste0(Y, ".t2e"); stopifnot(all(c(y, t, X, covs) %in% names(dat)))
	dat2 <- dat; x <- dat2[[X]]
	if(is.factor(x) || is.character(x)) { x <- factor(x, levels = if(is.null(ord_levels)) unique(x) else ord_levels); lv <- levels(x)
	} else { lv <- sort(unique(x[!is.na(x)])); x <- factor(x, levels = lv) }
	if(length(lv) < 2) stop("Exposure has <2 non-missing levels: ", X)
	# -------- (A) 分组：其余组 vs 参考组 --------
	dat2[[X]] <- relevel(x, ref = as.character(if(is.null(ref_level)) lv[1] else ref_level))
	f1 <- as.formula(paste0("Surv(", t, ", ", y, ") ~ ", paste(c(X, covs, if(!is.null(strata_var)) paste0("strata(", strata_var, ")")), collapse = " + ")))
	fit1 <- survival::coxph(f1, data = dat2)
	tab1 <- broom::tidy(fit1, exponentiate = TRUE, conf.int = TRUE) %>% filter(grepl(paste0("^", X), term)) %>%
		transmute(outcome = Y, exposure = X, type = "Cat", level = term, total = fit1$n, case = fit1$nevent, HR = estimate, LCI = conf.low, HCI = conf.high, P = p.value)
	# -------- (B) 趋势：按等级顺序当连续 --------
	dat2[[paste0(X, "_ord")]] <- as.numeric(dat2[[X]])
	f2 <- as.formula(paste0("Surv(", t, ", ", y, ") ~ ", paste(c(paste0(X, "_ord"), covs, if(!is.null(strata_var)) paste0("strata(", strata_var, ")")), collapse = " + ")))
	fit2 <- survival::coxph(f2, data = dat2)
	tab2 <- broom::tidy(fit2, exponentiate = TRUE, conf.int = TRUE) %>% filter(term == paste0(X, "_ord")) %>%
		transmute(outcome = Y, exposure = X, type = "Ptrend", level = paste0("trend(", paste(lv, collapse = "-"), ")"), total = fit2$n, case = fit2$nevent, HR = estimate, LCI = conf.low, HCI = conf.high, P = p.value)
	bind_rows(tab1, tab2)
}
run_models_diet <- function(dat, Ys, diets, FUN, model_fun = get_model_list, strata_var = NULL, exposure_suffix = "\\.std$", ...) {
	lapply(Ys, function(Y) lapply(diets, function(diet) {
		if(length(x <- grep(paste0("^diet\\.", diet, exposure_suffix), names(dat), value = TRUE)) != 1) stop("Cannot uniquely find exposure for diet: ", diet)
		lapply(names(ml <- model_fun(diet)), function(m) { out <- FUN(dat=dat, Y=Y, X=x, covs=ml[[m]], strata_var=strata_var, ...); out$model <- m; out$diet <- diet; out }) |> dplyr::bind_rows()
	}) |> dplyr::bind_rows()) |> dplyr::bind_rows()
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
# 🚩 Pattern
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
res_std <- run_models_diet(dat = dat, Ys = Ys, diets = diet_lst, FUN = run_cox, model_fun = get_model_list, strata_var = "region", exposure_suffix = "\\.std$")
	res_std %>% count(outcome, model)
	write.csv(res_std, "res_std.csv", row.names = FALSE)
res_q5 <- run_models_diet(dat = dat, Ys = Ys, diets = diet_lst, FUN = run_cox_cat, model_fun = get_model_list, strata_var = "region", exposure_suffix = "\\.q5$", ref_level = 0)
	write.csv(res_q5, "res_q5.csv", row.names = FALSE)
res_3c <- run_models_diet(dat = dat, Ys = Ys, diets = diet_lst, FUN = run_cox_cat, model_fun = get_model_list, strata_var = "region", exposure_suffix = "\\.3c$",ord_levels = c("low", "middle", "high"), ref_level = "low")
	write.csv(res_3c, "res_3c.csv", row.names = FALSE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
# 🚩 Food
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
field <- read_excel(paste0(dir0, "/files/foods.xlsx"), sheet = "Sheet2") %>% filter (doseFood > 0) ##🏮doseScreen
	field$Food.Group[field$name == "p103100"] <- "Organ_meats" # 🏮 修改特定食物的分组
dat.food <- readRDS(paste0(indir, "/Rdata/diet.fudan.rds")) %>% mutate(eid = as.character(eid)) %>% select(-any_of(c("energy.kJ", "enroll", "start_date")))
dat.grp <- dat.food %>% select(eid)
	grp.link <- split(field$name, field$Food.Group)
	for (i in names(grp.link)) {dat.grp[[i]] <- rowSums(as.matrix(dat.food[, grp.link[[i]], drop = FALSE]), na.rm = TRUE)}
	length(names(dat.grp)) 
dat <- dat %>% left_join(dat.food, by = "eid") %>% left_join(dat.grp, by = "eid")
# saveRDS(dat, "dat.grp.rds") 

covs.food <- setdiff(unique(c(covs.base, covs.add, "fhist_stroke", dxs)), "region") # strata_var = "region"
Y <- "cvd_stroke_i" # "qi_stroke"
y <- paste0(Y, ".Yt2e"); t <- paste0(Y, ".t2e")

dat1 <- dat %>% mutate(across(all_of(field$name), ~ifelse(. == 0, 0, 1))) # 把食物变量二分类化
output <- data.frame()
for (p in field$name) {
	temp <- dat1 %>% select(eid, all_of(c(y, t, p, covs.food, "region"))) %>% rename(pheno = all_of(p)) %>% na.omit()
	out <- run_cox(dat = temp, Y = Y, X = "pheno", covs = covs.food, strata_var = "region")
	total <- temp %>% group_by(pheno) %>% summarise(total = n(), .groups = "drop"); case <- temp %>% group_by(pheno) %>% summarise(case = sum(.data[[y]]), .groups = "drop")
	zph <- cox.zph(coxph(as.formula(paste0("Surv(", t, ", ", y, ") ~ pheno + ", paste(covs.food, collapse = " + "), " + strata(region)")), data = temp))
	output0 <- data.frame(exposure = rep(p, 2), class = 1:2, outcome = rep(Y, 2), total = total$total[match(c(0, 1), total$pheno)], case = case$case[match(c(0, 1), case$pheno)], HR = c(1, out$HR), LCI = c(1, out$LCI), UCI = c(1, out$HCI), P = c(1, out$P), res.zph = zph$table["pheno", "p"])
	output <- rbind(output, output0)
}
#FDR <- function(p, n) {row <- length(p); for (i in 1:n) p[seq(i, row, n)] <- p.adjust(p[seq(i, row, n)], method = "fdr"); p}
#output$FDR <- FDR(output$P, 2)
output_final <- output %>% filter(class == 2) %>% mutate(FDR = p.adjust(P, method = "fdr"))
output_ref <- output %>% filter(class == 1) %>% mutate(FDR = 1)
output <- bind_rows(output_ref, output_final) %>% arrange(exposure, class)
write.csv(output, paste0("Food_", Y, "_.csv"), row.names = FALSE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Group
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
func_quin <- function(x){ # 五分位
	breaks <- quantile(x, probs = seq(0, 1, 1/5), na.rm = TRUE)
	if(length(unique(breaks)) == 2){a <- factor(ifelse(x > 0, 1, 0))
	}else if(length(unique(breaks)) < 6){a <- cut(x, breaks = c(min(x, na.rm = TRUE) - 1, unique(breaks)), include.lowest = TRUE)
	}else{a <- cut(x, breaks = unique(breaks), include.lowest = TRUE)}
	return(a)
}
dat <- dat %>% mutate(
	Berries_Citrus = as.integer(func_quin(Berries + Citrus)),
	Citrus_ApplesPears_Driedfruit = as.integer(func_quin(Citrus + Apples_Pears + Dried_fruit))
#	Processed_red_meats = as.integer(func_quin(Processed_meats + Red_meats))
)
dat1 <- dat %>% mutate(across(all_of(unique(field$Food.Group)), func_quin))

# 统计（各分位的中位数 / Range / n / case）
stat <- function(x, data_q, data_num, Y){
	y <- paste0(Y, ".Yt2e")
	temp <- data_q %>% select(eid, all_of(y), quantile = all_of(x)) %>% left_join(data_num %>% select(eid, Food = all_of(x)), by = "eid")
	temp %>% group_by(quantile) %>% 
		summarise(Median = round(median(Food, na.rm = TRUE), 2), Range = paste(round(range(Food, na.rm = TRUE), 2), collapse = "~"), n = n(), case = sum(.data[[y]], na.rm = TRUE), .groups = "drop") %>% 
		mutate(pheno = x, .before = "quantile") %>% mutate(`Consumption level` = paste0(Median, " (", Range, ")"), .after = "quantile")
}
stat_data <- data.frame()
for (i in unique(field$Food.Group)) {
	stat0 <- stat(i, dat1, dat, Y)
	stat_data <- rbind(stat_data, stat0)
}
write.csv(stat_data, "Group_stat.csv", row.names = FALSE)

# cox
dat2 <- mutate(dat1, across(all_of(unique(field$Food.Group)), as.integer))  # 把 quintile 因子转成整数
output <- data.frame()
for (ref in 1:5) {
	for (component in unique(field$Food.Group)) {
		temp <- dat2 %>% select(eid, all_of(c(t, y, component, covs.food, "region"))) %>% rename(pheno = all_of(component))
		lv <- sort(unique(temp$pheno))
		ref_use <- min(ref, max(lv))
		lv <- c(ref_use, lv[lv != ref_use])
		temp$pheno <- factor(temp$pheno, levels = lv)
		fit <- coxph(as.formula(paste0("Surv(", t, ", ", y, ") ~ pheno + ", paste(covs.food, collapse = " + "), " + strata(region)")), data = temp)
		s <- summary(fit); z <- cox.zph(fit); coef <- as.data.frame(s$coefficients); ci <- as.data.frame(s$conf.int)
		k <- length(lv) - 1L
		output0 <- data.frame(Ref = rep(ref_use, length(lv)), exposure = rep(component, length(lv)), outcome = rep(Y, length(lv)), class = levels(temp$pheno), HR = c(1, coef$`exp(coef)`[1:k]), LCI = c(1, ci$`lower .95`[1:k]), HCI = c(1, ci$`upper .95`[1:k]), P = c(1, coef$`Pr(>|z|)`[1:k]), res.zph = z$table["pheno", "p"])
		output <- rbind(output, output0)
	}
}
write.csv(output, paste0("Group_", Y, "_.csv"), row.names = FALSE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩判断两个结果的方向一致性
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 稳健性检验
dat.sens <- dat %>% filter(cvd == "No", primary_ahtn == "No", t2dm == "No")  #172995 
# paste0("Group_", Y, "_sens.csv")
get_group_dir <- function(f) {
	read.csv(f) %>% mutate(Ref = as.integer(Ref), class = suppressWarnings(as.integer(as.character(class)))) %>%
		filter(!is.na(class), Ref == 1) %>% group_by(exposure) %>%
		summarise(
			n_level = n(), min_class = class[which.min(HR)], max_class = max(class, na.rm = TRUE),
			best_pos = ifelse(max_class == 1, NA_real_, (min_class - 1) / (max_class - 1)),
			shape = case_when(is.na(best_pos) ~ NA_character_, best_pos <= 1/3 ~ "low_intake_best", best_pos >= 2/3 ~ "high_intake_best", TRUE ~ "middle_best"),
			slope = if(n() >= 3) coef(lm(log(HR) ~ class))[2] else NA_real_,
			dir = case_when(is.na(slope) ~ NA_character_, slope > 0 ~ "higher_intake_higher_risk", slope < 0 ~ "higher_intake_lower_risk", TRUE ~ "flat"),
			p_min = min(P[class != 1], na.rm = TRUE), .groups = "drop"
		)
}
dir_cmp <- full_join(get_group_dir("Group_cvd_stroke_i_.csv"), get_group_dir("Group_cvd_stroke_i_sens.csv"), by = "exposure", suffix = c("_main", "_sens")) %>% mutate(
    same_dir = dir_main == dir_sens, same_shape = shape_main == shape_sens,
    broad_consistency = case_when(same_dir & same_shape ~ "high_consistency", same_dir ~ "same_direction_shape_shift", same_shape ~ "same_best_zone_but_slope_diff", TRUE ~ "direction_flip")
) %>% arrange(broad_consistency, p_min_main, p_min_sens)
write.csv(dir_cmp, "Group_direction_compare.csv", row.names = FALSE)
# 结果：21 个候选组里大致可分为三类
# 15 个属于 high_consistency，说明主分析和无病敏感性分析在“最低风险区间”和“整体方向”上都一致；
# 5 个属于 same_direction_shape_shift，说明方向一致，但最佳摄入区间位置有些移动；
# 1 个（Wine）属于 same_best_zone_but_slope_diff，说明最低风险区间一致，但线性趋势方向有差别。
# 也就是说，21 个候选组里没有一个是“主分析保护、敏感性分析有害”这种彻底翻向的组。


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 ML 数据准备💁
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
group_res <- read.csv(paste0("Group_", Y, "_.csv"))
group_sig_anyref <- group_res %>% mutate(class = suppressWarnings(as.integer(as.character(class)))) %>% filter(!is.na(class)) %>% group_by(exposure) %>%
	summarise(p_min = min(P, na.rm = TRUE), best_ref = Ref[which.min(P)], best_class = class[which.min(P)], n_sig = sum(P < 0.05, na.rm = TRUE), .groups = "drop") %>%
	filter(n_sig >= 1)  # Any Ref 筛显著组:不固定 Q1 为参考，只要任意 Ref 下任意一个 level 显著，就纳入候选
	group_sig <- group_sig_anyref$exposure; cat("候选显著组数量 =", length(group_sig), "\n"); print(group_sig)
	score.tbl <- group_sig_anyref %>% transmute(exposure, score = p_min)
# Spearman 相关聚类去冗余: 若 |rho| > 0.6，认为属于同一 cluster; 每个 cluster 保留 Cox 证据最强（P最小）的那个
reduce_correlated <- \(cor_mat, score_tbl, cutoff = .6, method = "average") {
	hc <- hclust(as.dist(1 - abs(cor_mat)), method = method); cl <- cutree(hc, h = 1 - cutoff)
	sel <- sapply(sort(unique(cl)), \(k) score_tbl %>% filter(exposure %in% names(cl)[cl == k]) %>% arrange(score) %>% slice(1) %>% pull(exposure))
	list(selected = unname(sel), clusters = cl, hc = hc)
}
cor_spear <- cor(dat %>% select(all_of(group_sig)) %>% as.matrix(), method = "spearman", use = "pairwise.complete.obs")
	red <- reduce_correlated(cor_spear, score.tbl, cutoff = .6)
	feat_use <- red$selected
	cat("去冗余后用于机器学习的食物组数量 =", length(feat_use), "\n"); print(feat_use)
	plot(hclust(as.dist(1 - abs(cor_spear)), method = "average"), main = "Spearman clustering of food groups")
dat_res <- dat %>% select(eid, all_of(covs.food), all_of(feat_use))
# 残差化（residualization）: 每个 group 对 covs.food 回归，取 residual
	get_residual <- \(y, cov_df) as.numeric(residuals(lm(y ~ ., data = data.frame(y = y, cov_df), na.action = na.exclude)))
	for (v in feat_use) dat_res[[paste0(v, "_res")]] <- get_residual(dat_res[[v]], dat_res[, covs.food, drop = FALSE])
	group_residual_vars <- paste0(feat_use, "_res")
dat_ml <- dat_res %>% left_join(dat %>% select(eid, all_of(paste0(Y, c(".t2e", ".Yt2e"))), region_code), by = "eid") %>% select(eid, all_of(group_residual_vars), all_of(paste0(Y, c(".t2e", ".Yt2e"))), region_code)
write.csv(dat_ml, "dat_ml.csv", row.names = FALSE)
saveRDS(dat_ml, "dat_ml.rds")

cat("候选组 =", length(group_sig), "; 去冗余后 =", length(feat_use), "; ML维度 =", nrow(dat_ml), "x", ncol(dat_ml), "; 残差化变量数 =", length(group_residual_vars), "\n")
dat_ml %>% group_by(region_code) %>% summarise(n = n(), case = sum(.data[[y]], na.rm = TRUE), case_pct = round(mean(.data[[y]], na.rm = TRUE) * 100, 3), .groups = "drop")
colSums(is.na(dat_ml))
summary(dat_ml[, group_residual_vars])


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 👉去跑机器学习的python代码
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

