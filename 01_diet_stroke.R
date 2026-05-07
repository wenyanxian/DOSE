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
# p6154	drug.aspirin1	all; p20107 fhist_father; p20110 fhist_mother; p20111 fhist_siblings
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
	age.b = factor(ifelse(age < 60, "<60", "≥60"), levels = c("<60", "≥60")),
	bmi.b = factor(ifelse(bmi < 25, "<25", ">=25"), levels = c("<25", ">=25")), 
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
table_s <- covs_miss %>% filter(n_missing > 0) %>% select(variable, n_missing) %>% rename(`Covariates` = variable, `No. of sample with missing covariates` = n_missing); writexl::write_xlsx(table_s, "Table_S_Missing_Covariates.xlsx")

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
tab_event <- tab_event %>% mutate(outcome = recode(outcome, !!!out_map)); write.csv(tab_event, "tab_event.csv", row.names = FALSE)

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
dat <- dat1 %>% select(-all_of(covs.all)) %>% left_join(dat.imp %>% select(eid, all_of(covs.all)), by = "eid") 
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

# 🚩 Food
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

# 🚩 Group
func_quin <- function(x){ # 五分位
	breaks <- quantile(x, probs = seq(0, 1, 1/5), na.rm = TRUE)
	if(length(unique(breaks)) == 2){a <- factor(ifelse(x > 0, 1, 0))
	}else if(length(unique(breaks)) < 6){a <- cut(x, breaks = c(min(x, na.rm = TRUE) - 1, unique(breaks)), include.lowest = TRUE)
	}else{a <- cut(x, breaks = unique(breaks), include.lowest = TRUE)}
	return(a)
}
dat <- dat %>% mutate(Citrus_ApplesPears_Driedfruit = Citrus + Apples_Pears + Dried_fruit)
food_groups <- unique(c(unique(field$Food.Group), "Citrus_ApplesPears_Driedfruit"))
dat1 <- dat %>% mutate(across(all_of(food_groups), func_quin))
stat <- function(x, data_q, data_num, Y){ 
	y <- paste0(Y, ".Yt2e")
	temp <- data_q %>% select(eid, all_of(y), quantile = all_of(x)) %>% left_join(data_num %>% select(eid, Food = all_of(x)), by = "eid")
	temp %>% group_by(quantile) %>% 
		summarise(Median = round(median(Food, na.rm = TRUE), 2), Range = paste(round(range(Food, na.rm = TRUE), 2), collapse = "~"), n = n(), case = sum(.data[[y]], na.rm = TRUE), .groups = "drop") %>% 
		mutate(pheno = x, .before = "quantile") %>% mutate(`Consumption level` = paste0(Median, " (", Range, ")"), .after = "quantile")
}
stat_data <- data.frame()
for (i in food_groups) {
	stat0 <- stat(i, dat1, dat, Y)
	stat_data <- rbind(stat_data, stat0)
}
write.csv(stat_data, "Group_stat.csv", row.names = FALSE)
dat2 <- mutate(dat1, across(all_of(food_groups), as.integer)) 
output <- data.frame()
for (ref in 1:5) {
	for (component in food_groups) {
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


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 🚩 RCS plot + FigS2 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pacman::p_load(rms, ggplot2, dplyr, patchwork)
fmt_p <- \(p) ifelse(is.na(p), "NA", ifelse(p < .001, "< 0.001", sprintf("= %.3f", p)))
get_k3 <- \(x) {
    x <- x[is.finite(x)]; k <- quantile(x, c(.1, .5, .9), names = F, type = 2)
    if(length(unique(k)) < 3 && length(x1 <- x[x > 0]) >= 20) k <- quantile(x1, c(.1, .5, .9), names = F, type = 2)
    if(length(unique(k)) < 3) NULL else k
}
run_rcs_one <- \(dat, t, e, x, covs, strat = "region") {
    d <- na.omit(dat[c(t, e, x, covs, strat)])
    if(nrow(d) < 100 || n_distinct(d[[x]]) < 5 || is.null(k <- get_k3(d[[x]]))) return(NULL)    
    dd <<- datadist(d); options(datadist = "dd")
    f <- as.formula(paste0("Surv(", t, ",", e, ")~rcs(", x, ",c(", paste(k, collapse=","), "))+", paste(c(covs, paste0("strat(", strat, ")")), collapse="+")))
    fit <- cph(f, data = d, x = T, y = T)
    an <- anova(fit); p1 <- an[x, "P"]; p2 <- an[grep("Nonlinear", rownames(an))[1], "P"]    
    xm <- quantile(d[[x]], .995)
    pr <- Predict(fit, name=x, ref.zero=T, fun=exp, np=200) %>% as.data.frame() %>% filter(!!sym(x) <= xm)
    h <- hist(d[[x]][d[[x]] <= xm], breaks=20, plot=F); ym <- max(pr$upper)*1.08; sc <- ym/max(h$counts)   
    p <- ggplot(pr, aes(.data[[x]])) +
        geom_histogram(data = filter(d, .data[[x]] <= xm), aes(y = after_stat(count)*sc), bins=20, fill="#A6CEE3", col="white", alpha=.7) +
        geom_hline(yintercept=1, lty=2, col="grey60") + geom_line(aes(y=yhat), col="#E41A1C", lwd=1.1) + geom_line(aes(y=lower), col="#E41A1C", lty=2, lwd=.8) + geom_line(aes(y=upper), col="#E41A1C", lty=2, lwd=.8) +
        scale_y_continuous("HR (95% CI)", sec.axis = sec_axis(~./sc, name="Frequency")) +
        annotate("text", x = min(pr[[x]]) + .02*diff(range(pr[[x]])), y = ym, label = paste0("P overall ", fmt_p(p1), "\nP non-linear ", fmt_p(p2)), hjust=0, vjust=1, size=3.5) +
        labs(title = ifelse(exists("nice_lab") && x %in% names(nice_lab), nice_lab[x], x), x = NULL) +
        coord_cartesian(xlim = c(min(pr[[x]]), xm)) + theme_bw(11) + theme(panel.grid = element_blank(), plot.title = element_text(hjust=.5, face="bold", size=11, margin=margin(b=10)), axis.title.y.left = element_text(size=9, margin=margin(r=5)), axis.title.y.right = element_text(size=9, margin=margin(l=5)), axis.text = element_text(size=8), plot.margin = margin(15, 15, 15, 15))    
    list(plot = p, st = data.frame(Exposure = x, P_Overall = p1, P_Nonlinear = p2, P_All_T = fmt_p(p1), P_NL_T = fmt_p(p2)))
}
out_dir <- "RCS"; plot_dir <- file.path(out_dir, "single_plots")
dir.create(out_dir, F, T); dir.create(plot_dir, F, T)
res <- Filter(Negate(is.null), lapply(setNames(food_groups, food_groups), \(x) tryCatch(run_rcs_one(dat, paste0(Y, ".t2e"), paste0(Y, ".Yt2e"), x, covs.food), error=\(e) NULL)))
plts <- lapply(res, `[[`, "plot")
rcs_p <- bind_rows(lapply(res, `[[`, "st")) %>% arrange(P_Nonlinear, P_Overall)
write.csv(rcs_p, file.path(out_dir, paste0("RCS_", Y, ".csv")), row.names=F)
invisible(Map(\(p, nm) ggsave(file.path(plot_dir, paste0(nm, "_", Y, ".pdf")), p, width=5.2, height=4.2), plts, names(plts)))
save_grid <- \(plist, suffix, title) {
    if(!length(plist)) return()
    p_grid <- wrap_plots(plist, ncol=4) + plot_annotation(title=title, theme=theme(plot.title=element_text(size=16, face="bold", hjust=.5)))
    ggsave(file.path(out_dir, paste0("fig_rcs_", suffix, "_", Y, ".pdf")), p_grid, width = 4*4.5, height = ceiling(length(plist)/4)*3.5)
}
save_grid(plts, "all", paste0("Restricted cubic spline associations of significant food groups with ", out_map[Y]))
save_grid(plts[rcs_p$Exposure[rcs_p$P_Nonlinear < .05]], "nonlinear", "Restricted cubic spline associations of food groups with significant nonlinearity")
rcs_p


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Fig 2
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pacman::p_load(dplyr, tidyr, ggplot2, stringr, scales, ggnewscale, cowplot)
norm_key <- function(x) str_to_lower(str_squish(as.character(x)))
excl <- "Citrus_ApplesPears_Driedfruit"; v_exp <- setdiff(unique(field$Food.Group), excl)
dat_mean <- dat %>% select(all_of(v_exp)) %>% summarise(across(everything(), ~mean(.x, na.rm=TRUE))) %>% pivot_longer(everything(), names_to="exposure", values_to="mean_intake")
cat_map <- field %>% distinct(exposure=Food.Group, Category) %>% 
	mutate(Category=case_when(
				Category %in% c("Meats","Fish_fish_dishes") ~ "Meats & Fish", 
				Category %in% c("Non_alcoholic_beverages","Alcoholic_beverages") ~ "Water & beverage",
				Category=="Fried_and_fast_foods" ~ "Fried & Fast food", Category=="Legumes_Nuts" ~ "Legumes & Nuts", Category=="Snack_Pastry" ~ "Snack & Pastry", Category=="Cereal_Bread" ~ "Cereal & Bread", Category=="Fats_Oil" ~ "Fats & Oil", Category=="Dairy_dairy_free_products" ~ "Dairy", TRUE ~ str_replace_all(as.character(Category), "_", " ")
			), Category=replace_na(Category, "Unknown")
)
raw_group <- read.csv(sprintf("Group_%s_.csv", Y)) %>% filter(!exposure %in% excl, exposure %in% v_exp)
rcs_dat <- read.csv(file.path("RCS", sprintf("RCS_%s.csv", Y))) %>% rename(exposure=Exposure) %>% filter(!exposure %in% excl)
sig_groups <- raw_group %>% transmute(exposure, class=as.integer(class), Ref=as.integer(Ref), P) %>% filter(!is.na(class), !is.na(P), class != Ref) %>% count(exposure, wt=P < .05, name="n_sig") %>% filter(n_sig > 0) %>% pull(exposure)
sig_key <- norm_key(sig_groups)
plot_dat <- raw_group %>% mutate(class=as.integer(class), Ref=as.integer(Ref)) %>% group_by(exposure) %>% mutate(best_ref=class[Ref==1][which.min(replace(HR[Ref==1], is.na(HR[Ref==1]), Inf))]) %>% ungroup() %>% filter(Ref==best_ref) %>% left_join(cat_map, by="exposure") %>% left_join(dat_mean, by="exposure") %>% left_join(rcs_dat, by="exposure") %>% group_by(exposure, Category) %>% mutate(logHR=log(HR), ref=class==best_ref, score=min(ifelse(ref|is.na(P), 1, P), na.rm=TRUE), eff=max(abs(logHR[!ref]), na.rm=TRUE), star=case_when(ref|is.na(P)~"", P<.001~"***", P<.01~"**", P<.05~"*", TRUE~""), nl_logp=-log10(pmax(P_Nonlinear, 1e-300)), nl_star=case_when(is.na(P_Nonlinear)~"", P_Nonlinear<.001~"***", P_Nonlinear<.01~"**", P_Nonlinear<.05~"*", TRUE~"")) %>% ungroup() %>% mutate(fw_sig=norm_key(exposure) %in% sig_key) %>% arrange(Category, desc(fw_sig), score, desc(eff), exposure) %>% mutate(exposure=factor(exposure, levels=unique(exposure)), id=as.integer(exposure))
n_id <- max(plot_dat$id); gap_l <- -0.35; gap_r <- n_id + 0.95; x_s <- gap_r - gap_l
y_c_min <- 1.02; y_c_max <- 1.18; q_t <- .38; q_p <- .41; y_q_b <- 1.35; y_q_t <- y_q_b + 4*q_p + q_t
y_n_min <- y_q_t + .18; y_n_max <- y_n_min + .20; y_b_b <- y_n_max + .18
b_ls <- data.frame(y=c(y_q_b-.08, y_q_t+.055, y_n_min-.055, y_n_max+.055, y_b_b-.055, y_b_b))
col_hr <- c("#5F8796","#FAF7F0","#BE665E"); col_nl <- c("#F2F5F5","#6F8790"); col_mean <- c("#E3EEE5","#6E9C7B")
pal_cat <- c("Cereal & Bread"="#C5A56A","Dairy"="#8FAF9A","Eggs"="#D8BF73","Fats & Oil"="#B4A37B","Fried & Fast food"="#CF9273","Fruits"="#D1A09A","Legumes & Nuts"="#A1B36F","Meats & Fish"="#91AFC1","Snack & Pastry"="#C3899D","Vegetables"="#7EA986","Water & beverage"="#89A9BF","Unknown"="#B8B8B8")
group_info <- distinct(plot_dat, exposure, Category, id, fw_sig, mean_intake, nl_logp, nl_star) %>% mutate(exp_key=norm_key(exposure), lab=str_replace_all(exposure,"_"," "), mean_h=rescale(mean_intake,to=c(.20,.82)), bar_t=y_b_b+mean_h, ang0=(90-360*(id-gap_l)/x_s+180)%%360-180, hj=ifelse(ang0 < -90 | ang0 > 90, 1, 0), ang=ifelse(ang0 < -90 | ang0 > 90, ang0+180, ang0), lab_y=max(y_b_b+mean_h,na.rm=TRUE)+.22, sig_lab=exp_key %in% sig_key, lab_col=ifelse(sig_lab, col_hr[3], "#2F3742"))
strip_dat <- plot_dat %>% mutate(ymin=y_q_b+(class-1)*q_p, ymax=ymin+q_t, ymid=(ymin+ymax)/2, star_col=ifelse(logHR>0,col_hr[3],col_hr[1]))
c_sep <- group_info %>% group_by(Category) %>% summarise(x=max(id)+.5,.groups="drop") %>% head(-1)
mx_q <- max(abs(strip_dat$logHR[is.finite(strip_dat$logHR)]),na.rm=TRUE); mx_nl <- max(group_info$nl_logp[is.finite(group_info$nl_logp)],na.rm=TRUE)
mean_min <- min(group_info$mean_intake,na.rm=TRUE); mean_max <- max(group_info$mean_intake,na.rm=TRUE)
cats <- sort(unique(group_info$Category)); cat_cols <- pal_cat[cats]; cat_cols[is.na(cat_cols)] <- "#B8B8B8"
gap_w <- (gap_r-(n_id+.5))+(.5-gap_l); x_ann <- n_id+.5+gap_w/2; if(x_ann > gap_r) x_ann <- gap_l+(x_ann-gap_r)
ann_df <- tibble::tibble(x=x_ann, y=c(y_b_b+.34,(y_n_min+y_n_max)/2,y_q_t+.075,y_q_b+4*q_p+q_t/2,y_q_b+q_t/2), lab=c("Average\nconsumption level\n(servings/day)","RCS\n(nonlinear P)","HR","Q5","Q1"), col=c(col_mean[2],col_nl[2],"#374151","#6B7280","#6B7280"), face=c("bold","bold","bold","plain","plain"), size=c(2.75,2.70,3.20,2.55,2.55))
p_main <- ggplot() +
	geom_segment(data=b_ls, aes(x=.5,xend=max(group_info$id)+.5,y=y,yend=y), color="#E7E7E7", linewidth=.25) + geom_segment(data=c_sep, aes(x=x,xend=x,y=y_c_min,yend=y_b_b+.80), color="#ECECEC", linewidth=.22, linetype="dashed") +
	geom_rect(data=group_info, aes(xmin=id-.47,xmax=id+.47,ymin=y_c_min,ymax=y_c_max,fill=Category), color="white", linewidth=.09) +
	scale_fill_manual(values=cat_cols, guide="none") + ggnewscale::new_scale_fill() +
	geom_rect(data=strip_dat, aes(xmin=id-.47,xmax=id+.47,ymin=ymin,ymax=ymax,fill=logHR), color="white", linewidth=.07) +
	geom_point(data=filter(strip_dat,ref), aes(x=id,y=ymid), shape=21, size=.48, stroke=.18, fill="white", color="#1F2937") +
	geom_text(data=filter(strip_dat,!ref & P<.05), aes(x=id,y=ymid,label=star,color=star_col), size=2.55, fontface="bold") +
	scale_fill_gradient2(low=col_hr[1],mid=col_hr[2],high=col_hr[3],midpoint=0,limits=c(-mx_q,mx_q),oob=squish,guide="none") + ggnewscale::new_scale_fill() +
	geom_rect(data=group_info, aes(xmin=id-.44,xmax=id+.44,ymin=y_n_min,ymax=y_n_max,fill=nl_logp), color="white", linewidth=.07, alpha=.90) +
	geom_text(data=filter(group_info,nl_star!=""), aes(x=id,y=(y_n_min+y_n_max)/2,label=nl_star), size=2.45, fontface="bold", color="#1F2937") +
	scale_fill_gradient(low=col_nl[1],high=col_nl[2],limits=c(0,max(2,mx_nl)),oob=squish,guide="none") + ggnewscale::new_scale_fill() +
	geom_rect(data=group_info, aes(xmin=id-.47,xmax=id+.47,ymin=y_b_b,ymax=bar_t,fill=mean_intake), color="white", linewidth=.045, alpha=.93) +
	scale_fill_gradient(low=col_mean[1],high=col_mean[2],guide="none") +
	geom_text(data=ann_df, aes(x=x,y=y,label=lab,color=col,fontface=face,size=size), angle=0, hjust=.5, vjust=.5, lineheight=.82) +
	geom_text(data=group_info, aes(x=id,y=lab_y,label=lab,angle=ang,hjust=hj,color=lab_col), size=3.00, lineheight=.78, fontface="bold") +
	scale_color_identity() + scale_size_identity() + coord_polar(clip="off") +
	scale_x_continuous(limits=c(gap_l,gap_r),breaks=NULL) + scale_y_continuous(limits=c(-.72,max(group_info$lab_y)+.88),breaks=NULL) +
	theme_void() + theme(plot.margin=margin(2,2,2,2))
mk_g <- function(y, cols, x0=.25, x1=.75, n=100) data.frame(x=seq(x0,x1,length.out=n), y=y, w=(x1-x0)/n*1.05, fill=colorRampPalette(cols)(n))
nr <- ceiling(length(cats)/2)
cat_df <- tibble::tibble(cat=cats, col=unname(cat_cols[cats]), i=seq_along(cats)) %>% mutate(col_id=(i-1)%/%nr+1, row_id=(i-1)%%nr+1, x=ifelse(col_id==1,.12,.56), y=.88-(row_id-1)*.04)
b_hr <- mk_g(.365,col_hr); b_nl <- mk_g(.225,col_nl); b_mn <- mk_g(.085,col_mean)
p_leg <- ggplot() + theme_void() + coord_cartesian(xlim=0:1,ylim=0:1,clip="off") +
	annotate("text",.5,.97,label="Food Category",fontface="bold",size=3.0,color="#1F2937") +
	geom_rect(data=cat_df,aes(xmin=x,xmax=x+.012,ymin=y-.014,ymax=y+.014,fill=col),color=NA) +
	geom_text(data=cat_df,aes(x=x+.03,y=y,label=cat),hjust=0,size=1.90,color="#4B5563") +
	annotate("text",.5,.425,label="Q1–Q5 log(HR)",size=2.15,fontface="bold",color="#374151") +
	geom_tile(data=b_hr,aes(x,y,fill=fill),width=b_hr$w,height=.021) +
	annotate("text",c(.25,.50,.75),.323,label=c(sprintf("%.2f",-mx_q),"0",sprintf("%.2f",mx_q)),size=1.70,color="#6B7280") +
	annotate("text",.5,.285,label="RCS (−log10 nonlinear P)",size=2.05,fontface="bold",color="#374151") +
	geom_tile(data=b_nl,aes(x,y,fill=fill),width=b_nl$w,height=.021) +
	annotate("text",c(.25,.75),.183,label=c("0",sprintf("%.1f",max(2,mx_nl))),size=1.70,color="#6B7280") +
	annotate("text",.5,.145,label="Average consumption level\n(servings/day)",size=1.95,fontface="bold",color="#374151") +
	geom_tile(data=b_mn,aes(x,y,fill=fill),width=b_mn$w,height=.021) +
	annotate("text",c(.25,.75),.043,label=c(sprintf("%.2f",mean_min),sprintf("%.2f",mean_max)),size=1.70,color="#6B7280") +
	scale_fill_identity()
p_f <- ggdraw() + draw_plot(p_main,0,0,1,1) + draw_plot(p_leg,.402,.402,.196,.196)
ggsave("Fig2.pdf", p_f, width=14, height=14)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
# 🚩 FigS1 食物组组间相关性
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pacman::p_load(dplyr, pheatmap, tibble)
cor.mat <- cor(select(dat, all_of(unique(field$Food.Group))), method = "spearman", use = "pairwise.complete.obs")
hc <- hclust(as.dist(1 - cor.mat), method = "average")
pheatmap::pheatmap(cor.mat, cluster_rows = hc, cluster_cols = hc, color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(200), breaks = seq(-1, 1, length.out = 201), border_color = NA, fontsize_row = 8, fontsize_col = 8, angle_col = 45, main = "Spearman correlation and hierarchical clustering of food groups", filename = "group_corr_cluster_heatmap.pdf", width = 13, height = 12)
cor.long <- as.data.frame(as.table(cor.mat)) %>% rename(group1 = Var1, group2 = Var2, rho = Freq) %>% filter(group1 != group2) %>% rowwise() %>% mutate(pair = paste(sort(c(group1, group2)), collapse = " | ")) %>% ungroup() %>% distinct(pair, .keep_all = TRUE) %>% arrange(desc(abs(rho)))
write.csv(cor.long, "group_corr_long.csv", row.names = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
# 🚩 Fig S3
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pacman::p_load(dplyr, ggplot2, ggrepel, stringr)
dat <- read.csv(paste0("Food_", Y, "_.csv")) %>% filter(class == 2) %>% mutate(logP = -log10(P), HR_type = ifelse(HR < 1, "HR < 1", "HR > 1")) %>% left_join(select(field, Category, Food.Group, Food.Item, exposure = name), by = "exposure") %>%
    mutate(Category=case_when(
				Category %in% c("Meats","Fish_fish_dishes") ~ "Meats & Fish", 
				Category %in% c("Non_alcoholic_beverages","Alcoholic_beverages") ~ "Water & beverage",
				Category=="Fried_and_fast_foods" ~ "Fried & Fast food", Category=="Legumes_Nuts" ~ "Legumes & Nuts", Category=="Snack_Pastry" ~ "Snack & Pastry", Category=="Cereal_Bread" ~ "Cereal & Bread", Category=="Fats_Oil" ~ "Fats & Oil", Category=="Dairy_dairy_free_products" ~ "Dairy", TRUE ~ str_replace_all(as.character(Category), "_", " ")
			), Category=replace_na(Category, "Unknown")
)
fdr_p <- if(any(dat$FDR < 0.05, na.rm = TRUE)) max(dat$P[dat$FDR < 0.05], na.rm = TRUE) else NULL
p_lines <- c(0.05, 0.01, 0.001)
grp <- dat %>% group_by(Category, Food.Group) %>% summarise(mP = min(P, na.rm=T), sig = sum(P < 0.05, na.rm=T), .groups="drop") %>% arrange(Category, mP, desc(sig)) %>% mutate(x = row_number())
ax <- grp %>% group_by(Category) %>% summarise(cx = mean(x), bnd = max(x) + 0.5, .groups="drop")
set.seed(12345)
df <- dat %>% left_join(select(grp, Food.Group, x), by="Food.Group") %>% group_by(Food.Group) %>% mutate(xj = x + runif(n(), -0.3, 0.3)) %>% ungroup()
p <- ggplot(df, aes(x = xj, y = logP)) + geom_vline(xintercept = head(ax$bnd, -1), linetype = "dashed", color = "grey85", linewidth = 0.5) + geom_hline(yintercept = -log10(p_lines), linetype = "dotted", color = "grey70", linewidth = 0.4) +
    annotate("text", x = 0.5, y = -log10(p_lines) + 0.05, label = paste0("P==", p_lines), parse = TRUE, size = 2.2, color = "grey50", hjust = 0) + 
    {if(!is.null(fdr_p)) geom_hline(yintercept = -log10(fdr_p), linetype = "dashed", color = "#BE665E", linewidth = 0.6)} +
    {if(!is.null(fdr_p)) annotate("text", x = 0.5, y = -log10(fdr_p) + 0.1, label = "FDR < 0.05", size = 2.8, color = "#BE665E", hjust = 0, fontface="bold")} +
    geom_point(aes(color = Category, shape = HR_type, size = total), alpha = 0.7) +
    geom_text_repel(data = filter(df, P < 0.05), aes(label = Food.Item), size = 2.8, max.overlaps = 30, box.padding = 0.3, segment.color = "grey50", show.legend = FALSE, color = "black") +
    scale_shape_manual(name = "Effect Direction", values = c("HR < 1" = 16, "HR > 1" = 17)) +
    scale_size_continuous(name = "Sample size", range = c(2, 6)) +
    scale_x_continuous(breaks = grp$x, labels = str_replace_all(grp$Food.Group, "_", " "), expand = expansion(mult = c(0.01, 0.01))) +
    labs(x = "Food Groups (ordered by Category)", y = expression(-log[10](italic(P)~value)), title = paste0("Manhattan Plot: Associations between individual foods and ", Y)) +  
    theme_bw(base_size = 12) + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 60, hjust = 1, size = 7), legend.position = "right", strip.background = element_blank())
ggsave(paste0("FigS3_", Y, ".pdf"), p, width = 18, height = 8)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 ML 数据准备
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
group_res <- read.csv(paste0("Group_", Y, "_.csv"))
group_sig_anyref <- group_res %>% mutate(class = suppressWarnings(as.integer(as.character(class)))) %>% filter(!is.na(class)) %>% group_by(exposure) %>%
	summarise(p_min = min(P, na.rm = TRUE), best_ref = Ref[which.min(P)], best_class = class[which.min(P)], n_sig = sum(P < 0.05, na.rm = TRUE), .groups = "drop") %>%
	filter(n_sig >= 1)  # Any Ref 筛显著组:不固定 Q1 为参考，只要任意 Ref 下任意一个 level 显著，就纳入候选
	group_sig <- group_sig_anyref$exposure; cat("候选显著组数量 =", length(group_sig), "\n"); print(group_sig)
	score.tbl <- group_sig_anyref %>% transmute(exposure, score = p_min)
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
	get_residual <- \(y, cov_df) as.numeric(residuals(lm(y ~ ., data = data.frame(y = y, cov_df), na.action = na.exclude)))
	for (v in feat_use) dat_res[[paste0(v, "_res")]] <- get_residual(dat_res[[v]], dat_res[, covs.food, drop = FALSE])
	group_residual_vars <- paste0(feat_use, "_res")
dat_ml <- dat_res %>% left_join(dat %>% select(eid, all_of(paste0(Y, c(".t2e", ".Yt2e"))), region_code), by = "eid") %>% select(eid, all_of(group_residual_vars), all_of(paste0(Y, c(".t2e", ".Yt2e"))), region_code)
write.csv(dat_ml, "dat_ml.csv", row.names = FALSE)
saveRDS(dat_ml, "dat_ml.rds")


