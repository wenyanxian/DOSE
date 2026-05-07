pacman::p_load(readxl, dplyr, tidyr, data.table, tidyverse, ggrepel, lubridate, purrr, survival, ggplot2, forcats, stringr)

dir0 <- ifelse(Sys.info()[["sysname"]] == "Windows", "D:/", "/work/sph-huangj")
indir = paste0(dir0, "/data/ukb/phe")
invisible(lapply(c("phe.f.R", "assoc.f.R", "plot.f.R"), function(f) source(file.path(dir0, "scripts", "f", f))))

setwd("D:/analysis/stroke")

## ------------------------------------------------------------------
## DOSE score
## ------------------------------------------------------------------
dat <- readRDS("dat.grp.rds")
cols  <- c("Yogurt","Allium_vegetables","Snacks_and_pastries","Legumes","Citrus","Apples_Pears","Dried_fruit")
dat.ml <- dat %>% mutate( # 🚩 DOSE-5
#	Citrus.dose = ifelse(Citrus > 0, 1, 0),
#	Apples_Pears.dose = ifelse(Apples_Pears > 0 & Apples_Pears <= 2, 1, 0),
#	Dried_fruit.dose = ifelse(Dried_fruit > 0, 1, 0),	
    Yogurt.dose = ifelse(Yogurt > 0 & Yogurt <= 1, 1, 0),
    Allium_vegetables.dose = ifelse(Allium_vegetables > 0 & Allium_vegetables <= 0.5, 1, 0),
    Snacks_and_pastries.dose = ifelse(Snacks_and_pastries > 2.25 & Snacks_and_pastries <= 5.31, 1, 0),
    Legumes.dose = ifelse(Legumes > 0, 1, 0),
	Citrus_ApplesPears_Driedfruit = Citrus + Apples_Pears + Dried_fruit,
    Fruit.dose = ifelse(Citrus_ApplesPears_Driedfruit > 0, 1, 0),
    diet.dose5.sum = rowSums1(across(c(Yogurt.dose, Allium_vegetables.dose, Snacks_and_pastries.dose, Legumes.dose, Fruit.dose)))
) %>% select(eid, diet.dose5.sum, ends_with(".dose"))  
summary(dat.ml$diet.dose5.sum)
saveRDS(dat.ml, "dose_score.rds")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Fig3c
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C <- list(hi="#2F5D8A", lo="#D9E1EA"); t_sz <- 11
thm_dose <- theme_bw() + theme(
	panel.border=element_blank(), axis.line=element_blank(), axis.ticks=element_blank(), panel.grid=element_blank(), legend.position="bottom", 
	axis.text.y=element_text(size=t_sz, face="bold", color="black"), axis.text.x=element_text(size=t_sz, face="bold", color="grey40"), axis.title.x=element_text(size=t_sz+1, face="bold", margin=margin(t=10)), 
	plot.margin=margin(10,20,10,20)
)
df_p <- tibble(
	f = factor(c("Fruit","Legumes","Yogurt","Allium vegetables","Snacks and pastries"), levels = rev(c("Fruit","Legumes","Yogurt","Allium vegetables","Snacks and pastries"))),
	l_l = c("0", "0", "0", "0", "0–2.25"), m_l = c(">0", ">0", ">0–1", ">0–0.5", ">2.25–5.31"), h_l = c("", "", ">1", ">0.5", ">5.31"),
	l_s = "0", m_s = "1", h_s = c("1","1","0","0","0")
)%>% pivot_longer(-f, names_to=c("lvl", ".value"), names_sep="_") %>% mutate(lvl=factor(lvl, levels=c("l","m","h")), x=as.numeric(lvl), y=as.numeric(f), x1=case_when(f %in% c("Fruit","Legumes") & lvl=="m" ~ 1.5, TRUE ~ x-0.5), x2=case_when(f %in% c("Fruit","Legumes") & lvl=="m" ~ 3.5, TRUE ~ x+0.5))
p <- ggplot(df_p, aes(x, f)) +
  geom_rect(aes(xmin=0.5, xmax=3.5, ymin=y-0.45, ymax=y+0.45), fill="#F2F2F2", inherit.aes=F) +
  geom_rect(data=filter(df_p, s=="1"), aes(xmin=x1, xmax=x2, ymin=y-0.45, ymax=y+0.45), fill=C$hi) +
  geom_segment(aes(x=1.5, xend=1.5, y=y-0.45, yend=y+0.45), color="white", linewidth=1.5) +
  geom_segment(data=filter(df_p, !f %in% c("Fruit","Legumes")), aes(x=2.5, xend=2.5, y=y-0.45, yend=y+0.45), color="white", linewidth=1.5) +
  geom_text(aes(label=l, color=s), size=3.8, fontface="bold") +
  scale_color_manual(values=c("0"="grey40", "1"="white"), guide="none") +
  scale_x_continuous(breaks=1:3, labels=c("Low","Medium","High"), position="bottom", expand=c(0,0)) +
  labs(x="Consumption level", y=NULL) + thm_dose
ggsave("Fig3c.pdf", p, width=5, height=3.5, device=cairo_pdf)


## ------------------------------------------------------------------
## 🚩 Pattern
## ------------------------------------------------------------------
std_10_90 <- function(x){x <- as.numeric(x); q <- quantile(x, c(0.10, 0.90), na.rm = TRUE, names = FALSE); if(!all(is.finite(q)) || q[2] == q[1]) rep(NA_real_, length(x)) else (x - q[1]) / (q[2] - q[1])}
score_q5 <- function(x){x <- as.numeric(x); out <- rep(NA_real_, length(x)); ok <- is.finite(x); out[ok] <- c(0, 25, 50, 75, 100)[dplyr::ntile(x[ok], 5)]; out}

dat1 <- readRDS("dat1.rds")
dat.imp <- readRDS("dat.imp.rds")
dat.ml <- readRDS("dose_score.rds")

dxs <- c("primary_ahtn", "t2dm", "cvd")
drugs <- c("drug.lipid", "drug.htn", "drug.dm", "drug.aspirin")
covs.base <- c("age", "sex.b", paste0("PC", 1:10))
covs.add <- c("ethnic.c3", "edu.c", "tdi", "job.b", "living_alone.b", "enroll", "energy.kJ", "smoke.c", "days_pa_mod.c", "sleep_duration.b", "bmi") # strata(region)
covs.all <- unique(c(covs.base, covs.add, dxs, drugs, "fhist_stroke", "alcohol.c")) 
diet_with_alcohol <- c("medi24", "mind", "ahei", "redii", "fds", "upf"); diet_without_alcohol <- c("dash", "hpdi", "phdi", "digm", "hlcd", "hlfd", "modern", "dose5") 
diet_lst <- c(diet_with_alcohol, diet_without_alcohol)
get_model_list <- function(diet) { #调整酒精
	alc_in_score <- diet %in% diet_with_alcohol
	m1 <- c(covs.base, covs.add, if (!alc_in_score) "alcohol.c", "fhist_stroke"); m2 <- c(m1, dxs); m3 <- c(m2, drugs) 
	list(m1 = unique(m1), m2 = unique(m2), m3 = unique(m3))
}
get_covs <- function(diet, model = c("m1", "m2", "m3")) {model <- match.arg(model); get_model_list(diet)[[model]]}

dat <- dat1 %>% select(-all_of(covs.all)) %>% left_join(dat.imp %>% select(eid, all_of(covs.all)), by = "eid") %>% left_join(dat.ml, by = "eid")
sum_cols <- names(dat)[grepl("^diet\\..+\\.sum$", names(dat))]
dat <- dat %>% mutate(
	across(all_of(sum_cols), std_10_90, .names = "{sub('\\\\.sum$', '', .col)}.std"),
	across(all_of(sum_cols), score_q5, .names = "{sub('\\\\.sum$', '', .col)}.q5")
)


score_q4_cut <- function(x){
	x <- as.numeric(x); q <- unique(quantile(x, c(0, .25, .5, .75, 1), na.rm = TRUE, names = FALSE))
	if(length(q) < 3) return(factor(rep(NA_character_, length(x)))); cut(x, breaks = q, include.lowest = TRUE, right = TRUE)
}
dat <- dat %>% mutate(
	diet.dose5.4g = case_when(diet.dose5.sum <= 1 ~ "0-1", diet.dose5.sum == 2 ~ "2", diet.dose5.sum == 3 ~ "3", diet.dose5.sum >= 4 ~ "4-5", TRUE ~ NA_character_),
	diet.dose5.4g = factor(diet.dose5.4g, levels = c("0-1","2","3","4-5")), diet.dash.4g = score_q4_cut(diet.dash.sum)
)
table(dat$diet.dose5.4g, useNA = "always")
table(dat$diet.dash.4g, useNA = "always")
table(dat$diet.dose5.sum, dat$diet.dose5.4g, useNA = "always")
table(dat$diet.dash.sum, dat$diet.dash.4g, useNA = "always")

covs.food <- c("age","sex.b",paste0("PC",1:10),"ethnic.c3","edu.c","tdi","job.b","living_alone.b","enroll","energy.kJ","smoke.c","days_pa_mod.c","sleep_duration.b","bmi","fhist_stroke","primary_ahtn","t2dm","cvd")
covs.dose <- c(covs.food, "alcohol.c")
Ys <- c("qi_stroke", "cvd_stroke_i", "cvd_stroke_ih", "cvd_stroke_sh")
Y <- "cvd_stroke_i" # "qi_stroke"
dose_lst <- c("dose5") # dose7

tab_event_incident <- function(dat, Ys, exp_var = NULL, labels = NULL) { 
	out <- bind_rows(lapply(Ys, function(Y, y=paste0(Y,".Yt2e"), t=paste0(Y,".t2e"), p=paste0("prior_",Y)) {
		if (!all(c(y, t) %in% names(dat))) return(NULL) 
		d <- dat %>% filter(if(p %in% names(.)) !coalesce(.data[[p]], FALSE) else TRUE, !is.na(.data[[y]]))
		f <- function(df, e, l) df %>% summarise(total=n(), case=sum(.data[[y]]==1, na.rm=T), person_years=round(sum(.data[[t]], na.rm=T),1), .groups="drop") %>% mutate(outcome=Y, exposure=e, level=l, case_pct=round(case/total*100,2), ir_1e5=round(case/person_years*1e5,1))
		bind_rows(f(d, "Overall", "Overall"), if (!is.null(exp_var) && exp_var %in% names(d)) d %>% filter(!is.na(.data[[exp_var]])) %>% group_by(v=.data[[exp_var]]) %>% f(exp_var,"") %>% arrange(v) %>% mutate(level=paste0("Q",row_number())) %>% select(-v))
	}))
	if (!is.null(labels)) mutate(out, label=labels[outcome], .after=outcome) else out
}
tab_event <- bind_rows(lapply(diet_lst, function(d) {tab_event_incident(dat = dat, Ys = Ys, exp_var = paste0("diet.", d, ".4g"))})) %>% distinct() 
write.csv(tab_event, "tab_event_4g.csv", row.names = FALSE)	

# cox
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

res_std <- run_models_diet(dat = dat, Ys = Ys, diets = diet_lst, FUN = run_cox, model_fun = get_model_list, strata_var = "region", exposure_suffix = "\\.std$")
	res_std %>% count(outcome, model)
write.csv(res_std, "res_std.csv", row.names = FALSE)
res_4g_dose <- run_models_diet(dat, Ys, "dose5", run_cox_cat, get_model_list, "region", "\\.4g$", ref_level = "0-1", ord_levels = c("0-1","2","3","4-5")) 
res_4g_dash <- run_models_diet(dat, Ys, "dash", run_cox_cat, get_model_list, "region", "\\.4g$", ref_level = levels(dat$diet.dash.4g)[1], ord_levels = levels(dat$diet.dash.4g))
res_4g <- bind_rows(res_4g_dash, res_4g_dose) %>%
	group_by(outcome, exposure, model) %>% mutate(lvl = ifelse(type == "Cat", paste0("Q", cumsum(type == "Cat") + 1), NA)) %>% ungroup() %>%
	left_join(select(filter(tab_event, level != "Overall"), outcome, exposure, level, t=total, c=case), by=c("outcome", "exposure", "lvl"="level")) %>%
	mutate(total = coalesce(t, total), case = coalesce(c, case)) %>% select(-lvl, -t, -c)
write.csv(res_4g, "res_4g_dash_dose5.csv", row.names = FALSE)


# RCS
res_rcs <- run_rcs_one(dat, paste0(Y, ".t2e"), paste0(Y, ".Yt2e"), paste0("diet.", dose_lst[1], ".sum"), covs.dose)
ggsave(paste0("RCS_", dose_lst[1], "_", Y, ".pdf"), res_rcs$plot, width = 5.2, height = 4.2); res_rcs$st; res_rcs$plot


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 表1
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pacman::p_load(gtsummary, flextable, tidyverse)
v_cont <- c("age", "tdi", "enroll", "energy.kJ", "bmi")
v_cat <- c("sex.b", "ethnic.c3", "edu.c", "job.b", "living_alone.b", "smoke.c", "alcohol.c", "days_pa_mod.c", "sleep_duration.b", "fhist_stroke", "primary_ahtn", "t2dm", "cvd", "drug.lipid", "drug.htn", "drug.dm", "drug.aspirin")
dat %>% select(Group = diet.dose5.4g, `DOSE mean` = diet.dose5.sum, all_of(v_cont), all_of(v_cat)) %>% filter(!is.na(Group)) %>%
  tbl_summary(by = Group, type = list(`DOSE mean` ~ "continuous"), statistic = list(all_continuous() ~ "{mean} ± {sd}", all_categorical() ~ "{n} ({p}%)"), digits = list(all_continuous() ~ 1, all_categorical() ~ c(0, 1)), missing = "no",
    label = list(
      age ~ "Age (years)", sex.b ~ "Sex", tdi ~ "Townsend deprivation index", enroll ~ "Number of dietary assessments", energy.kJ ~ "Energy intake (kJ)", bmi ~ "BMI (kg/m²)", ethnic.c3 ~ "Ethnicity", edu.c ~ "Education level", job.b ~ "Employment status", living_alone.b ~ "Living alone", smoke.c ~ "Smoking status", alcohol.c ~ "Alcohol status", days_pa_mod.c ~ "Physical activity (days/week)", sleep_duration.b ~ "Healthy sleep duration",
      fhist_stroke ~ "Family history of stroke",
      primary_ahtn ~ "History of Hypertension", t2dm ~ "History of T2DM", cvd ~ "History of CVD",
      drug.lipid ~ "Lipid-lowering drugs", drug.htn ~ "Antihypertensive drugs", drug.dm ~ "Antidiabetic drugs", drug.aspirin ~ "Aspirin use"
)) %>% add_overall(last = FALSE, col_label = "**Overall**") %>% add_p(test = list(all_continuous() ~ "kruskal.test", all_categorical() ~ "chisq.test"), pvalue_fun = ~ style_pvalue(.x, digits = 3)
) %>% bold_labels() %>% as_flex_table() %>% flextable::save_as_docx(path = "Table1_Dose5_4g.docx")

pacman::p_load(dplyr, tidyr, stringr, purrr, forestploter, grid, ggplot2, writexl)
dir0 <- ifelse(Sys.info()[["sysname"]] == "Windows", "D:/", "/work/sph-huangj"); indir <- paste0(dir0, "/analysis/stroke"); outdir <- paste0(indir, "/result"); if(!dir.exists(outdir)) dir.create(outdir, recursive = TRUE); setwd(indir)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Fig3d
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
out_map <- c(qi_stroke="Overall stroke", cvd_stroke_i="Ischemic stroke", cvd_stroke_ih="Intracerebral hemorrhage", cvd_stroke_sh="Subarachnoid hemorrhage")
diet_map <- c(medi24="MedDiet", dash="DASH", mind="MIND", hpdi="hPDI", ahei="AHEI", phdi="PHDI", digm="DI-GM", modern="Modern", upf="UPF", hlcd="HLCD", hlfd="HLFD", redii="rEDII", fds="FDS", dose5="DOSE")
keep_diets <- c("DOSE", "DASH", "MIND", "AHEI", "MedDiet", "hPDI", "PHDI", "DI-GM", "HLFD", "HLCD")
C_S <- list(dose="#2F5D8A", sig="#4A7EBB", nonsig="#BFC7D1")
fp <- function(p, f) paste0(ifelse(p<0.001, "<0.001", sprintf("%.3f", p)), ifelse(f<0.05, "†", ""))
res_all <- read.csv("res_std.csv") %>% mutate(diet=recode(str_remove_all(exposure, "^diet\\.|\\.(std|q5)$"), !!!diet_map), out_lab=recode(outcome, !!!out_map), model_lab=recode(model, m1="Model 1", m2="Model 2", m3="Model 3")) %>% filter(outcome %in% c("qi_stroke", "cvd_stroke_i"), diet %in% keep_diets, model %in% c("m1","m2","m3")) %>% group_by(outcome, model) %>% mutate(FDR=p.adjust(P, "fdr"), `HR (95% CI)`=sprintf("%.2f (%.2f-%.2f)", HR, LCI, HCI), `P value`=fmt_p(P), FDR_fmt=fmt_p(FDR)) %>% ungroup()
	diet_lvls <- res_all %>% filter(outcome=="cvd_stroke_i", model=="m2") %>% arrange(HR) %>% pull(diet)
	tab_s10 <- res_all %>% mutate(diet=factor(diet, diet_lvls)) %>% arrange(outcome, diet, model) %>% select(out_lab, diet, total, case, model_lab, `HR (95% CI)`, `P value`, FDR_fmt) %>% pivot_wider(names_from=model_lab, values_from=c(`HR (95% CI)`, `P value`, FDR_fmt), names_glue="{model_lab} {.value}") %>% transmute(Outcome=out_lab, `Diet score`=as.character(diet), Total=total, Event=case, `Model 1 HR (95% CI)`, `Model 1 P value`, `Model 1 FDR`=`Model 1 FDR_fmt`, `Model 2 HR (95% CI)`, `Model 2 P value`, `Model 2 FDR`=`Model 2 FDR_fmt`, `Model 3 HR (95% CI)`, `Model 3 P value`, `Model 3 FDR`=`Model 3 FDR_fmt`)
	write.csv(tab_s10, paste0(outdir, "/Table_S10_DietScores_Stroke_IS.csv"), row.names=FALSE)	
res <- res_all %>% filter(model=="m2") %>% group_by(outcome) %>% mutate(est=sprintf("%.2f (%.2f-%.2f)", HR, LCI, HCI), p_lab=fp(P, FDR), p_col=case_when(diet=="DOSE"~C_S["dose"], P<0.05~C_S["sig"], TRUE~C_S["nonsig"]), evt_str=paste0(format(case, big.mark=","), " / ", format(total, big.mark=","))) %>% ungroup() %>% mutate(diet=factor(diet, diet_lvls)) %>% arrange(outcome, diet)
	dp <- res %>% group_split(outcome) %>% map_dfr(\(x) bind_rows(tibble(diet=paste0(x$out_lab[1], " (", x$evt_str[1], ")"), is_sum=TRUE, HR=NA, LCI=NA, HCI=NA, est=NA, p_lab=NA, p_col=NA), x %>% mutate(diet=paste0("   ", diet), is_sum=FALSE)))
p <- forest(dp %>% transmute(`Dietary Pattern`=diet, ` `=strrep(" ", 38), `Hazard ratio (95% CI)`=est, `P value`=p_lab) %>% 
	mutate(across(everything(), ~replace_na(.x, ""))), est=dp$HR, lower=dp$LCI, upper=dp$HCI, ci_column=2, ref_line=1, xlim=c(0.6, 1.3), ticks_at=seq(0.6, 1.3, 0.1), theme=forest_theme(base_size=7.5, ci_lwd=1.8, ci_Theight=0.08, refline_gp=gpar(lwd=0.8, lty="dashed", col="grey30"), vert_line=c(0.7,0.8,0.9,1.1,1.2), vert_line_gp=list(gpar(lty="dotted", col="grey90")), core=list(bg_params=list(fill="white"), padding=unit(c(0.1,4), "mm")), header=list(fg_params=list(fontface="bold", size=8), bg_params=list(fill="white")))) %>%
	add_border(part="header", where="bottom") %>% edit_plot(row=which(dp$is_sum), gp=gpar(fontface="bold")) %>% edit_plot(row=0, col=4, which="text", gp=gpar(fontface="bold.italic"))
	for(i in which(!dp$is_sum)) p <- p %>% edit_plot(row=i, col=2, which="ci", gp=gpar(col=dp$p_col[i], fill=dp$p_col[i]))
	pdf(paste0(outdir, "/Fig_Stroke_Forest_m2.pdf"), 7.0, 6.0)
	grid.draw(p); grid.text("† FDR < 0.05", .5, .02, gp=gpar(fontsize=7, fontface="italic"))
	dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Fig3e RCS + 折线图
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pacman::p_load(rms, ggplot2, dplyr, patchwork, stringr)
res_4g <- read.csv("res_4g_dash_dose5.csv")
extract_rcs_full <- \(dat, t, e, x, covs, strat = "region", outcome_label) {
    d <- na.omit(dat[c(t, e, x, covs, strat)])
    if(nrow(d) < 100 || n_distinct(d[[x]]) < 5 || is.null(k <- get_k3(d[[x]]))) return(NULL)    
    dd <<- datadist(d); options(datadist = "dd")
    f <- as.formula(paste0("Surv(", t, ",", e, ")~rcs(", x, ",c(", paste(k, collapse=","), "))+", paste(c(covs, paste0("strat(", strat, ")")), collapse="+")))
    fit <- cph(f, data = d, x = T, y = T)
    an <- anova(fit); p_all <- an[x, "P"]; p_non <- an[grep("Nonlinear", rownames(an))[1], "P"]; xm <- quantile(d[[x]], .995)
    pr <- Predict(fit, name = x, ref.zero = T, fun = exp, np = 200) %>% as.data.frame() %>% filter(!!sym(x) <= xm) %>% mutate(outcome = outcome_label, p_overall = p_all, p_nonlinear = p_non)
    return(pr)
}
make_final_plot <- function(outcome_display_name, plot_title, rcs_data) {
    trend_data <- df_trend %>% filter(outcome == outcome_display_name)
    raw_name <- ifelse(outcome_display_name == "Ischemic stroke", "cvd_stroke_i", "qi_stroke")    
    p_trend_val <- res_4g %>% filter(outcome == raw_name, model == "m2", type == "Ptrend") %>% pull(P) %>% .[1]   
    p_trend_lab <- if(is.na(p_trend_val)) "" else if(p_trend_val < 0.001) "P-trend < 0.001" else paste0("P-trend = ", sprintf("%.3f", p_trend_val))
    p_main <- ggplot(trend_data, aes(grp, HR)) + geom_hline(yintercept = 1, linetype = "dashed", color = "grey55", linewidth = .55) + geom_line(aes(group = 1), color = "#2F5D8A", linewidth = 1) + geom_errorbar(aes(ymin = LCI, ymax = HCI), color = "#2F5D8A", width = .12, linewidth = .8) + geom_point(color = "#2F5D8A", size = 3) + geom_text(aes(label = ifelse(grp == "0–1", "Ref.", sprintf("%.2f", HR))),nudge_x = 0.22, vjust = -0.8, size = 3.2, fontface = "bold", hjust = 0) + annotate("text", x = 4.4, y = 1.05, label = p_trend_lab, hjust = 1, size = 3.8, fontface = "bold", color = "black") + scale_y_continuous("HR (95% CI)", limits = c(0.3, 1.1), breaks = seq(0.6, 1.2, 0.2)) + labs(x = "DOSE score", title = plot_title) + theme_bw(base_size = 11) + theme(panel.grid = element_blank(), plot.title = element_text(size = 13, face = "bold"), axis.title = element_text(face = "bold"), axis.text = element_text(color = "black"))
    rcs_p_txt <- paste0("P-overall ", fmt_p(rcs_data$p_overall[1]), "\nP-non-linear ", fmt_p(rcs_data$p_nonlinear[1]))    
    p_inset <- ggplot(rcs_data, aes(x = .data[[x_var]], y = yhat)) + geom_hline(yintercept = 1, linetype = "dashed", color = "grey55", linewidth = 0.4) + geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#2F5D8A", alpha = 0.15) + geom_line(color = "#2F5D8A", linewidth = 0.8) + annotate("text", x = -Inf, y = Inf, label = rcs_p_txt, hjust = -0.1, vjust = 1.2, size = 2.4, fontface = "bold") + labs(x = "Score", y = "HR") + theme_bw(base_size = 7) + theme(panel.grid = element_blank(), axis.title = element_text(size = 7, face = "bold"), plot.background = element_rect(fill = "white", color = "grey20", linewidth = 0.4))
    p_main + inset_element(p_inset, left = 0.01, bottom = 0.03, right = 0.38, top = 0.45, align_to = "panel")
}
pr_isc_full <- extract_rcs_full(dat, "cvd_stroke_i.t2e", "cvd_stroke_i.Yt2e", x_var, covs.dose, "region", "Ischemic stroke")
pr_all_full <- extract_rcs_full(dat, "qi_stroke.t2e", "qi_stroke.Yt2e", x_var, covs.dose, "region", "Overall stroke")
p_final <- make_final_plot("Ischemic stroke", "Ischemic Stroke", pr_isc_full) / 
           make_final_plot("Overall stroke", "Overall Stroke", pr_all_full)
ggsave("Figure_3e.pdf", p_final, width = 6, height = 10.5, device = cairo_pdf)


## ------------------------------------------------------------------
##  🚩 other health-related outcomes
## ------------------------------------------------------------------
## 先提取 icd10_other.lst
code0 <- readRDS(paste0(indir, "/Rdata/icd10.code0.rds"))
date0 <- readRDS(paste0(indir, "/Rdata/icd10.date0.rds"))
fn <- fread(paste0(indir, "/common/icd10_other.lst"), sep = "\t", header = FALSE, fill = TRUE) %>% as.data.frame()
colnames(fn)[1:3] <- c("pattern", "name", "label")
dat.other <- data.frame(eid = code0$eid, stringsAsFactors = FALSE)
for(i in 1:nrow(fn)) {print(i); dat.other <- get_icd("icd10", code0, date0, dat.other, fn$pattern[i], fn$name[i])}
saveRDS(dat.other, paste0(indir, "/Rdata/fod.icd10.other.rds"))

primary_stroke <- c("qi_stroke", "cvd_stroke_i", "cvd_stroke_ih", "cvd_stroke_sh")
fn_other <- fread(paste0(indir, "/common/icd10_other.lst"), sep="\t", header=FALSE, fill=TRUE) %>% as.data.frame()
colnames(fn_other)[1:3] <- c("pattern","name","label")
fn_other <- fn_other %>% filter(name != "stroke")  # remove duplicate I6[01]|I6[34]
fn_stroke <- tibble(pattern = NA_character_, name = primary_stroke, label = c("Overall stroke", "Ischemic stroke", "Intracerebral hemorrhage", "Subarachnoid hemorrhage"))
fn_other2 <- bind_rows(fn_stroke, fn_other) %>% distinct(name, .keep_all=TRUE)
Ys_other2 <- fn_other2$name
write.csv(fn_other2, "outcome_wide_label_list.csv", row.names=FALSE)

fod.other <- readRDS(paste0(indir, "/Rdata/fod.icd10.other.rds")) %>% mutate(eid = as.character(eid))
dat_ow <- dat %>% dplyr::select(eid, all_of(covs.all), region, names(dat)[grepl("^diet\\.dose", names(dat))], birth_date, start_date, date_lost, date_death, any_of(paste0("fod_icd10_", primary_stroke))) %>% left_join(fod.other, by="eid")
for(Y in Ys_other2){
	dat_ow <- dat_ow %>% select(-any_of(c(paste0(Y,".Yt2e"), paste0(Y,".t2e"), paste0("prior_",Y))))
	dat_ow <- t2e(dat_ow, NA, paste0("fod_icd10_",Y), "birth_date", "start_date", "date_lost", "date_death", date_follow_end, Y, "year")
}
dat_ow <- dat_ow %>% mutate(across(ends_with("Yt2e"), as.integer))
for(Y in Ys_other2){
	fv <- paste0("fod_icd10_", Y)
	dat_ow[[paste0("prior_",Y)]] <- if(fv %in% names(dat_ow)){!is.na(dat_ow[[fv]]) & !is.na(dat_ow$start_date) & dat_ow[[fv]] <= dat_ow$start_date
	} else FALSE
}

tab_event_incident_by_exp <- function(dat, Ys, exp_var, labels=NULL){ #  事件数统计（按每个 outcome 排除基线既往病）
	out <- bind_rows(lapply(Ys, \(Y){
		y <- paste0(Y,".Yt2e"); t <- paste0(Y,".t2e"); p <- paste0("prior_",Y)
		if(!all(c(y,t,p,exp_var) %in% names(dat))) return(NULL)
		d <- dat %>% filter(!coalesce(.data[[p]], FALSE), !is.na(.data[[y]]), !is.na(.data[[exp_var]]))
		bind_rows(
			d %>% summarise(total=n(), case=sum(.data[[y]]==1, na.rm=TRUE), .groups="drop") %>% mutate(outcome=Y, exposure="Overall", level="Overall"),
			d %>% group_by(level=.data[[exp_var]]) %>% summarise(total=n(), case=sum(.data[[y]]==1, na.rm=TRUE), .groups="drop") %>% mutate(outcome=Y, exposure=exp_var)
		)
	}))
	if(!is.null(labels)) out <- out %>% left_join(labels %>% select(name,label), by=c("outcome"="name"))
	out
}
tab_event_other_4g <- tab_event_incident_by_exp(dat_ow, Ys_other2, "diet.dose5.4g", fn_other2)
write.csv(tab_event_other_4g, "tab_event_other_4g_include_primary_stroke.csv", row.names=FALSE)

# Outcome-wide association
covs.Outcomewide <- c("age","sex.b",paste0("PC",1:10),"ethnic.c3","edu.c","tdi","job.b","living_alone.b","enroll","energy.kJ","smoke.c","days_pa_mod.c","sleep_duration.b","bmi","alcohol.c")
dose_lst <- "dose5"
run_models_dose_other <- function(dat, Ys, diet="dose5", FUN, covs=covs.Outcomewide, strata_var="region", exposure_suffix="\\.std$", ...){
	bind_rows(lapply(Ys, \(Y){
		p <- paste0("prior_", Y)
		dY <- if(p %in% names(dat)) dat %>% filter(!coalesce(.data[[p]], FALSE)) else dat
		x <- grep(paste0("^diet\\.", diet, exposure_suffix), names(dY), value=TRUE)
		if(length(x) != 1) stop("Cannot uniquely find exposure for diet: ", diet)
		FUN(dat=dY, Y=Y, X=x, covs=covs, strata_var=strata_var, ...) %>% mutate(diet=diet, model="outcome_wide")
	}))
}
res_std_other <- run_models_dose_other(dat_ow, Ys_other2, dose_lst, run_cox, covs.Outcomewide, "region", "\\.std$") %>% left_join(fn_other2 %>% select(name,label), by=c("outcome"="name")) %>% group_by(exposure) %>% mutate(FDR=p.adjust(P, "fdr")) %>% ungroup()
write.csv(res_std_other, "res_std_other_include_primary_stroke.csv", row.names=FALSE)
res_4g_other <- run_models_dose_other(dat_ow, Ys_other2, dose_lst, run_cox_cat, covs.Outcomewide, "region", "\\.4g$", ref_level="0-1", ord_levels=c("0-1","2","3","4-5")) %>% left_join(fn_other2 %>% select(name,label), by=c("outcome"="name")) %>% group_by(exposure, type, level) %>% mutate(FDR=p.adjust(P, "fdr")) %>% ungroup()
write.csv(res_4g_other,  "res_4g_other_include_primary_stroke.csv", row.names=FALSE)

# 表格 Table S23: outcome-wide associations excluding primary stroke outcomes
pacman::p_load(dplyr, data.table, stringr, tidyr, ggplot2, writexl, rlang)
fmt_p  <- \(p) ifelse(is.na(p), "", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
fmt_hr <- \(hr,l,h) ifelse(is.na(hr), "", sprintf("%.2f (%.2f-%.2f)", hr,l,h))
get_label <- function(d){lab <- intersect(c("label","label.x","label.y"), names(d)); if(length(lab)==0) return(mutate(d, label=NA_character_)); mutate(d, label=coalesce(!!!syms(lab))) %>% select(-any_of(setdiff(lab, "label")))}
ord <- fn_other2$name
res_std <- read.csv("res_std_other_include_primary_stroke.csv") %>% left_join(fn_other2 %>% select(name,label), by=c("outcome"="name")) %>% get_label()
res_4g <- read.csv("res_4g_other_include_primary_stroke.csv") %>% left_join(fn_other2 %>% select(name,label), by=c("outcome"="name")) %>% get_label()
ev4g <- read.csv("tab_event_other_4g_include_primary_stroke.csv") %>% filter(exposure=="diet.dose5.4g") %>% left_join(fn_other2 %>% select(name,label), by=c("outcome"="name")) %>% get_label() %>% mutate(q_level=as.character(level))
std_row <- res_std %>% transmute(outcome, label, q_level=NA_character_, row_order=0, `DOSE diet score`="Per 10th-to-90th percentile increment", Total=total, Events=case, `HR (95% CI)`=fmt_hr(HR,LCI,HCI), `P value`=fmt_p(P), FDR=fmt_p(FDR))
ref_row <- ev4g %>% filter(q_level=="0-1") %>% transmute(outcome, label, q_level, row_order=1, `DOSE diet score`="Q1 [0-1]", Total=total, Events=case, `HR (95% CI)`="Ref.", `P value`="", FDR="")
cat_row <- res_4g %>% filter(type=="Cat") %>% mutate(
	q_level=case_when(str_detect(level,"4-5$")~"4-5", str_detect(level,"3$")~"3", str_detect(level,"2$")~"2", TRUE~NA_character_),
	row_order=match(q_level, c("2","3","4-5"))+1,
	`DOSE diet score`=recode(q_level, `2`="Q2 (1-2]", `3`="Q3 (2-3]", `4-5`="Q4 (3-5]")
) %>% left_join(ev4g %>% select(outcome, q_level, Total=total, Events=case), by=c("outcome","q_level")) %>% transmute(outcome, label, q_level, row_order, `DOSE diet score`, Total, Events, `HR (95% CI)`=fmt_hr(HR,LCI,HCI), `P value`=fmt_p(P), FDR=fmt_p(FDR))
trend_row <- res_4g %>% filter(type=="Ptrend") %>% transmute(outcome, label, q_level=NA_character_, row_order=5, `DOSE diet score`="P for trend", Total=NA_integer_, Events=NA_integer_, `HR (95% CI)`="", `P value`=fmt_p(P), FDR=fmt_p(FDR))
tab_s23 <- bind_rows(std_row, ref_row, cat_row, trend_row) %>% mutate(outcome=factor(outcome, levels=ord)) %>% arrange(outcome, row_order) %>% group_by(outcome) %>% mutate(Outcome=ifelse(row_number()==1, as.character(label), "")) %>% ungroup() %>%
	transmute(Outcome, `DOSE diet score`, Total=ifelse(is.na(Total), "", format(Total, big.mark=",")), Events=ifelse(is.na(Events), "", format(Events, big.mark=",")),`HR (95% CI)`, `P value`, FDR)
write.csv(tab_s23, "Table_S23_Outcome_wide_DOSE_include_primary_stroke.csv", row.names=FALSE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩Fig5 Outcome-wide
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pacman::p_load(dplyr, ggplot2, stringr, tidyr, rlang)
dat_fig5 <- read.csv("res_std_other_include_primary_stroke.csv") %>% left_join(fn_other2 %>% select(name, label), by = c("outcome" = "name")) %>% get_label() %>%
	mutate(label = ifelse(is.na(label) | label == "", outcome, label), FDR = p.adjust(P, "fdr"), sig  = FDR < 0.05,
		Category = case_when(
			outcome %in% c(primary_stroke, "cvd","htn","ihd","arrhythmia","cvd_afib","hf","cvd_vte","cerebro","pad") ~ "Circulatory system disorders",
			outcome %in% c("infect","infect_bact","infect_viral") ~ "Infections",
			outcome %in% c("cancer","ca_crc","ca_esophagus","ca_stomach","ca_lung","ca_melanoma","ca_breast","ca_prostate","ca_leukemia") ~ "Cancer",
			outcome %in% c("blood_immune","anemia") ~ "Blood & immune",
			outcome %in% c("endo","dm","obesity") ~ "Endocrine",
			outcome %in% c("mental","dementia","substance_abuse","psychotic","mood","neurotic") ~ "Mental & behavioural",
			outcome %in% c("neuro","parkinson","epilepsy","tia","sleep_disorder") ~ "Nervous system",
			outcome %in% c("eye","cataract","glaucoma","ear") ~ "Eye & ear",
			outcome %in% c("resp","flu_pneumonia","copd","asthma") ~ "Respiratory",
			outcome %in% c("digestive","ibd","liver","pancreatitis") ~ "Digestive",
			outcome %in% c("skin","skin_infect_eczema") ~ "Skin",
			outcome %in% c("musculo","ra_disorder","oa","sciatica","soft_tissue","oporosis") ~ "Musculoskeletal",
			outcome %in% c("gu","renal_failure","ckd") ~ "Genitourinary",
			TRUE ~ "Other")
)
pal_disease <- c("Circulatory system disorders" = "#D73027", "Infections" = "#CF9273", "Cancer" = "#D1A09A", "Blood & immune" = "#D8BF73", "Endocrine" = "#A1B36F", "Mental & behavioural" = "#7EA986", "Nervous system" = "#91AFC1", "Eye & ear" = "#89A9BF", "Respiratory" = "#8FAF9A", "Digestive" = "#C5A56A", "Skin" = "#C3899D", "Musculoskeletal" = "#B4A37B", "Genitourinary" = "#BE665E", "Other" = "#B8B8B8")
dat_fig5 <- dat_fig5 %>% mutate(Category = factor(Category, levels = names(pal_disease))) %>% arrange(Category, P, abs(HR - 1))
base_inner <- 1
plot0 <- dat_fig5 %>% transmute(individual = label, group = Category, HR, LCI, HCI, P, sig, plot_value = base_inner - HR, plot_LCI = base_inner - HCI, plot_HCI = base_inner - LCI, lab_col = ifelse(sig, "#2F3742", "#9CA3AF"), line_col = ifelse(sig, "#374151", "#BDBDBD"))
gap_cat <- plot0 %>% distinct(group) %>% slice(1:(n() - 1)) %>% uncount(1) %>% mutate(individual = paste0("cat_gap_", group, "_", row_number()), P = Inf)
gap_n <- 6
gap0 <- tibble(individual = paste0("gap_", 1:gap_n), group = factor(NA, levels = names(pal_disease)), HR = NA_real_, LCI = NA_real_, HCI = NA_real_, P = Inf, plot_value = NA_real_, plot_LCI = NA_real_, plot_HCI = NA_real_, lab_col = NA_character_, line_col = NA_character_)
plotdat <- bind_rows(gap0, bind_rows(plot0, gap_cat) %>% arrange(group, P))
total_n <- nrow(plotdat); axis_x <- (gap_n + 1) / 2
start_offset_rad <- -2 * pi * (axis_x - .5) / total_n
plotdat <- plotdat %>% mutate(xid = row_number(), theta_deg = (360 * (xid - .5) / total_n + start_offset_rad * 180 / pi) %% 360, angle = 90 - theta_deg, hjust = ifelse(theta_deg > 180, 1, 0), angle = ifelse(theta_deg > 180, angle + 180, angle))
ticks <- tibble(tick_val = c(1, .9, .8, .7, .6), plot_loc = base_inner - tick_val, txt = sprintf("%.1f", tick_val))
axis_line <- tibble(x = axis_x, xend = axis_x, y = min(ticks$plot_loc), yend = max(ticks$plot_loc))
axis_tick <- ticks %>% mutate(x = axis_x - .35, xend = axis_x + .35, xlab = axis_x + 1.2)
p5 <- ggplot(plotdat, aes(xid, plot_value)) +
	geom_hline(yintercept = ticks$plot_loc, colour = "#ECECEC", linewidth = .3) +
	geom_col(aes(fill = group), width = .85, na.rm = TRUE) +
	geom_errorbar(aes(ymin = plot_LCI, ymax = plot_HCI, colour = line_col), width = .3, linewidth = .4, na.rm = TRUE) +
	geom_text(data = filter(plotdat, !str_detect(individual, "^gap_|^cat_gap_")), aes(x = xid, y = plot_HCI + .05, label = individual, angle = angle, hjust = hjust, colour = lab_col), size = 2.4, inherit.aes = FALSE, na.rm = TRUE) +
	geom_segment(data = axis_line, aes(x = x, xend = xend, y = y, yend = yend), inherit.aes = FALSE, linewidth = .5, colour = "#6B7280") + geom_segment(data = axis_tick, aes(x = x, xend = xend, y = plot_loc, yend = plot_loc), inherit.aes = FALSE, linewidth = .5, colour = "#6B7280") +
	geom_text(data = axis_tick, aes(x = xlab, y = plot_loc, label = txt), inherit.aes = FALSE, hjust = 0, size = 3, colour = "#4B5563", fontface = "bold") +
	annotate("text", x = axis_x, y = max(ticks$plot_loc) + .06, label = "HR", size = 3.3, fontface = "bold", colour = "#2F3742") +
	scale_fill_manual(values = pal_disease, na.translate = FALSE, name = "Disease Category") + scale_colour_identity() + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(limits = c(min(plotdat$plot_LCI, na.rm = TRUE) - .1, max(plotdat$plot_HCI, na.rm = TRUE) + .28), expand = c(0, 0)) +
	coord_polar(start = start_offset_rad, clip = "off") +
	theme_void() + theme(legend.position = "right", legend.title = element_text(face = "bold", size = 11), legend.text = element_text(size = 10), plot.margin = margin(20, 20, 20, 20))
ggsave("Fig5.pdf", p5, width = 13, height = 10, device = cairo_pdf)
