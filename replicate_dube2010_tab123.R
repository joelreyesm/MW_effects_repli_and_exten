# =============================================================
# Replication: Dube, Lester, Reich (2010)
# "Minimum Wages and Low-Wage Jobs: Evidence from
#  Cross-State Comparisons Across Time"
# Replication of Tables 1, 2, and 3
#
# Data files (must be in working directory):
#   QCEWindustry_minwage_all.dta
#   QCEWindustry_minwage_contig.dta
#
# Corrections applied (relative to original):
#   1. add_leads_contig: group_by(countyreal) instead of pair_id_county
#      (leads must follow a single county's MW through time)
#   2. Table 3: use differenced leads per equation (8) of the paper
#      DlnMW_t0_t4  = lnMW_f4  - lnMW       [g4:  1-yr-ahead level]
#      DlnMW_t4_t12 = lnMW_f4  - lnMW_f12   [g12: 3-yr-ahead level, sign per paper]
#   3. Spec 3/4: renamed absorb_censdiv / absorb_msa to avoid swap;
#      Spec 3 = censdiv x period FE + state trends (collinearity warning is benign)
#      Spec 4 = MSA x period FE
#   4. get_lincom: correct delta-method SE via matrix quadratic form
# =============================================================

library(haven)
library(dplyr)
library(fixest)

# ---- Load and clean data ------------------------------------
dat_all    <- read_dta("QCEWindustry_minwage_all.dta")    |> zap_labels()
dat_contig <- read_dta("QCEWindustry_minwage_contig.dta") |> zap_labels()


# =============================================================
# TABLE 1: Descriptive Statistics
# =============================================================

cat("\n")
cat("=============================================================\n")
cat("TABLE 1: DESCRIPTIVE STATISTICS\n")
cat("=============================================================\n")

print_desc <- function(df, panel_label, has_pairs = FALSE) {
  cat(sprintf("\n--- %s ---\n", panel_label))
  
  d <- df |> filter(nonmissing_rest_both == 66)
  
  # -- Main county-level variables --
  main_vars <- c("countypop2000", "cntypopdens", "cntyarea",
                 "emp_rest_both", "AWW_rest_both", "minwage")
  cat(sprintf("\n%-22s  %13s  %13s\n", "Variable", "Mean", "Std. Dev."))
  cat(strrep("-", 52), "\n")
  for (v in main_vars) {
    cat(sprintf("%-22s  %13.3f  %13.3f\n", v,
                mean(d[[v]], na.rm = TRUE),
                sd(d[[v]],   na.rm = TRUE)))
  }
  
  # -- Industry employment and wages --
  # For all-county dataset: use ACFS (not ACFSRETAIL which is contig-only)
  # For contig dataset:     use ACFSRETAIL
  cat(sprintf("\n%-14s  %12s  %12s  %12s  %12s\n",
              "Industry", "Mean Emp", "SD Emp", "Mean AWW", "SD AWW"))
  cat(strrep("-", 66), "\n")
  
  industry_keys <- if ("empACFSRETAIL" %in% names(df)) {
    c("ACFSRETAIL", "RETAIL", "MFG", "TOT")
  } else {
    c("ACFS", "RETAIL", "MFG", "TOT")
  }
  
  for (k in industry_keys) {
    nm_col  <- paste0("nonmissing_", k)
    emp_col <- paste0("emp", k)
    aww_col <- paste0("AWW", k)
    if (!emp_col %in% names(df) || !aww_col %in% names(df)) {
      cat(sprintf("%-14s  %12s  %12s  %12s  %12s  (not available)\n",
                  k, "", "", "", ""))
      next
    }
    dk <- if (nm_col %in% names(df)) df |> filter(.data[[nm_col]] == 66) else d
    cat(sprintf("%-14s  %12.3f  %12.3f  %12.3f  %12.3f\n", k,
                mean(dk[[emp_col]], na.rm = TRUE),
                sd(dk[[emp_col]],   na.rm = TRUE),
                mean(dk[[aww_col]], na.rm = TRUE),
                sd(dk[[aww_col]],   na.rm = TRUE)))
  }
  
  # -- Sample sizes --
  n_cty <- d |> distinct(county) |> nrow()
  cat(sprintf("\nNumber of counties:     %d\n", n_cty))
  
  if (has_pairs) {
    n_pairs <- df |> filter(nonmissing_both_pair == 66) |>
      distinct(pair_id) |> nrow()
    cat(sprintf("Number of county pairs: %d\n", n_pairs))
  }
  
  # -- Minimum-wage events per county --
  divisor <- if (has_pairs) 2 else 1
  if ("event_type" %in% names(d)) {
    fed_ev <- d |> filter(event_type == 1) |>
      count(county) |>
      summarise(m = mean(n / divisor)) |> pull(m)
    sta_ev <- d |> filter(event_type == 2) |>
      count(county) |>
      summarise(m = mean(n / divisor)) |> pull(m)
    cat(sprintf("Mean federal events / county: %.3f\n", fed_ev))
    cat(sprintf("Mean state   events / county: %.3f\n", sta_ev))
  }
}

print_desc(dat_all,    "All-County Sample",  has_pairs = FALSE)
print_desc(dat_contig, "County-Pair Sample", has_pairs = TRUE)


# =============================================================
# TABLE 2: Minimum Wage Effects on Earnings and Employment
#
# Spec 1: All counties,            county + period FE
# Spec 2: All counties,            county + period x CensDiv FE
# Spec 3: All counties,            county + period x CensDiv FE
#                                  + state-specific linear trends
# Spec 4: All counties,            county + period x MSA FE
# Spec 5: Contiguous county pairs, county + period FE
# Spec 6: Contiguous county pairs, county + pair x period FE
#
# Control sets:
#   c2  AWW: no extra controls
#   c2  EMP: lnpop
#   c3  AWW: lnAWW_TOT
#   c3  EMP: lnpop + lnemp_TOT
#
# SE clustering:
#   Specs 1-4: clustered by state_fips
#   Specs 5-6: two-way CGM by state_fips x bordersegment
# =============================================================

cat("\n")
cat("=============================================================\n")
cat("TABLE 2: MINIMUM WAGE EFFECTS ON EARNINGS AND EMPLOYMENT\n")
cat("=============================================================\n")
cat("Coefficient on log(minimum wage) | SEs in parentheses\n")
cat("c2 = no extra / pop control; c3 = pop + total emp/wage control\n")

# ---- Prepare all-county analysis dataset --------------------
d_all <- dat_all |>
  filter(nonmissing_rest_both == 66) |>
  mutate(
    absorb_censdiv = as.integer(interaction(period, censusdiv, drop = TRUE)),
    absorb_msa     = as.integer(interaction(period, cbmsa,     drop = TRUE)),
    state_fips_f   = factor(state_fips),
    censdiv_f      = factor(censusdiv)
  )

# ---- Helper: pull MW coefficient and SE ---------------------
mw_coef <- function(mod) {
  b  <- coef(mod)["lnMW"]
  se <- se(mod)["lnMW"]
  c(est = unname(b), se = unname(se))
}

# ---- Print one specification row ----------------------------
print_spec_row <- function(label, aww_c2, aww_c3, emp_c2, emp_c3, n) {
  cat(sprintf("\n  [%s]  N = %d\n", label, as.integer(n)))
  cat(sprintf("    AWW  c2: %8.4f  (%6.4f)   c3: %8.4f  (%6.4f)\n",
              aww_c2[1], aww_c2[2], aww_c3[1], aww_c3[2]))
  cat(sprintf("    EMP  c2: %8.4f  (%6.4f)   c3: %8.4f  (%6.4f)\n",
              emp_c2[1], emp_c2[2], emp_c3[1], emp_c3[2]))
}

# ---- Spec 1: All counties, common period FE -----------------
m_aww_c2_s1 <- feols(lnAWW_rest_both ~ lnMW             | county + period, data = d_all, cluster = ~state_fips)
m_aww_c3_s1 <- feols(lnAWW_rest_both ~ lnMW + lnAWW_TOT | county + period, data = d_all, cluster = ~state_fips)
m_emp_c2_s1 <- feols(lnemp_rest_both ~ lnMW + lnpop     | county + period, data = d_all, cluster = ~state_fips)
m_emp_c3_s1 <- feols(lnemp_rest_both ~ lnMW + lnpop + lnemp_TOT | county + period, data = d_all, cluster = ~state_fips)
print_spec_row("Spec 1: All counties, period FE",
               mw_coef(m_aww_c2_s1), mw_coef(m_aww_c3_s1),
               mw_coef(m_emp_c2_s1), mw_coef(m_emp_c3_s1),
               nobs(m_aww_c2_s1))

# ---- Spec 2: All counties, census division x period FE ------
m_aww_c2_s2 <- feols(lnAWW_rest_both ~ lnMW             | county + absorb_censdiv, data = d_all, cluster = ~state_fips)
m_aww_c3_s2 <- feols(lnAWW_rest_both ~ lnMW + lnAWW_TOT | county + absorb_censdiv, data = d_all, cluster = ~state_fips)
m_emp_c2_s2 <- feols(lnemp_rest_both ~ lnMW + lnpop     | county + absorb_censdiv, data = d_all, cluster = ~state_fips)
m_emp_c3_s2 <- feols(lnemp_rest_both ~ lnMW + lnpop + lnemp_TOT | county + absorb_censdiv, data = d_all, cluster = ~state_fips)
print_spec_row("Spec 2: All counties, censdiv x period FE",
               mw_coef(m_aww_c2_s2), mw_coef(m_aww_c3_s2),
               mw_coef(m_emp_c2_s2), mw_coef(m_emp_c3_s2),
               nobs(m_aww_c2_s2))

# ---- Spec 3: Census division x period FE + state trends -----
# censdiv x period dummies and state x period slopes are included as
# explicit regressors under a common period FE. feols drops the subset of
# censdiv x period cells that are fully spanned by the state slopes;
# the collinearity warning is expected and benign.
m_aww_c2_s3 <- feols(lnAWW_rest_both ~ lnMW + censdiv_f:factor(period) + state_fips_f:period | county + period,
                     data = d_all, cluster = ~state_fips)
m_aww_c3_s3 <- feols(lnAWW_rest_both ~ lnMW + lnAWW_TOT + censdiv_f:factor(period) + state_fips_f:period | county + period,
                     data = d_all, cluster = ~state_fips)
m_emp_c2_s3 <- feols(lnemp_rest_both ~ lnMW + lnpop + censdiv_f:factor(period) + state_fips_f:period | county + period,
                     data = d_all, cluster = ~state_fips)
m_emp_c3_s3 <- feols(lnemp_rest_both ~ lnMW + lnpop + lnemp_TOT + censdiv_f:factor(period) + state_fips_f:period | county + period,
                     data = d_all, cluster = ~state_fips)
print_spec_row("Spec 3: All counties, censdiv x period FE + state trends",
               mw_coef(m_aww_c2_s3), mw_coef(m_aww_c3_s3),
               mw_coef(m_emp_c2_s3), mw_coef(m_emp_c3_s3),
               nobs(m_aww_c2_s3))

# ---- Spec 4: All counties, MSA x period FE ------------------
m_aww_c2_s4 <- feols(lnAWW_rest_both ~ lnMW             | county + absorb_msa, data = d_all, cluster = ~state_fips)
m_aww_c3_s4 <- feols(lnAWW_rest_both ~ lnMW + lnAWW_TOT | county + absorb_msa, data = d_all, cluster = ~state_fips)
m_emp_c2_s4 <- feols(lnemp_rest_both ~ lnMW + lnpop     | county + absorb_msa, data = d_all, cluster = ~state_fips)
m_emp_c3_s4 <- feols(lnemp_rest_both ~ lnMW + lnpop + lnemp_TOT | county + absorb_msa, data = d_all, cluster = ~state_fips)
print_spec_row("Spec 4: All counties, MSA x period FE",
               mw_coef(m_aww_c2_s4), mw_coef(m_aww_c3_s4),
               mw_coef(m_emp_c2_s4), mw_coef(m_emp_c3_s4),
               nobs(m_aww_c2_s4))

# ---- Prepare contiguous-county analysis dataset -------------
# pair_id format: "SSCCC-SSCCC" e.g. "01003-12033"
# positions: 1-2=state A, 3-5=county A, 6=hyphen, 7-8=state B, 9-11=county B
d_ct <- dat_contig |>
  filter(nonmissing_rest_both == 66) |>
  mutate(
    state_a       = as.integer(substr(pair_id, 1, 2)),
    state_b       = as.integer(substr(pair_id, 7, 8)),
    st_min        = pmin(state_a, state_b),
    st_max        = pmax(state_a, state_b),
    bordersegment = as.integer(interaction(st_min, st_max, drop = TRUE)),
    pair_period   = as.integer(interaction(pair_id, period, drop = TRUE))
  )

# Quick check — print unique state_b values; should all be valid 2-digit FIPS
cat(sprintf("\n  state_b range: %d to %d (should be valid FIPS, no NAs)\n",
            min(d_ct$state_b, na.rm = TRUE),
            max(d_ct$state_b, na.rm = TRUE)))
cat(sprintf("  NAs in state_b: %d\n", sum(is.na(d_ct$state_b))))
cat(sprintf("  Unique bordersegments: %d\n", n_distinct(d_ct$bordersegment)))

# ---- Specs 5-6: contiguous pairs, two-way clustered SEs -----
for (s in 5:6) {
  ab  <- if (s == 5) "period" else "pair_period"
  lbl <- if (s == 5) "Spec 5: Contig. pairs, period FE" else
    "Spec 6: Contig. pairs, pair x period FE"
  
  f_aww_c2 <- as.formula(paste0("lnAWW_rest_both ~ lnMW             | county + ", ab))
  f_aww_c3 <- as.formula(paste0("lnAWW_rest_both ~ lnMW + lnAWW_TOT | county + ", ab))
  f_emp_c2 <- as.formula(paste0("lnemp_rest_both ~ lnMW + lnpop     | county + ", ab))
  f_emp_c3 <- as.formula(paste0("lnemp_rest_both ~ lnMW + lnpop + lnemp_TOT | county + ", ab))
  
  clust <- ~state_fips + bordersegment
  
  m_aww_c2 <- feols(f_aww_c2, data = d_ct, cluster = clust)
  m_aww_c3 <- feols(f_aww_c3, data = d_ct, cluster = clust)
  m_emp_c2 <- feols(f_emp_c2, data = d_ct, cluster = clust)
  m_emp_c3 <- feols(f_emp_c3, data = d_ct, cluster = clust)
  
  print_spec_row(lbl,
                 mw_coef(m_aww_c2), mw_coef(m_aww_c3),
                 mw_coef(m_emp_c2), mw_coef(m_emp_c3),
                 nobs(m_aww_c2))
}


# =============================================================
# TABLE 3: Pre-Existing Trends in Employment/Earnings
#
# Uses differenced leads per equation (8) of the paper.
#   DlnMW_t0_t4  = ln(MW_{t+4})  - ln(MW_t)      [g4  in paper: 1-yr-ahead level]
#   DlnMW_t4_t12 = ln(MW_{t+4})  - ln(MW_{t+12}) [g12 in paper: 3-yr-ahead level,
#                                                   sign reversed so negative = pre-trend]
# Regression: depvar ~ lnMW + DlnMW_t0_t4 + DlnMW_t4_t12 + lnpop | county + absorb
#
# Reported (matching paper's Table 3 notation):
#   (a) DlnMW_t4_t12 alone           -> gt12 coefficient
#   (b) DlnMW_t4_t12 + DlnMW_t0_t4  -> trend = gt4 - gt12
#   (c) DlnMW_t0_t4 alone            -> gt4  coefficient
#
# Three specifications:
#   Spec 1: All counties,     county + period FE,        cluster state
#   Spec 2: Cross-state MSA,  county + period x MSA FE,  cluster state
#   Spec 3: Contig. pairs,    county + pair x period FE,  two-way cluster
# =============================================================

cat("\n")
cat("=============================================================\n")
cat("TABLE 3: PRE-EXISTING TRENDS (LEADS OF LOG MINIMUM WAGE)\n")
cat("=============================================================\n")
cat("Dep vars: lnemp_rest_both, lnemp_TOT, lnAWW_rest_both, lnAWW_TOT\n")
cat("(a) = gt12 coeff; (b) = trend (gt4-gt12); (c) = gt4 coeff\n")

# ---- Generate differenced leads -----------------------------
# Computed on the FULL panel before sample-period restriction so that
# observations near the end of the sample window can borrow MW values
# from periods outside the regression sample (e.g. 2007-2009).

add_leads_all <- function(df) {
  df |>
    filter(nonmissing_rest_both == 66) |>
    arrange(countyreal, period) |>
    group_by(countyreal) |>
    mutate(
      lnMW_f4      = dplyr::lead(lnMW,  4),
      lnMW_f12     = dplyr::lead(lnMW, 12),
      DlnMW_t0_t4  = lnMW_f4  - lnMW,      # g4:  ln(MW_{t+4})  - ln(MW_t)
      DlnMW_t4_t12 = lnMW_f12 - lnMW_f4   # g12: ln(MW_{t+12}) - ln(MW_{t+4})
      #   negative when MW rises between t+4 and t+12
      #   so negative coeff = employment falling ahead of hike
    ) |>
    ungroup()
}

add_leads_contig <- function(df) {
  df |>
    filter(nonmissing_rest_both == 66, !is.na(pair_id_county)) |>
    arrange(countyreal, period) |>
    group_by(countyreal) |>
    mutate(
      lnMW_f4      = dplyr::lead(lnMW,  4),
      lnMW_f12     = dplyr::lead(lnMW, 12),
      DlnMW_t0_t4  = lnMW_f4  - lnMW,      # g4:  ln(MW_{t+4})  - ln(MW_t)
      DlnMW_t4_t12 = lnMW_f12 - lnMW_f4    # g12: ln(MW_{t+12}) - ln(MW_{t+4})
    ) |>
    ungroup()
}

# Generate leads on full panel, then restrict to sample period
samp_all_full <- add_leads_all(dat_all)

samp_all <- samp_all_full |>
  filter(year >= 1990,
         !(year > 2006 | (year == 2006 & quarter >= 3))) |>
  mutate(
    absorb_msa = as.integer(interaction(period, cbmsa, drop = TRUE))
  )

samp_contig_full <- add_leads_contig(dat_contig)

samp_contig <- samp_contig_full |>
  filter(year >= 1990,
         !(year > 2006 | (year == 2006 & quarter >= 3))) |>
  mutate(
    state_a       = as.integer(substr(pair_id, 1, 2)),
    state_b       = as.integer(substr(pair_id, 7, 8)),
    st_min        = pmin(state_a, state_b),
    st_max        = pmax(state_a, state_b),
    bordersegment = as.integer(interaction(st_min, st_max, drop = TRUE)),
    pair_period   = as.integer(interaction(pair_id, period, drop = TRUE))
  )

# ---- Pre-trends regression helper --------------------------
# Regressors: lnMW (level) + DlnMW_t0_t4 + DlnMW_t4_t12
# Delta-method SE via matrix quadratic form: sqrt(w' V w)
run_pretrend_reg <- function(df, depvar, absorb_var, clust_fml) {
  fml <- as.formula(
    paste0(depvar,
           " ~ lnMW + DlnMW_t0_t4 + DlnMW_t4_t12 + lnpop | county + ", absorb_var)
  )
  
  mod <- tryCatch(
    feols(fml, data = df, vcov = clust_fml),
    error = function(e) NULL
  )
  if (is.null(mod)) return(list(a = c(NA, NA), b = c(NA, NA), c = c(NA, NA), N = NA))
  
  co <- coef(mod)
  vc <- vcov(mod)
  
  # delta-method SE = sqrt(w' V w)
  get_lincom <- function(vars) {
    w  <- rep(1, length(vars))
    b  <- sum(co[vars])
    se <- sqrt(as.numeric(t(w) %*% vc[vars, vars, drop = FALSE] %*% w))
    c(est = b, se = se)
  }
  
  list(
    a = get_lincom("DlnMW_t4_t12"),                        # gt12 alone
    b = get_lincom(c("DlnMW_t4_t12", "DlnMW_t0_t4")),     # trend = gt4 - gt12
    c = get_lincom("DlnMW_t0_t4"),                         # gt4 alone
    N = nobs(mod)
  )
}

# ---- Print one Table 3 block --------------------------------
depvars_t3 <- c("lnemp_rest_both", "lnemp_TOT", "lnAWW_rest_both", "lnAWW_TOT")
col_lbls   <- c("Emp-Rest", "Emp-TOT", "AWW-Rest", "AWW-TOT")

print_t3_block <- function(df, spec_label, absorb_var, clust_fml) {
  cat(sprintf("\n  Specification: %s\n", spec_label))
  
  res <- setNames(
    lapply(depvars_t3, run_pretrend_reg,
           df = df, absorb_var = absorb_var, clust_fml = clust_fml),
    depvars_t3
  )
  
  N_val <- res[[1]]$N
  
  cat(sprintf("  %-13s", ""))
  for (lbl in col_lbls) cat(sprintf("  %16s", lbl))
  cat(sprintf("  %6s\n", "N"))
  
  row_labels <- c("(a) gt12 alone", "(b) trend     ", "(c) gt4  alone")
  row_keys   <- c("a", "b", "c")
  
  for (rk in seq_along(row_keys)) {
    cat(sprintf("  %-13s", row_labels[rk]))
    for (dv in depvars_t3) {
      r <- res[[dv]][[row_keys[rk]]]
      if (anyNA(r)) {
        cat(sprintf("  %7s (%6s)", "NA", "NA"))
      } else {
        cat(sprintf("  %7.4f (%6.4f)", r[1], r[2]))
      }
    }
    if (rk == 1) cat(sprintf("  %6d", as.integer(N_val))) else cat(sprintf("  %6s", ""))
    cat("\n")
  }
}

# ---- Spec 1: All counties, period FE ------------------------
print_t3_block(samp_all,
               "1: All counties, period FE",
               "period",
               ~state_fips)

# ---- Spec 2: Cross-state MSA, period x MSA FE ---------------
print_t3_block(samp_all |> filter(sample_flag == 2),
               "2: Cross-state MSA, period x MSA FE",
               "absorb_msa",
               ~state_fips)

# ---- Spec 3: Contiguous pairs, pair x period FE -------------
print_t3_block(samp_contig,
               "3: Contiguous pairs, pair x period FE",
               "pair_period",
               ~state_fips + bordersegment)

cat("\n")
cat("=============================================================\n")
cat("Replication complete.\n")
cat("=============================================================\n")