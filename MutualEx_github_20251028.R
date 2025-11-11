# =========================================================
# TTE Simulation (Step1–4 ) — GitHub-ready
#  - Stabilized IP weights (pooled trimming), Modified Poisson (weighted)
#  - PS(A) / PS(B) learned separately (where applicable); multinomial GPS for 4-category model
#  - Single seed at top; no EPS; drop non-finite
#  - Unified plotting via dictionaries + 3 helpers
#  - Offsets: only where strictly needed for readability
# =========================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(purrr)
  library(ggplot2); library(patchwork)
  library(nnet);    library(splines)
})

## ---- Global seed ----
set.seed(123)

## ---- Truth & global configs ----
RR_A <- 0.3
RR_B <- 0.9
TRUE_MARGINAL_RR <- RR_B / RR_A

TRIM_LEVELS    <- c(0, 0.001, 0.005, 0.01, 0.02, 0.03, 0.05)
N_SIM_DEFAULT  <- 1000             # increase for the final run if needed
N_GRID         <- c(100, 500, 1000, 3000, 10000)
TRIMS_FOR_SIZE <- c(0.00, 0.01, 0.03)

## ---- Unified dictionaries (labels/colors/linetypes) ----
STRATEGY_LABELS <- c(
  "cond_poisson"          = "Conditional RR (reference)",
  "tte_ipw_unrestricted"  = "Strategy 0-1: Unrestricted TTE",
  "tte_ipw_subset"        = "Strategy 0-2: Subset-Weighted TTE",
  "tte_ipw_strict"        = "Strategy 0-3: Contamination-Aware TTE",
  "exclusive"             = "Strategy 1: Exclusive",
  "ipw"                   = "Strategy 2: IPW (A_only vs B_only)",
  "gps_q"                 = "Strategy 3: GPS (quadratic)",
  "gps_spline"            = "Strategy 4: GPS (splines)",
  "dual"                  = "Strategy 5: Dual-IPW",
  "newdual"               = "Strategy 6: Sequential Dual-IPW"
)
STRATEGY_COLORS <- c(
  "cond_poisson"="#7f7f7f","tte_ipw_unrestricted"="#1f77b4","tte_ipw_subset"="#ff7f0e",
  "tte_ipw_strict"="#2ca02c","exclusive"="#1f77b4","ipw"="#ff7f0e","gps_q"="#2ca02c",
  "gps_spline"="#d62728","dual"="#9467bd","newdual"="#8c564b"
)
STRATEGY_LINETYPE <- c(
  "cond_poisson"="dashed","tte_ipw_unrestricted"="solid","tte_ipw_subset"="solid",
  "tte_ipw_strict"="solid","exclusive"="dashed","ipw"="dashed","gps_q"="dashed",
  "gps_spline"="dashed","dual"="dashed","newdual"="dashed"
)

# =========================================================
# Utilities: trimming, ESS, modified Poisson
# =========================================================
apply_trimming <- function(weights, level) {
  if (level == 0) return(weights)
  wpos <- weights[is.finite(weights) & weights > 0]
  if (!length(wpos)) return(rep(NA_real_, length(weights)))
  q <- stats::quantile(wpos, probs = c(level, 1 - level), na.rm = TRUE)
  out <- weights
  out[weights > 0 & (weights < q[1] | weights > q[2])] <- NA_real_
  out
}

compute_ess <- function(w, grp = NULL) {
  ess_all <- if (all(is.na(w))) NA_real_ else (sum(w, na.rm = TRUE)^2) / sum(w^2, na.rm = TRUE)
  if (is.null(grp)) return(list(ESS_total = ess_all, ESS_A = NA_real_, ESS_B = NA_real_))
  ess_A <- if (all(is.na(w[grp == "A"]))) NA_real_ else (sum(w[grp == "A"], na.rm = TRUE)^2) / sum(w[grp == "A"]^2, na.rm = TRUE)
  ess_B <- if (all(is.na(w[grp == "B"]))) NA_real_ else (sum(w[grp == "B"], na.rm = TRUE)^2) / sum(w[grp == "B"]^2, na.rm = TRUE)
  list(ESS_total = ess_all, ESS_A = ess_A, ESS_B = ess_B)
}

estimate_rr_poisson <- function(long) {
  dat <- long %>% dplyr::filter(is.finite(w), w > 0, group %in% c("A","B")) %>%
    dplyr::mutate(B = as.integer(group == "B"))
  if (!nrow(dat)) return(list(rr_pois = NA_real_))
  fit <- tryCatch(glm(O ~ B, data = dat, weights = w, family = poisson()), error = function(e) NULL)
  if (is.null(fit) || is.na(coef(fit)["B"])) return(list(rr_pois = NA_real_))
  list(rr_pois = exp(coef(fit)["B"]))
}

# =========================================================
# Data generators
# =========================================================
generate_data <- function(n, scenario = 1, type = c("mutual", "nonmutual")) {
  type <- match.arg(type)
  L <- rnorm(n)
  if (scenario == 1) {
    p_A <- plogis(-0.5 + 1.2 * L)
    p_B <- plogis(-1.5 + 0.7 * L)
  } else if (scenario == 2) {
    p_A <- plogis(-0.5 + 1.2 * L)
    p_B <- plogis(-0.5 - 1.2 * L)
  } else stop("scenario must be 1 or 2")
  if (type == "mutual") { A <- rbinom(n,1,p_A); B <- 1 - A } else { A <- rbinom(n,1,p_A); B <- rbinom(n,1,p_B) }
  log_risk <- log(0.1) + 0.3 * L + log(RR_A) * A + log(RR_B) * B
  p_O <- pmin(exp(log_risk), 1)
  O   <- rbinom(n, 1, p_O)
  tibble(L, A, B, O)
}

generate_data_step2 <- function(n, a, b) {
  L <- rnorm(n)
  pA <- plogis(a + 1.2*L)
  pB <- plogis(b + 0.7*L)
  A <- rbinom(n,1,pA); B <- rbinom(n,1,pB)
  pO <- pmin(exp(log(0.1)+0.3*L + log(RR_A)*A + log(RR_B)*B), 1)
  O  <- rbinom(n,1,pO)
  tibble(L,A,B,O)
}

generate_data_step4 <- function(n, scenario) {
  L <- rnorm(n)
  if (scenario == 1) {
    p_A <- plogis(-0.5 + 1.2 * L)
    p_B <- plogis(-1.5 + 0.7 * L)
  } else if (scenario == 2) {
    p_A <- plogis(-0.5 + 1.2 * L)
    p_B <- plogis(-0.5 - 1.2 * L)
  } else stop("scenario must be 1 or 2")
  A <- rbinom(n, 1, p_A)
  B <- rbinom(n, 1, p_B)
  p_O <- pmin(exp(log(0.1) + 0.3*L + log(RR_A)*A + log(RR_B)*B), 1)
  O   <- rbinom(n, 1, p_O)
  tibble(L, A, B, O)
}

# =========================================================
# Estimators
# =========================================================
estimate_conditional_rr <- function(dat) {
  fit <- tryCatch(glm(O ~ factor(B) + L, data = dat, family = poisson()), error = function(e) NULL)
  if (is.null(fit)) return(list(est = NA_real_, ESS = list(ESS_total=NA_real_,ESS_A=NA_real_,ESS_B=NA_real_)))
  cf <- coef(fit)
  b_idx <- which(names(cf) %in% c("factor(B)1","B"))
  if (length(b_idx) == 0L || is.na(cf[b_idx[1]])) {
    return(list(est = NA_real_, ESS = list(ESS_total=NA_real_,ESS_A=NA_real_,ESS_B=NA_real_)))
  }
  list(est = exp(cf[b_idx[1]]), ESS = list(ESS_total=NA_real_,ESS_A=NA_real_,ESS_B=NA_real_))
}

# Strategy 0-1 / 0-2 / 0-3
estimate_rr_01 <- function(dat, trim_level = 0) {
  psA <- glm(A ~ L, data = dat, family = binomial())
  psB <- glm(B ~ L, data = dat, family = binomial())
  d2  <- dat %>% mutate(psA = predict(psA, type="response"), psB = predict(psB, type="response"))
  pA <- mean(d2$A==1); pB <- mean(d2$B==1)
  long <- bind_rows(
    d2 %>% transmute(group="A", O=O, w = ifelse(A==1, pA/psA, 0)),
    d2 %>% transmute(group="B", O=O, w = ifelse(B==1, pB/psB, 0))
  )
  long$w <- apply_trimming(long$w, trim_level)
  long <- long %>% filter(!is.na(w))
  list(est = estimate_rr_poisson(long)$rr_pois, ESS = compute_ess(long$w, grp = long$group))
}

estimate_rr_02 <- function(dat, trim_level = 0) {
  psA <- glm(A ~ L, data = dat, family = binomial())
  psB <- glm(B ~ L, data = dat, family = binomial())
  d2  <- dat %>% mutate(psA = predict(psA, type="response"), psB = predict(psB, type="response"))
  subset_ab <- d2 %>% filter((A==1 & B==0) | (A==0 & B==1))
  pAonly <- mean(subset_ab$A==1); pBonly <- 1 - pAonly
  long <- bind_rows(
    subset_ab %>% filter(A==1 & B==0) %>% transmute(group="A", O=O, w = pAonly/psA),
    subset_ab %>% filter(A==0 & B==1) %>% transmute(group="B", O=O, w = pBonly/psB)
  )
  long$w <- apply_trimming(long$w, trim_level)
  long <- long %>% filter(!is.na(w))
  list(est = estimate_rr_poisson(long)$rr_pois, ESS = compute_ess(long$w, grp = long$group))
}

estimate_rr_03 <- function(dat, trim_level = 0) {
  psA <- glm(A ~ L, data = dat, family = binomial())
  psB <- glm(B ~ L, data = dat, family = binomial())
  dat2 <- dat %>% mutate(pA = predict(psA, type="response"), pB = predict(psB, type="response"))
  ps_B0_given_A1 <- glm(I(B==0) ~ L, data = dat2 %>% filter(A==1), family = binomial())
  ps_A0_given_B1 <- glm(I(A==0) ~ L, data = dat2 %>% filter(B==1), family = binomial())
  dat2 <- dat2 %>%
    mutate(
      pb0_a1 = ifelse(A==1, predict(ps_B0_given_A1, newdata = dat2, type="response"), NA_real_),
      pa0_b1 = ifelse(B==1, predict(ps_A0_given_B1, newdata = dat2, type="response"), NA_real_)
    )
  subset_ab <- dat2 %>% filter((A==1 & B==0) | (A==0 & B==1))
  pAonly <- mean(subset_ab$A==1); pBonly <- 1 - pAonly
  long <- bind_rows(
    subset_ab %>% filter(A==1 & B==0) %>% transmute(group="A", O=O, w = pAonly / (pA * pb0_a1)),
    subset_ab %>% filter(A==0 & B==1) %>% transmute(group="B", O=O, w = pBonly / (pB * pa0_b1))
  )
  long$w <- apply_trimming(long$w, trim_level)
  long <- long %>% filter(is.finite(w), !is.na(w), w > 0)
  list(est = estimate_rr_poisson(long)$rr_pois, ESS = compute_ess(long$w, grp = long$group))
}

# Step4 methods (A_only vs B_only comparisons)
compute_weights_analysis2 <- function(method, data) {
  if (method == "exclusive") {
    df <- data %>% filter((A==1 & B==0) | (A==0 & B==1))
    psA <- glm(A ~ L, data = df, family = binomial())
    psB <- glm(B ~ L, data = df, family = binomial())
    pAhat <- predict(psA, type="response"); pBhat <- predict(psB, type="response")
    pA <- mean(df$A==1); pB <- mean(df$B==1)
    df$weight <- ifelse(df$A==1, pA/pAhat, pB/pBhat)
    df$weight[!is.finite(df$weight)] <- NA_real_
    return(df %>% select(L,A,B,O,weight))
  }
  if (method == "ipw") {
    df <- data %>% filter((A==1 & B==0) | (A==0 & B==1))
    ps  <- glm(I(A==1) ~ L, data = df, family = binomial())
    phat <- predict(ps, type="response")
    pA <- mean(df$A==1)
    df$weight <- ifelse(df$A==1, pA/phat, (1-pA)/(1-phat))
    df$weight[!is.finite(df$weight)] <- NA_real_
    return(df %>% select(L,A,B,O,weight))
  }
  if (method %in% c("gps_q","gps_spline")) {
    df0 <- data %>%
      mutate(cat = case_when(
        A==1 & B==0 ~ "A_only",
        A==0 & B==1 ~ "B_only",
        A==1 & B==1 ~ "Both",
        TRUE        ~ "Neither"
      )) %>% mutate(cat = factor(cat, levels = c("A_only","B_only","Both","Neither")))
    gps <- if (method=="gps_q") multinom(cat ~ L + I(L^2), data = df0, trace = FALSE)
    else                 multinom(cat ~ splines::ns(L, df = 4), data = df0, trace = FALSE)
    probs <- as.data.frame(predict(gps, type="probs"))
    need_cols <- levels(df0$cat); miss <- setdiff(need_cols, colnames(probs))
    for (m in miss) probs[[m]] <- NA_real_
    probs <- probs[, need_cols, drop = FALSE]
    colnames(probs) <- paste0("gps_", need_cols)
    df0 <- bind_cols(df0, probs)
    pmarg <- prop.table(table(df0$cat))
    df <- df0 %>% filter(cat %in% c("A_only","B_only")) %>%
      mutate(
        weight = case_when(
          cat=="A_only" ~ as.numeric(pmarg["A_only"]) / gps_A_only,
          cat=="B_only" ~ as.numeric(pmarg["B_only"]) / gps_B_only
        ),
        A = as.integer(cat=="A_only"),
        B = as.integer(cat=="B_only")
      )
    df$weight[!is.finite(df$weight)] <- NA_real_
    return(df %>% select(L,A,B,O,weight))
  }
  if (method == "dual") {
    psA <- glm(A ~ L, data = data, family = binomial())
    psB <- glm(B ~ L, data = data, family = binomial())
    df <- data %>%
      mutate(
        pAh = predict(psA, type="response"),
        pBh = predict(psB, type="response"),
        Pr_Aonly = pAh * (1 - pBh),
        Pr_Bonly = pBh * (1 - pAh)
      ) %>%
      filter((A==1 & B==0) | (A==0 & B==1)) %>%
      mutate(treated = ifelse(A==1, 1L, 0L))
    pAonly <- mean(df$treated==1); pBonly <- 1 - pAonly
    df$weight <- ifelse(df$treated==1, pAonly/df$Pr_Aonly, pBonly/df$Pr_Bonly)
    df$weight[!is.finite(df$weight)] <- NA_real_
    return(df %>% select(L,A,B,O,weight))
  }
  if (method == "newdual") {
    psA <- glm(A ~ L, data = data, family = binomial())
    psB <- glm(B ~ L, data = data, family = binomial())
    df <- data %>% mutate(pAh = predict(psA, type="response"),
                          pBh = predict(psB, type="response"))
    ps_B0_given_A1 <- glm(I(B==0) ~ L, data = df %>% filter(A==1), family = binomial())
    ps_A0_given_B1 <- glm(I(A==0) ~ L, data = df %>% filter(B==1), family = binomial())
    df$pb0_a1 <- NA_real_; df$pa0_b1 <- NA_real_
    df$pb0_a1[df$A==1] <- predict(ps_B0_given_A1, newdata = df[df$A==1,,drop=FALSE], type="response")
    df$pa0_b1[df$B==1] <- predict(ps_A0_given_B1, newdata = df[df$B==1,,drop=FALSE], type="response")
    df <- df %>%
      mutate(Pr_Aonly = pAh * pb0_a1,
             Pr_Bonly = pBh * pa0_b1) %>%
      filter((A==1 & B==0) | (A==0 & B==1)) %>%
      mutate(treated = ifelse(A==1, 1L, 0L))
    pA <- mean(df$treated==1); pB <- 1 - pA
    df$weight <- ifelse(df$treated==1, pA/df$Pr_Aonly, pB/df$Pr_Bonly)
    df$weight[!is.finite(df$weight)] <- NA_real_
    return(df %>% select(L,A,B,O,weight))
  }
  stop("Unknown method")
}

estimate_rr_analysis2 <- function(dat, method, trim_level = 0) {
  dfw <- compute_weights_analysis2(method, dat)
  dfw$weight <- apply_trimming(dfw$weight, trim_level)
  dfw <- dfw %>% filter(is.finite(weight)) %>%
    mutate(group = if_else(B==1 & A==0, "B", "A")) %>%
    transmute(group, O, w = weight)
  # fallback if one arm vanishes after trimming
  nA <- sum(dfw$group=="A" & is.finite(dfw$w) & dfw$w>0)
  nB <- sum(dfw$group=="B" & is.finite(dfw$w) & dfw$w>0)
  if (nA==0 || nB==0) {
    dfw <- compute_weights_analysis2(method, dat) %>%
      filter(is.finite(weight)) %>%
      mutate(group = if_else(B==1 & A==0, "B", "A")) %>%
      transmute(group, O, w = weight)
  }
  ests <- estimate_rr_poisson(dfw)
  ess  <- compute_ess(dfw$w, grp = dfw$group)
  list(est = ests$rr_pois, ESS = ess)
}

# =========================================================
# Monte Carlo runners
# =========================================================
run_mc_strategy0 <- function(n, scenario, type = c("mutual","nonmutual"),
                             strategies = c("cond_poisson","tte_ipw_unrestricted","tte_ipw_subset","tte_ipw_strict"),
                             n_sim = N_SIM_DEFAULT, trim_levels = TRIM_LEVELS) {
  type <- match.arg(type)
  out <- expand.grid(sim=1:n_sim, trim=trim_levels, strategy=strategies) %>%
    mutate(est=NA_real_, bias=NA_real_, ESS_total=NA_real_, ESS_A=NA_real_, ESS_B=NA_real_)
  for (i in 1:n_sim) {
    dat <- generate_data(n, scenario = scenario, type = type)
    for (s in strategies) for (t in trim_levels) {
      res <- tryCatch({
        if (s=="cond_poisson") estimate_conditional_rr(dat)
        else if (s=="tte_ipw_unrestricted") estimate_rr_01(dat, t)
        else if (s=="tte_ipw_subset")       estimate_rr_02(dat, t)
        else                                estimate_rr_03(dat, t)
      }, error=function(e) list(est=NA_real_, ESS=list(ESS_total=NA,ESS_A=NA,ESS_B=NA)))
      idx <- out$sim==i & out$trim==t & out$strategy==s
      out$est[idx] <- res$est
      out$bias[idx] <- res$est - TRUE_MARGINAL_RR
      out$ESS_total[idx] <- res$ESS$ESS_total
      out$ESS_A[idx]     <- res$ESS$ESS_A
      out$ESS_B[idx]     <- res$ESS$ESS_B
    }
  }
  summary <- out %>%
    group_by(strategy, trim) %>%
    summarise(
      MC_bias = mean(bias, na.rm=TRUE),
      MC_sd   = sd(bias, na.rm=TRUE),
      MSE     = mean(bias^2, na.rm=TRUE),
      ESS_total_mean = mean(ESS_total, na.rm=TRUE),
      ESS_A_mean = mean(ESS_A, na.rm=TRUE),
      ESS_B_mean = mean(ESS_B, na.rm=TRUE),
      valid_n = sum(is.finite(est)),
      .groups="drop"
    ) %>% mutate(scenario=paste0("Scenario ", scenario), type=type)
  list(raw=out, summary=summary)
}

run_mc_step3 <- function(n=10000, n_sim=N_SIM_DEFAULT, trim_levels=TRIM_LEVELS) {
  combos <- tidyr::expand_grid(sim=1:n_sim, trim=trim_levels,
                               type=c("mutual","nonmutual"),
                               strat=c("cond_poisson","tte_ipw_unrestricted","tte_ipw_subset","tte_ipw_strict"))
  out <- combos %>% mutate(est=NA_real_, bias=NA_real_, ESS_total=NA_real_, ESS_A=NA_real_, ESS_B=NA_real_)
  for (i in 1:n_sim) for (tp in c("mutual","nonmutual")) {
    dat <- generate_data(n, scenario=2, type=tp)
    for (t in trim_levels) {
      # cond
      r <- estimate_conditional_rr(dat)
      idx <- with(out, sim==i & type==tp & trim==t & strat=="cond_poisson")
      out$est[idx] <- r$est; out$bias[idx] <- r$est - TRUE_MARGINAL_RR
      # 0-1
      r <- estimate_rr_01(dat, t)
      idx <- with(out, sim==i & type==tp & trim==t & strat=="tte_ipw_unrestricted")
      out$est[idx] <- r$est; out$bias[idx] <- r$est - TRUE_MARGINAL_RR
      out$ESS_total[idx] <- r$ESS$ESS_total; out$ESS_A[idx] <- r$ESS$ESS_A; out$ESS_B[idx] <- r$ESS$ESS_B
      if (tp=="nonmutual") {
        # 0-2
        r <- estimate_rr_02(dat, t)
        idx <- with(out, sim==i & type==tp & trim==t & strat=="tte_ipw_subset")
        out$est[idx] <- r$est; out$bias[idx] <- r$est - TRUE_MARGINAL_RR
        out$ESS_total[idx] <- r$ESS$ESS_total; out$ESS_A[idx] <- r$ESS$ESS_A; out$ESS_B[idx] <- r$ESS$ESS_B
        # 0-3
        r <- estimate_rr_03(dat, t)
        idx <- with(out, sim==i & type==tp & trim==t & strat=="tte_ipw_strict")
        out$est[idx] <- r$est; out$bias[idx] <- r$est - TRUE_MARGINAL_RR
        out$ESS_total[idx] <- r$ESS$ESS_total; out$ESS_A[idx] <- r$ESS$ESS_A; out$ESS_B[idx] <- r$ESS$ESS_B
      }
    }
  }
  summ <- out %>% group_by(type, strat, trim) %>%
    summarise(MC_bias=mean(bias,na.rm=TRUE), MC_sd=sd(bias,na.rm=TRUE), MSE=mean(bias^2,na.rm=TRUE),
              ESS_total_mean=mean(ESS_total,na.rm=TRUE), ESS_A_mean=mean(ESS_A,na.rm=TRUE),
              ESS_B_mean=mean(ESS_B,na.rm=TRUE), valid_n=sum(is.finite(est)), .groups="drop")
  list(raw=out, summary=summ)
}

run_mc_step4 <- function(n, scenario, methods, n_sim=N_SIM_DEFAULT, trim_levels=TRIM_LEVELS) {
  out <- expand.grid(sim=1:n_sim, trim=trim_levels, method=methods) %>%
    mutate(est=NA_real_, bias=NA_real_, ESS_total=NA_real_, ESS_A=NA_real_, ESS_B=NA_real_)
  for (i in 1:n_sim) {
    dat <- generate_data_step4(n, scenario)
    for (m in methods) for (t in trim_levels) {
      res <- tryCatch(estimate_rr_analysis2(dat, m, t),
                      error=function(e) list(est=NA_real_, ESS=list(ESS_total=NA,ESS_A=NA,ESS_B=NA)))
      idx <- out$sim==i & out$trim==t & out$method==m
      out$est[idx] <- res$est; out$bias[idx] <- res$est - TRUE_MARGINAL_RR
      out$ESS_total[idx] <- res$ESS$ESS_total; out$ESS_A[idx] <- res$ESS$ESS_A; out$ESS_B[idx] <- res$ESS$ESS_B
    }
  }
  summary <- out %>% group_by(method, trim) %>%
    summarise(MC_bias=mean(bias,na.rm=TRUE), MC_sd=sd(bias,na.rm=TRUE), MSE=mean(bias^2,na.rm=TRUE),
              ESS_total_mean=mean(ESS_total,na.rm=TRUE), ESS_A_mean=mean(ESS_A,na.rm=TRUE),
              ESS_B_mean=mean(ESS_B,na.rm=TRUE), valid_n=sum(is.finite(est)), .groups="drop") %>%
    mutate(scenario=paste0("Scenario ", scenario))
  list(raw=out, summary=summary)
}

# =========================================================
# Plot helpers (keep these for consistency across all figures)
# =========================================================
.add_strategy_label <- function(df, key_col) {
  key_col <- rlang::ensym(key_col)
  key_chr <- as.character(df[[as.character(rlang::as_name(key_col))]])
  lab <- unname(STRATEGY_LABELS[key_chr])
  lab_levels <- unique(unname(STRATEGY_LABELS[names(STRATEGY_LABELS) %in% key_chr]))
  df %>% dplyr::mutate(.key = key_chr, .label = factor(lab, levels = lab_levels))
}
.make_scales_for_labels <- function(labels_factor) {
  labs <- levels(labels_factor)
  inv  <- setNames(names(STRATEGY_LABELS), STRATEGY_LABELS)
  keys <- unname(inv[labs])
  col_vec <- unname(STRATEGY_COLORS[keys]); names(col_vec) <- labs
  lt_vec  <- unname(STRATEGY_LINETYPE[keys]); names(lt_vec) <- labs
  list(colors = col_vec, linetypes = lt_vec, labels = labs)
}

plot_bias_by_trim <- function(df, key_col=c("strategy","method"),
                              title_txt="", ylab="Mean Bias") {
  nm  <- match.arg(key_col)
  df2 <- .add_strategy_label(df, nm)
  sc  <- .make_scales_for_labels(df2$.label)
  ggplot(df2, aes(x=trim*100, y=MC_bias, color=.label, fill=.label, linetype=.label, group=.label)) +
    geom_ribbon(aes(ymin=MC_bias - MC_sd, ymax=MC_bias + MC_sd), alpha=.18, color=NA) +
    geom_line(size=0.5) +
    geom_hline(yintercept=0, linetype="dashed", color="black", linewidth=0.5) +
    scale_color_manual(values=sc$colors, breaks=sc$labels, drop=FALSE) +
    scale_fill_manual(values=sc$colors,  breaks=sc$labels, drop=FALSE) +
    scale_linetype_manual(values=sc$linetypes, breaks=sc$labels, drop=FALSE) +
    labs(x="Trimming (%)", y=ylab, color="Strategy", fill="Strategy", linetype="Strategy", title=title_txt) +
    theme_minimal(base_size=14) +
    theme(legend.position="bottom", legend.box="horizontal")
}

plot_mse_by_trim <- function(df, key_col=c("strategy","method"), title_txt="MSE by trimming") {
  nm  <- match.arg(key_col)
  df2 <- .add_strategy_label(df, nm)
  sc  <- .make_scales_for_labels(df2$.label)
  ggplot(df2, aes(x=trim*100, y=MSE, color=.label, linetype=.label, group=.label)) +
    geom_line(size=0.5) + geom_point(size=1.2) +
    scale_color_manual(values=sc$colors, breaks=sc$labels, drop=FALSE) +
    scale_linetype_manual(values=sc$linetypes, breaks=sc$labels, drop=FALSE) +
    labs(x="Trimming (%)", y="MSE", color="Strategy", linetype="Strategy", title=title_txt) +
    theme_minimal(base_size=14) + theme(legend.position="bottom")
}

plot_ess_by_trim <- function(df, key_col=c("strategy","method"), title_txt="ESS(total) by trimming") {
  nm <- match.arg(key_col)
  df2 <- df %>% group_by(.data[[nm]]) %>% filter(any(is.finite(ESS_total_mean))) %>% ungroup()
  df2 <- .add_strategy_label(df2, nm)
  sc  <- .make_scales_for_labels(df2$.label)
  ggplot(df2, aes(x=trim*100, y=ESS_total_mean, color=.label, linetype=.label, group=.label)) +
    geom_line(size=0.5) + geom_point(size=1.2) +
    scale_color_manual(values=sc$colors, breaks=sc$labels, drop=FALSE) +
    scale_linetype_manual(values=sc$linetypes, breaks=sc$labels, drop=FALSE) +
    labs(x="Trimming (%)", y="ESS (total mean)", color="Strategy", linetype="Strategy", title=title_txt) +
    theme_minimal(base_size=14) + theme(legend.position="bottom")
}

# =========================================================
# Step2: both-rate sensitivity (parameters + runner)
# =========================================================
param_list <- list(
  both70 = list(a= 1.9, b= 1.9, label="Both: 70.2%"),
  both50 = list(a= 1.0, b= 1.0, label="Both: 50.5%"),
  both30 = list(a= 0.1, b= 0.1, label="Both: 30.8%"),
  both10 = list(a=-1.1, b=-1.1, label="Both: 10.7%"),
  both03 = list(a=-2.2, b=-2.2, label="Both: 3.1%")
)

run_mc_step2_oneparam <- function(a, b, n=10000, n_sim=N_SIM_DEFAULT, trim_levels=TRIM_LEVELS) {
  res <- expand.grid(sim=1:n_sim, trim=trim_levels, strategy=c("cond","0-1","0-2","0-3")) %>%
    mutate(est=NA_real_, bias=NA_real_, ESS_total=NA_real_)
  for (i in 1:n_sim) {
    dat <- generate_data_step2(n, a, b)
    for (t in trim_levels) {
      rC <- estimate_conditional_rr(dat)
      res$est[res$sim==i & res$trim==t & res$strategy=="cond"]  <- rC$est
      res$bias[res$sim==i & res$trim==t & res$strategy=="cond"] <- rC$est - TRUE_MARGINAL_RR
      r1 <- estimate_rr_01(dat, t); r2 <- estimate_rr_02(dat, t); r3 <- estimate_rr_03(dat, t)
      for (lab in c("0-1","0-2","0-3")) {
        rr  <- list(`0-1`=r1,`0-2`=r2,`0-3`=r3)[[lab]]
        idx <- res$sim==i & res$trim==t & res$strategy==lab
        res$est[idx]       <- rr$est
        res$bias[idx]      <- rr$est - TRUE_MARGINAL_RR
        res$ESS_total[idx] <- rr$ESS$ESS_total
      }
    }
  }
  res %>%
    group_by(strategy, trim) %>%
    summarise(
      MC_bias        = mean(bias, na.rm=TRUE),
      MC_sd          = sd(bias, na.rm=TRUE),
      MSE            = mean(bias^2, na.rm=TRUE),
      ESS_total_mean = mean(ESS_total, na.rm=TRUE),
      valid_n        = sum(is.finite(est)),
      .groups = "drop"
    )
}

# normalizer to global keys
normalize_step2_keys <- function(df_step2) {
  keymap <- c("cond"="cond_poisson","0-1"="tte_ipw_unrestricted","0-2"="tte_ipw_subset","0-3"="tte_ipw_strict")
  df_step2 %>% mutate(strategy = unname(keymap[strategy]))
}


# =========================================================
# Example runs (toggle printing as needed)
# =========================================================

## Step1 (Scenario 1): mutual vs nonmutual
res_step1_mutual <- run_mc_strategy0(
  n = 10000, scenario = 1, type = "mutual",
  strategies = c("cond_poisson","tte_ipw_unrestricted"),
  n_sim = N_SIM_DEFAULT
)
res_step1_nonmut <- run_mc_strategy0(
  n = 10000, scenario = 1, type = "nonmutual",
  strategies = c("cond_poisson","tte_ipw_unrestricted","tte_ipw_subset","tte_ipw_strict"),
  n_sim = N_SIM_DEFAULT
)
step1_summary <- bind_rows(res_step1_mutual$summary, res_step1_nonmut$summary)

# Figures
# print(plot_bias_by_trim(step1_summary, key_col="strategy",
#                         title_txt="Step 1 (Scenario 1): Mean Bias (±SD) — mutual & nonmutual"))
# print(plot_mse_by_trim(step1_summary, key_col="strategy",
#                        title_txt="Step 1 (Scenario 1): MSE — mutual & nonmutual"))
# print(plot_ess_by_trim(step1_summary, key_col="strategy",
#                        title_txt="Step 1 (Scenario 1): ESS — mutual & nonmutual"))

## Step2: both-rate sensitivity panels
res_step2_list <- purrr::imap(param_list, ~{
  sm <- run_mc_step2_oneparam(.x$a, .x$b)
  sm$both_label <- .x$label
  sm
})
step2_summary <- bind_rows(res_step2_list) %>%
  mutate(both_order = factor(both_label,
                             levels=c("Both: 3.1%","Both: 10.7%","Both: 30.8%","Both: 50.5%","Both: 70.2%")))
step2_summary_norm <- normalize_step2_keys(step2_summary)

# Example: MSE / ESS with unified helpers (facets by both proportion)
# p2_mse <- plot_mse_by_trim(step2_summary_norm, key_col="strategy",
#                            title_txt="Step 2: MSE across both-treated proportions") + facet_wrap(~ both_order)
# p2_ess <- plot_ess_by_trim(step2_summary_norm, key_col="strategy",
#                            title_txt="Step 2: ESS across both-treated proportions") + facet_wrap(~ both_order)
# print(p2_mse); print(p2_ess)

## Step3 (Scenario 2): mutual & nonmutual in one object
res_step3 <- run_mc_step3(n = 10000, n_sim = N_SIM_DEFAULT)
# print(plot_bias_by_trim(res_step3$summary %>% rename(strategy=strat), key_col="strategy",
#                         title_txt="Step 3 (Scenario 2): Mean Bias (±SD) — mutual & nonmutual"))
# print(plot_mse_by_trim(res_step3$summary %>% rename(strategy=strat), key_col="strategy",
#                        title_txt="Step 3 (Scenario 2): MSE — mutual & nonmutual"))
# print(plot_ess_by_trim(res_step3$summary %>% rename(strategy=strat), key_col="strategy",
#                        title_txt="Step 3 (Scenario 2): ESS — mutual & nonmutual"))

## Step4 (Scenarios 1 & 2): six methods on A_only vs B_only
methods_4  <- c("exclusive","ipw","gps_q","gps_spline","dual","newdual")
res_step4_s1 <- run_mc_step4(n = 10000, scenario = 1, methods = methods_4, n_sim = N_SIM_DEFAULT)
res_step4_s2 <- run_mc_step4(n = 10000, scenario = 2, methods = methods_4, n_sim = N_SIM_DEFAULT)

# print(plot_bias_by_trim(res_step4_s1$summary %>% rename(strategy=method), key_col="strategy",
#                         title_txt="Step 4-1 (Scenario 1): Bias (±SD) by trimming"))
# print(plot_mse_by_trim(res_step4_s1$summary %>% rename(strategy=method), key_col="strategy",
#                        title_txt="Step 4-1 (Scenario 1): MSE by trimming"))
# print(plot_ess_by_trim(res_step4_s1$summary %>% rename(strategy=method), key_col="strategy",
#                        title_txt="Step 4-1 (Scenario 1): ESS by trimming"))
# print(plot_bias_by_trim(res_step4_s2$summary %>% rename(strategy=method), key_col="strategy",
#                         title_txt="Step 4-2 (Scenario 2): Bias (±SD) by trimming"))
# print(plot_mse_by_trim(res_step4_s2$summary %>% rename(strategy=method), key_col="strategy",
#                        title_txt="Step 4-2 (Scenario 2): MSE by trimming"))
# print(plot_ess_by_trim(res_step4_s2$summary %>% rename(strategy=method), key_col="strategy",
#                        title_txt="Step 4-2 (Scenario 2): ESS by trimming"))

## =========================================================
## Sample-size sensitivity (Bias only) — minimal utilities
## =========================================================

## Step1 (Scenario 1): mutual / nonmutual
##  - mutual: cond_poisson, tte_ipw_unrestricted
##  - nonmutual: cond_poisson, tte_ipw_unrestricted, tte_ipw_subset, tte_ipw_strict
run_size_bias_step1 <- function(n_vec = c(100, 500, 1000, 3000, 10000),
                                type_vec = c("mutual","nonmutual"),
                                n_sim = 200,
                                trims   = c(0.00, 0.01, 0.03)) {
  out_list <- list()
  for (tp in type_vec) for (n_ in n_vec) {
    strats <- if (tp == "mutual") c("cond_poisson","tte_ipw_unrestricted") else
      c("cond_poisson","tte_ipw_unrestricted","tte_ipw_subset","tte_ipw_strict")
    sm <- run_mc_strategy0(
      n = n_, scenario = 1, type = tp,
      strategies = strats, n_sim = n_sim, trim_levels = trims
    )$summary %>%
      dplyr::transmute(
        n = n_, type = tp, strategy,
        trim, MC_bias, MC_sd, valid_n
      )
    out_list[[paste(tp, n_, sep = "_")]] <- sm
  }
  dplyr::bind_rows(out_list)
}

## Step4: Scenario 1 & 2, six methods
run_size_bias_step4 <- function(n_vec = c(100, 500, 1000, 3000, 10000),
                                scenario_vec = c(1, 2),
                                methods = c("exclusive","ipw","gps_q","gps_spline","dual","newdual"),
                                n_sim = 200,
                                trims = c(0.00, 0.01, 0.03)) {
  out_list <- list()
  for (sc in scenario_vec) for (n_ in n_vec) {
    sm <- run_mc_step4(
      n = n_, scenario = sc, methods = methods,
      n_sim = n_sim, trim_levels = trims
    )$summary %>%
      dplyr::transmute(
        n = n_, scenario = paste0("Scenario ", sc),
        method, trim, MC_bias, MC_sd, valid_n
      )
    out_list[[paste(sc, n_, sep = "_")]] <- sm
  }
  dplyr::bind_rows(out_list)
}

## Optional: very small plotting helper (Bias±SD vs n; log-x)
plot_bias_vs_n_min <- function(df, key_col = c("strategy","method"),
                               facet_by = NULL,
                               title_txt = "Bias ± SD vs n") {
  key_col <- match.arg(key_col)
  nm <- rlang::sym(key_col)
  df %>%
    dplyr::mutate(key = as.character(!!nm),
                  label = STRATEGY_LABELS[key]) %>%
    ggplot2::ggplot(
      ggplot2::aes(x = n, y = MC_bias, color = label,
                   group = interaction(label, trim),
                   linetype = as.factor(trim))
    ) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = MC_bias - MC_sd, ymax = MC_bias + MC_sd, fill = label),
                         alpha = .15, color = NA) +
    ggplot2::geom_line(size = 0.6) +
    ggplot2::geom_point(size = 1.2) +
    ggplot2::scale_x_log10(breaks = c(100, 500, 1000, 3000, 10000)) +
    ggplot2::labs(x="Sample size (log10)", y="Mean Bias (±SD)",
                  color="Strategy", fill="Strategy", linetype="Trimming",
                  title = title_txt) +
    { if (!is.null(facet_by)) ggplot2::facet_wrap(stats::as.formula(paste("~", facet_by))) else ggplot2::geom_blank() } +
    ggplot2::theme_minimal(base_size=14) +
    ggplot2::theme(legend.position="bottom", legend.box="vertical")
}


## ===== Step1: Bias-only size sensitivity (Scenario 1) =====
size_bias_step1 <- run_size_bias_step1(
  n_vec = c(100, 500, 1000, 3000, 10000),
  type_vec = c("mutual","nonmutual"),
  n_sim = 200,
  trims = c(0.00, 0.01, 0.03)
)
# For plotting/saving (optional):
# print(plot_bias_vs_n_min(size_bias_step1, key_col="strategy",
#                          facet_by = "type",
#                          title_txt = "Step 1: Bias ± SD vs n (Scenario 1)"))

## ===== Step4: Bias-only size sensitivity (Scenarios 1 & 2) =====
size_bias_step4 <- run_size_bias_step4(
  n_vec = c(100, 500, 1000, 3000, 10000),
  scenario_vec = c(1, 2),
  methods = c("exclusive","ipw","gps_q","gps_spline","dual","newdual"),
  n_sim = 200,
  trims = c(0.00, 0.01, 0.03)
)
# print(plot_bias_vs_n_min(size_bias_step4 %>% dplyr::rename(strategy = method),
#                          key_col="strategy", facet_by = "scenario",
#                          title_txt = "Step 4: Bias ± SD vs n (Scenarios 1 & 2)"))


# =========================================================
# Styling note:
#  - To harmonize final colors/linetypes across all figures,
#    edit STRATEGY_COLORS and STRATEGY_LINETYPE above.
#  - If tiny x-offsets are desired for overplotted markers,
#    add `x = trim*100 + <small offset>` in the ggplot aes().
# =========================================================
