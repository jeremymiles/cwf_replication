

knitr::opts_chunk$set(echo = FALSE)
suppressMessages({
  library(dplyr)
  library(ggplot2)
  library(glue)
  library(lavaan)
  library(moments)
  library(parallel)
  library(psych)
  library(rlecuyer)
  library(semptools)
  library(semPlot)
  library(simsem)
  library(snow)
})


# Set up 3 SimDataDist objects, each with 9 variables, and with skew and 
# kurtosis of 0, 0 (normal), 2, 7 and 3, 21, representing normal, moderately
# non-normal and severely non-normal 
dist_1 <- simsem::bindDist(skewness = 0, kurtosis = 0)
dist_2 <- simsem::bindDist(skewness = 2, kurtosis = 7)
dist_3 <- simsem::bindDist(skewness = 3, kurtosis = 21)

# Data generating process model without cross loadings and 
# without variances (because these vary across the two models)
dgp_model_no_cross_loadings <- '
  # data used for models 1 and 2
  # loadings
  f1 =~ 0.7 * y1 + 0.7 * y2 + 0.7 * y3
  f2 =~ 0.7 * y4 + 0.7 * y5 + 0.7 * y6
  f3 =~ 0.7 * y7 + 0.7 * y8 + 0.7 * y9
  # factor variances
  f1 ~~ 1 * f1
  f2 ~~ 1 * f2
  f3 ~~ 1 * f3
  # factor covariances
  f1 ~~ 0.3 * f2
  f1 ~~ 0.3 * f3
  f2 ~~ 0.3 * f3
  # measured variable error variances
  y1 ~~ 0.51 * y1
  y2 ~~ 0.51 * y2
  y3 ~~ 0.51 * y3
  y4 ~~ 0.51 * y4
  y5 ~~ 0.51 * y5
  y6 ~~ 0.51 * y6
  y7 ~~ 0.51 * y7
  y8 ~~ 0.51 * y8
  y9 ~~ 0.51 * y9
  '



dgp_model_cross_loadings <- 
  glue::glue(
    dgp_model_no_cross_loadings,
    '# cross loadings
    ',
    'f2 =~ 0.35 * y7
    ',
    'f3 =~ 0.35 * y6
    '
  )


# Model specification 1
model_1 <-
  '# No cross loadings
  f1 =~ y1 + y2 + y3
  f2 =~ y4 + y5 + y6
  f3 =~ y7 + y8 + y9
  f1 ~~ f2
  f1 ~~ f3
  f2 ~~ f3
  '

# Model specification 2
model_2 <-
  # misspecification of inclusion
  glue::glue(
    model_1,
    # add params that are zero in pop
    'f2 =~ y8
    ',
    'f3 =~ y5'
  )

# Model specification 3
model_3 <-
  # misspecification of exclusion
  # same as model 1 but fit to dgp_cross_loadings
  model_1

# Model specification 4
model_4 <-
  # misspecification of exclusion AND inclusion
  # same as model 2 but fit to dgp_cross_loadings
  model_2






GenerateDFs <- function(n) {
  # Function to generate raw data 
  # Returns list of data frames.
  # 2 datasets (data_1 without cross loadings, data_2 with)
  # with 3 levels of non-normality (1:none, 2:moderate, 3:extreme)
  data_1_dist_1 <-
    # normal
    simsem::generate(
      model = dgp_model_no_cross_loadings, n = n,  indDist=dist_1
    )
  data_1_dist_2 <-
    # low skew and kurtosis
    simsem::generate(
      model = dgp_model_no_cross_loadings, n = n,  indDist=dist_2)
  data_1_dist_3 <-
    # low skew and kurtosis
    simsem::generate(
      model = dgp_model_no_cross_loadings, n = n,  indDist=dist_3)
  data_2_dist_1 <-
    # normal
    simsem::generate(
      model = dgp_model_cross_loadings, n = n,  indDist=dist_1)
  data_2_dist_2 <-
    # low skew and kurtosis
    simsem::generate(
      model = dgp_model_cross_loadings, n = n,  indDist=dist_2
    )
  data_2_dist_3 <-
    # low skew and kurtosis
    simsem::generate(
      model = dgp_model_cross_loadings, n = n,  indDist=dist_3
    )
  return(
    list(
      data_1_dist_1 = data_1_dist_1,
      data_1_dist_2 = data_1_dist_2,
      data_1_dist_3 = data_1_dist_3,
      data_2_dist_1 = data_2_dist_1,
      data_2_dist_2 = data_2_dist_2,
      data_2_dist_3 = data_2_dist_3
    )
  )
}

RunOneSimulation <- function(n) {
  # RunOneSimulation runs one simulation at the specified sample size for all 
  # 4 models, 3 distributions, and 3 estimators
  dfs <- GenerateDFs(n)
  
  # fit models
  # Model 1 - correct / correct
  
  list_res <- list()
  
  temp_fit <- lavaan::cfa(
    model = model_1, data = dfs$data_1_dist_1, estimator = "MLM"
  ) 
  list_res$fit_m1_dist1_ml <- GetChi(temp_fit, scaled = FALSE)
  list_res$fit_m1_dist1_sb <- GetChi(temp_fit, scaled = TRUE)
  
  temp_fit <- lavaan::cfa(
    model = model_1, data = dfs$data_1_dist_2, estimator = "MLM"
  )
  list_res$fit_m1_dist2_ml <- GetChi(temp_fit, scaled = FALSE)
  list_res$fit_m1_dist2_sb <- GetChi(temp_fit, scaled = TRUE)
  
  
  temp_fit <- lavaan::cfa(
    model = model_1, data = dfs$data_1_dist_3, estimator = "MLM"
  )
  list_res$fit_m1_dist3_ml <- GetChi(temp_fit, scaled = FALSE)
  list_res$fit_m1_dist3_sb <- GetChi(temp_fit, scaled = TRUE)
  
  
  list_res$fit_m1_dist1_wls <- lavaan::cfa(
    model = model_1, data = dfs$data_1_dist_1, estimator = "WLS"
  )%>% GetChi()
  
  list_res$fit_m1_dist2_wls <- lavaan::cfa(
    model = model_1, data = dfs$data_1_dist_2, estimator = "WLS"
  ) %>% GetChi()
  
  list_res$fit_m1_dist3_wls <- lavaan::cfa(
    model = model_1, data = dfs$data_1_dist_3, estimator = "WLS"
  ) %>% GetChi()
  
  # model 2
  temp_fit <- lavaan::cfa(
    model = model_2, data = dfs$data_1_dist_1, estimator = "MLM"
  ) 
  list_res$fit_m2_dist1_ml <- GetChi(temp_fit, scaled = FALSE)
  list_res$fit_m2_dist1_sb <- GetChi(temp_fit, scaled = TRUE)
  
  
  
  temp_fit <- lavaan::cfa(
    model = model_2, data = dfs$data_1_dist_2, estimator = "MLM"
  )
  list_res$fit_m2_dist2_ml <- GetChi(temp_fit, scaled = FALSE)
  list_res$fit_m2_dist2_sb <- GetChi(temp_fit, scaled = TRUE)
  
  
  temp_fit <- lavaan::cfa(
    model = model_2, data = dfs$data_1_dist_3, estimator = "MLM"
  )
  list_res$fit_m2_dist3_ml <- GetChi(temp_fit, scaled = FALSE)
  list_res$fit_m2_dist3_sb <- GetChi(temp_fit, scaled = TRUE)
  
  
  list_res$fit_m2_dist1_wls <- lavaan::cfa(
    model = model_2, data = dfs$data_1_dist_1, estimator = "WLS"
  ) %>% GetChi()
  
  
  list_res$fit_m2_dist2_wls <- lavaan::cfa(
    model = model_2, data = dfs$data_1_dist_2, estimator = "WLS"
  ) %>% GetChi()
  
  list_res$fit_m2_dist3_wls <- lavaan::cfa(
    model = model_2, data = dfs$data_1_dist_3, estimator = "WLS"
  ) %>% GetChi()
  
  # model 3
  temp_fit <- lavaan::cfa(
    model = model_3, data = dfs$data_2_dist_1, estimator = "MLM"
  )
  list_res$fit_m3_dist1_ml <- GetChi(temp_fit, scaled = FALSE)
  list_res$fit_m3_dist1_sb <- GetChi(temp_fit, scaled = TRUE)
  
  
  temp_fit <- lavaan::cfa(
    model = model_3, data = dfs$data_2_dist_2, estimator = "MLM"
  ) 
  list_res$fit_m3_dist2_ml <- GetChi(temp_fit, scaled = FALSE)
  list_res$fit_m3_dist2_sb <- GetChi(temp_fit, scaled = TRUE)
  
  temp_fit <- lavaan::cfa(
    model = model_3, data = dfs$data_2_dist_3, estimator = "MLM"
  ) 
  list_res$fit_m3_dist3_ml <- GetChi(temp_fit, scaled = FALSE)
  list_res$fit_m3_dist3_sb <- GetChi(temp_fit, scaled = TRUE)
  
  
  
  list_res$fit_m3_dist1_wls <- lavaan::cfa(
    model = model_3, data = dfs$data_2_dist_1, estimator = "WLS"
  ) %>% GetChi()
  
  list_res$fit_m3_dist2_wls <- lavaan::cfa(
    model = model_3, data = dfs$data_2_dist_2, estimator = "WLS"
  ) %>% GetChi()
  
  list_res$fit_m3_dist3_wls <- lavaan::cfa(
    model = model_3, data = dfs$data_2_dist_3, estimator = "WLS"
  ) %>% GetChi()
  
  # model 4
  temp_fit <- lavaan::cfa(
    model = model_4, data = dfs$data_2_dist_1, estimator = "MLM"
  ) 
  list_res$fit_m4_dist1_ml <- GetChi(temp_fit, scaled = FALSE)
  list_res$fit_m4_dist1_sb <- GetChi(temp_fit, scaled = TRUE)
  
  temp_fit <- lavaan::cfa(
    model = model_4, data = dfs$data_2_dist_2, estimator = "MLM"
  )  
  list_res$fit_m4_dist2_ml <- GetChi(temp_fit, scaled = FALSE)
  list_res$fit_m4_dist2_sb <- GetChi(temp_fit, scaled = TRUE)
  
  list_res$fit_m4_dist3_ml <- lavaan::cfa(
    model = model_4, data = dfs$data_2_dist_3, estimator = "MLM"
  ) 
  list_res$fit_m4_dist3_ml <- GetChi(temp_fit, scaled = FALSE)
  list_res$fit_m4_dist3_sb <- GetChi(temp_fit, scaled = TRUE)
  
  
  list_res$fit_m4_dist1_wls <- lavaan::cfa(
    model = model_4, data = dfs$data_2_dist_1, estimator = "WLS"
  ) %>% GetChi()
  
  list_res$fit_m4_dist2_wls <- lavaan::cfa(
    model = model_4, data = dfs$data_2_dist_2, estimator = "WLS"
  ) %>% GetChi()
  
  list_res$fit_m4_dist3_wls <- lavaan::cfa(
    model = model_4, data = dfs$data_2_dist_3, estimator = "WLS"
  ) %>% GetChi()
  
  return(list_res)  
}



GetChi <- function(fitted_model, scaled = FALSE) {
  # Convenience function to extract fit indices from lavaan object
  # if model is invalid (negative variances, non-convergence)
  # returns NA
  if (!lavTech(fitted_model, "converged")) {
    # Check for convergence
    return(NA)
  }
  if (
    # check observed variables have positive residual variance
    (lavTech(fitted_model, "coef")$theta |> 
     diag() |> 
     min()) <= 0 
  ) {
    return(NA)
  }
  if (!scaled) {
    fit_stats <- fitMeasures(fitted_model)[c("chisq", "df", "pvalue")]
    return(fit_stats)
    
  }
  if (scaled) {
    fit_stats <- fitMeasures(fitted_model)[c("chisq.scaled", "df.scaled", "pvalue.scaled")]
    names(fit_stats) <- c("chisq", "df", "pvalue")
    return(fit_stats)
  }
}


GetMeanChi <- function(res, model) {
  # Function to take list of objects and return mean chi-square, proportion  of 
  # converged models, and rejection rate (based on p < 0.05)
  df_all_res <- lapply(res, function(x) {
    if(is.na(x[[model]][[1]])) { return(c(chisq = NA, pvalue = NA)) }
    return(c(chisq = x[[model]][["chisq"]], pvalue = x[[model]][["pvalue"]] ))
  }) %>% dplyr::bind_rows()
  df_all_res$reject <- df_all_res$pvalue < 0.05
  
  return(
    list(
      mean_chi = mean(df_all_res$chisq, na.rm = TRUE),
      prop_converged = sum(!is.na(df_all_res$pvalue)),
      reject_rate = mean(df_all_res$reject, na.rm = TRUE)
    )
  )
}
