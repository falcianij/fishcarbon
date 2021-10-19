#' Calculate Optimal & Suboptimal Spectra Within Interval
#'
#' This function is used to find optimal and suboptimal average spectra within the provided interval of the biomass vs yield plot.
#'
#' @param strategies_df The dataframe of all strategy metrics computed
#'
#' @param sim0 The simulation object of the unfished scenario
#'
#' @param lower.break The lower bound of the interval in yield
#'
#' @param upper.break The upper bound of the interval in yield
#'
#' @return A dataframe of weights and the normalized size spectrum for the optimal & suboptimal strategies
#'
#' @export

BvY_intervals <- function(strategies_df, sim0, lower.break, upper.break)
{
  percentage.taken <- 0.2
  last.years <- 100
  params <- sim0@params

  # Creating interval
  strategies_df.subset <- filter(strategies_df,
                                 (Y > lower.break) & (Y < upper.break))

  # Getting top and bottom 20% of strategies in interval
  strategies_df.subset_sub <- strategies_df.subset %>%
    arrange(residuals) %>%
    filter(row_number() < nrow(strategies_df.subset) * percentage.taken)

  strategies_df.subset_opt <- strategies_df.subset %>%
    arrange(desc(residuals)) %>%
    filter(row_number() < nrow(strategies_df.subset) * percentage.taken)

  # Getting average stats to calculate average strategies
  strategies_df.subset_sub <- strategies_df.subset_sub %>%
    summarize(F_small = mean(F_small),
              F_med = mean(F_med),
              F_large = mean(F_large),
              S_small = mean(S_small),
              S_med = mean(S_med),
              S_large = mean(S_large))

  strategies_df.subset_opt<- strategies_df.subset_opt %>%
    summarize(F_small = mean(F_small),
              F_med = mean(F_med),
              F_large = mean(F_large),
              S_small = mean(S_small),
              S_med = mean(S_med),
              S_large = mean(S_large))

  # Run suboptimal strategy
  gear_params(params)$sel <- rep(as.numeric(strategies_df.subset_sub[1, 4:6]), each = 3)
  gear_params(params)$knife_edge_size <- species_params(params)$w_inf * gear_params(params)$sel
  sim1_sub <- project(params, t_max = 300, effort = c("small" = strategies_df.subset_sub$F_small,
                                                      "medium" = strategies_df.subset_sub$F_med,
                                                      "large" = strategies_df.subset_sub$F_large))

  # Run optimal strategy
  gear_params(params)$sel <- rep(as.numeric(strategies_df.subset_opt[1, 4:6]), each = 3)
  gear_params(params)$knife_edge_size <- species_params(params)$w_inf * gear_params(params)$sel
  sim1_opt <- project(params, t_max = 300, effort = c("small" = strategies_df.subset_opt$F_small,
                                                      "medium" = strategies_df.subset_opt$F_med,
                                                      "large" = strategies_df.subset_opt$F_large))

  # Relative size spectra
  NvsW0 <- colSums(apply(tail(N(sim0), last.years), c(2,3), mean))
  NvsW1_sub <- colSums(apply(tail(N(sim1_sub), last.years), c(2,3), mean))
  NvsW1_opt <- colSums(apply(tail(N(sim1_opt), last.years), c(2,3), mean))

  relative.NvsW_sub <- NvsW1_sub / NvsW0
  relative.NvsW_opt <- NvsW1_opt / NvsW0

  relative.NvsW_sub_df <- data.frame(w = w(params), relative.N = relative.NvsW_sub) %>%
    mutate(type = "suboptimal")
  relative.NvsW_opt_df <- data.frame(w = w(params), relative.N = relative.NvsW_opt) %>%
    mutate(type = "optimal")
  relative.NvsW_df <- bind_rows(relative.NvsW_sub_df, relative.NvsW_opt_df)

  return(relative.NvsW_df)
}
