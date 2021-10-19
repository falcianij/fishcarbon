#' Calculate Simulation Metrics
#'
#' These functions are used to calculate size spectrum metrics.
#'
#' @param sim0 The simulation object of the unfished scenario
#'
#' @param sim1 The simulation object of the fished scenario
#'
#' @param stratID A unique strategy identification number to be attached to the calculated metrics
#'
#' @param metrics0 The metrics as computed by `unharvestedMetrics()` associated with the unfished scenario
#'
#' @return A matrix of the input controls and the calculated community properties
#'
#' @export

unharvestedMetrics <- function(sim0)
{
  # Last years that are averaged
  last.years <- 100

  # Unharvested strategy identification number
  stratID <- 1

  # Misc. calculations that get used again below
  biomass <- tail(sweep(sim0@n, 3, sim0@params@w * sim0@params@dw, "*"), last.years)
  Btotal <- mean(apply(biomass, 1, sum))
  slope_b <- total_n <- apply(tail(sweep(apply(sim0@n, c(1,3), sum), 2, sim0@params@w, "*"), 100), 2, mean)
  summary_fit <- summary(lm(log(slope_b) ~ log(sim0@params@w)))

  # Input controls
  F_vec <- getEffort(sim0)[1,]
  S_vec <- gear_params(sim0@params)$sel[seq(1, 9, 3)]
  #%>% group_by(gear) %>% summarize(sel = first(sel)) %>% .$sel

  # Species-specific metrics
  SSB_vec <- apply(tail(getSSB(sim0), last.years), 2, mean)

  # Class-specific metrics
  E_vec <- rowsum(as.numeric((SSB_vec / SSB_vec) < 0.1), rep(1:3, each = 3))
  Y_vec <- apply(apply(sweep(tail(getFMortGear(sim0), 100), c(1,3,4), biomass, "*"), c(1,2), sum), 2, mean) / Btotal

  # Community metrics
  B <- Btotal / Btotal
  Bflux <- mean(apply(sweep(biomass, c(2,3), tail(sim0@params@mu_b, last.years), "*"), 1, sum)) / Btotal
  Y <- sum(Y_vec)
  E <- sum(E_vec)
  Wmean <- mean(tail(getMeanWeight(sim0), last.years))
  Winfmax <- mean(tail(getMeanMaxWeight(sim0), last.years)[,2])
  propLarge <- mean(tail(getProportionOfLargeFish(sim0, threshold_w = species_params(sim0@params)$w_inf[6]), 100))
  slope <- summary_fit$coefficients[2,1]
  intercept <- summary_fit$coefficients[1,1]
  r2 <- summary_fit$r.squared

  # Metrics compilation
  metrics <- c(stratID,
               F_vec, S_vec,
               B, Bflux,
               Y, Y_vec,
               E, E_vec,
               Wmean, Winfmax, propLarge, slope, intercept, r2,
               SSB_vec, Btotal)
  metrics <- matrix(metrics, nrow = 1)
  colnames(metrics) <- c("stratID",
                         "F_small", "F_med", "F_large", "S_small", "S_med", "S_large",
                         "B", "Bflux",
                         "Y", "Y_small", "Y_med", "Y_large",
                         "E", "E_small", "E_med", "E_large",
                         "Wmean", "Winfmax", "propLarge", "slope", "intercept", "r2",
                         "SSB_1", "SSB_2", "SSB_3", "SSB_4", "SSB_5", "SSB_6", "SSB_7", "SSB_8", "SSB_9", "Btotal")

  return(metrics)
}


strategyMetrics <- function(sim1, stratID, metrics0)
{
  # Last years that are averaged
  last.years <- 100

  # Misc. calculations that get used again below
  biomass <- tail(sweep(sim1@n, 3, sim1@params@w * sim1@params@dw, "*"), last.years)
  Btotal <- mean(apply(biomass, 1, sum))
  slope_b <- total_n <- apply(tail(sweep(apply(sim1@n, c(1,3), sum), 2, sim1@params@w, "*"), 100), 2, mean)
  summary_fit <- summary(lm(log(slope_b) ~ log(sim1@params@w)))

  # Input controls
  F_vec <- getEffort(sim1)[1,]
  S_vec <- gear_params(sim1@params)$sel[seq(1, 9, 3)]

  # Species-specific metrics
  SSB_vec <- apply(tail(getSSB(sim1), last.years), 2, mean)

  # Class-specific metrics
  E_vec <- rowsum(as.numeric((SSB_vec / metrics0[24:32]) < 0.1), rep(1:3, each = 3))
  Y_vec <- apply(apply(sweep(tail(getFMortGear(sim1), 100), c(1,3,4), biomass, "*"), c(1,2), sum), 2, mean) / metrics0[,"Btotal"]

  # Community metrics
  B <- Btotal / metrics0[,"Btotal"]
  Bflux <- mean(apply(sweep(biomass, c(2,3), tail(sim1@params@mu_b, last.years), "*"), 1, sum)) / metrics0[,"Btotal"]
  Y <- sum(Y_vec)
  E <- sum(E_vec)
  Wmean <- mean(tail(getMeanWeight(sim1), last.years))
  Winfmax <- mean(tail(getMeanMaxWeight(sim1), last.years)[,2])
  propLarge <- mean(tail(getProportionOfLargeFish(sim1, threshold_w = species_params(sim1@params)$w_inf[6]), 100))
  slope <- summary_fit$coefficients[2,1]
  intercept <- summary_fit$coefficients[1,1]
  r2 <- summary_fit$r.squared

  # Metrics compilation
  metrics <- c(stratID,
               F_vec, S_vec,
               B, Bflux,
               Y, Y_vec,
               E, E_vec,
               Wmean, Winfmax, propLarge, slope, intercept, r2)
  metrics <- matrix(metrics, nrow = 1)
  colnames(metrics) <- c("stratID",
                         "F_small", "F_med", "F_large", "S_small", "S_med", "S_large",
                         "B", "Bflux",
                         "Y", "Y_small", "Y_med", "Y_large",
                         "E", "E_small", "E_med", "E_large",
                         "Wmean", "Winfmax", "propLarge", "slope", "intercept", "r2")

  return(metrics)
}
