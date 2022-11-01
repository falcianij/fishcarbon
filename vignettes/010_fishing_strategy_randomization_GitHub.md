Fishing Strategy Randomization
================
Jonathan E. Falciani[^1]
01 November 2022

# Aim and Scope

Here we develop the code that generates the data we want to analyze. We
will create a fish community consisting of 9 species that each differ by
their asymptotic weight $W_\infty$. These 9 species will be divided into
3 different classes (these can be changed in the `initialize simulation`
chunk in the Simulations section of this document):

-   3 *small* species ($W_\infty < 100$ g),
-   3 *medium* species ($100$g $< W_\infty < 1.5$ kg), and
-   3 *large* species ($1.5$ kg $< W_\infty$).

Each of the 3 classes will be harvested on a species-specific basis,
however each species within a class will share the same harvesting
parameters $S$ and $F$. Each species will be selected by a knife-edge
selectivity function harvesting all fish at or above a proportion of
$W_\infty$, called *selectivity* $S$ ($0 < S < 1$). The selected fish of
each species will be fished at a *fishing mortality* $F$ ($0 < F < 3$).
Since $S$ and $F$ are shared by the species within a class, there are 6
independent variables. Each of these 6 variables will be uniformly
sampled from their predesignated ranges shown above.

-   $S_{S}$ and $F_{S}$
-   $S_{M}$ and $F_{M}$
-   $S_{L}$ and $F_{L}$.

Each run will allow the simulation to run for 300 years to guarantee it
reaches a steady state. Each run will collect information on biomass
$B$, yield $Y$, extirpations $E$, biomass flux $B_{flux}$, and average
size $W_{avg}$.

First we need to load in the `mizer`
[package](https://sizespectrum.org/mizer).

``` r
library(mizer)
```

The other R packages needed are:

``` r
library(tidyverse)
library(here)
```

Run `devtools::load_all` in order to load all of the necessary functions
inside the `fishcarbon` R project. If you are not, then see the
functions vignette for the functions needed.

``` r
devtools::load_all("./", export_all = TRUE)
```

$~$

# Simulations

Each strategy will be given a unique `stratID`. Since we aim to generate
upwards of 10,000 strategies which each take \~30 seconds to run, we
will chunk this code into parts where we only generate a small subset of
the 10,000 strategies as defined by `no.simulations`. This will be added
to the `fishcarbon_metrics_9sp_workingcopy.csv` file every time it is
run.

If you are running these files for the first time, there will be only 1
simulation in `fishcarbon_metrics_9sp_workingcopy.csv` for
initialization. If you have run this code before, there should be
additional simulations.

``` r
strategies_df <- read.csv(file.path(here(), 'data_out',
                                    'fishcarbon_metrics_9sp_workingcopy.csv'),header = TRUE)
```

This is where we specify the initial parameters for each strategy. We
specify the number of fish in the ecosystem with `no.fish`, the number
of fisheries (groups) with `no.fishery`, and the number of fish in each
fishery with `no.fish.fishery`. We also specify how many strategies we
want to generate during this run with `no.simulations`.

``` r
no.fish = 9
no.fishery = 3
no.fish.fishery = no.fish / no.fishery

no.simulations = 9
```

## Unharvested Simulation

This generates the unharvested control strategy. It then plots the size
spectrum, feeding level vs size, mortality vs size, biomass vs time, and
reproduction.

``` r
# Trait-based model parameters
params <- newTraitParams(no_sp = no.fish)

# Experimental variables
S.list <- rep(.01, no.fishery)
F.list <- rep(0, no.fishery)
  
gear_params(params) <- data.frame(species = species_params(params)$species,
                                  gear = c(rep("small", no.fish.fishery), rep("medium", no.fish.fishery), rep("large", no.fish.fishery)),
                                  sel_func = "knife_edge",
                                  sel = rep(S.list, each = no.fish.fishery),
                                  knife_edge_size = species_params(params)$w_inf * rep(S.list, each = no.fish.fishery))

# Simulation run
sim0 <- project(params, t_max = 300, effort = c("small" = F.list[1],
                                                "medium" = F.list[2],
                                                "large" = F.list[3]))
plot(sim0)
```

![](C:/Users/jonat/AppData/Local/Temp/Rtmpc1pAhR/preview-7314ffed85.dir/010_fishing_strategy_randomization_GitHub_files/figure-gfm/unharvested%20simulation-1.png)<!-- -->
![mizer
plots](/figures/010_fishing_strategy_randomization/unharvested-simulation.png)

## Unharvested Metrics

We will save the unharvested metrics for comparison to the experimental
strategies. For details, see `metrics.R`.

``` r
metrics0 <- unharvestedMetrics(sim0, no.fish, no.fishery)
metrics_full <- metrics0[, -24:-length(metrics0)]
```

## Simulation Loop

Now we will loop through randomized strategies that are compared to the
metrics of the unharvested strategy `metrics0`. This will compute
`counter` strategies that will be assigned a unique `stratID` starting
from wherever `strategies_df` left off.

``` r
counter <- no.simulations
strat.start <- max(strategies_df$stratID)

metrics_full_init <- matrix(0,
                            nrow = counter,
                            ncol = length(metrics_full))
metrics_full <- rbind(metrics_full, metrics_full_init)

for (count in 2:(1+counter))
{
  strat <- strat.start + count - 1
  # Experimental variables
  F.list <- runif(no.fishery, 0, 3)
  S.list <- runif(no.fishery, 0, 1)

  # Trait-based model parameters
  gear_params(params)$sel <- rep(S.list, each = no.fish.fishery)
  gear_params(params)$knife_edge_size <- species_params(params)$w_inf * gear_params(params)$sel

  # Simulation run
  sim1 <- project(params, t_max = 300, effort = c("small" = F.list[1],
                                                  "medium" = F.list[2],
                                                  "large" = F.list[3]))
  
  # Add metrics
  metrics_full[count,] <- strategyMetrics(sim1, strat, metrics0, no.fish, no.fishery)
  print(paste0(count, ", ", strat))
  #fprintf("%d, %d\n", count, strat)
}
```

    ## [1] "2, 1"
    ## [1] "3, 2"
    ## [1] "4, 3"
    ## [1] "5, 4"
    ## [1] "6, 5"
    ## [1] "7, 6"
    ## [1] "8, 7"
    ## [1] "9, 8"
    ## [1] "10, 9"

To write to the csv file, we convert the matrix to a dataframe. We then
append it to the `strategies_df` so we can add more in the future.

``` r
metrics_df <- as.data.frame(metrics_full[-1,])
strategies_df <- bind_rows(strategies_df, metrics_df)

rownames(strategies_df) <- NULL
```

| stratID |   F_small |     F_med |   F_large |   S_small |     S_med |   S_large |         B |     Bflux |         Y |   Y_small |     Y_med |   Y_large |   E | E_small | E_med | E_large |     Wmean |   Winfmax | propLarge |      slope |  intercept |        r2 |
|--------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----:|--------:|------:|--------:|----------:|----------:|----------:|-----------:|-----------:|----------:|
|       0 | 1.1884765 | 0.7387810 | 0.5656855 | 0.8227803 | 0.2800604 | 0.4049423 | 0.8614603 | 0.1965383 | 0.1209249 | 0.0656374 | 0.0000034 | 0.0552841 |   4 |       0 |     3 |       1 | 0.6766237 |  766.1184 | 0.1547663 | -0.8362179 |  -7.899007 | 0.8019967 |
|       1 | 1.7905948 | 0.7144968 | 0.0978281 | 0.3756304 | 0.7830636 | 0.7307587 | 1.0564039 | 0.2924925 | 0.1340808 | 0.1009541 | 0.0296594 | 0.0034673 |   2 |       2 |     0 |       0 | 1.2229067 | 2105.5374 | 0.2402056 | -0.6237987 |  -7.646785 | 0.9026944 |
|       2 | 1.4944723 | 0.7206439 | 2.5495914 | 0.7741052 | 0.0251570 | 0.9655723 | 0.9352614 | 0.2384810 | 0.0716186 | 0.0679541 | 0.0000000 | 0.0036645 |   4 |       0 |     3 |       1 | 0.7398750 | 1239.7095 | 0.3778454 | -0.8799255 |  -8.027288 | 0.7267015 |
|       3 | 2.1480937 | 1.6004738 | 1.3065030 | 0.3490952 | 0.4350357 | 0.1893396 | 0.8792182 | 0.2276708 | 0.5870231 | 0.1016397 | 0.4853833 | 0.0000000 |   3 |       0 |     0 |       3 | 1.3976975 |  247.3081 | 0.0000000 | -2.8557578 | -12.560386 | 0.5207438 |
|       4 | 0.9604125 | 2.8031479 | 0.1393625 | 0.3912340 | 0.3629795 | 0.9333549 | 0.9414770 | 0.2534735 | 0.2141240 | 0.2129327 | 0.0000000 | 0.0011913 |   5 |       1 |     3 |       1 | 0.8168442 | 1254.6054 | 0.4249110 | -0.9759692 |  -8.280261 | 0.6321033 |
|       5 | 2.2132553 | 1.6928814 | 2.5413713 | 0.0151682 | 0.4963109 | 0.1032635 | 0.9506646 | 0.2487359 | 0.4819394 | 0.0000000 | 0.4819394 | 0.0000000 |   6 |       3 |     0 |       3 | 3.0058695 |  257.7873 | 0.0000000 | -7.2041594 | -22.191685 | 0.4483129 |

$~$

# Plots

``` r
ggplot(strategies_df, aes(x=Y, y=B, color=Wmean )) +
  geom_point() +
  geom_hline(yintercept = 1, linetype = "dotted") +
  ggtitle("Average annual yield versus average annual biomass at equilibrium") +
  ylab("Average annual biomass (g)") + xlab("Average annual yield (g)") +
  scale_color_gradient(low = "blue", high = "red")
```

![](C:/Users/jonat/AppData/Local/Temp/Rtmpc1pAhR/preview-7314ffed85.dir/010_fishing_strategy_randomization_GitHub_files/figure-gfm/biomass%20vs%20yield%20plot-1.png)<!-- -->
![BvY
plot](/figures/010_fishing_strategy_randomization/biomass-vs-yield-plot.png)

``` r
library(corrplot)
M = cor(strategies_df)
testRes = cor.mtest(strategies_df, conf.level = 0.95)
corrplot(M, method = 'square', type = 'lower', p.mat = testRes$p, sig.level = 0.10)
```

![](C:/Users/jonat/AppData/Local/Temp/Rtmpc1pAhR/preview-7314ffed85.dir/010_fishing_strategy_randomization_GitHub_files/figure-gfm/correlation-1.png)<!-- -->
![corr
plots](/figures/010_fishing_strategy_randomization/correlation.png)

[^1]: Gulf of Maine Research Institute, Temple University. Currently:
    Technical University of Denmark. <jonathan.e.falciani@gmail.com>
