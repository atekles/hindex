context("Simulating h-index")
library(foreach)

test_that("initialization", {

  set.seed(12345)

  n <- 20
  dcitationsSpeed <- 2  # in common log log notation: beta
  dcitationsPeak <- 3 # t_peak, i.e. the age at which expected citations are max
  dcitationsMean <- 2 # max expected citations in a certain year
  # (expected citations in the best year)
  dcitationsDispersion <- 1.1

  dcitationsAlpha <-  # the alpha in common log log notation
    dcitationsPeak /
    (((dcitationsSpeed - 1) / (dcitationsSpeed + 1)) ^ (1 / dcitationsSpeed))
  dcitationsLoglogFactor <-
    dcitationsMean / (((dcitationsSpeed / dcitationsAlpha) *
                         ((dcitationsPeak / dcitationsAlpha) ^
                            (dcitationsSpeed - 1))) /
                        ((1 + (dcitationsPeak / dcitationsAlpha) ^
                            dcitationsSpeed) ^ 2))

  simulationData <- setup_simulation(n = n, boost = FALSE,
                                     distr_initial_papers = 'poisson',
                                     distr_citations = 'poisson',
                                     dcitations_alpha = dcitationsAlpha,
                                     dcitations_dispersion = dcitationsDispersion,
                                     dcitations_loglog_factor = dcitationsLoglogFactor,
                                     dcitations_speed = dcitationsSpeed,
                                     dpapers_pois_lambda = 2, alpha_share = .33)

  expect_equal(simulationData$scientists$scientist, 1:n)  # test if scientists ids are euqal to their position in data
  expect_equal(length(simulationData), 2)
  expect_equal(names(simulationData), c('papers', 'scientists'))
  expect_equal(nrow(simulationData$scientists), n)
  expect_equal(names(simulationData$scientists), c('scientist', 'h0', 'hAlpha0'))
  expect_equal(length(as.list(simulationData$scientists$h0)), n)
  expect_equal(length(as.list(simulationData$scientists$hAlpha0)), n)
  expect_equal(nrow(simulationData$papers), max(simulationData$papers$paper))

  # other combinations
  simulate_hindex(runs = 2, n = n, periods = 3)
  simdata <- simulate_hindex(runs = 2, n = n, periods = 3,
                  dpapers_pois_lambda = 4, dcitations_speed = 3,
                  dcitations_peak = 2, dcitations_mean = 5)
  simulate_hindex(runs = 2, n = n, periods = 3, distr_initial_papers = 'nbinomial',
                  distr_citations = 'nbinomial', coauthors = 3, strategic_teams = TRUE,
                  diligence_share = .8, diligence_corr = .5, selfcitations = TRUE,
                  update_alpha_authors = TRUE, boost = TRUE, boost_size = .3,
                  alpha_share = .4)
  simulate_hindex(runs = 2, n = n, periods = 3, distr_initial_papers = 'nbinomial',
                  dpapers_nbinom_dispersion = , dpapers_nbinom_mean = ,
                  distr_citations = 'nbinomial', dcitations_speed = 3,
                  dcitations_peak = 3, dcitations_mean = 2, dcitations_dispersion = 1.3)

  plot_hsim(simdata, plot_hindex = TRUE, plot_halpha = TRUE,
            group_boundaries = c(3, 5))
  plot_hsim(simdata, plot_hindex = FALSE, plot_halpha = TRUE,
            group_boundaries = list(c(1, 2), c(3, 4), c(5, Inf)))
  plot_hsim(simdata, plot_halpha = TRUE,
            group_boundaries = c(5), plot_group_diffs = TRUE)
  plot_hsim(simdata, plot_hindex = TRUE,
            group_boundaries = list(c(0, 5), c(4, 9)), plot_group_diffs = TRUE)
  plot_hsim(simdata, plot_halpha = TRUE,
            group_boundaries = 'median', exclude_group_boundaries = TRUE,
            plot_group_diffs = TRUE)

})
