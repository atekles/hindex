context("Simulating h-index")
library(foreach)

test_that("h_sim", {

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
                                     subgroups_distr = 1,
                                     subgroup_advantage = 1,
                                     init_type = 'fixage',
                                     distr_initial_papers = 'poisson',
                                     max_age_scientists = 5,
                                     productivity = exp(2.466973) * (80 / 100) ^ 2.47832,
                                     distr_citations = 'poisson',
                                     dcitations_alpha = dcitationsAlpha,
                                     dcitations_dispersion = dcitationsDispersion,
                                     dcitations_loglog_factor = dcitationsLoglogFactor,
                                     dcitations_speed = dcitationsSpeed,
                                     dpapers_pois_lambda = 2, alpha_share = .33)

  expect_equal(simulationData$scientists$scientist, 1:n)  # test if scientists ids are euqal to their position in data
  expect_equal(length(simulationData), 3)
  expect_equal(names(simulationData), c('papers', 'scientists', 'scientistsAgeInit'))
  expect_equal(nrow(simulationData$scientists), n)
  expect_equal(length(simulationData$scientistsAgeInit), n)
  expect_equal(names(simulationData$scientists), c('scientist', 'subgroup', 'h0', 'toppapers0', 'hAlpha0'))
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
  simulate_hindex(runs = 2, n = n, periods = 3, strategic_teams = TRUE,
                  subgroups_distr = .3)

  plot_hsim(simdata, plot_hindex = TRUE, plot_halpha = TRUE,
            group_boundaries = c(3, 5))
  plot_hsim(simdata, plot_hindex = TRUE, plot_halpha = TRUE,
            plot_mindex = TRUE, group_boundaries = c(3, 5))
  plot_hsim(simdata, plot_hindex = TRUE, plot_halpha = TRUE,
            plot_toppapers = TRUE, group_boundaries = c(3, 5))
  plot_hsim(simdata, plot_hindex = TRUE, plot_halpha = FALSE,
            plot_toppapers = TRUE, group_boundaries = c(3, 5))
  plot_hsim(simdata, plot_hindex = FALSE, plot_halpha = TRUE,
            group_boundaries = list(c(1, 2), c(3, 4), c(5, Inf)))
  plot_hsim(simdata, plot_halpha = TRUE,
            group_boundaries = c(5), plot_group_diffs = TRUE)
  plot_hsim(simdata, plot_toppapers = TRUE,
            group_boundaries = c(5), plot_group_diffs = TRUE)
  plot_hsim(simdata, plot_mindex = TRUE,
            group_boundaries = c(5), plot_group_diffs = TRUE)
  plot_hsim(simdata, plot_hindex = TRUE,
            group_boundaries = list(c(0, 5), c(4, 9)), plot_group_diffs = TRUE)
  plot_hsim(simdata, plot_halpha = TRUE,
            group_boundaries = 'median', exclude_group_boundaries = TRUE,
            plot_group_diffs = TRUE)


  simdata <- simulate_hindex(runs = 2, n = n, periods = 3,
                             subgroups_distr = .5,
                             dpapers_pois_lambda = 4, dcitations_speed = 3,
                             dcitations_peak = 2, dcitations_mean = 5)
  plot_hsim(simdata, plot_halpha = TRUE,
            subgroups = TRUE,
            plot_group_diffs = TRUE)

  skip_on_cran()

  expect_warning(simulate_hindex(runs = 2, n = n, periods = 3, strategic_teams = TRUE,
                                 subgroups_distr = .1))
  expect_warning(simulate_hindex(runs = 2, n = n, periods = 3, strategic_teams = TRUE,
                                 subgroups_distr = .05))
  expect_warning(simulate_hindex(runs = 2, n = n, periods = 3, strategic_teams = TRUE,
                                 subgroups_distr = .9))
  expect_warning(simulate_hindex(runs = 2, n = n, periods = 3, strategic_teams = TRUE,
                                 subgroups_distr = .95))

  init_types <- c('fixage', 'varage')
  ns <- c(20, 30, 100, 100, 100)
  periodss <- c(1:5)
  subgroups_distrs <- c(.1, .2, .5, .9, 1)
  subgroup_advantages <- c(-10, -.3, 0, .1, 100)
  subgroup_exchanges <- c(0, .3, .5, .8, .1)
  coauthorss <- c(2, 3, 5, 10, 100)
  max_age_scientistss <- c(1, 2, 4, 10, 100)
  dpapers_pois_lambdas <- c(1:3, 10.321, 100)
  productivitys <- c(0, 5, 40.3, 80, 100)
  dcitations_speeds <- c(1.1, 3, 2, 10, 100)
  dcitations_peaks <- c(1, 3.1, 5, 10, 100)
  dcitations_means <- c(1, 2.1, 5, 10, 100)
  strategic_teamss <- c(T, F)
  diligence_shares <- c(0, .2, .5, .9, 1)
  diligence_corrs <- c(0, .3, .5, .8, 1)
  selfcitationss <- c(T, F)
  update_alpha_authorss <- c(T, F)
  boosts <- c(T, F)
  boost_sizes <- c(-100, 0, .01, 10, 100)
  alpha_shares <- c(.01, .3, .5, .99, 1)
  dpapers_nbinom_dispersions <- c(.01, 1, 5, 10.123, 100)
  dpapers_nbinom_means <- c(.01, 1, 5, 10.123, 100)
  dcitations_dispersions <- c(1, 2, 5, 10.123, 100)

  message('tests for poisson distrs...')

  for (i in 1:5) {
    message('test run ', i, '...')
    simulate_hindex(runs = 2, n = ns[i], periods = periodss[i],
                    subgroups_distr = subgroups_distrs[i],
                    subgroup_advantage = subgroup_advantages[i],
                    subgroup_exchange = subgroup_exchanges[i],
                    coauthors = coauthorss[i],
                    max_age_scientists = max_age_scientistss[i],
                    init_type = init_types[i %% 2 + 1],
                    distr_initial_papers = 'poisson',
                    dpapers_pois_lambda = dpapers_pois_lambdas[i],
                    productivity = productivitys[i],
                    distr_citations = 'poisson',
                    dcitations_speed = dcitations_speeds[i],
                    dcitations_peak = dcitations_peaks[i],
                    dcitations_mean = dcitations_means[i],
                    strategic_teams = strategic_teamss[i %% 2 + 1],
                    diligence_share = diligence_shares[i],
                    diligence_corr = diligence_corrs[i],
                    selfcitations = selfcitationss[i %% 2 + 1],
                    update_alpha_authors = update_alpha_authorss[i %% 2 + 1],
                    boost = boosts[i %% 2 + 1],
                    boost_size = boost_sizes[i],
                    alpha_share = alpha_shares[i]
                    )
  }

  message('tests for nbinomial distrs...')

  for (i in 1:5) {
    message('test run ', i, '...')
    simulate_hindex(runs = 2, n = ns[i], periods = periodss[i],
                    subgroups_distr = subgroups_distrs[i],
                    subgroup_advantage = subgroup_advantages[i],
                    subgroup_exchange = subgroup_exchanges[i],
                    coauthors = coauthorss[i],
                    max_age_scientists = max_age_scientistss[i],
                    init_type = init_types[i %% 2 + 1],
                    distr_initial_papers = 'nbinomial',
                    dpapers_nbinom_dispersion = dpapers_nbinom_dispersions[i],
                    dpapers_nbinom_mean = dpapers_nbinom_means[i],
                    productivity = productivitys[i],
                    distr_citations = 'nbinomial',
                    dcitations_speed = dcitations_speeds[i],
                    dcitations_peak = dcitations_peaks[i],
                    dcitations_mean = dcitations_means[i],
                    dcitations_dispersion = dcitations_dispersions[i],
                    strategic_teams = strategic_teamss[i %% 2 + 1],
                    diligence_share = diligence_shares[i],
                    diligence_corr = diligence_corrs[i],
                    selfcitations = selfcitationss[i %% 2 + 1],
                    update_alpha_authors = update_alpha_authorss[i %% 2 + 1],
                    boost = boosts[i %% 2 + 1],
                    boost_size = boost_sizes[i],
                    alpha_share = alpha_shares[i]
    )
  }

})
