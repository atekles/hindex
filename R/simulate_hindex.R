#' Simulate h-index and h-alpha values
#'
#' Simulate the effect of publishing, being cited, and (strategic)
#' collaborating on the development of h-index and h-alpha values for a
#' specified set of agents.
#'
#' @param runs Number of times the simulation is repeated.
#' @param n Number of agents acting in each simulation.
#' @param periods Number of periods the agents collaborate across in each period.
#' @param distr_initial_papers Distribution of the papers the scientists have
#' already published at the start of the simulation. Currently, the poisson
#' distribution ("poisson") and the negative binomial distribution ("nbinomial")
#' are supported.
#' @param dpapers_pois_lambda The distribution parameter for a poisson
#' distribution of initial papers.
#' @param dpapers_nbinom_dispersion Dispersion parameter of a negative binomial
#' distribution of initial papers.
#' @param dpapers_nbinom_mean Expected value of a negative binomial
#' distribution of initial papers.
#' @param distr_citations Distribution of citations the papers get. The expected
#' value of this distribution follows a log-logistic function of time.
#' Currently, the poisson distribution ("poisson") and the negative binomial
#' distribution ("nbinomial") are supported.
#' @param dcitations_speed The steepness (shape parameter) of the log-logistic
#' time function of the expected citation values.
#' @param dcitations_peak The period after publishing when the expected value
#' of the citation distribution reaches its maximum.
#' @param dcitations_mean The maximum expected value of the citation
#' distribution (at period dcitations_peak after publishing, the citation
#' distribution has dcitations_mean).
#' @param dcitations_dispersion For a negative binomial citation distribution,
#' dcitations_dispersion is a factor by which the variance exceeds the expected
#' value.
#' @param coauthors Average number of coauthors publishing papers.
#' @param strategic_teams If this parameter is set to TRUE, agents with high
#' h-index avoid co-authorships with agents who have equal or higher h-index
#' values (they strategically select co-authors to improve their h-alpha index).
#' This is implemented by assigning the agents with the highest h-index values
#' to separate teams and randomly assigning the other agents to the teams.
#' Otherwise, the collaborating agents are assigned to co-authorships at random.
#' @param diligence_share The share of agents publishing in each period.
#' @param diligence_corr The correlation between the initial h-index value and
#' the probability to publish in a given period. This parameter only has an
#' effect if diligence_share < 1.
#' @param selfcitations If this parameter is set to TRUE, a paper gets one
#' additional citation if at least one of its authors has a h-index value
#' that exceeds the number of previous citations of the paper by one or two.
#' This reflects agents strategically citing their own papers with citations
#' just below their h-index to accelerate the growth of their h-index.
#' @param update_alpha_authors If this parameter is set to TRUE, the alpha
#' author of newly written papers is determined every period based on the
#' current h-index values of its authors. Without this option, the alpha
#' author is determined when the paper is written and held constant from then on.
#' @param boost If this parameter is set to TRUE, papers of agents with a higher
#' h-index are cited more frequently than papers of agents with lower h-index.
#' For each team, this effect is based on the team's co-author with the
#' highest h-index within this team.
#' @param boost_size Magnitude of the boost effect. For every additional h point
#' of a paper's co-author who has the highest h-index among all of the paper's
#' co-authors, citations of the paper are increased by boost_size, rounded to the next
#' integer.
#' @param alpha_share The share of previously published papers where the
#' corresponding agent is alpha author.
#'
#' @return For each run, the h-index values and the h-alpha values for each
#' period are stored in a list of lists.
#'
#' @export
#' @importFrom foreach "%do%"
#'
#' @examples
#' set.seed(123)
#' simdata <- simulate_hindex(runs = 2, n = 20, periods = 3)
#' plot_hsim(simdata, plot_hindex = TRUE)
simulate_hindex <- function(runs = 1, n = 100, periods = 20,
                            distr_initial_papers = 'poisson',
                            dpapers_pois_lambda = 2,
                            dpapers_nbinom_dispersion = 1.1, dpapers_nbinom_mean = 2,
                            distr_citations = 'poisson', dcitations_speed = 2,
                            dcitations_peak = 3, dcitations_mean = 2,
                            dcitations_dispersion = 1.1,
                            coauthors = 5, strategic_teams = FALSE,
                            diligence_share = 1, diligence_corr = 0,
                            selfcitations = FALSE, update_alpha_authors = FALSE,
                            boost = FALSE, boost_size = .1,
                            alpha_share = .33) {

  # dcitations_speed: in common log log notation: beta

  if (coauthors <= 1) {
    stop('average teamsize has to be greater than 1')
  }

  # merton effect
  # papers get additional citations proportional to the highest h value
  # among its authors (in previous period); currently the merton effect only
  # comes into effect in the first period of the paper

  dcitations_alpha <-  # the alpha in common log log notation
    dcitations_peak /
    (((dcitations_speed - 1) / (dcitations_speed + 1)) ^ (1 / dcitations_speed))
  dcitations_loglog_factor <-
    dcitations_mean / (((dcitations_speed / dcitations_alpha) *
                         ((dcitations_peak / dcitations_alpha) ^
                            (dcitations_speed - 1))) /
                        ((1 + (dcitations_peak / dcitations_alpha) ^
                            dcitations_speed) ^ 2))

  hValuesRuns <- list()
  hAlphaValuesRuns <- list()

  for (currentRun in 1:runs) {

    message(paste('run ', currentRun, '...', sep = ''))

    simulationData <- setup_simulation(n = n, boost = boost,
       boost_size = boost_size, distr_initial_papers = distr_initial_papers,
       distr_citations = distr_citations, dcitations_alpha = dcitations_alpha,
       dcitations_dispersion = dcitations_dispersion,
       dcitations_loglog_factor = dcitations_loglog_factor,
       dcitations_speed = dcitations_speed, dpapers_pois_lambda = dpapers_pois_lambda,
       dpapers_nbinom_dispersion = dpapers_nbinom_dispersion, dpapers_nbinom_mean = dpapers_nbinom_mean,
       alpha_share = alpha_share)

    hAlphaValues <- list(simulationData$scientists$hAlpha0)
    names(hAlphaValues) <- c('period-0')
    hValues <- list(simulationData$scientists$h0)
    names(hValues) <- c('period-0')
    nextPaperId <- nrow(simulationData$papers) + 1

    ## iterate over periods

    for (currentPeriod in 1:periods) {

      simulationData$papers$age <- simulationData$papers$age + 1

      # determine author teams
      if (diligence_share != 1) {

        diligence <- diligence_corr * hValues[['period-0']] +
          sqrt(1 - diligence_corr ^ 2) * stats::rnorm(length(hValues[['period-0']]))
        activeScientists <-
          which(diligence > stats::quantile(diligence,
                                     prob = c(1 - diligence_share)))

        nTeams <- floor(length(activeScientists) / coauthors)

        if (strategic_teams) {

          scientistsHOrder <- order(-hValues[[length(hValues)]])
          authorsTeams <- vector(mode = 'numeric',
                                 length = nrow(simulationData$scientists))
          authorsTeams[scientistsHOrder[activeScientists][1:nTeams]] <- 1:nTeams
          authorsTeams[scientistsHOrder[activeScientists][(nTeams + 1):length(activeScientists)]] <-
            sample(nTeams, length(activeScientists) - nTeams, replace = TRUE)

        } else {

          authorsTeams <- vector(mode = 'integer',
                                 length = nrow(simulationData$scientists))
          authorsTeams[activeScientists] <-
            sample(nTeams, length(activeScientists), replace = TRUE)

        }

        # 0-elements in authorsTeams correspond to authors not active in this period

      } else {

        nTeams <- floor(nrow(simulationData$scientists) / coauthors)

        if (strategic_teams) {

          scientistsHOrder <- order(-hValues[[length(hValues)]])
          authorsTeams <- vector(mode = 'numeric',
                                 length = nrow(simulationData$scientists))
          authorsTeams[scientistsHOrder[1:nTeams]] <- 1:nTeams
          authorsTeams[scientistsHOrder[(nTeams + 1):nrow(simulationData$scientists)]] <-
            sample(nTeams, nrow(simulationData$scientists) - nTeams, replace = TRUE)

        } else {

          authorsTeams <-
            sample(nTeams, nrow(simulationData$scientists), replace = TRUE)

        }

      }

      # add paper for each team, intial age = 1
      newPaperIds <- nextPaperId:(nextPaperId + nTeams)
      nextPaperId <- nextPaperId + nTeams + 1
      currentTeam <- 1
      newPapers <- foreach::foreach(currentTeam = 1:nTeams, .combine = 'rbind') %do% {
        # get indices of scientists in this team
        currentAuthors <- which(authorsTeams == currentTeam)
        if (length(currentAuthors) == 0) {
          return(NULL)
        }
        maxH <- max(hValues[[length(hValues)]][currentAuthors])
        alphaAuthors <-
          hValues[[length(hValues)]][currentAuthors] == maxH
        #     -> more than one alpha author possible?

        if (boost) {
          mertonBonus <- round(boost_size * maxH)
          return(cbind(newPaperIds[currentTeam], currentAuthors, 1, 0, alphaAuthors, mertonBonus))
        } else {
          return(cbind(newPaperIds[currentTeam], currentAuthors, 1, 0, alphaAuthors))
        }
      }

      if (boost) {
        colnames(newPapers) <- c( 'paper', 'scientist','age', 'citations', 'alpha', 'merton')
      } else {
        colnames(newPapers) <- c( 'paper', 'scientist','age', 'citations', 'alpha')
      }
      simulationData$papers <- rbind(simulationData$papers, newPapers)

      # update alpha authors if specified
      cfun <- function(a, b) {return(NULL)}
      if (update_alpha_authors) {
        currentPaper <- 0
        foreach::foreach(currentPaper = unique(simulationData$papers$paper),
                                         .combine = 'cfun') %do% {
          currentScientists <- simulationData$papers$scientist[
            simulationData$papers$paper == currentPaper]
          maxH <- max(hValues[[length(hValues)]][currentScientists])
          alphaAuthors <- hValues[[length(hValues)]][currentScientists] == maxH

          paperIndices <- which(simulationData$papers$paper == currentPaper)
          scientistsMatches <-
            match(currentScientists, simulationData$papers$scientist[paperIndices])
          simulationData$papers$alpha[paperIndices][scientistsMatches] <- alphaAuthors

          return(NULL)
        }
      }

      # get citations for all papers (not just the new ones...), considering their age
      if (distr_citations == 'uniform') {
        stop('uniform citation distributino not supported any more')
        # newPapersCitations <-
          # sample(dcitationsUnifMin:dcitationsUnifMax, length(teams), replace = TRUE)
      } else if (distr_citations == 'poisson') {

        simulationData$papers$citations <- apply(simulationData$papers, MARGIN = 1,
                                          FUN = function(currentPaper) {
          lambda <- dcitations_loglog_factor *
            (((dcitations_speed / dcitations_alpha) * ((currentPaper['age'] / dcitations_alpha) ^
                                                       (dcitations_speed - 1))) /
               ((1 + (currentPaper['age'] / dcitations_alpha) ^ dcitations_speed) ^ 2))
          return(currentPaper['citations'] + stats::rpois(1, lambda))
        })

      } else if (distr_citations == 'nbinomial') {

        simulationData$papers$citations <- apply(simulationData$papers, MARGIN = 1,
                                          FUN = function(currentPaper) {
          currentExp <- dcitations_loglog_factor *
            (((dcitations_speed / dcitations_alpha) * ((currentPaper['age'] / dcitations_alpha) ^
                                                       (dcitations_speed - 1))) /
               ((1 + (currentPaper['age'] / dcitations_alpha) ^ dcitations_speed) ^ 2))
          currentP <- currentExp / (currentExp * dcitations_dispersion)
          currentN <- (currentExp * currentP) / (1 - currentP)
          return(currentPaper['citations'] + stats::rnbinom(1, size = currentN, prob = currentP))
        })

      } else {
        stop('citation distribution not supported')
      }

      if (boost) {
        simulationData$papers$citations <-
          simulationData$papers$citations + simulationData$papers$merton
      }

      # selfcitations
      if (selfcitations) {

        # all scientists with a paper written in current period
        papersActiveScientists <- which(simulationData$papers$scientist %in% newPapers[ , 'scientist'])

        paperIndex <- 0
        selfcitedActivePapers <- foreach::foreach(paperIndex = papersActiveScientists,
                                            .combine = 'c') %do% {
          scientistCurrentH <- hValues[[length(hValues)]][
            simulationData$papers$scientist[paperIndex]]
          if (simulationData$papers$citations[paperIndex] == scientistCurrentH - 1
              || simulationData$papers$citations[paperIndex] == scientistCurrentH - 2) {
            return(TRUE)
          } else {
            return(FALSE)
          }
        }

        # get all rows of the selfcited papers
        # (not just the rows corresponding to the selfcitING authors)
        selfcitedPapersIndices <- simulationData$papers$paper %in%
          simulationData$papers$paper[papersActiveScientists[selfcitedActivePapers]]

        simulationData$papers$citations[selfcitedPapersIndices] <-
          simulationData$papers$citations[selfcitedPapersIndices] + 1

      }

      # compute h, hAlpha
      newHs <- stats::aggregate(simulationData$papers[ , 'citations'],
                  by = list(simulationData$papers[ , 'scientist']),
                  FUN = function(citations) {
                    length(which(citations >= rank(-citations, ties.method = 'first')))
                  }
      )
      # store new h-index values
      periodLabel <- paste('period-', currentPeriod, sep = '')
      hValues[[periodLabel]] <- vector(mode = 'numeric', length = n)
      hValues[[periodLabel]][newHs$Group.1] <- newHs$x
      rm(newHs)

      newHAlphas <- stats::aggregate(1:nrow(simulationData$papers),
                  by = list(simulationData$papers[ , 'scientist']),
                  FUN = function(x) {
                    length(which(simulationData$papers[x, 'citations'] >=
                                   rank(-simulationData$papers[x, 'citations'],
                                        ties.method = 'first') &  # core papers
                                 simulationData$papers[x, 'alpha']))  # alpha papers
                  }
      )
      hAlphaValues[[periodLabel]] <- vector(mode = 'numeric', length = n)
      hAlphaValues[[periodLabel]][newHAlphas$Group.1] <- newHAlphas$x
      rm(newHAlphas)

    }

    hValuesRuns[[paste('run-', currentRun, sep = '')]] <- hValues
    hAlphaValuesRuns[[paste('run-', currentRun, sep = '')]] <- hAlphaValues

  }

  res <- list(hValuesRuns, hAlphaValuesRuns)
  names(res) <- c('h', 'h_alpha')
  return(res)

}

setup_simulation <- function(n, boost, boost_size = 0,
                             distr_initial_papers, distr_citations,
                             dcitations_loglog_factor, dcitations_alpha,
                             dcitations_speed, dcitations_dispersion,
                             dpapers_pois_lambda = NULL,
                             dpapers_nbinom_dispersion = NULL, dpapers_nbinom_mean = NULL,
                             alpha_share) {

  ## initialization

  # create n scientists and their papers

  if (distr_initial_papers == 'uniform') {
    stop('uniform paper distribution not supported any more')
    # noPapers <- sample(dpapersUnifMin:dpapersUnifMax, n, replace = TRUE)
  } else if (distr_initial_papers == 'poisson') {
    noPapers <- stats::rpois(n = n, lambda = dpapers_pois_lambda)
  } else if (distr_initial_papers == 'nbinomial') {
    noPapers <-
      stats::rnbinom(n = n, size = dpapers_nbinom_dispersion, mu = dpapers_nbinom_mean)
  } else {
    stop('paper distribution not supported')
  }

  zeroPaperScientists <- list() # scientists may have zero papers at start;
  # necessary to collect these scientists in order to consider them later on
  initialScientist <- 1

  # calculate distribution parameters
  maxPaperAge <- 5
  currentPaperAge <- 1
  if (distr_citations == 'uniform') {
    stop('uniform citation distribution not supported')
  } else if (distr_citations == 'poisson') {

    lambdas <- vapply(1:maxPaperAge, function(currentPaperAge) {
      dcitations_loglog_factor *
        (((dcitations_speed / dcitations_alpha) * ((currentPaperAge / dcitations_alpha) ^
                                                     (dcitations_speed - 1))) /
           ((1 + (currentPaperAge / dcitations_alpha) ^ dcitations_speed) ^ 2))
    }, double(1))

  } else if (distr_citations == 'nbinomial') {

    nbinomExps <- vapply(1:maxPaperAge, function(currentPaperAge) {
      dcitations_loglog_factor *
        (((dcitations_speed / dcitations_alpha) * ((currentPaperAge / dcitations_alpha) ^
                                                     (dcitations_speed - 1))) /
           ((1 + (currentPaperAge / dcitations_alpha) ^ dcitations_speed) ^ 2))
    }, double(1))
    # TODO review STATA
    nbinomP <- 1 / dcitations_dispersion
    nbinomNs <- (nbinomExps * nbinomP) / (1 - nbinomP)

  } else {
    stop('citation distribution not supported')
  }

  papers <- cbind(rep(1:n, noPapers),
                  sample(maxPaperAge, sum(noPapers), replace = TRUE))
  zeroPaperScientists <- which(noPapers == 0)

  # each element of papersOlderEq: indices of papers that are at least
  # the age of the index in papersOlderEq
  papersOlderEq <- lapply(1:maxPaperAge, function(currentPaperAge) {
    which(papers[, 2] >= currentPaperAge)
  })

  if (distr_citations == 'poisson') {

    # add citations to papers
    # for each possible age: draw as many poisson distributed values as papers
    # with at least this age are in the data
    papers <- cbind(papers, foreach::foreach(currentPaperAge = 1:maxPaperAge,
                                             .combine = `+`) %do% {

      currentAgeCitations <- vector(mode = 'integer', length = nrow(papers))
      currentAgeCitations[papersOlderEq[[currentPaperAge]]] <-
        stats::rpois(length(papersOlderEq[[currentPaperAge]]),
                     lambdas[currentPaperAge])
      return(currentAgeCitations)

    })

  } else if (distr_citations == 'nbinomial') {

    # add citations to papers
    # for each possible age: draw as many nbinomial distributed values as papers
    # with at least this age are in the data
    papers <- cbind(papers, foreach::foreach(currentPaperAge = 1:maxPaperAge,
                                             .combine = `+`) %do% {

      currentAgeCitations <- vector(mode = 'integer', length = nrow(papers))
      currentAgeCitations[papersOlderEq[[currentPaperAge]]] <-
        stats::rnbinom(length(papersOlderEq[[currentPaperAge]]),
                       size = nbinomNs[currentPaperAge],
                       prob = nbinomP)
      return(currentAgeCitations)

    })

  }

  # define alpha papers

  if (boost) {
    papers <- cbind(1:nrow(papers), papers, stats::runif(nrow(papers)) < alpha_share, 0)
    colnames(papers) <- c( 'paper', 'scientist','age', 'citations', 'alpha', 'merton')
  } else {
    papers <- cbind(1:nrow(papers), papers, stats::runif(nrow(papers)) < alpha_share)
    colnames(papers) <- c( 'paper', 'scientist','age', 'citations', 'alpha')
  }

  # calculate h

  scientists <- stats::aggregate(papers[ , 'citations'],
                          by = list(papers[ , 'scientist']),
                          FUN = function(citations) {
                            # rank(-citations, ties.method = 'first')
                            #   --> paper order by desc nr of citations
                            length(which(citations >=
                                           rank(-citations, ties.method = 'first')))
                          }
  )
  colnames(scientists) <- c('scientist', 'h0')

  if (length(zeroPaperScientists) > 0) {
    # add scientists with no paper (wouldn't be considered in previous step)
    newScientists <- cbind(zeroPaperScientists, 0)
    colnames(newScientists) <- c('scientist', 'h0')
    scientists <- rbind(scientists, newScientists)
  }

  # make sure scientists ids correspond to their indices in scientists data:
  scientists <- scientists[order(scientists$scientist), ]

  if (boost) {
    # determine merton bonus for papers
    # all alpha papers: merton bonus based on their h
    # non-alpha papers: randomly define the max h for the team -> a bit higher
    #     than their own h
    papers[ , 'merton'] <- scientists$h0[match(papers[ , 'scientist'], scientists$scientist)]
    papers[ , 'merton'][!papers[ , 'alpha']] <-
      papers[ , 'merton'][!papers[ , 'alpha']] +
      sample(5, length(which(!papers[ , 'alpha'])), replace = TRUE)
    papers[ , 'merton'] <- round(papers[ , 'merton'] * boost_size)
  }

  # calculate h alpha

  hAlphas <- stats::aggregate(1:nrow(papers),
                      by = list(papers[ , 'scientist']),
                      FUN = function(scientistPapers) {
                        length(which(papers[scientistPapers, 'citations'] >=
                                       rank(-papers[scientistPapers, 'citations'],
                                            ties.method = 'first') &  # core papers
                                       papers[scientistPapers, 'alpha'])) # alpha papers
                      })
  colnames(hAlphas) <- c('scientist', 'hAlpha0')
  scientists$hAlpha0 <- 0
  scientists$hAlpha0[hAlphas$scientist] <- hAlphas$hAlpha0
  #   -> possible because scientists id corresponds to their index in scientists data
  rm(hAlphas)

  return(list(papers = data.frame(papers), scientists = data.frame(scientists)))

}

