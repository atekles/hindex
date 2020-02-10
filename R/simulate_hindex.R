#' Simulate h-index and h-alpha values
#'
#' Simulate the effect of publishing, being cited, and (strategic)
#' collaborating on the development of h-index and h-alpha values for a
#' specified set of agents.
#'
#' @param runs Number of times the simulation is repeated.
#' @param n Number of agents acting in each simulation.
#' @param periods Number of periods the agents collaborate across in each period.
#' @param subgroups_distr Share of scientists in the first subgroup among all
#' scientists
#' @param subgroup_advantage Factor by which citations of papers
#' published by agents of subgroup 2 exceed those of papers published by
#' subgroup 1. This option is intended to reflect subdisciplines with
#' different citation levels.
#' @param subgroup_exchange Share of agents publishing (alone or in collaboration)
#' with the other subgroup in each period. For example, when specifying
#' subgroup_exchange = .1, 10\% of each subgroup join the other subgroup each period.
#' @param init_type Type of the initial setup. May be 'fixage' or 'varage'.
#' For init_type = 'fixage', all initial papers have the same age (specified
#' by max_age_scientists). For init_type = 'varage', papers get a random age
#' which is less than or equal to max_age_scientists.
#' @param distr_initial_papers Distribution of the papers the scientists have
#' already published at the start of the simulation. Currently, the poisson
#' distribution ("poisson") and the negative binomial distribution ("nbinomial")
#' are supported.
#' @param max_age_scientists Maximum age of scientists at the start of the
#' simulation. For init_type = varage, a random age less than or equal to
#' max_age_scientists is assigned to the initial papers. For init_type = fixage,
#' all papers are max_age_scientists old.
#' @param dpapers_pois_lambda The distribution parameter for a poisson
#' distribution of initial papers.
#' @param dpapers_nbinom_dispersion Dispersion parameter of a negative binomial
#' distribution of initial papers.
#' @param dpapers_nbinom_mean Expected value of a negative binomial
#' distribution of initial papers.
#' @param productivity The share of papers published by the 20\% most
#' productive agents in percentage. This parameter is only used for init_type = 'varage'.
#' For init_type = 'fixage', diligence_share and diligence_corr can be used to
#' control the productivity of scientists.
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
#' @param diligence_share The share of agents publishing in each period. Only
#' used for init_type = 'fixage'.
#' @param diligence_corr The correlation between the initial h-index value and
#' the probability to publish in a given period. This parameter only has an
#' effect if diligence_share < 1. Only used for init_type = 'fixage'.
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
                            subgroups_distr = 1, subgroup_advantage = 1,
                            subgroup_exchange = 0, init_type = 'fixage',
                            distr_initial_papers = 'poisson',
                            max_age_scientists = 5,
                            dpapers_pois_lambda = 2,
                            dpapers_nbinom_dispersion = 1.1, dpapers_nbinom_mean = 2,
                            productivity = 80,
                            distr_citations = 'poisson', dcitations_speed = 2,
                            dcitations_peak = 3, dcitations_mean = 2,
                            dcitations_dispersion = 1.1,
                            coauthors = 5, strategic_teams = FALSE,
                            diligence_share = 1, diligence_corr = 0,
                            selfcitations = FALSE, update_alpha_authors = FALSE,
                            boost = FALSE, boost_size = .1,
                            alpha_share = .33) {

  # check parameters
  if (runs <= 0 || runs %% 1 != 0) {
    stop('runs must be an integer > 0')
  }

  if (n <= 1 || n %% 1 != 0) {
    stop('n must be an integer > 1')
  }

  if (periods <= 0 || periods %% 1 != 0) {
    stop('periods must be an integer > 0')
  }

  if (subgroups_distr <= 0 || subgroups_distr > 1) {
    stop('subgroups_distr must be > 0 and <= 1')
  }

  if (subgroup_exchange < 0 || subgroup_exchange > 1) {
    stop('subgroup_exchange must be between 0 and 1')
  }

  if (coauthors <= 1 || coauthors %% 1 != 0) {
    stop('average teamsize has to be an integer greater than 1')
  }

  if (init_type == 'fixage') {
    if (!missing(productivity)) {
      warning('productivity is not used for init_type fixage;
              specify init_type varage in order to use productivity')
    }
  } else if (init_type == 'varage') {
    if (!missing(diligence_share) || !missing(diligence_corr)) {
      warning('diligence_share and diligence_corr are not used for init_type varage;
              specify init_type fixage in order to use diligence_share or diligence_corr')
    }
  } else {
    stop('init_type not supported')
  }

  if (distr_initial_papers != 'poisson' && distr_initial_papers != 'nbinomial') {
    stop('distr_initial_papers must be either \'poisson\' or \'binomial\'')
  }

  if (max_age_scientists < 1 || max_age_scientists %% 1 != 0) {
    stop('max_age_scientists must be an integer >= 1')
  }

  if (dpapers_pois_lambda <= 0) {
    stop('dpapers_pois_lambda must be > 0')
  }

  if (dpapers_nbinom_dispersion <= 0) {
    stop('dpapers_nbinom_dispersion must be > 0')
  }

  if (dpapers_nbinom_mean <= 0) {
    stop('dpapers_nbinom_mean must be > 0')
  }

  if (productivity < 0 || productivity > 100) {
    stop('productivity must be between 0 an 100')
  }

  if (distr_citations != 'poisson' && distr_citations != 'nbinomial') {
    stop('distr_citations must be either \'poisson\' or \'binomial\'')
  }

  if (dcitations_speed <= 1) {
    stop('dcitations_speed must be > 1')
  }

  if (dcitations_peak <= 0) {
    stop('dcitations_peak must be > 0')
  }

  if (dcitations_mean <= 0) {
    stop('dcitations_mean must be > 0')
  }

  if (dcitations_dispersion < 1) {
    stop('dcitations_dispersion must be >= 1')
  }

  if (diligence_share < 0 || diligence_share > 1) {
    stop('diligence_share must be between 0 and 1')
  }

  if (diligence_corr < 0 || diligence_corr > 1) {
    stop('diligence_corr must be between 0 and 1')
  }

  if (alpha_share < 0 || alpha_share > 1) {
    stop('alpha_share must be between 0 and 1')
  }

  # some constant distribution/productivity parameters
  # dcitations_speed: in common log log notation: beta
  dcitations_alpha <-  # the alpha in common log log notation
    dcitations_peak /
    (((dcitations_speed - 1) / (dcitations_speed + 1)) ^ (1 / dcitations_speed))
  dcitations_loglog_factor <-
    dcitations_mean / (((dcitations_speed / dcitations_alpha) *
                         ((dcitations_peak / dcitations_alpha) ^
                            (dcitations_speed - 1))) /
                        ((1 + (dcitations_peak / dcitations_alpha) ^
                            dcitations_speed) ^ 2))
  productivity_param <- exp(2.466973) * (productivity / 100) ^ 2.47832

  hValuesRuns <- list()
  hAlphaValuesRuns <- list()
  toppaperValuesRuns <- list()
  mindexValuesRuns <- list()

  for (currentRun in 1:runs) {

    message(paste('run ', currentRun, '...', sep = ''))

    simulationData <- setup_simulation(n = n, boost = boost,
       boost_size = boost_size, subgroups_distr = subgroups_distr,
       subgroup_advantage = subgroup_advantage, init_type = init_type,
       distr_initial_papers = distr_initial_papers,
       max_age_scientists = max_age_scientists, productivity_param = productivity_param,
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
    toppaperValues <- list(simulationData$scientists$toppapers0)
    names(toppaperValues) <- c('period-0')
    mindexValues <- list(simulationData$scientists$h0[simulationData$scientistsAgeInit > 0] /
                           simulationData$scientistsAgeInit[simulationData$scientistsAgeInit > 0])
    names(mindexValues) <- c('period-0')
    nextPaperId <- nrow(simulationData$papers) + 1

    initial_papers <- simulationData$papers$paper

    ## iterate over periods

    for (currentPeriod in 1:periods) {

      simulationData$papers$age <- simulationData$papers$age + 1

      # update alpha authors if specified
      # this is done before creating new papers, because the alpha authors
      # are determined for new papers anyway --> avoid doing this twice
      # for the new papers
      cfun <- function(a, b) {return(NULL)}
      if (update_alpha_authors) {
        currentPaper <- 0
        # do not update alpha and merton for the initial papers because
        # the max h-index for these author teams are (partly) assigned randomly
        # (assuming that there are co-authors not included in the simulation as agents)
        period_papers <- setdiff(unique(simulationData$papers$paper), initial_papers)
        foreach::foreach(currentPaper = period_papers,
                         .combine = 'cfun') %do% {

                           paperIndices <- which(simulationData$papers$paper == currentPaper)

                           currentScientists <- simulationData$papers$scientist[paperIndices]
                           maxH <- max(hValues[[length(hValues)]][currentScientists])
                           alphaAuthors <- hValues[[length(hValues)]][currentScientists] == maxH

                           # scientistsMatches <-
                           #   match(currentScientists, simulationData$papers$scientist[paperIndices])
                           simulationData$papers$alpha[paperIndices] <- alphaAuthors
                           if (boost) {
                             simulationData$papers$merton[paperIndices] <-
                               round(boost_size * maxH)
                           }

                           return(NULL)

                         }
      }

      # determine author teams

      if (init_type == 'fixage') {

        if (diligence_share != 1) {
          diligence <- diligence_corr * hValues[['period-0']] +
            sqrt(1 - diligence_corr ^ 2) * stats::rnorm(length(hValues[['period-0']]))
          activeScientists <-
            which(diligence > stats::quantile(diligence,
                                              prob = (1 - diligence_share)))
        } else {
          activeScientists <- 1:nrow(simulationData$scientists)
        }

      } else if (init_type == 'varage') {

        activeScientists <-
          which(simulationData$scientists$productivity >
                  stats::runif(nrow(simulationData$scientists)))

      }

      nTeams <- length(activeScientists) / coauthors
      nTeamsGroup1 <- round(nTeams * subgroups_distr)
      nTeamsGroup2 <- round(nTeams * (1-subgroups_distr))

      if (nTeamsGroup1 == 0) {
        warning('no teams in subgroup one (before subgroup exchange); increase subgroups_distr or n')
      }
      if (subgroups_distr < 1 && nTeamsGroup2 == 0) {
        warning('no teams in subgroup two (before subgroup exchange); decrease subgroups_distr or increase n')
      }

      # indices of active scientists in subgroup 1 in simulationData$scientists
      activeScientistsGroup1 <-
        activeScientists[
          which(simulationData$scientists$subgroup[activeScientists] == 1)]
      activeScientistsGroup2 <-
        activeScientists[
          which(simulationData$scientists$subgroup[activeScientists] == 2)]

      # exchange between groups
      if (subgroup_exchange > 0) {
        if (length(activeScientistsGroup1) == 0) {
          # no scientists can exchange from group 1 to group 2
          fromOneToTwo <- vector(mode = 'logical', length = 0)
        } else {
          fromOneToTwo <- stats::runif(length(activeScientistsGroup1)) < subgroup_exchange
        }
        if (subgroups_distr < 1 && nTeamsGroup2 == 0) {
          # no scientists can exchange from group 2 to group 1
          fromOneToTwo <- vector(mode = 'logical', length = 0)
        } else {
          fromTwoToOne <- stats::runif(length(activeScientistsGroup2)) < subgroup_exchange
        }
        fromOneToTwoIds <- activeScientistsGroup1[fromOneToTwo]
        fromTwoToOneIds <- activeScientistsGroup2[fromTwoToOne]
        activeScientistsGroup1 <-
          c(activeScientistsGroup1[!fromOneToTwo], fromTwoToOneIds)
        activeScientistsGroup2 <-
          c(activeScientistsGroup2[!fromTwoToOne], fromOneToTwoIds)
      }

      # due to exchange, less scientists may be active in a group than teams are specified
      # limit the number of teams to the number of scientists at max
      nTeamsGroup1 <- min(nTeamsGroup1, length(activeScientistsGroup1))
      nTeamsGroup2 <- min(nTeamsGroup2, length(activeScientistsGroup2))
      nTeams <- nTeamsGroup1 + nTeamsGroup2   # correct for rounding errors

      # if nTeams == 0: no new papers to create
      if (nTeams > 0) {

        # for each team one entry (will be filled with author ids)
        teamsAuthors <- list()

        if (strategic_teams) {

          if (nTeamsGroup1 > 0) {

            scientistsHOrder <- order(-hValues[[length(hValues)]][activeScientistsGroup1])

            # fill entries in teamsAuthors for first subgroup with best authors in group 1
            foreach::foreach(currentTeam = 1:nTeamsGroup1) %do% {
              teamsAuthors[[currentTeam]] <- activeScientistsGroup1[scientistsHOrder[currentTeam]]
              return(NULL)
            }

            # for each author in group 1, assign team id randomly
            authorsTeams1 <- sample(nTeamsGroup1, length(activeScientistsGroup1) - nTeamsGroup1, replace = TRUE)
            # fill teams based on these random team ids
            foreach::foreach(currentTeam = 1:nTeamsGroup1) %do% {
              teamsAuthors[[currentTeam]] <- c(teamsAuthors[[currentTeam]],
                                               activeScientistsGroup1[
                                                 scientistsHOrder[(nTeamsGroup1 + 1):length(activeScientistsGroup1)][
                                                   authorsTeams1 == currentTeam]]
              )
              return(NULL)
            }

          }

          # same for scientists in subgroup 2; if no subgroup2,
          # the following lines don't change anything
          if (nTeamsGroup2 > 0) {
            scientistsHOrder <- order(-hValues[[length(hValues)]][activeScientistsGroup2])
            foreach::foreach(currentTeam = 1:nTeamsGroup2) %do% {
              teamsAuthors[[currentTeam + nTeamsGroup1]] <- activeScientistsGroup2[scientistsHOrder[currentTeam]]
              return(NULL)
            }
            authorsTeams2 <- sample(nTeamsGroup2, length(activeScientistsGroup2) - nTeamsGroup2, replace = TRUE)
            foreach::foreach(currentTeam = 1:nTeamsGroup2) %do% {
              teamsAuthors[[currentTeam + nTeamsGroup1]] <- c(teamsAuthors[[currentTeam + nTeamsGroup1]],
                                                              activeScientistsGroup2[
                                                                scientistsHOrder[(nTeamsGroup2 + 1):length(activeScientistsGroup2)][
                                                                  authorsTeams2 == currentTeam]]
              )
              return(NULL)
            }
          }

        } else {

          if (nTeamsGroup1 > 0) {
            authorsTeams1 <- sample(nTeamsGroup1, length(activeScientistsGroup1), replace = TRUE)
            foreach::foreach(currentTeam = 1:nTeamsGroup1) %do% {
              teamsAuthors[[currentTeam]] <- activeScientistsGroup1[authorsTeams1 == currentTeam]
              return(NULL)
            }
          }

          if (nTeamsGroup2 > 0) {
            authorsTeams2 <- sample(nTeamsGroup2, length(activeScientistsGroup2), replace = TRUE)
            foreach::foreach(currentTeam = 1:nTeamsGroup2) %do% {
              teamsAuthors[[currentTeam + nTeamsGroup1]] <- activeScientistsGroup2[authorsTeams2 == currentTeam]
              return(NULL)
            }
          }

        }

        # TODO ver
        # if nTeams1 > 0: all activeScientists1 in unlist(teamsAuthors)
        if (nTeamsGroup1 > 0 && !all(activeScientistsGroup1 %in% unlist(teamsAuthors))) {
          stop('not all active scientists assigned to teams in subgroup 1')
        }
        # if nTeams2 > 0: all activeScientists2 in unlist(teamsAuthors)
        if (nTeamsGroup2 > 0 && !all(activeScientistsGroup2 %in% unlist(teamsAuthors))) {
          stop('not all active scientists assigned to teams in subgroup 2')
        }
        # length(teamsAuthors) == nTeams == nTeams1 + nTeams2
        if (length(teamsAuthors) != nTeams) {
          # if a team has no authors assigned (can happen because team memberships are assigned randomly with replacement)
          stop('no teams does not equal the no of expected teams')
        }
        # no author is assigned to both subgroups
        if (nTeamsGroup1 > 0 && nTeamsGroup2 > 0 && length(intersect(unlist(teamsAuthors[1:nTeamsGroup1]), unlist(teamsAuthors[(nTeamsGroup1 + 1):(nTeams)]))) != 0) {
          stop('some authors are assigned to multiple subgroups')
        }

        # add paper for each team, intial age = 1
        newPaperIds <- nextPaperId:(nextPaperId + nTeams - 1)
        nextPaperId <- nextPaperId + nTeams

        # TODO ver
        if (nTeams > 0 && length(newPaperIds) != nTeams) {
          stop('length(newPaperIds) != length(newPapers')
        }

        # create new papers
        currentTeam <- 1
        newPapers <- foreach::foreach(currentTeam = 1:nTeams, .combine = 'rbind') %do% {
          # get indices of scientists in this team
          currentAuthors <- teamsAuthors[[currentTeam]]
          if (length(currentAuthors) == 0) {
            return(NULL)
          }
          maxH <- max(hValues[[length(hValues)]][currentAuthors])
          alphaAuthors <-
            hValues[[length(hValues)]][currentAuthors] == maxH
          #     -> more than one alpha author possible?

          if (currentTeam <= nTeamsGroup1) {
            currentSubgroup <- 1
          } else {
            currentSubgroup <- 2
          }

          if (boost) {
            mertonBonus <- round(boost_size * maxH)
            return(cbind(newPaperIds[currentTeam], currentAuthors, 1, 0, alphaAuthors, mertonBonus, currentSubgroup))
          } else {
            return(cbind(newPaperIds[currentTeam], currentAuthors, 1, 0, alphaAuthors, currentSubgroup))
          }
        }

        # TODO ver each paper exactly assigned to one subgroup
        if (nrow(unique(newPapers[ , c(1, ncol(newPapers))])) != length(unique(newPapers[ , 1]))) {
          stop('papers assigned to more than one subgroup')
        }

        if (boost) {
          colnames(newPapers) <- c( 'paper', 'scientist','age', 'citations', 'alpha', 'merton', 'subgroup')
        } else {
          colnames(newPapers) <- c( 'paper', 'scientist','age', 'citations', 'alpha', 'subgroup')
        }
        simulationData$papers <- rbind(simulationData$papers, newPapers)

      }

      # get citations for all papers (not just the new ones...), considering their age

      papers_age_subgroup <- unique(simulationData$papers[ , c('paper', 'age', 'subgroup')])
      if (distr_citations == 'uniform') {
        stop('uniform citation distributino not supported any more')
        # newPapersCitations <-
          # sample(dcitationsUnifMin:dcitationsUnifMax, length(teams), replace = TRUE)
      } else if (distr_citations == 'poisson') {

        currentLambdas <- dcitations_loglog_factor *
          (((dcitations_speed / dcitations_alpha) * ((papers_age_subgroup$age / dcitations_alpha) ^
                                                       (dcitations_speed - 1))) /
             ((1 + (papers_age_subgroup$age / dcitations_alpha) ^ dcitations_speed) ^ 2))
        newCitations <- stats::rpois(length(currentLambdas), currentLambdas)

      } else if (distr_citations == 'nbinomial') {

        currentExps <- dcitations_loglog_factor *
          (((dcitations_speed / dcitations_alpha) * ((papers_age_subgroup$age / dcitations_alpha) ^
                                                       (dcitations_speed - 1))) /
             ((1 + (papers_age_subgroup$age / dcitations_alpha) ^ dcitations_speed) ^ 2))
        currentPs <- currentExps / (currentExps * dcitations_dispersion)
        currentNs <- (currentExps * currentPs) / (1 - currentPs)
        newCitations <- stats::rnbinom(length(currentExps), size = currentNs, prob = currentPs)

      } else {
        stop('citation distribution not supported')
      }

      # add subgroup advantage
      group2Papers <- papers_age_subgroup$subgroup == 2
      newCitations[group2Papers] <- newCitations[group2Papers] * subgroup_advantage

      # add new citations to the papers
      unique_paper_matches <- match(simulationData$papers$paper, papers_age_subgroup$paper)
      simulationData$papers$citations <-
        simulationData$papers$citations + newCitations[unique_paper_matches]

      # add merton bonus if specified
      if (boost) {
        simulationData$papers$citations <-
          simulationData$papers$citations + simulationData$papers$merton
      }

      # selfcitations
      if (selfcitations) {

        # all authorships (rows in simulationData$papers) of active scientists
        papersActiveScientists <-
          which(simulationData$papers$scientist %in% newPapers[ , 'scientist'])

        # for each authorship of active scientists, get h index
        papersActiveScientistsH <-
          hValues[[length(hValues)]][
            # scientist id for each authorship of active scientists
            simulationData$papers$scientist[papersActiveScientists]]
        # rows in simulationData$papers[papersActiveScientists, ] for authorships causing a selfcitation
        selfcitedActivePapers <-
          simulationData$papers$citations[papersActiveScientists] ==
          papersActiveScientistsH - 1 |
          simulationData$papers$citations[papersActiveScientists] ==
          papersActiveScientistsH - 2

        # get all rows of the selfcited papers (rows in simulationData$papers)
        # (not just the rows corresponding to the selfciting authors)
        selfcitedPapersIndices <- simulationData$papers$paper %in%
          simulationData$papers$paper[papersActiveScientists[selfcitedActivePapers]]

        simulationData$papers$citations[selfcitedPapersIndices] <-
          simulationData$papers$citations[selfcitedPapersIndices] + 1

      }

      # determine toppapers

      if (subgroups_distr == 1) {
        papersUnique <- unique(simulationData$papers[ , c('paper', 'citations')])
        #   -> this step is necessary, because simulationData$papers$citations > p90
        #       would be on the level of authorships
        # TODO ver
        if (nrow(papersUnique) != length(unique(simulationData$papers[ , c('paper')]))) {
          stop('nrow(papersUnique) != no of distinct papers 3')
        }
        p90 <- stats::quantile(papersUnique[ , 'citations'], prob = .9, type = 1)
        top10Papers <- simulationData$papers$paper[simulationData$papers$citations > p90]
        toppapersIndicesOld <- simulationData$papers$paper %in% top10Papers
        toppapersIndices <- simulationData$papers$citations > p90
        # TODO ver
        if (!all(toppapersIndices == toppapersIndicesOld)) {
          stop('toppapersIndices != toppapersIndicesOld')
        }
      } else {

        papersUniqueGroup1 <- unique(simulationData$papers[simulationData$papers$subgroup == 1,
                                                           c('paper', 'citations')])
        papersUniqueGroup2 <- unique(simulationData$papers[simulationData$papers$subgroup == 2,
                                                           c('paper', 'citations')])
        # TODO ver
        if ((nrow(papersUniqueGroup1) + nrow(papersUniqueGroup2)) !=
            length(unique(simulationData$papers[ , c('paper')]))) {
          # each paper must be assigned to exactly one subgroup
          stop('nrow(papersUniqueGroup1) + nrow(papersUniqueGroup2) != no of distinct papers')
        }
        p90Group1 <- stats::quantile(
          papersUniqueGroup1[ , 'citations'], prob = .9, type = 1
        )
        p90Group2 <- stats::quantile(
          papersUniqueGroup2[ , 'citations'], prob = .9, type = 1
        )
        top10PapersGroup1 <-
          intersect(
            simulationData$papers$paper[simulationData$papers$citations > p90Group1],
            papersUniqueGroup1[ , 'paper']
          )
        top10PapersGroup2 <-
          intersect(
            simulationData$papers$paper[simulationData$papers$citations > p90Group2],
            papersUniqueGroup2[ , 'paper']
          )
        top10Papers <- union(top10PapersGroup1, top10PapersGroup2)
        #   -> all papers which are a top paper in at least one of the two groups
        #       AND have at least one author in this group
        toppapersIndicesOld <- simulationData$papers$paper %in% top10Papers
        toppapersIndices <-
          (simulationData$papers$citations > p90Group1 & simulationData$papers$subgroup == 1) |
          (simulationData$papers$citations > p90Group2 & simulationData$papers$subgroup == 2)
        # TODO ver
        if (!all(toppapersIndices == toppapersIndicesOld)) {
          stop('toppapersIndices != toppapersIndicesOld')
        }
      }
      # count top papers
      newToppapers <- stats::aggregate(toppapersIndices,
                                        by = list(simulationData$papers$scientist),
                                        FUN = function(top10s) {
                                          length(which(top10s))
                                        }
      )
      # store new toppaper values
      periodLabel <- paste('period-', currentPeriod, sep = '')
      toppaperValues[[periodLabel]] <- vector(mode = 'numeric', length = n)
      toppaperValues[[periodLabel]][newToppapers$Group.1] <- newToppapers$x

      # compute h, hAlpha
      newHs <- stats::aggregate(simulationData$papers[ , 'citations'],
                  by = list(simulationData$papers[ , 'scientist']),
                  FUN = function(citations) {
                    length(which(citations >= rank(-citations, ties.method = 'first')))
                  }
      )
      # store new h-index values
      hValues[[periodLabel]] <- vector(mode = 'numeric', length = n)
      hValues[[periodLabel]][newHs$Group.1] <- newHs$x

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

      # calculate mindex values
      mindexValues[[periodLabel]] <- vector(mode = 'numeric', length = n)
      mindexValues[[periodLabel]][newHs$Group.1] <-
        newHs$x / (simulationData$scientistsAgeInit[newHs$Group.1] + currentPeriod)

    }

    runLabel <- paste('run-', currentRun, sep = '')
    hValuesRuns[[runLabel]] <- hValues
    hAlphaValuesRuns[[runLabel]] <- hAlphaValues
    toppaperValuesRuns[[runLabel]] <- toppaperValues
    mindexValuesRuns[[runLabel]] <- mindexValues

  }

  res <- list(hValuesRuns, hAlphaValuesRuns, toppaperValuesRuns, mindexValuesRuns, simulationData$scientists$subgroup)
  names(res) <- c('h', 'h_alpha', 'top10_papers', 'mindex', 'subgroup')
  return(res)

}

setup_simulation <- function(n, boost, boost_size = 0,
                             subgroups_distr, subgroup_advantage,
                             init_type, distr_initial_papers,
                             max_age_scientists,
                             productivity_param, distr_citations,
                             dcitations_loglog_factor, dcitations_alpha,
                             dcitations_speed, dcitations_dispersion,
                             dpapers_pois_lambda = NULL,
                             dpapers_nbinom_dispersion = NULL, dpapers_nbinom_mean = NULL,
                             alpha_share) {

  if (subgroups_distr < 0 || subgroups_distr > 1) {
    stop('invalid argument for subgroups_distr; must be >0 and <=1')
  }

  ## initialization

  # create n scientists and their papers

  if (init_type == 'fixage') {

    # in stata: init_type 1

    if (distr_initial_papers == 'uniform') {
      stop('uniform paper distribution not supported any more')
      # noPapers <- sample(dpapersUnifMin:dpapersUnifMax, n, replace = TRUE)
    } else if (distr_initial_papers == 'poisson') {
      noPapers <- stats::rpois(n = n, lambda = dpapers_pois_lambda)
    } else if (distr_initial_papers == 'nbinomial') {
      # in R:
      #   -> size is the number of successful trials (dispersion parameter)
      #   -> mu is the mean
      #   -> no need to calculate the probability for single trials or number
      #       of failures (as in Stata)
      # in Stata: rnbinom(n, p)
      #   -> if n is an integer, this is the number of failures before the nth success
      #   -> if n is not an integer?
      #   -> p is the probability of success for a single trial
      noPapers <-
        stats::rnbinom(n = n, size = dpapers_nbinom_dispersion, mu = dpapers_nbinom_mean)
    } else {
      stop('paper distribution not supported')
    }

    papers_age <- sample(max_age_scientists, sum(noPapers), replace = TRUE)

  } else if (init_type == 'varage') {

    # create productivity for each scientist
    scientists_prod <- stats::runif(n) ^ productivity_param
    # assign random age to each scientist
    scientists_age <- sample(max_age_scientists, n, replace = TRUE)
    # decide in which periods the scientists have published --> paper + age
    papers_age_list <- purrr::map2(scientists_age, scientists_prod, function(current_age, current_prod) {
      which(stats::runif(current_age) <= current_prod)
    })
    noPapers <- vapply(papers_age_list, length, integer(1))
    current_scientist_papers <- NULL
    papers_age <- foreach::foreach(current_scientist_papers = papers_age_list,
                                   .combine = c) %do% {
      return(current_scientist_papers)
    }
    rm(current_scientist_papers)

  } else {
    stop('init_type not supported')
  }

  zeroPaperScientists <- list() # scientists may have zero papers at start;
  # necessary to collect these scientists in order to consider them later on
  initialScientist <- 1

  # calculate distribution parameters
  maxPaperAge <- max_age_scientists
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
    nbinomP <- 1 / dcitations_dispersion
    nbinomNs <- (nbinomExps * nbinomP) / (1 - nbinomP)

  } else {
    stop('citation distribution not supported')
  }

  # determine subgroups
  subgroups <- vector(mode = 'integer', length = n)
  subgroups[] <- 1
  if (subgroups_distr < 1) {
    subgroup_break <- round(subgroups_distr * n)
    if (subgroup_break >= 1 && subgroup_break < n) {
      subgroups[(subgroup_break + 1):n] <- 2
    }
  }

  papersSubgroup <- rep(subgroups, noPapers)
  papers <- cbind(rep(1:n, noPapers), papers_age)
  zeroPapers <- which(noPapers == 0)
  zeroPaperScientists <- cbind(zeroPapers, subgroups[zeroPapers])

  # each element of papersOlderEq: indices of papers that are at least
  # the age of the index in papersOlderEq
  papersOlderEq <- lapply(1:maxPaperAge, function(currentPaperAge) {
    which(papers[, 2] >= currentPaperAge)
  })

  # assign citations to papers

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
      group2Papers <- papersSubgroup == 2
      currentAgeCitations[group2Papers] <-
        currentAgeCitations[group2Papers] * subgroup_advantage
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
      group2Papers <- papersSubgroup == 2
      currentAgeCitations[group2Papers] <-
        currentAgeCitations[group2Papers] * subgroup_advantage
      return(currentAgeCitations)

    })

  }

  # combine, name columns and define alpha papers
  if (boost) {
    papers <- cbind(1:nrow(papers), papers, stats::runif(nrow(papers)) < alpha_share, 0, papersSubgroup)
    colnames(papers) <- c('paper', 'scientist', 'age', 'citations', 'alpha', 'merton', 'subgroup')
  } else {
    papers <- cbind(1:nrow(papers), papers, stats::runif(nrow(papers)) < alpha_share, papersSubgroup)
    colnames(papers) <- c('paper', 'scientist', 'age', 'citations', 'alpha', 'subgroup')
  }

  # determine age of scientists
  scientistsAgeInit <- stats::aggregate(papers[ , 'age'],
                              by = list(papers[ , 'scientist']),
                              FUN = function(ages) {
                                return(max(ages))
                              })
  names(scientistsAgeInit) <- c('scientist', 'age')

  # determine top 10% papers

  if (subgroups_distr == 1) {
    # unique stuff: not absolutely necessary here, but consistent with
    # top paper detection in subsequent simulation runs
    papersUnique <- unique(papers[ , c('paper', 'citations')])
    # TODO ver
    if (nrow(papersUnique) != length(unique(papers[ , c('paper')]))) {
      stop('length(papersUnique) != no of distinct papers 1')
    }
    p90 <- stats::quantile(papersUnique[ , 'citations'], prob = .9, type = 1)
    top10Papers <- papers[ , 'citations'] > p90
  } else {
    papersUnique <- unique(papers[ , c('paper', 'subgroup', 'citations')])
    # TODO ver
    if (nrow(papersUnique) != length(unique(papers[ , c('paper')]))) {
      # each paper must be assigned to exactly one subgroup
      stop('length(papersUnique) != no of distinct papers 2')
    }
    p90Group1 <- stats::quantile(
      papersUnique[papersUnique[ , 'subgroup'] == 1, 'citations'], prob = .9, type = 1
    )
    p90Group2 <- stats::quantile(
      papersUnique[papersUnique[ , 'subgroup'] == 2, 'citations'], prob = .9, type = 1
    )
    top10PapersGroup1 <-
      papers[ , 'citations'] > p90Group1 & papers[ , 'subgroup'] == 1
    top10PapersGroup2 <-
      papers[ , 'citations'] > p90Group2 & papers[ , 'subgroup'] == 2
    top10Papers <- top10PapersGroup1 | top10PapersGroup2
  }

  # calculate h

  # assumption: all of a scientist's papers published in one subgroup
  # this is true because paper subgroups are determined only based on
  # scientist subgroup for initial setup
  # TODO ver
  if (nrow(unique(papers[ , c('scientist', 'subgroup')])) != length(unique(papers[ , 'scientist']))) {
    stop('scientists with papers in multiple subgroups in initial setup')
  }
  scientists <- stats::aggregate(papers[ , 'citations'],
                          by = list(papers[ , 'scientist'], papers[ , 'subgroup']),
                          FUN = function(citations) {
                            # rank(-citations, ties.method = 'first')
                            #   --> paper order by desc nr of citations
                            length(which(citations >=
                                           rank(-citations, ties.method = 'first')))
                          }
  )
  colnames(scientists) <- c('scientist', 'subgroup', 'h0')

  # count top papers
  initialTop10s <- stats::aggregate(top10Papers,
                   by = list(papers[ , 'scientist']),
                   FUN = function(top10s) {
                     length(which(top10s == 1))
                   }
  )
  colnames(initialTop10s) <- c('scientist', 'topPapers')
  top10Matches <- match(scientists[ , 'scientist'], initialTop10s[ , 'scientist'])
  #   -> for each row in scientists: matching row in initialTop10s
  scientists <- cbind(scientists, initialTop10s[top10Matches, 'topPapers'])
  colnames(scientists) <- c('scientist', 'subgroup', 'h0', 'toppapers0')

  if (length(zeroPaperScientists) > 0) {
    # add scientists with no paper (wouldn't be considered in previous step)
    newScientists <- cbind(zeroPaperScientists, 0, 0)
    colnames(newScientists) <- c('scientist', 'subgroup', 'h0', 'toppapers0')
    scientists <- rbind(scientists, newScientists)
    newAges <- cbind(zeroPaperScientists[ , 1], 0)
    colnames(newAges) <- c('scientist', 'age')
    scientistsAgeInit <- rbind(scientistsAgeInit, newAges)
  }

  # make sure scientists' ids correspond to their indices in scientists data:
  scientists <- scientists[order(scientists$scientist), ]
  # age of scientists ordered by scientist id
  # TODO ver
  if (!all(scientists$scientist == scientistsAgeInit$scientist[order(scientistsAgeInit$scientist)])) {
    stop('scientist ids in scientists and scientistsAgeInit do not match')
  }
  scientistsAgeInit <- scientistsAgeInit$age[order(scientistsAgeInit$scientist)]
  # TODO ver
  if (!all(scientists$scientist == 1:n)) {
    stop('not exactly one row in scientists for each agent')
  }

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

  # add productivity of scientists
  if (init_type == 'varage') {
    scientists$productivity <- scientists_prod
  }

  return(list(papers = data.frame(papers), scientists = data.frame(scientists),
         scientistsAgeInit = scientistsAgeInit))

}

