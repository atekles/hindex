#' Plot the result of simulate_hindex
#'
#' Plot the result of a simulation computed by simulate_hindex.
#'
#' @param simdata The result of a simulation returned
#' by \code{\link{simulate_hindex}}.
#' @param plot_hindex If this parameter is set to TRUE, the h-index values are
#' plotted.
#' @param plot_halpha If this parameter is set to TRUE, the h-alpha values are
#' plotted.
#' @param plot_toppapers If this parameter is set to TRUE, the numbers of
#' top-10\% papers are plotted.
#' @param plot_mindex If this parameter is set to TRUE, the mindex values are
#' plotted.
#' @param subgroups If this parameter is set to TRUE, the subgroups in simdata
#' are considered in the visualization.
#' @param group_boundaries A list of vectors or a vector of integers specifying
#' the groups for plotting the h-index/h-alpha values separately for each of
#' these groups. The groups are defined based on the initial h-index of the
#' agents. If a list is specified, each element must be a vector of length 2
#' representing the lower and the upper bound for the initial h-index (if the
#' boundaries are included in the corresponding intervals is specified by the
#' exclude_group_boundaries parameter).
#' If a vector of integers is specified, each element in group_boundaries
#' separates two groups such that all agents with an initial h-index below
#' this boundary (and equal to or above any lower boundary; if
#' exclude_group_boundaries is set to TRUE, the initial h-index has to be
#' above any lower boundary) are in the first group, and all agents with
#' an initial h-index equal to or above this boundary (and below any higher
#' boundary) are in the second group.
#' @param plot_group_diffs If this parameter is specified, the difference
#' between the groups that are specified by group_boundaries is plotted.
#' @param exclude_group_boundaries If this parameter is set to TRUE, the
#' scientists are grouped such that those scientists whose initial h-index
#' is equal to a boundary are not included.
#'
#' @return A ggplot object (\code{\link[ggplot2]{ggplot}}).
#'
#' @export
#' @importFrom foreach "%do%"
#'
#' @examples
#' set.seed(123)
#' simdata <- simulate_hindex(runs = 2, n = 20, periods = 3)
#' plot_hsim(simdata, plot_hindex = TRUE, plot_halpha = TRUE)
plot_hsim <- function(simdata, plot_hindex = FALSE, plot_halpha = FALSE,
                      plot_toppapers = FALSE,
                      plot_mindex = FALSE,
                      subgroups = FALSE,
                      group_boundaries = NULL,
                      exclude_group_boundaries = FALSE,
                      plot_group_diffs = FALSE) {

  if (!plot_hindex & !plot_halpha & !plot_toppapers & !plot_mindex) {
    stop('at least one of the parameters plot_hindex, plot_halpha,
         plot_toppapers, plot_mindex must be TRUE')
  }

  if (subgroups && (!is.null(group_boundaries) || exclude_group_boundaries)) {
    warning('if subgroups is set to TRUE, group_boundaries and exclude_group_boundaries
            have no effect')
  }

  if (typeof(simdata) != 'list' ||
      length(unique(lapply(simdata[1:4], length))) != 1 ||
      length(unique(unlist(
        lapply(simdata[1:4], function(currentRun) {lapply(currentRun, length)})
      ))) != 1 ||
      length(simdata[[1]][[1]][[1]]) != length(simdata[[5]])) {
    # check if each run has same no of periods, and each period has
    # same no of returned lists (hindex values, h alpha values, toppapers, etc.)
    stop('structure of simdata not correct; make sure to use the result
         returned by simulate_hindex')
  }

  # groups -> initial h values (e.g. low vs high...)
  if (!is.null(group_boundaries)) {
    if (typeof(group_boundaries) == 'list') {

      # use vector for each group indicating lower and upper boundary

      groups <- length(group_boundaries)

      groupBoundariesOrdered <- lapply(group_boundaries, function(currentEl) {
        if (!typeof(currentEl) %in% c('double', 'integer') ||
            length(currentEl) != 2) {
          stop('if group_boundaries are given as a list, each element
               has to be a vector of length 2')
        }
        if (currentEl[1] <= currentEl[2]) {
          return(currentEl)
        } else {
          return(c(currentEl[2], currentEl[1]))
        }
      })

      groupOrder <- order(unlist(lapply(groupBoundariesOrdered, `[`, 1)))
      groupBoundariesOrdered <- groupBoundariesOrdered[groupOrder]

    } else if (typeof(group_boundaries) %in% c('double', 'integer')) {

      # use numeric to indicate the threshold for determining groups

      group_boundaries_ordered <-
        unique(group_boundaries[order(group_boundaries)])
      groups <- length(group_boundaries_ordered) + 1

      groupBoundariesOrdered <- lapply(1:groups, FUN = function(currentGroup) {
        if (currentGroup == 1) {
          return(c(0, group_boundaries_ordered[1] - .Machine$double.eps))
        } else if (currentGroup == groups) {
          return(c(group_boundaries_ordered[currentGroup - 1], Inf))
        } else {
          (c(group_boundaries_ordered[currentGroup - 1],
             group_boundaries_ordered[currentGroup] - .Machine$double.eps))
        }
      })

    } else if (typeof(group_boundaries) == 'character') {

      # use function to determine the threshold between groups

      if (group_boundaries == 'median') {
        threshold <- stats::median(simdata$h[[1]][[1]])
        groups <- 2
        groupBoundariesOrdered <- list(c(0, threshold), c(threshold, Inf))
      } else {
        stop('function for determining groups not supported')
      }

    } else {
      stop('if group_boundaries is specified,
           it must be either a list or a vector')
    }
  } else {
    groups <- 1
  }

  if (subgroups) {
    groups <- length(unique(simdata$subgroup))
  }

  # plot_group_diffs: plot the difference of index/indices for two groups
  # two groups have to be specified for this
  if (plot_group_diffs) {
    if (is.null(group_boundaries) && !subgroups) {
      stop('in order to plot differences between groups of scientists with
           different initial h-index, group_boundaries must be specified')
    }
    if (groups != 2) {
      stop('differences between groups of scientists is only possible for two groups')
    }
  }

  if (plot_hindex) {

    nPeriods <- length(simdata$h[[1]])

    hRun <- NULL
    hMeansRuns <- foreach::foreach(hRun = simdata$h,
                                   .combine = '+') %do% {

       if (groups > 1) {
         if (subgroups) {
           groupsScientists <- lapply(1:groups, function(currentGroup) {which(simdata$subgroup == currentGroup)})
         } else {
           groupsScientists <- lapply(1:groups, function(currentGroup) {
             # TODO
             if (exclude_group_boundaries) {
               currentScientists <- which(hRun[[1]] >
                                            groupBoundariesOrdered[[currentGroup]][1] &
                                            hRun[[1]] < groupBoundariesOrdered[[currentGroup]][2])
               if (length(currentScientists) <= 0) {
                 warning(paste('no scientist in group', currentGroup, sep = ' '))
               }
             } else {
               currentScientists <- which(hRun[[1]] >=
                                            groupBoundariesOrdered[[currentGroup]][1] &
                                            hRun[[1]] <= groupBoundariesOrdered[[currentGroup]][2])
               if (length(currentScientists) <= 0) {
                 warning(paste('no scientist in group', currentGroup, sep = ' '))
               }
             }
             return(currentScientists)
           })
         }
       } else {
         groupsScientists <- list(1:length(hRun[[1]]))
       }

       hPeriod <- NULL
       runMeans <- foreach::foreach(hPeriod = hRun, .combine = 'cbind') %do% {

         # get mean for each group
         runPeriodGroupMeans <- vapply(1:groups, FUN = function(currentGroup) {
           mean(hPeriod[groupsScientists[[currentGroup]]])
         }, FUN.VALUE = double(1))
         return(t(runPeriodGroupMeans))

       }

       return(runMeans)

    }

    hMeansRuns <- hMeansRuns / length(simdata$h)
    vals <- t(hMeansRuns)
    period <-rep(1:nPeriods, each = groups)
    hInitGroup <- as.integer(round(rep(1:groups, nPeriods)))
    indexType <- 'h-index'

    plotData <- data.frame(vals = vals, period = period,
                           hInitGroup = hInitGroup, indexType = indexType)

    if (plot_group_diffs) {
      diffs <- vals[c(FALSE, TRUE)] - vals[c(TRUE, FALSE)]
      plotData <- rbind(plotData,
                        data.frame(vals = diffs, period = 1:nPeriods,
                                   hInitGroup = groups + 1, indexType = indexType))
    }

  }

  if (plot_halpha) {

    nPeriods <- length(simdata$h[[1]])

    hAlphaRunIndex <- 0
    hAlphaMeansRuns <- foreach::foreach(hAlphaRunIndex = 1:length(simdata$h_alpha),
      .combine = '+') %do% {

       # hAlphaRunIndex <- hAlphaRunIndex + 1

       if (groups > 1) {
         if (subgroups) {
           groupsScientists <- lapply(1:groups, function(currentGroup) {which(simdata$subgroup == currentGroup)})
         } else {
           groupsScientists <- lapply(1:groups, function(currentGroup) {
             # TODO
             if (exclude_group_boundaries) {
               currentScientists <- which(simdata$h[[hAlphaRunIndex]][[1]] >
                                            groupBoundariesOrdered[[currentGroup]][1] &
                                            simdata$h[[hAlphaRunIndex]][[1]] <
                                            groupBoundariesOrdered[[currentGroup]][2])
               if (length(currentScientists) <= 0) {
                 warning(paste('no scientist in group', currentGroup, sep = ' '))
               }
             } else {
               currentScientists <- which(simdata$h[[hAlphaRunIndex]][[1]] >=
                                            groupBoundariesOrdered[[currentGroup]][1] &
                                            simdata$h[[hAlphaRunIndex]][[1]] <=
                                            groupBoundariesOrdered[[currentGroup]][2])
               if (length(currentScientists) <= 0) {
                 warning(paste('no scientist in group', currentGroup, sep = ' '))
               }
             }
             return(currentScientists)
           })
         }
       } else {
         groupsScientists <- list(1:length(simdata$h[[1]][[1]]))
       }

       hAlphaPeriod <- NULL
       runMeans <- foreach::foreach(hAlphaPeriod = simdata$h_alpha[[hAlphaRunIndex]], .combine = 'cbind') %do% {

         # get mean for each group
         runPeriodGroupMeans <- vapply(1:groups, FUN = function(currentGroup) {
           mean(hAlphaPeriod[groupsScientists[[currentGroup]]])
         }, FUN.VALUE = double(1))
         return(t(runPeriodGroupMeans))

       }

       return(runMeans)

    }

    hAlphaMeansRuns <- hAlphaMeansRuns / length(simdata$h_alpha)
    vals <- t(hAlphaMeansRuns)
    period <-rep(1:nPeriods, each = groups)
    hInitGroup <- as.integer(round(rep(1:groups, nPeriods)))
    indexType <- 'h-alpha'

    if (plot_hindex) {
      plotData <- rbind(plotData, data.frame(vals = vals, period = period,
                                             hInitGroup = hInitGroup, indexType = indexType))
    } else {
      plotData <- data.frame(vals = vals, period = period,
                             hInitGroup = hInitGroup, indexType = indexType)
    }

    if (plot_group_diffs) {
      diffs <- vals[c(FALSE, TRUE)] - vals[c(TRUE, FALSE)]
      plotData <- rbind(plotData,
                        data.frame(vals = diffs, period = 1:nPeriods,
                                   hInitGroup = groups + 1, indexType = indexType))
    }

  }

  if (plot_toppapers) {

    nPeriods <- length(simdata$h[[1]])

    toppapersRunIndex <- 0
    toppapersMeansRuns <- foreach::foreach(toppapersRunIndex = 1:length(simdata$top10_papers),
                                        .combine = '+') %do% {

                                          # toppapersRunIndex <- toppapersRunIndex + 1

                                          if (groups > 1) {
                                            if (subgroups) {
                                              groupsScientists <- lapply(1:groups, function(currentGroup) {which(simdata$subgroup == currentGroup)})
                                            } else {
                                              groupsScientists <- lapply(1:groups, function(currentGroup) {
                                                # TODO
                                                if (exclude_group_boundaries) {
                                                  currentScientists <- which(simdata$h[[toppapersRunIndex]][[1]] >
                                                                               groupBoundariesOrdered[[currentGroup]][1] &
                                                                               simdata$h[[toppapersRunIndex]][[1]] <
                                                                               groupBoundariesOrdered[[currentGroup]][2])
                                                  if (length(currentScientists) <= 0) {
                                                    warning(paste('no scientist in group', currentGroup, sep = ' '))
                                                  }
                                                } else {
                                                  currentScientists <- which(simdata$h[[toppapersRunIndex]][[1]] >=
                                                                               groupBoundariesOrdered[[currentGroup]][1] &
                                                                               simdata$h[[toppapersRunIndex]][[1]] <=
                                                                               groupBoundariesOrdered[[currentGroup]][2])
                                                  if (length(currentScientists) <= 0) {
                                                    warning(paste('no scientist in group', currentGroup, sep = ' '))
                                                  }
                                                }
                                                return(currentScientists)
                                              })
                                            }
                                          } else {
                                            groupsScientists <- list(1:length(simdata$h[[1]][[1]]))
                                          }

                                          toppapersPeriod <- NULL
                                          runMeans <- foreach::foreach(toppapersPeriod = simdata$top10_papers[[toppapersRunIndex]], .combine = 'cbind') %do% {

                                            # get mean for each group
                                            runPeriodGroupMeans <- vapply(1:groups, FUN = function(currentGroup) {
                                              mean(toppapersPeriod[groupsScientists[[currentGroup]]])
                                            }, FUN.VALUE = double(1))
                                            return(t(runPeriodGroupMeans))

                                          }

                                          return(runMeans)

                                        }

    toppapersMeansRuns <- toppapersMeansRuns / length(simdata$top10_papers)
    vals <- t(toppapersMeansRuns)
    period <-rep(1:nPeriods, each = groups)
    hInitGroup <- as.integer(round(rep(1:groups, nPeriods)))
    indexType <- 'toppapers'

    if (plot_hindex || plot_halpha) {
      plotData <- rbind(plotData, data.frame(vals = vals, period = period,
                                             hInitGroup = hInitGroup, indexType = indexType))
    } else {
      plotData <- data.frame(vals = vals, period = period,
                             hInitGroup = hInitGroup, indexType = indexType)
    }

    if (plot_group_diffs) {
      diffs <- vals[c(FALSE, TRUE)] - vals[c(TRUE, FALSE)]
      plotData <- rbind(plotData,
                        data.frame(vals = diffs, period = 1:nPeriods,
                                   hInitGroup = groups + 1, indexType = indexType))
    }

  }

  if (plot_mindex) {

    nPeriods <- length(simdata$h[[1]])

    mindexRunIndex <- 0
    mindexMeansRuns <- foreach::foreach(mindexRunIndex = 1:length(simdata$mindex),
                                           .combine = '+') %do% {

                                             # mindexRunIndex <- mindexRunIndex + 1

                                             if (groups > 1) {
                                               if (subgroups) {
                                                 groupsScientists <- lapply(1:groups, function(currentGroup) {which(simdata$subgroup == currentGroup)})
                                               } else {
                                                 groupsScientists <- lapply(1:groups, function(currentGroup) {
                                                   # TODO
                                                   if (exclude_group_boundaries) {
                                                     currentScientists <- which(simdata$h[[mindexRunIndex]][[1]] >
                                                                                  groupBoundariesOrdered[[currentGroup]][1] &
                                                                                  simdata$h[[mindexRunIndex]][[1]] <
                                                                                  groupBoundariesOrdered[[currentGroup]][2])
                                                     if (length(currentScientists) <= 0) {
                                                       warning(paste('no scientist in group', currentGroup, sep = ' '))
                                                     }
                                                   } else {
                                                     currentScientists <- which(simdata$h[[mindexRunIndex]][[1]] >=
                                                                                  groupBoundariesOrdered[[currentGroup]][1] &
                                                                                  simdata$h[[mindexRunIndex]][[1]] <=
                                                                                  groupBoundariesOrdered[[currentGroup]][2])
                                                     if (length(currentScientists) <= 0) {
                                                       warning(paste('no scientist in group', currentGroup, sep = ' '))
                                                     }
                                                   }
                                                   return(currentScientists)
                                                 })
                                               }
                                             } else {
                                               groupsScientists <- list(1:length(simdata$h[[1]][[1]]))
                                             }

                                             mindexPeriod <- NULL
                                             runMeans <- foreach::foreach(mindexPeriod = simdata$mindex[[mindexRunIndex]], .combine = 'cbind') %do% {

                                               # get mean for each group
                                               runPeriodGroupMeans <- vapply(1:groups, FUN = function(currentGroup) {
                                                 mean(mindexPeriod[groupsScientists[[currentGroup]]])
                                               }, FUN.VALUE = double(1))
                                               return(t(runPeriodGroupMeans))

                                             }

                                             return(runMeans)

                                           }

    mindexMeansRuns <- mindexMeansRuns / length(simdata$mindex)
    vals <- t(mindexMeansRuns)
    period <-rep(1:nPeriods, each = groups)
    hInitGroup <- as.integer(round(rep(1:groups, nPeriods)))
    indexType <- 'mindex'

    if (plot_hindex || plot_halpha || plot_toppapers) {
      plotData <- rbind(plotData, data.frame(vals = vals, period = period,
                                             hInitGroup = hInitGroup, indexType = indexType))
    } else {
      plotData <- data.frame(vals = vals, period = period,
                             hInitGroup = hInitGroup, indexType = indexType)
    }

    if (plot_group_diffs) {
      diffs <- vals[c(FALSE, TRUE)] - vals[c(TRUE, FALSE)]
      plotData <- rbind(plotData,
                        data.frame(vals = diffs, period = 1:nPeriods,
                                   hInitGroup = groups + 1, indexType = indexType))
    }

  }

  # initialize ggplot
  sim_plot <- ggplot2::ggplot(plotData,
                              ggplot2::aes(x = period, y = vals,
                                           group = interaction(hInitGroup, indexType),
                                           color = factor(hInitGroup),
                                           linetype = indexType)) +
    ggplot2::geom_line() +
    ggplot2::xlab('Period') +
    ggplot2::theme(axis.title.y = ggplot2::element_blank())

  if (groups > 1) {

    if (subgroups) {
      labels <- as.character(1:groups)
    } else {
      labels <- vector(mode = 'character', length = groups)
      for (currentGroup in 1:groups) {
        if (exclude_group_boundaries) {
          labels[currentGroup] <- eval(bquote(expression(
            h[init] %in% group("(", list(
              .(groupBoundariesOrdered[[currentGroup]][1]),
              .(groupBoundariesOrdered[[currentGroup]][2])
            ), ")")
          )))
        } else {
          labels[currentGroup] <- eval(bquote(expression(
            h[init] %in% group("[", list(
              .(groupBoundariesOrdered[[currentGroup]][1]),
              .(groupBoundariesOrdered[[currentGroup]][2])
            ), ")")
          )))
        }
      }
      if (typeof(group_boundaries) %in% c('double', 'integer')) {
        labels[1] <- eval(bquote(expression(
          h[init] < .(groupBoundariesOrdered[[1]][2])
        )))
        if (exclude_group_boundaries) {
          labels[groups] <- eval(bquote(expression(
            h[init] > .(groupBoundariesOrdered[[groups]][1])
          )))
        } else {
          labels[groups] <- eval(bquote(expression(
            h[init] >= .(groupBoundariesOrdered[[groups]][1])
          )))
        }
      }
    }

    if (plot_group_diffs) {
      labels <- c(labels, 'Difference')
    }

    sim_plot <- sim_plot +
      ggplot2::scale_color_discrete(name = 'Initial h-index',
                                    labels = labels)
  } else {
    sim_plot <- sim_plot +
      ggplot2::scale_color_discrete(guide = FALSE)
  }

  if (plot_hindex + plot_halpha + plot_toppapers == 1) {
    sim_plot <- sim_plot +
      ggplot2::scale_linetype_discrete(guide = FALSE)
  } else {
    sim_plot <- sim_plot +
      ggplot2::scale_linetype_discrete(name = 'Index')
  }

  return(sim_plot)

}
