#'@export
#'
#'@title Extract number of marked individuals
#'
#'@description Extracts numbers first marked at each occasion and in each group
#'    from a model.
#'
#'@param x (required) a MARK model object as returned by
#'    \code{\link[RMark]{mark}()}.
#'
#'@return
#'   Returns a 2D matrix of numbers first marked in each group (row) by occasion
#'   (column).
#'
#'@author Greg Robertson
extract.n.marked <- function(x) {
  n.groups <-
    ifelse(is.null(x$data$group.covariates),
           1,
           nlevels(x$data$data$group))


  if (n.groups == 1) {
    marked <- rep(0, x$nocc)
    marked[as.numeric(dimnames(table(regexpr('1', x$data$data$ch)))[[1]])] <- t(table(regexpr('1', x$data$data$ch)))
    marked <- t(marked[-x$nocc])
  } else {
      marked <- matrix(0, nrow = n.groups, ncol = x$nocc)
      for(g in 1:n.groups){
        marked[g,as.numeric(dimnames(table(regexpr('1', x$data$data$ch[x$data$data$group == g])))[[1]])] <- t(table(regexpr('1', x$data$data$ch[x$data$data$group == g])))
      }
      marked <- marked[,-x$nocc]
  }

  return(marked)
}

#'
#'@title Extract values from a Mark object
#'
#'@description This function takes a fitted CJS Mark model, and extracts relevant details
#'  (numbers marked and recaptured in each group(s) and estimated \eqn{\phi} and *p* values)
#'   needed to simulate capture histories.
#'
#'
#'@param x (required) a MARK model object as returned by
#'    \code{\link[RMark]{mark}()}.
#'
#'@details
#'  This function is designed for relatively straight forward models, although it
#'  can handle age effects, as estimates are extracted directly from their PIMs.
#'  Any individual covariates are ignored, however covariates in the design
#'  matrix would be included, as they are reflected in the PIMs.
#'
#'@return
#'  A list with needed values to simulate a capture history. Generally sent to
#'  \code{\link[CMRhelper]{simul.boot}()} as the next step.
#'
#'
#'@author
#'  Greg Robertson
#'
extract.model <- function(x) {
  n.occasions <- x$nocc
  n.groups <- x$number.of.groups
  marked <- extract.n.marked(x)

  phi <- matrix(as.vector(t(t(
    sapply(1:n.groups, function (g) {
      as.vector(t(summary.mark(x)$reals$Phi[[g]]$pim))
    })
  ))), ncol = n.occasions - 1, byrow = TRUE)

  if (!is.null(x$parameters$Phi$fixed$value)  |
      (!is.null(x$fixed$index) &
       sum(x$fixed$index < (max(
         x$design.data$Phi$model.index
       ))) > 0)) {
    phi.index <- as.vector(sapply(1:n.groups, function (g)
      x$pims$Phi[[]]$pim))

    phi.vec <-
      as.vector(sapply(1:n.groups, function (g)
        (phi[(1 + ((g - 1) * (n.occasions - 1))):((n.occasions - 1) * (g)), 1:(n.occasions - 1)])))

    phi.vec[which(phi.index %in% x$fixed$index)] <-
      x$fixed$value[which(x$fixed$index %in% phi.index)]

    phi.list <- lapply(1:n.groups, function (g) {
      matrix(phi.vec[(1 + ((g - 1) * ((
        n.occasions - 1
      ) ^ 2))):((n.occasions - 1) ^ 2 * (g))],
      ncol = n.occasions - 1, byrow = FALSE)
    })

    phi <- do.call(rbind, phi.list)

  }

  p <- matrix(as.vector(t(t(
    sapply(1:n.groups, function (g) {
      as.vector(t(summary.mark(x)$reals$p[[g]]$pim))
    })
  ))), ncol = n.occasions - 1, byrow = TRUE)

  if (!is.null(x$parameters$p$fixed$value)  |
      (!is.null(x$fixed$index) &
       sum(x$fixed$index > (max(
         x$design.data$Phi$model.index
       ))) > 0)) {
    p.index <-
      as.vector(sapply(1:n.groups, function (g)
        x$pims$p[[g]]$pim))

    p.vec <-
      as.vector(sapply(1:n.groups, function (g)
        (p[(1 + ((g - 1) * (n.occasions - 1))):((n.occasions - 1) * (g)), 1:(n.occasions - 1)])))

    p.vec[which(p.index %in% x$fixed$index)] <-
      x$fixed$value[which(x$fixed$index %in% p.index)]

    p.list <- lapply(1:n.groups, function (g) {
      p.hold <-
        p.vec[(1 + ((g - 1) * ((
          n.occasions - 1
        ) ^ 2))):((n.occasions - 1) ^ 2 * (g))]

      matrix(p.hold, ncol = n.occasions - 1, byrow = FALSE)

    })

    p <- do.call(rbind, p.list)

  }

  Phi <- phi[rep(1:nrow(phi), times = as.vector(t(marked))), ]

  P <- p[rep(1:nrow(p), times = as.vector(t(marked))), ]

  CH <- matrix(0, ncol = n.occasions, nrow = sum(marked))

  #define a vector with marking occasion
  mark.occ <- NULL

  for (i in 1:n.groups) {
    for (j in 1:(n.occasions - 1)) {
      add <- rep(j, each = marked[i, j])
      mark.occ <- c(mark.occ, add)
    }
  }

  res <- list(
    marked = marked,
    Phi = Phi,
    P = P,
    n.occasions = n.occasions,
    n.groups = n.groups,
    mark.occ = mark.occ,
    CH = CH
  )

  return(res)

}

#'
#'@title Simulate capture histories from a fitted CJS model
#'
#'@description Simulates a new capture history with the same structure and rates.
#'
#'@param extract (required) a list as returned from \code{\link[CMRhelper]{extract.model}()}.
#'
#'@details
#'  This function is designed for relatively straight forward models, although it
#'  can handle age effects, as estimates are extracted directly from their PIMs.
#'  Any individual covariates are ignored, however covariates in the design
#'  matrix would be included, as they are reflected in the PIMs.
#'
#'@return
#'  A simulated capture history formatted as a data.frame ready for analysis with
#'  usual RMARK functions. Generally sent to \code{\link[RMark]{process.data}()} as the
#'  next step.
#'
#'
#'@author
#'  Simulation piece largely taken from Kéry and Shaub's (2012) simul.cjs function
#'  (see \code{\link[CMRhelper]{simul.cjs}()}) -
#'  further modified by Sarah Gutowsky and Greg Robertson.
#'
simul.boot <- function(extract) {
  a <- within(extract, {
    #fill in CH
    for (i in 1:sum(marked))
    {
      CH[i, mark.occ[i]] <- 1
      if (mark.occ[i] == n.occasions)
        next
      for (t in (mark.occ[i] + 1):n.occasions)
      {
        #survive?
        sur <- rbinom(1, 1, Phi[i, t - 1])
        if (sur == 0)
          break #move to next
        #recaptured?
        rp <- rbinom(1, 1, P[i, t - 1])
        if (rp == 1)
          CH[i, t] <- 1
      } #t
    } #i
  })

  chout <-
    data.frame(ch = CMRhelper::pasty(a$CH),
               group = rep(1:a$n.groups, times = rowSums(a$marked)))
  return(chout)
}

#'
#'@title Perform repeated mark analyses on simulated capture histories
#'
#'@description Use RMark to setup and run multiple analyses
#'
#'@param x (required) a MARK model object as returned by
#'    \code{\link[RMark]{mark}()}.
#'
#'@param reps (required) number of repetitions
#'
#'@param tsm (optional, default \code{FALSE}) whether to include a time since
#'  marking effect in the Mark models fit to the simulated capture histories.
#'
#'@details
#'  Uses \code{\link[CMRhelper]{simul.boot}()} to generate simulated capture histories
#'  from a fitted model. Then uses \code{RMark} to fit models to these simulated
#'  histories.
#'
#'@return
#'  A list containing mean deviance along with 95\% confidence limits.
#'
#'@author
#'  Sarah Gutowksy, with further options by Greg Robertson
#'
sims <- function(x, reps, tsm = FALSE)
{
  # Get model elements needed for simul.boot()
  extract.list <- extract.model(x)

  deviance <- dim(reps)
  for (i in 1:reps)
  {
    cat("iteration = ", i, "\n")

    if (x$number.of.groups == 1) {
      sim.processed <-
        RMark::process.data(simul.boot(extract.list), model = "CJS")
      sim.ddl = RMark::make.design.data(sim.processed)
      global.sim <- RMark::mark(
        sim.processed,
        sim.ddl,
        model.parameters = list(
          Phi = list(formula =  ~ time),
          p = list(formula =  ~ time)
        ),
        output = FALSE,
        silent = TRUE
      )
    } else {
      sim.processed <-
        RMark::process.data(simul.boot(extract.list),
                            model = "CJS",
                            groups = "group")
      sim.ddl = RMark::make.design.data(sim.processed)
      if (tsm == TRUE) {
        sim.ddl <- RMark::add.design.data(
          sim.processed,
          sim.ddl,
          parameter = "Phi",
          type = "age",
          bins = c(0, 1, 17),
          name = "tsm",
          right = FALSE,
          replace = TRUE
        )
        global.sim <- RMark::mark(
          sim.processed,
          sim.ddl,
          model.parameters = list(
            Phi = list(formula =  ~ tsm + time * group),
            p = list(formula =  ~ time * group)
          ),
          output = FALSE,
          silent = TRUE
        )
      }  else {
        global.sim <- RMark::mark(
          sim.processed,
          sim.ddl,
          model.parameters = list(
            Phi = list(formula =  ~ time * group),
            p = list(formula =  ~ time * group)
          ),
          output = FALSE,
          silent = TRUE
        )
      }
    }

    deviance[i] <- global.sim$results$lnl
  }
  out <-
    list(
      deviance.mean = mean(deviance),
      deviance.025 = quantile(deviance, 0.025),
      deviance.975 = quantile(deviance, 0.975)
    )
}

#'
#'@export
#'
#'@title Calculate a bootstrapped c-hat from a fitted CJS Mark model
#'
#'@description This function takes a fitted mark model and computes bootstrapped
#'    mean (and 95% CI) c-hat.
#'
#'@param x (required) a MARK model object as returned by
#'    \code{\link[RMark]{mark}()}.
#'
#'@param reps (required) number of repetitions
#'
#'@param tsm (optional, default \code{FALSE})  whether to include a time since
#'  marking effect in the Mark models fit to the simulated capture histories.
#'
#'@details
#'    Makes use of \code{\link[CMRhelper]{sims}()} to simulate and
#'    analyse a set capture histories based on the model. This process is
#'    repeated \code{reps} times and a bootstrapped c-hat value is calculated.
#'
#'@return
#'    Returns a 1x3 matrix with columns for: bootstrapped c.hat, upper (97.5%), and lower(2.5%)
#'    confidence limits.
#'
#'@author
#'  Sarah Gutowsky, Greg Robertson
#'
bootstrap.deviance <- function(x, reps, tsm = FALSE) {
  sim.out <- sims(x, reps, tsm)
  data.deviance <- x$results$lnl
  sim.ci <- c(sim.out$deviance.025, sim.out$deviance.975)
  cat(
    "data deviance = ",
    data.deviance,
    ", simulation mean = ",
    sim.out$deviance.mean,
    ", simulation 95%CI = ",
    sim.ci,
    "\n"
  )

  #modified c.hat
  c.hat <- data.deviance / sim.out$deviance.mean
  cat("modified c.hat =", c.hat,  "\n")

  cbind(c.hat,
        data.deviance / sim.out$deviance.975,
        data.deviance / sim.out$deviance.025)
}
