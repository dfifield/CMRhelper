#'@export
#'
#'@title Extract marked individuals
#'
#'@description Extracts numbers first marked at each occasion and in each group
#'    from a model.
#'
#'@param x (required) a MARK model object as returned by
#'    \link[RMark]{mark}().
#'
#'@return
#'   Returns a 2D matrix of capture histories with one row per individual and one
#'   column per occasion.
#'
#'@author modified from XXXX by Greg Robertson - this was taken from the interweb ...

# use this version everywhere (remove extra args from calls to the other one)

Marked<-function(x) {

  n.groups <- ifelse(is.null(x$data$group.covariates), 1, nlevels(x$data$data$group))

  marked<-matrix(nrow=n.groups,ncol=x$data$nocc-1)

  for(g in 1:n.groups)
  {

    if(n.groups == 1){
      ch <- x$data$data$ch
    } else {
      ch <- x$data$data$ch[x$data$data$group == g]
    }

    for(i in 1:x$data$nocc-1)
    {
      ch1<-ch[(as.numeric(substr(ch,1,i)))==1]
      marked[g,i]<-length(ch1)
    }
  }
  return(marked)
}

#'@title Extract marked individuals - old way
#'
#'@description Extracts numbers first marked at each occasion and in each group.
#'
#'@param data (required) the data.
#'
#'@param n.occasions (required) number of occasions.
#'
#'@param groups (required) a vector of group labels.
#'
#'@details This function is deprecated and may go away.
#'
#'    The default param values (which were specific to the dipper data)
#'    have been removed. Also fixed a bug where "groups" was referred to as
#'    "Groups" in the function.
#'@return
#'   Returns a 2D matrix of capture histories with one row per individual and one
#'   column per occasion.
#'
#'@author ?? - from interweb - need to find source

Marked.3 <- function(data, n.occasions, groups)
{
  group<-data[,2]
  marked<-matrix(nrow=length(groups),ncol=n.occasions)
  for(g in 1:length(groups))
  {
    data_<-subset(data,group==groups[g])
    data_
    ch<-data_$ch
    for(i in 1:n.occasions)
    {
      ch1<-ch[(as.numeric(substr(ch,1,i)))==1]
      marked[g,i]<-length(ch1)
    }
  }
  return(marked)
}



#'
#'@title Simulate capture histories from a fitted CJS model
#'
#'@description This function takes a fitted CJS mark model, extracts relevant details
#'  (numbers marked and recaptured in each group(s) and estimated \phi and *p* values)
#'  and simulates a new capture history with the same structure and rates.
#'
#'@param x (required) a MARK model object as returned by
#'    \link[RMark]{mark}().
#'
#'@details
#'  This function is designed for relatively straight forward models, although it
#'  can handle age effects, as estimates are extracted directly from their PIMs.
#'  Any individual covariates are ignored, however covariates #'  in the design
#'  matrix would be included, as they are reflected in the PIMs.
#'
#'@return
#'  A simulated capture history formatted as a data.frame ready for analysis with
#'  usual RMARK functions. Generally sent to \link[RMark]{process.data}() as the
#'  next step.
#'
#'
#'@author
#'  Simulation piece taken from Kéry and Shaub's (2012) simul.cjs function -
#'  further modified by Sarah Gutowsky and Greg Robertson
#'

# add a new function to split simul.boot into two steps
# 1) extract needed stuff from model
# 2) do the simulation
# right now both are done with every call to simul.boot - but don't need
# to run 1) every time


extract.model<-function(x) {

  n.occasions <- x$nocc
  n.groups <- x$number.of.groups
  marked <- Marked(x)

  phi <- matrix(
    as.vector(
      t(
        t(
          sapply(1:n.groups, function (g) {
            as.vector(t(summary.mark(x)$reals$Phi[[g]]$pim))
          }
          )))), ncol = n.occasions - 1, byrow = TRUE)

  if(!is.null(x$parameters$Phi$fixed$value)  | (!is.null(x$fixed$index) & sum(x$fixed$index < (max(x$design.data$Phi$model.index)))> 0)){

    phi.index <- as.vector(sapply(1:n.groups, function (g) x$pims$Phi[[

    ]]$pim))

    phi.vec <- as.vector(sapply(1:n.groups, function (g) (phi[(1+((g-1)*(n.occasions-1))):((n.occasions - 1) * (g)), 1:(n.occasions - 1)])))

    phi.vec[which(phi.index %in% x$fixed$index)] <- x$fixed$value[which(x$fixed$index %in% phi.index)]

    phi.list <- lapply(1:n.groups, function (g) {
      matrix(phi.vec[(1+((g-1)*((n.occasions-1)^2))):((n.occasions - 1)^2 * (g))],
             ncol = n.occasions - 1, byrow = FALSE)
    })

    phi <- do.call(rbind, phi.list)

  }

  p <- matrix(
    as.vector(
      t(
        t(
          sapply(1:n.groups, function (g) {
            as.vector(t(summary.mark(x)$reals$p[[g]]$pim))
          }
          )))), ncol = n.occasions - 1, byrow = TRUE)

  if(!is.null(x$parameters$p$fixed$value)  | (!is.null(x$fixed$index) & sum(x$fixed$index > (max(x$design.data$Phi$model.index)))> 0)){

    p.index <- as.vector(sapply(1:n.groups, function (g) x$pims$p[[g]]$pim))

    p.vec <- as.vector(sapply(1:n.groups, function (g) (p[(1+((g-1)*(n.occasions-1))):((n.occasions - 1) * (g)), 1:(n.occasions - 1)])))

    p.vec[which(p.index %in% x$fixed$index)] <- x$fixed$value[which(x$fixed$index %in% p.index)]

    p.list <- lapply(1:n.groups, function (g) {

      p.hold <- p.vec[(1+((g-1)*((n.occasions-1)^2))):((n.occasions - 1)^2 * (g))]

      matrix(p.hold, ncol = n.occasions - 1, byrow = FALSE)

    })

    p <- do.call(rbind, p.list)

  }

  Phi <- phi[rep(1:nrow(phi), times = as.vector(t(marked))), ]

  P <- p[rep(1:nrow(p), times = as.vector(t(marked))), ]

  CH<-matrix(0,ncol=n.occasions,nrow=sum(marked))

  #define a vector with marking occasion
  mark.occ <- NULL

  for(i in 1:n.groups) {
    for (j in 1:(n.occasions -1)) {
      add <- rep(j, each = marked[i, j])
      mark.occ <- c(mark.occ,add)
    }
  }

  res <- list(marked = marked,
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
#'@description This function takes a fitted CJS mark model, extracts relevant details
#'  (numbers marked and recaptured in each group(s) and estimated \phi and *p* values)
#'  and simulates a new capture history with the same structure and rates.
#'
#'@param x (required) a MARK model object as returned by
#'    \link[RMark]{mark}().
#'
#'@details
#'  This function is designed for relatively straight forward models, although it
#'  can handle age effects, as estimates are extracted directly from their PIMs.
#'  Any individual covariates are ignored, however covariates #'  in the design
#'  matrix would be included, as they are reflected in the PIMs.
#'
#'@return
#'  A simulated capture history formatted as a data.frame ready for analysis with
#'  usual RMARK functions. Generally sent to \link[RMark]{process.data}() as the
#'  next step.
#'
#'
#'@author
#'  Simulation piece taken from Kéry and Shaub's (2012) simul.cjs function -
#'  further modified by Sarah Gutowsky and Greg Robertson
#'

# add a new function to split simul.boot into two steps
# 1) extract needed stuff from model
# 2) do the simulation
# right now both are done with every call to simul.boot - but don't need
# to run 1) every time

simul.boot <- function(extract) {

  a <- within(extract,{
    #fill in CH
    for (i in 1:sum(marked))
    {
      CH[i,mark.occ[i]]<-1
      if (mark.occ[i]==n.occasions) next
      for(t in (mark.occ[i]+1):n.occasions)
      {
        #survive?
        sur<-rbinom(1,1,Phi[i,t-1])
        if(sur==0) break #move to next
        #recaptured?
        rp<-rbinom(1,1,P[i,t-1])
        if(rp==1) CH[i,t]<-1
      } #t
    } #i
  })

  chout <- data.frame(ch = CMRhelper::pasty(a$CH), group = rep(1:a$n.groups, times = rowSums(a$marked)))
  return(chout)
}

#'
#'@title Perform a mark analysis on simulated capture histories
#'
#'@description ??
#'
#'@param x (required) a MARK model object as returned by
#'    \link[RMark]{mark}().
#'
#'@details
#'  Any pertinent details....
#'
#'@return
#'
#'  A data.frame "out" of mean deviance values along with 95% confidence limits
#'
#'@author
#'  Sarah Gutowksy, with further options by Greg Robertson
#'
sims<-function(x, reps, tsm = FALSE, ...)
{
  deviance<-dim(reps)
  for (i in 1:reps)
  {
    cat("iteration = ", iter <- i, "\n")

    if(x$number.of.groups == 1){
      sim.processed <- RMark::process.data(simul.boot(...), model="CJS")
      sim.ddl=RMark::make.design.data(sim.processed)
      global.sim<-RMark::mark(sim.processed,sim.ddl,
                              model.parameters=list(Phi=list(formula=~time),p=list(formula=~time)),
                              output=FALSE,silent=TRUE)
    } else {
      sim.processed <- RMark::process.data(simul.boot(...), model="CJS",groups = "group")
      sim.ddl=RMark::make.design.data(sim.processed)
      if(tsm == TRUE){
        sim.ddl <- RMark::add.design.data(sim.processed, sim.ddl,
                                          parameter="Phi", type="age", bins=c(0,1, 17),name="tsm",
                                          right = FALSE, replace = TRUE)
        global.sim<-RMark::mark(sim.processed,sim.ddl,
                                model.parameters=list(Phi=list(formula=~tsm+time*group),p=list(formula=~time*group)),
                                output=FALSE,silent=TRUE)
      }  else {
        global.sim<-RMark::mark(sim.processed,sim.ddl,
                                model.parameters=list(Phi=list(formula=~time*group),p=list(formula=~time*group)),
                                output=FALSE,silent=TRUE)
      }
    }

    deviance[i]<-global.sim$results$lnl
  }
  out<-list(deviance.mean=mean(deviance),deviance.025=quantile(deviance,0.025),deviance.975=quantile(deviance,0.975))
}

#'
#'@export
#'
#'@title Calculate a bootstrapped c-hat value from a fitted CJS mark model
#'
#'@description This function takes a fitted mark model and passes it to
#'    \link[CMRhelper]{simul.boot}() to extract relevant information from the
#'    model simulates a capture history based on the model. This process is
#'    repeated @param reps times and a bootstrapped c-hat value is calculated.
#'
#'@param x (required) a MARK model object as returned by
#'    \link[RMark]{mark}().
#'
#'@param reps (required)
#'
#'@param tsm (optional)
#'
#'@details
#'  Any pertinent details....
#'
#'@return
#' Returns a bootstrapped c.hat.
#'
#'@author
#'  Sarah Gutowsky, Greg Robertson
#'
bootstrap.deviance <- function(x, reps, tsm = FALSE) {

  extract.list <- extract.model(x)
  sim.out <- sims(x, reps, tsm, extract.list)
  fetch.deviance <- function(y) y$results$lnl
  data.deviance <- fetch.deviance(x)
  sim.ci<-c(sim.out$deviance.025,sim.out$deviance.975)
  cat("data deviance = ",data.deviance, ", simulation mean = ", sim.out$deviance.mean, ", simulation 95%CI = ", sim.ci,  "\n")

  #modified c.hat
  c.hat<-data.deviance/sim.out$deviance.mean
  cat("modified c.hat =",c.hat,  "\n")

  cbind(c.hat, data.deviance/sim.out$deviance.975, data.deviance/sim.out$deviance.025)
}
