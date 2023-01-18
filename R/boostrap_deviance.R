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
#'@author Greg Robertson based on 3 argument version from ????

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

#'
#'@title ??
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
#'@author
#'  ???
#'
simul.boot<-function(x) {

  n.occasions <- x$data$nocc
  n.groups <- ifelse(is.null(x$data$group.covariates), 1, nlevels(x$data$data$group))
  marked <- Marked(x)

  phi <- matrix(
    as.vector(
      t(
        t(
          sapply(1:n.groups, function (g) {
            as.vector(t(summary.mark(x)$reals$Phi[[g]]$pim))
            }
            )))), ncol = n.occasions - 1, byrow = TRUE)

  if(!is.null(x$parameters$Phi$fixed$value)){

    phi.index <- as.vector(sapply(1:n.groups, function (g) x$pims$Phi[[g]]$pim))

    phi.vec <- as.vector(sapply(1:n.groups, function (g) (phi[(1+((g-1)*(n.occasions-1))):((n.occasions - 1) * (g)), 1:(n.occasions - 1)])))

    phi.vec[which(phi.index %in% x$fixed$index)] <- x$fixed$value

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

  if(!is.null(x$parameters$p$fixed$value)){

    p.index <- as.vector(sapply(1:n.groups, function (g) x$pims$p[[g]]$pim))

    p.vec <- as.vector(sapply(1:n.groups, function (g) (p[(1+((g-1)*(n.occasions-1))):((n.occasions - 1) * (g)), 1:(n.occasions - 1)])))

    p.vec[which(p.index %in% x$fixed$index)] <- x$fixed$value

    # p <- matrix(p.vec, ncol = n.occasions - 1, byrow = FALSE)

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
  mark.occ <- rep(rep(1:dim(marked)[2], n.groups), times = as.vector(t(marked)))


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

  chout <- data.frame(ch = pasty(CH), group = rep(1:n.groups, times = rowSums(marked)))
  return(chout)
}


#'
#'@title ??
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
#'@author
#'  ???
#'
sims<-function(x, reps, tsm = FALSE)
{
  deviance<-dim(reps)
  for (i in 1:reps)
  {
    cat("iteration = ", iter <- i, "\n")

    if(is.null(x$data$group.covariates)){
      sim.processed <- process.data(simul.boot(x), model="CJS")
      sim.ddl=make.design.data(sim.processed)
      global.sim<-mark(sim.processed,sim.ddl,
                       model.parameters=list(Phi=list(formula=~time),p=list(formula=~time)),
                       output=FALSE,silent=TRUE)
    } else {
      sim.processed <- process.data(simul.boot(x), model="CJS",groups = "group")
      sim.ddl=make.design.data(sim.processed)
      if(tsm == TRUE){
        sim.ddl <- add.design.data(sim.processed, sim.ddl,
                                   parameter="Phi", type="age", bins=c(0,1, 17),name="tsm",
                                   right = FALSE, replace = TRUE)
        global.sim<-mark(sim.processed,sim.ddl,
                         model.parameters=list(Phi=list(formula=~tsm+time*group),p=list(formula=~time*group)),
                         output=FALSE,silent=TRUE)
      }  else {
        global.sim<-mark(sim.processed,sim.ddl,
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
#'@title ??
#'
#'@description ??
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
#'@return ??
#'
#'@author
#'  Greg Robertson?
#'
bootstrap.deviance <- function(x, reps, tsm = FALSE) {

  sim.out <- sims(x, reps, tsm)
  fetch.deviance <- function(y) y$results$lnl
  data.deviance <- fetch.deviance(x)
  sim.ci<-c(sim.out$deviance.025,sim.out$deviance.975)
  cat("data deviance = ",data.deviance, ", simulation mean = ", sim.out$deviance.mean, ", simulation 95%CI = ", sim.ci,  "\n")

  #modified c.hat
  c.hat<-data.deviance/sim.out$deviance.mean
  cat("modified c.hat =",c.hat,  "\n")

}
