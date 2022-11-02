#'@export
#'
#'@title Simulate CJS data for 1 group at a time
#'
#'@description This function blah blah blah...
#'
#'@param phi (required, numeric matrix) matrix of survival probabilities
#'    (as a number between 0 and 1) with one row per group, one column per
#'    number of occasions - 1.
#'
#'@param p (required, numeric matrix) matrix of recapture probabilities
#'    (as a number between 0 and 1) with one row per group, one column per
#'    number of occasions - 1.
#'
#'@param marked (required, integer, length >= 1) matrix of number of animals
#'    marked with one row per group, and one column per number of occasions - 1.
#'
#'@param tsm (optional). Boolean (default \code{FALSE}) indicating whether a
#'    time since marking effect should be included.
#'
#'@details
#'  Any pertinent details....
#'
#'@return
#'   Returns a 2D matrix of capture histories with one row per individual and 1
#'   column per occasion.
#'
#'@author
#'  Who Do?
#'
simul.cjs<-function(phi,p,marked, tsm = FALSE)
{
  n.occasions<-length(p)+1

  if(tsm == 0){
    Phi<-matrix(phi,n.occasions-1,nrow=sum(marked),byrow=T)
    P<-matrix(p,n.occasions-1,nrow=sum(marked),byrow=T)
  }

  if(tsm == 1){
    Phi <- matrix(0, ncol = n.occasions-1,nrow=sum(marked))
    for (i in 1:length(marked)){
      Phi[(sum(marked[1:i])-marked[i]+1):sum(marked[1:i]), i:(n.occasions-1)] <-
        matrix(rep(phi[c(i+n.occasions-2,seq(i, (n.occasions-2),
            length.out = n.occasions - i -1))], marked[i]),
            ncol = n.occasions-i, byrow = TRUE)
    }

    P<-matrix(p,n.occasions-1,nrow=sum(marked),byrow=T)
  }

  #n.occasions<-dim(Phi)[2]+1
  CH<-matrix(0,ncol=n.occasions,nrow=sum(marked))
  #define a vector with marking occasion
  mark.occ<-rep(1:length(marked),marked[1:length(marked)])
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
  return(CH)
}

#'@export
#'
#'@title Create capture history character strings
#'
#'@description Converts the input capture history matrix into a vector of
#'    capture history strings.
#'
#'@param x (required, numeric) matrix of capture histories (0/1) with one row
#'    per group and one column per number of occasions - 1.
#'
#'@return
#'   A vector of capture history strings with length equal to \code{nrow(x)}.
#'
#'@author
#'  Who Do?
#'
pasty<-function(x)
{
  k<-ncol(x) # XXX Not used?
  n<-nrow(x)
  out<-array(dim=n)
  for (i in 1:n)
  {
    out[i]<-paste(x[i,],collapse="")
  }
  return(out)
}


# XXX from "set up bootstrap.R"
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


# xxx from "CMR functions.R"
# # XXX need to remove reference to global variable "Groups" and
# replace with reference to param "groups"
Marked<-function(data=dipper,n.occasions=7,groups=Groups)
{
  group<-data[,2]
  marked<-matrix(nrow=length(Groups),ncol=n.occasions)
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

#'@export
#'
#'@title Create trap dependence
#'
#' @description Creates trap dependence covariates for a capture history
#'
#' @param ch (required, numeric vector) capture history vector (0/1 values only)
#'
#' @param varname (required, character string, default \code{td}) prefix for
#'     variable name
#'
#' @param begin.time (required, integer, default 1) time for first occasion
#' @return Returns a dataframe with trap-dependent variables named
#'     \code{varnamet+1}, ... , \code{varnamet+nocc-1} where \code{t} is
#'     begin.time and \code{nocc} is the number of occasions.
#'
#' @details Assumes that the time between occasions is 1.
#'
#' @author From RMARK chapter of gentle introduction:
#'  http://www.phidot.org/software/mark/docs/book/pdf/app_3.pdf
#'
create.td = function(ch,
                     varname = "td",
                     begin.time = 1)
{
  # turn vector of capture history strings into a vector of characters
  char.vec = unlist(strsplit(ch, ""))

  # test to make sure they only contain 0 or 1
  if (!all(char.vec %in% c(0, 1)))
    stop("Function only valid for CJS model without missing values")
  else
  {
    # get number of occasions (nocc) and change it into a matrix of numbers
    nocc = nchar(ch[1])
    tdmat = matrix(as.numeric(char.vec), ncol = nocc, byrow = TRUE)
    # remove the last column which is not used
    tdmat = tdmat[, 1:(nocc - 1)]
    # turn it into a dataframe and assign the field (column) names
    tdmat = as.data.frame(tdmat)
    names(tdmat) = paste(varname, (begin.time + 1):(begin.time + nocc -
                                                      1), sep = "")
    return(tdmat)
  }
}

#'@export
#'
#'@title Select model results
#'
#' @description Select a particular model of interest from a list of models.
#'
#' @param model_num (required, integer) the model number of interest
#'
#' @details Converts a numeric \code{model_num} offset  to a symbolic model
#'     name and extracts the named model from global variable \code{all.models} using
#'     global variable \code{model.list}. This requires the following syntax in
#'     the calling script:
#'
#'     \code{all.models = mark.wrapper(model.list, data = data.processed,
#'        ddl = data.ddl, threads = 2)}
#'
#'     It would be a good idea to replace these global vars with passed params
#'     to make this more generic.
#'
#' @return
#'    The chosen model.
#'
#' @author Who do?
#'
select_model_results <- function(model_num) {
  all.models[[paste0(model.list$Phi[model_num], '.', model.list$p[model_num])]]
}


