library(RMark)
library(magrittr)



# Sat Mar 19 09:49:53 2022 ------------------------------

# attempt to set up a general bootstrapping procedure based solely on an existing model

### first - step up a dummy ch-array with known rates

### use Kery and Shaub approach

n.occs <- 10
n.groups <- 3

marked <- matrix(c(rep(50, n.occs - 1),
                   rep(40, n.occs - 1),
                   rep(45, n.occs - 1)),
                 nrow = n.groups, byrow = TRUE)
phi <- matrix(c(rep(0.9, (n.occs - 1) ^ 2),
                rep(0.9,  (n.occs - 1) ^ 2),
                rep(0.2,  (n.occs - 1) ^ 2)),
              nrow = n.groups * (n.occs -1), ncol = n.occs -1, byrow = TRUE)
p <- matrix(c(rep(0.9, (n.occs - 1) ^ 2),
              rep(0.2,  (n.occs - 1) ^ 2),
              rep(0.2,  (n.occs - 1) ^ 2)),
            nrow = n.groups * (n.occs -1), ncol = n.occs -1, byrow = TRUE)

for(i in 1:n.groups){
    diag(phi[(1:(n.occs-1)) + ((i-1) * (n.occs - 1)),1:(n.occs-1)]) <- diag(phi[(1:(n.occs-1)) + ((i-1) * (n.occs - 1)),1:(n.occs-1)]) / 2
}

phi

Phi <- phi[rep(1:nrow(phi), times = as.vector(t(marked))), ]

P <- p[rep(1:nrow(p), times = as.vector(t(marked))), ]


# something with no releases in one year
# marked[,5] <- 0

# single group
# marked <- matrix(rep(50, n.occs - 1), nrow = n.groups, byrow = TRUE)
# something with no releases in one year
# marked <- matrix(c(rep(50, 4), 0, rep(50,4)), nrow = n.groups, byrow = TRUE)
# phi <- matrix(rep(0.8, n.occs - 1), nrow = n.groups, byrow = TRUE)
# p <- matrix(rep(0.8, n.occs - 1), nrow = n.groups, byrow = TRUE)

## this chunk for a model with group
mydata <- data.frame(ch = pasty(simul.cjs(phi,p,marked)), group = rep(1:n.groups, times = rowSums(marked)))
mydata <- data.frame(ch = pasty(CH), group = rep(1:n.groups, times = rowSums(marked)))

data.proc <- process.data(mydata,model="CJS", groups = "group")
data.ddl <- make.design.data(data.proc)
mymodel <- mark(data.proc, data.ddl,
                model.parameters = list(Phi = list(formula = ~time*group), p = list(formula = ~time*group)))

# add this in for tsm
data.ddl <- add.design.data(data.proc, data.ddl,
                            parameter="Phi", type="age", bins=c(0,1, 17),name="tsm",
                            right = FALSE, replace = TRUE)
mymodel <- mark(data.proc, data.ddl,
                model.parameters = list(Phi = list(formula = ~tsm+time*group), p = list(formula = ~time*group)))


## this chunk for a model with just time
mydata <- data.frame(ch = pasty(simul.cjs(phi,p,marked)))
mydata <- data.frame(ch = pasty(CH))
data.proc <- process.data(mydata,model="CJS")
data.ddl <- make.design.data(data.proc)
mymodel <- mark(data.proc, data.ddl,
                model.parameters = list(Phi = list(formula = ~time), p = list(formula = ~time)))

## and this chunk for fixed values with just time
mydata <- data.frame(ch = pasty(simul.cjs(phi,p,marked)))
data.proc <- process.data(mydata,model="CJS")
data.ddl <- make.design.data(data.proc)
mymodel <- mark(data.proc, data.ddl,
                model.parameters = list(Phi = list(formula = ~time, fixed=list(time=1, value = 0.8)), p = list(formula = ~time, fixed=list(time=3, value = 0.8))))

# alternate method of fixing values in the ddl directly


mydata <- data.frame(ch = pasty(simul.cjs(phi,p,marked)), group = rep(1:n.groups, times = rowSums(marked)))
data.proc <- process.data(mydata,model="CJS", groups = 'group')
data.ddl <- make.design.data(data.proc)
data.ddl$p$fix <-  NA
data.ddl$p$fix[data.ddl$p$time == 3 & data.ddl$p$group == 2] <-0.8

mymodel <- mark(data.proc, data.ddl,
                model.parameters = list(Phi = list(formula = ~time*group), p = list(formula = ~time*group)))


str(mymodel)


## try it with the dipper data set

data("dipper")

maledipper <- data.frame(ch = dipper[dipper$sex == "Male", 1])

data.proc <- process.data(maledipper,model="CJS")
data.ddl <- make.design.data(data.proc)
mymodel <- mark(data.proc, data.ddl,
                model.parameters = list(Phi = list(formula = ~time), p = list(formula = ~time)))


CJS.GOF.testing(data.proc)

str(mymodel)

bootstrap.deviance(mymodel, 10, tsm = TRUE)
