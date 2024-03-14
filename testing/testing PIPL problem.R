library(RMark)

sum(extract.n.marked(x))
test <- extract.model(x)

test$P[100,]

simul.boot(test)

see <- sims(x, 1)

bootstrap.deviance(x, 2, tsm = TRUE)
