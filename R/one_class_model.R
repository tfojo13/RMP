# N = number of observations; indexed j
# n = number of independent subjects; indexed i
# D = number of dimensions; indexed d
# P = number of points in transformation of time; indexed p
#
# M[j,1:P] = the time transformation of the time at which measurement y[j] was taken

# y[j] = the jth observation
# subj[j] = the individual on whom observation y[j] is made
# dimension[j] = the dimension (subscale) about which observation y[j] is made

# FIXED.BETA.MEANS[d,p]
# FIXED.BETA.SDS[d,p]
# IDENTITY[1:(P*D), 1:(P*D)]
# ZERO.VECTOR[1:(P*D)]

write.one.class.model <- function(model.file) {
cat("model{

##---------------------##
##-- THE FIXED BETAS --##
##---------------------##

for (d in 1:D)
{
    for (p in 1:P)
    {
        for (f in 1:n.fixed)
        {
#            fixed.beta[P*(d-1)+p,f,1] ~ dnorm(0,100)
            fixed.beta[p,d,f,1] ~ dnorm(0,100)
        }
    }
}


##----------------------------##
##-- THE INDIVIDUAL EFFECTS --##
##----------------------------##


individual.omega[1:(P*D), 1:(P*D),1] ~ dwish(IDENTITY[1:(P*D), 1:(P*D)], (P*D))
individual.sigma[1:(P*D), 1:(P*D),1] <- inverse(individual.omega[1:(P*D), 1:(P*D),1])

#pick a set of betas for everybody
for (i in 1:n)
{
    individual.beta[i,1:(P*D)] ~ dmnorm(ZERO.VECTOR[1:(P*D)], individual.omega[1:(P*D), 1:(P*D),1])
}


##---------------------##
##-- THE ERROR TERMS --##
##---------------------##

for (d in 1:D)
{
    tau.epsilon[d] ~ dgamma(0.001,0.001)
    sigma.epsilon[d] <- tau.epsilon[d]^-0.5
}

##------------------##
##-- THE OUTCOMES --##
##------------------##

for (j in 1:N)
{
    for (p in 1:P)
    {
        fixed.beta.sum[j,p] <- sum(fixed.beta[p,dimension[j],1:n.fixed,1] * fixed.covariates[j,1:n.fixed])
    }

    s[j] <- sum(M[j,1:P] * (fixed.beta.sum[j,1:P] + individual.beta[subj[j],P*(dimension[j]-1)+1:P]))

    y[j] ~ dnorm(s[j], tau.epsilon[dimension[j]])
}


}", fill=TRUE, file=model.file)
}
