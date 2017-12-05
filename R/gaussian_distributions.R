##---------------------------------------------------##
##-- HELPER CODE FOR FUNCTIONS INVOLVING GAUSSIANS --##
##
##   Major components are (1) functions for multivariate
##    gaussian distributions and (2) functions for
##    mixtures of (univariate) gaussians
##
##---------------------------------------------------##

##------------------------------------------------------------------------##
##-- GENERAL FUNCTIONS OPERATING ON MULTIVARIATE GAUSSIAN DISTRIBUTIONS --##
##------------------------------------------------------------------------##

create.MVG <- function(mu, sigma)
{
    if (is.null(dim(mu)))
        mu = matrix(mu, ncol=1)

    list(mu=mu, sigma=sigma, length=nrow(mu))
}

set.MVG.mean <- function(mvg, mu)
{
    if (is.null(dim(mu)))
        mu = matrix(mu, ncol=1)

    mvg$mu = mu
    mvg
}

marginalize.MVG <- function(mvg, indices, keep.selected.indices=T)
{
    if (keep.selected.indices)
        create.MVG(mvg$mu[indices], mvg$sigma[indices,indices])
    else
        create.MVG(mvg$mu[-indices], mvg$sigma[-indices,-indices])
}

#given should be a vector of the same length as the MVG, with
# numeric values for the variables on which to condition,
# and NA for the other variables
#if keep.givens.in.dist is true, they are kept in the distribution (with zero variance),
# otherwise, the returned distribution involves only the remaining unknowns
conditionalize.MVG <- function(mvg, given, keep.givens.in.dist=F)
{
    if (all(is.na(given)))
        return (mvg)

    unknown.indices = sapply(given, is.na)
    known.indices = !unknown.indices

    sigma11 = mvg$sigma[known.indices, known.indices]
    sigma11.inv = solve(sigma11)

    sigma22 = mvg$sigma[unknown.indices, unknown.indices]
    sigma12 = mvg$sigma[known.indices, unknown.indices]
    sigma21 = mvg$sigma[unknown.indices, known.indices]

    mu1 = mvg$mu[known.indices]
    mu2 = mvg$mu[unknown.indices]
    y1 = matrix(given[known.indices], ncol=1)

    sigma21.times.11inv = sigma21 %*% sigma11.inv

    new.mu = mu2 + sigma21.times.11inv %*% (y1-mu1)
    new.sigma = sigma22 - sigma21.times.11inv %*% sigma12

    if (keep.givens.in.dist)
    {
        known = !is.na(given)

        new.mu.with.given = numeric(mvg$length)
        new.mu.with.given[known] = given[known]
        new.mu.with.given[!known] = new.mu

        new.sigma.with.given = matrix(0, nrow=mvg$length, ncol=mvg$length)
        new.sigma.with.given[!known, !known] = new.sigma

        create.MVG(mu=new.mu.with.given, sigma=new.sigma.with.given)
    }
    else
    {
        create.MVG(mu=new.mu, sigma=new.sigma)
    }
}

#B represents the constant matrix by which the distribution is multiplied
multiply.MVG <- function(mvg, B)
{
    create.MVG(mu = B %*% mvg$mu,
               sigma = B %*% mvg$sigma %*% t(B))
}

combine.MVGs <- function(mvgs, probabilities=rep(1/length(mvgs),length(mvgs)))
{
    probabilities = probabilities / sum(probabilities)

    combined.mu = matrix(0,nrow=mvgs[[1]]$length)
    combined.expected.sq = expected.var = matrix(0, ncol=ncol(mvgs[[1]]$sigma), nrow=nrow(mvgs[[1]]$sigma))

    for (i in 1:length(mvgs))
    {
        mu = mvgs[[i]]$mu
        combined.mu =  combined.mu + mu * probabilities[i]
        combined.expected.sq = combined.expected.sq + mu %*% t(mu) * probabilities[i]
        expected.var = expected.var + mvgs[[i]]$sigma * probabilities[i]
    }

    create.MVG(mu=combined.mu,
               sigma=expected.var + combined.expected.sq - combined.mu %*% t(combined.mu))
}

r.mvg <- function(n, mvg)
{
    rmvnorm(n, mean=as.vector(mvg$mu), sigma=mvg$sigma, method='chol')
}


##----------------------------##
##-- MIXTURES  OF GAUSSIANS --##
##----------------------------##

#returns a list
#rv[[i]] corresponds to a mixture of (univariate) gaussians representing the ith dimension of the given mvgs
get.gaussian.mixtures.from.mvgs <- function(mvgs, weights=rep(1/length(mvgs), length(mvgs)))
{
    lapply(1:mvgs[[1]]$length, function(i){
        create.gaussian.mixture(means=sapply(mvgs, function(mvg){mvg$mu[i,]}),
                                sds=sapply(mvgs, function(mvg){sqrt(mvg$sigma[i,i])}),
                                weights=weights)
    })
}

create.gaussian.mixture <- function(means, sds=rep(1,length(means)), weights=rep(1/length(means), length(means)))
{
    rv = list(means=means, sds=sds, weights=weights/sum(weights), n=length(means))

    rv
}

dnorm.mix <- function(x, mix)
{
    sum(dnorm(x, mix$means, mix$sds) * mix$weights)
}

pnorm.mix <- function(q, mix, lower.tail = TRUE)
{
    sum(pnorm(q, mix$means, mix$sds, lower.tail=lower.tail) * mix$weights)
}

qnorm.mix <- function(p, mix, lower.tail = TRUE, log.p = F)
{
    rv = numeric()

    for (one.p in p)
    {
        fn = function(x){
            (pnorm.mix(x, mix, lower.tail=lower.tail) - one.p)^2
        }

        individual.qs = qnorm(one.p, mix$means, mix$sds)
        min.q = min(individual.qs)
        max.q = max(individual.qs)
        mean.q = mean(individual.qs)
if (is.na(min.q))
    browser()
        one.q = optimize(fn, c(min.q, max.q), tol=0.0001)$minimum

        rv = c(rv, one.q)
    }

    rv
}

rnorm.mix <- function(n, mix)
{
    picks = resample.likelihoods(mix$weights, n)
    sapply(picks, function(pick){
        rnorm(1, mix$means[pick], mix$sds[pick])
    })
}

mean.mix <- function(mix)
{
    sum(mix$means*mix$weights)
}
