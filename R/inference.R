##----------------------------------------------------------------------------------------##
##-- A NOTE ON INDEXING --
##------------------------
##
## The coefficient vectors (the beta's) are treated as column vector, indexed:
## intercept_dimension1, first_slope_dimension1, second_slope_addon_dimension1,
##  intercept_dimension2, first_slope_dimension2, second_slope_addon_dimension2,
##  etc
##
## When the set of observations is treated as single vector, it is a column vector indexed:
## obs_on_time1_dimension1, obs_on_time1_dimension2, ..., obs_on_time1_dimension_d,
##  obs_on_time1_dimension2, obs_on_time2_dimension2, ..., obs_on_time2_dimension_d,
##
##----------------------------------------------------------------------------------------##


##----------------------------------##
##-- THE MAIN PREDICTION FUNCTION --##
##----------------------------------##

#PARAMETERS:
#observations - a matrix, each row corresponding to an individual set of observations
#                         each column corresponsing to an observation on a dimension and time, ordered t1/d1, t1/d2, t1/d3, ... , t2/d1, t2/d2, ..., tn/d1, tn/d2, ...
#observed times - a vector corresponding to the times for the columns of observations
#predict times -  taken to be all times not observed
#
#returns a 3d array indexed [i, d/t, statistic]
#  where i denotes an individual
#  d/t denotes an observation for a dimension and time
#  statistic is one of: mean, ci.lower, ci.upper


#'@title Make predictions of Repeated Measurements
#'@param observations an array indexed [i,time,test] - where i denotes an individual
#'@param observed.times a numeric vector such that observed.times[t] is the time at which observations[,t,] were taken
#'@param covariates A row for each individual with the covariates used to fit the model
#'
#'@return A 4d array indexed[individual, time, test, statistic], where statistic is one of: mean, ci.lower, ci.upper, iq.lower, iq.upper
#'
#'@export
make.measurement.predictions <- function(dm,
                                         observations, observed.times,
                                         predict.times,
                                         covariates=NULL,
                                         ci.coverage=0.95,
                                         iq.coverage=0.5,
                                         sum=F)
{
    nT = length(observed.times)
    D = dm$D

    covariates = cbind(rep(1,dim(observations)[1]), covariates[,dm$fixed.covariate.names])

    #make sure we're a 3d array
    if (length(observed.times)==1)
        observations = array(observations, dim=c(prod(dim(observations))/nT/D,nT,D))

    #get alphas ready for intervals
    alpha.ci = 1-ci.coverage
    alpha.iq = 1-iq.coverage

    #set up observations into a matrix indexed [i, D*(t-1) + d
    n.predict = dim(observations)[1]
    obs = matrix(NaN, nrow=n.predict, ncol=length(observed.times)*D)
    for (i in 1:n.predict)
    {
        for (t in 1:length(observed.times))
        {
            obs[i, D*(t-1) + 1:D] = observations[i,t,]
        }
    }

    #set up times and masks
    all.times = c(observed.times, predict.times)
    observed.mask = c(rep(T, dm$D * length(observed.times)),
                      rep(F, dm$D * length(predict.times)))

    if (sum)
        num.predictions.per = length(predict.times)
    else
        num.predictions.per = length(predict.times) * dm$D

    if (sum)
        sum.matrix = create.sum.dimensions.matrix(dm, predict.times)

    num.iterations = get.num.iterations(dm)

    obs.dists = get.observation.distributions(dm, all.times)
    flattened.dists = flatten.list.first.level(obs.dists)

    M.all = create.time.design.matrix(dm, all.times)
    M.observed = create.time.design.matrix(dm, observed.times)

    predictions = sapply(1:n.predict, function(j)
    {
        print(paste0('Subj ', j))
        observed = obs[j,]

        keep.mask = !observed.mask
        keep.mask[observed.mask] = !is.na(observed)

        condition.on = rep(NaN, dm$D * length(predict.times) + length(observed))
        condition.on[observed.mask] = observed
        condition.on = condition.on[keep.mask]

        flattened.means = flatten.list.first.level(get.means.for.covariates(dm, M=M.all, covariates=covariates[j,]))

        dists = lapply(1:length(flattened.dists), function(index){
            dist = set.MVG.mean(flattened.dists[[index]], flattened.means[[index]])
            iter.dist = marginalize.MVG(dist, keep.mask)
            iter.dist = conditionalize.MVG(iter.dist, condition.on)
            if (sum)
                iter.dist = multiply.MVG(iter.dist, sum.matrix)
            iter.dist
        })
        weights = as.numeric(sapply(1:num.iterations, get.class.probabilities, dm, observed, observed.times, M.observed, covariates=covariates[j,], obs.dists))

        #some distributions have non-sensical sigmas, due to a very high individual.sigma
        throw.out = sapply(dists, function(dist){any(diag(dist$sigma)<0)})
        dists = dists[!throw.out]
        weights = weights[!throw.out]

        inner.predictions = array(0, dim=c(num.predictions.per, 5), dimnames=list(NULL, c('mean', 'ci.lower', 'iq.lower', 'iq.upper', 'ci.upper')))
        mixes = get.gaussian.mixtures.from.mvgs(dists, weights=weights)
        for (i in 1:num.predictions.per)
        {
            mix = mixes[[i]]
            inner.predictions[i,] = c(mean=mean.mix(mix), qnorm.mix(c(alpha.ci/2, alpha.iq/2, 1-alpha.iq/2, 1-alpha.ci/2), mix))
        }
        inner.predictions
    })

    predictions = t(predictions)
    dim(predictions) = c(n.predict, num.predictions.per, 5)

    preds.arr = array(0, dim=c(n.predict, length(predict.times), dm$D, 5),
                      dimnames=list(NULL, as.character(predict.times), as.character(dm$dimensions),  c('mean', 'ci.lower', 'iq.lower', 'iq.upper', 'ci.upper')))

    if (sum)
        D = 1

    for (i in 1:n.predict)
    {
        for (t in 1:length(predict.times))
        {
            preds.arr[i,t,,] = predictions[i, D*(t-1) + 1:D, ]
        }
    }

    preds.arr
}



##--------------------------------------------------------------------##
##-- FUNCTIONS FOR DETERMINING GENERAL DISTRIBUTION OF OBSERVATIONS --##
##--------------------------------------------------------------------##

##rv[[iteration]][[k]] is the mvg for class k at iteration
get.observation.distributions <- function(dm, times)
{
    num.iterations = get.num.iterations(dm)

    M = create.time.design.matrix(dm, times)
    iteration.mvgs = lapply(1:num.iterations, get.obs.dists.for.iteration, dm, times, M)

    iteration.mvgs
}

get.observation.distributions.and.weights <- function(dm, times)
{
    num.iterations = get.num.iterations(dm)
    dists = get.observation.distributions(dm, times)
    if (dm$K == 1)
        weights = rep(1, num.iterations)
    else
        weights = dm$BUGSoutput$sims.list$pi / num.iterations

    list(dists = flatten.list.first.level(dists), weights=as.numeric(weights))
}

#if aggregate is true, returns a single mvg object
#if false, returns a list of mvg objects
get.obs.dists.for.iteration <- function(iteration, dm, times, M)
{
    rv = list()
    for (k in 1:dm$K)
    {
        sigma.epsilon = get.sigma.epsilon(dm, iteration, times)

#        mu = get.means.for.iteration(dm, iteration, k, M)
        mu = rep(0, length(times)*dm$D)

        individual.sigma = get.individual.sigma(dm, iteration, k)

        rv[[k]] = create.MVG(mu = mu,
                             sigma = sigma.epsilon + M %*% individual.sigma %*% t(M))
    }

    rv
}


##---------------------------------------------------##
##-- DEALING WITH POSTERIOR DISTRIBUTIONS OF BETAS --##
##---------------------------------------------------##

## observations - an array indexed [i,time,test] - where i denotes an individual
##
## rv[[subj]][[iteration]][[k]] represents the posterior MVG for the individual betas for
##   individual subj based on iteration if they belong in class k
get.posterior.beta.distributions <- function(dm, observations, observed.times, covariates)
{
    subj = 1:dim(observations)[1]
    num.iterations = get.num.iterations(dm)

    covariates = cbind(rep(1,dim(observations)[1]), covariates[,dm$fixed.covariate.names])

    M = create.time.design.matrix(dm, observed.times)

    lapply(subj, function(j){
        lapply(1:num.iterations, function(iteration){
            lapply(1:dm$K, function(k){
#print(paste0('j=',j, ', iteration=',iteration))
                get.posterior.beta.dist.for.iteration.and.class(iteration, k, dm, observations[j,,], observed.times, M, covariates[j,])
            })
        })
    })
}

#observations here indexed [time,test]
get.posterior.beta.dist.for.iteration.and.class <- function(iteration, k, dm, observations, observed.times, M, covariates)
{
    if (all(is.na(observations)))
    {
        sigma = get.individual.sigma(dm, iteration, k)
        return (create.MVG(mu=rep(0, dim(sigma)[1]), sigma=sigma))
    }

    observations = as.numeric(t(observations))
    keep = !is.na(observations)
    M = M[keep,]

    omega.b = get.individual.omega(dm, iteration, k)
    omega.e = get.omega.epsilon(dm, iteration, observed.times)[keep,keep]

    y = matrix(observations[keep] - as.numeric(get.means.for.iteration.class.and.covariates(dm, iteration, k, M, covariates)), ncol=1)

    M.t.omega.e = t(M) %*% omega.e
    A = M.t.omega.e %*% M + omega.b
    b = M.t.omega.e %*% y
    A.inv = solve(A)

    create.MVG(A.inv %*% b, A.inv)
}

#'@title test function
#'@export
test.from.inside <- function()
{
    browser()
}


##------------------------------------------------##
##-- FUNCTIONS DEALING WITH CLASS PROBABILITIES --##
##------------------------------------------------##

#rv[subj, iteration, k] is the probability that individual subj is in class k at iteration
#'@title Get probabilites of being in each latent class
#'@export
get.all.class.probabilities <- function(dm, observations, observed.times, covariates)
{
    n.subj = dim(observations)[1]
    D=dm$D

    covariates = cbind(rep(1,dim(observations)[1]), covariates[,dm$fixed.covariate.names])

    obs = matrix(NaN, nrow=n.subj, ncol=length(observed.times)*D)
    for (i in 1:n.subj)
    {
        for (t in 1:length(observed.times))
        {
            obs[i, D*(t-1) + 1:D] = observations[i,t,]
        }
    }

    num.iterations = get.num.iterations(dm)
    obs.dists = get.observation.distributions(dm, observed.times)
    M = create.time.design.matrix(dm, observed.times)

    rv=sapply(1:n.subj, function(j){
        t(sapply(1:num.iterations, get.class.probabilities,
                 dm, obs[j,], observed.times, M, covariates=covariates[j,], all.obs.dists=obs.dists))
    })

    rv = t(rv)

    dim(rv) = c(n.subj, num.iterations, dm$K)

    rv
}

#if the observations are fewer than the num obs in the obs.dist,
#the given observations are presumed to be the first length(obsevations) observations on the distribution
get.class.probabilities <- function(iteration, dm, observations, times, M, covariates, all.obs.dists=NULL)
{
    if (dm$K==1)
        return (1)

    if (is.null(all.obs.dists))
        obs.dists = get.obs.dists.for.iteration(iteration, dm, times, M)
    else
        obs.dists = all.obs.dists[[iteration]]

    keep.indices = rep(F, obs.dists[[1]]$length)
    keep.indices[1:length(observations)] = !is.na(observations)

    p.class = dm$BUGSoutput$sims.list$pi[iteration,]
    if (all(is.na(observations)))
        return (p.class)

    observations = observations[!is.na(observations)]

    p.obs.given.class = sapply(1:dm$K, function(k){
        dist = set.MVG.mean(obs.dists[[k]], get.means.for.iteration.class.and.covariates(dm, iteration, k, M, covariates))
        dist = marginalize.MVG(dist, keep.indices)

        if (length(observations)==1)
            dnorm(observations, mean=as.numeric(dist$mu), sd=sqrt(as.numeric(dist$sigma)))
        else
            dmvnorm(observations, mean=dist$mu, sigma=dist$sigma)
    })

    p.obs.and.class = p.obs.given.class * p.class

    p.obs.and.class / sum(p.obs.and.class)
}

simulate.class <- function(class.probs)
{
    cum.probs = cumsum(class.probs)
    (1:length(class.probs))[runif(1,0,sum(class.probs))<=cum.probs][1]
}


##---------------------------------------------------##
##-- EXTRACT THE MODEL-ESTIMATED COVARIANCE MATRIX --##
##---------------------------------------------------##

#'@title Return the (average) model-estimated covariance matrix
#'
#'@param rmp An rmp object as returned by create.measurement.predictor
#'@param times A vector of times for which to generate a covariance matrix
#'
#'@return A T*D x T*D covariance matrix, where T=length(times) and D is the number of tests (dimensions).
#'The return value is indexed such that element [D*(t1-1)+d1, D*(t2-1)+d2] represents the covariance between a measurment on test d1 at time times[t1] with a measurement on test d2 at time times[t2]
#'
#'@export
get.overall.covariance.matrix <- function(rmp, times)
{
    dists.and.weights = get.observation.distributions.and.weights(rmp, times)

    combined.sigma = combine.MVGs(dists.and.weights$dists, dists.and.weights$weights)$sigma

    combined.sigma
}

##------------------------------------------------------------##
##-- HELPERS FOR SETTING UP DESIGN MATRICES (based on time) --##
##------------------------------------------------------------##

#rv is a matrix such that rv %*% col.vector.of.betas is the ordered
# set of observations for each time, dimension

# Maps a vector of betas [beta_d1_p1, beta_d1_p2,...,beta_dD_p1,...beta_dD_pP]
#  to a vector of observations [y_t1_d1, y_t1_d2, ..., y_t1_dD, ..., y_tT_d1, ..., y_tT_dD]
create.time.design.matrix <- function(dm, times)
{
    nT = length(times)
    D = dm$D
    P = dm$P

    design.matrix = matrix(0, nrow=nT*D, ncol=P*D)

    for (t in 1:nT)
    {
        for (d in 1:D)
        {
            design.matrix[(t-1)*D + d, P*(d-1) + 1:P] = dm$tt.func(times[t])
        }
    }

    design.matrix
}

create.slope.design.matrix <- function(dm, times)
{
    nT = length(times)
    D = dm$D
    P = dm$P

    design.matrix = matrix(0, nrow=nT*D, ncol=P*D)

    for (t in 1:nT)
    {
        for (d in 1:D)
        {
            design.matrix[(t-1)*D + d, P*(d-1) + 1:P] = dm$slope.t.func(times[t])
        }
    }

    design.matrix
}

create.sum.dimensions.matrix <- function(dm, times)
{
    nT = length(times)

    do.create.sum.dimensions.matrix(dm$D, nT)
}

do.create.sum.dimensions.matrix <- function(D, nT)
{
    rv = matrix(0, nrow=nT, ncol=nT*D)

    for (t in 1:nT)
    {
        rv[t,D*(t-1)+1:D] = rep(1,D)
    }

    rv
}

##----------------------------------------------------##
##-- GETTING THE MEANS OF DISTS BASED ON COVARIATES --##
##----------------------------------------------------##

get.means.for.covariates <- function(dm, M, covariates)
{
    num.iterations = get.num.iterations(dm)
    lapply(1:num.iterations, function(iter){
        lapply(1:dm$K, function(k){
            get.means.for.iteration.class.and.covariates(dm, iter, k, M, covariates)
        })
    })
}

get.means.for.iteration.class.and.covariates <- function(dm, iteration, k, M, covariates)
{
    #fixed.betas are indexed p,d,f,k

#    fixed.betas = matrix(dm$BUGSoutput$sims.list$fixed.beta[iteration,,,,k],
   #                         nrow=dm$P*dm$D, ncol=dm$n.fixed)
 #   fixed.beta.sums = colMeans(t(fixed.betas) * as.numeric(covariates))

  #  M %*% as.matrix(fixed.beta.sums, ncol=1)
    M %*% get.fixed.beta.sums(dm, iteration, k, covariates)
}

##------------------------------------------------------------##
##-- HELPERS FOR EXTRACTING PARAMETERS FROM THE MCMC OUTPUT --##
##------------------------------------------------------------##


get.fixed.beta.sums <- function(dm, iteration, k, covariates)
{
    if (is.null(dm$n.fixed))
    {
        fixed.beta.sums = dm$BUGSoutput$sims.list$fixed.beta[iteration,,k]
    }
    else
    {
        fixed.betas = matrix(dm$BUGSoutput$sims.list$fixed.beta[iteration,,,,k],
                             nrow=dm$P*dm$D, ncol=dm$n.fixed)
        fixed.beta.sums = colMeans(t(fixed.betas) * as.numeric(covariates))
    }

    matrix(fixed.beta.sums, ncol=1)
}

get.individual.sigma <- function(dm, iteration, k)
{
    dm$BUGSoutput$sims.list$individual.sigma[iteration,,,k]
}

get.individual.omega <- function(dm, iteration, k)
{
    solve(get.individual.sigma(dm, iteration, k))
}

get.individual.beta <- function(dm, iteration, subj)
{
    if (dm$dimensions.covary)
        dm$BUGSoutput$sims.list$individual.beta[iteration,subj,]
    else
    {
        subscale.dimension = dim(dm$BUGSoutput$sims.list$individual.beta.subscales)[3]
        rv = matrix(0, ncol=1, nrow=subscale.dimension*dm$num.dimensions)
        for (d in 1:dm$num.dimensions)
            rv[3*d-2:0] = dm$BUGSoutput$sims.list$individual.beta.subscales[iteration,subj,,d]
        rv
    }
}

get.num.iterations <- function(dm)
{
    dm$BUGSoutput$n.keep * dm$BUGSoutput$n.chains
}

get.sigma.epsilon <- function(dm, iteration, times)
{
    diag(rep(dm$BUGSoutput$sims.list$sigma.epsilon[iteration,]^2, length(times)))
}

get.omega.epsilon <- function(dm, iteration, times)
{
    diag(diag(get.sigma.epsilon(dm, iteration, times))^-1)
}

##-----------------------##
##-- LOW-LEVEL HELPERS --##
##-----------------------##

flatten.list.first.level <- function(l)
{
    rv = list()
    for (elem in l)
        rv = c(rv, elem)

    rv
}
