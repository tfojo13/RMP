
#'@title Set up a Predictor of Binary Events Based on Repeated Measurements
#'
#'@param rmp The results of a call to make.measurement.predictor
#'@param dataset A data frame with rows $time, $event (1 or 0), $id
#'@param event.name Specifying a name allows for more than one event to be predicted
#'
#'@export
#' @importFrom stats as.formula binomial cov dnorm glm
#' @importFrom stats optimize pnorm qnorm rnorm runif sd
create.event.predictor <- function(rmp,
                                   dataset,
                                   covariates,
                                   event.name='event',
                                   verbose=FALSE)
{
    num.iterations = get.num.iterations(rmp)
    N = dim(dataset)[1]

    if (is.null(rmp$event.predictor.coefficients))
        rmp$event.predictor.coefficients = list()

    if (is.null(rmp$id.to.subj))
        rmp$id.to.subj = sapply(rmp$subj.to.id, index.of, rmp$subj.to.id)
    rmp$class.coefficient.names = paste0('class', 1:rmp$K)
    rmp$level.coefficient.names = paste0('level_', rmp$dimensions)
    rmp$slope.coefficient.names = paste0('slope_', rmp$dimensions)
    covariate.names = c(rmp$class.coefficient.names, rmp$level.coefficient.names, rmp$slope.coefficient.names, rmp$fixed.covariate.names)

    rmp$event.predictor.coefficients[[event.name]] = t(sapply(1:num.iterations, function(iteration){
        if (verbose)
            print(paste0('Fitting logistic model for iteration ', iteration, ' (of ', num.iterations, ')'))
        reg.covariates = t(sapply(1:N, function(j){
            id = dataset[j,'id']
            subj = rmp$id.to.subj[as.character(id)]

            time = dataset[j, 'time']
            M = create.time.design.matrix(rmp, time)
            M.slope = create.slope.design.matrix(rmp, time)

            covariates.for.j = covariates[covariates$id==id,rmp$fixed.covariate.names]

            if (rmp$K == 1)
            {
                intercepts = 1
                eta = 1
            }
            else
            {
                eta = rmp$BUGSoutput$sims.list$eta[iteration,subj]
                intercepts = sapply(1:rmp$K, function(k){if (eta==k) 1 else 0})
            }

            betas = matrix(rmp$BUGSoutput$sims.list$individual.beta[iteration, subj,] +
                               get.fixed.beta.sums(rmp, iteration, eta, c(1,covariates.for.j)),
                           ncol=1)

            covs = as.numeric(c(intercepts, M %*% betas, M.slope %*% betas, covariates.for.j))
            names(covs) = covariate.names
            covs
        }))
        dataset.for.iter = cbind(dataset, reg.covariates)
        formula = as.formula(paste0('event~0+', paste0(covariate.names, collapse=' + ')))
        fit=glm(formula, family=binomial(link='logit'), data=dataset.for.iter)
        fit$coefficients
    }))

    rmp
}

#'@title Make Predictions About Binary Events Based on Repeated Measurements
#'
#'@param observations an array indexed [i,time,test] - where i denotes an individual
#'@param observed.times a numeric vector such that observed.times[t] is the time at which observations[,t,] were taken
#'
#'@export
make.event.predictions <- function(rmp,
                                   observations,
                                   observed.times,
                                   predict.times,
                                   covariates,
                                   event.name='event',
                                   verbose=F)
{
    if (verbose)
        print('Calculating posterior distributions on betas...')
    posterior.beta.distributions = get.posterior.beta.distributions(rmp, observations, observed.times, covariates)
    class.probabilities = get.all.class.probabilities(rmp, observations, observed.times, covariates)

    n.subj = dim(observations)[1]
    num.iterations = get.num.iterations(rmp)

    logit.means=sapply(1:length(predict.times), function(t){
        M = create.time.design.matrix(rmp, predict.times[t])
        M.slope = create.slope.design.matrix(rmp, predict.times[t])

        t(sapply(1:n.subj, function(j){
            if (verbose)
                print(paste0('making predictions for individual ', j))

            vals.per.iter = sapply(1:num.iterations, function(iteration){
                sum(sapply(1:rmp$K, function(k){
                    coefficients = rmp$event.predictor.coefficients[[event.name]][iteration,]
                    coefficients.for.levels = matrix(coefficients[rmp$slope.coefficient.names], nrow=1)
                    coefficients.for.slopes = matrix(coefficients[rmp$level.coefficient.names], nrow=1)
                    class.coefficient = coefficients[rmp$class.coefficient.names][k]
                    coefficients.for.covariates = coefficients[length(coefficients)+1-(rmp$n.fixed):1]

                    beta.dist = posterior.beta.distributions[[j]][[iteration]][[k]]
                    beta.dist$mu = beta.dist$mu + get.fixed.beta.sums(rmp, iteration, k, covariates[j,])

                    norm.dist = multiply.MVG(beta.dist, coefficients.for.levels %*% M + coefficients.for.slopes %*% M.slope)
                    norm.dist$mu = norm.dist$mu + class.coefficient + sum(coefficients.for.covariates * covariates[j,])

                    mean.logitnorm(as.numeric(norm.dist$mu), as.numeric(norm.dist$sigma)) * class.probabilities[j,iteration,k]
                }))
            })

            mean(vals.per.iter, na.rm=T)
        }))
    })

#    rowMeans(logit.means, dims=2)
    matrix(logit.means, ncol=length(predict.times))
}

# library(logitnorm)
# 
#' @importFrom logitnorm momentsLogitnorm
mean.logitnorm <- function(mu, sigma.sq)
{
    #    log.odds = saved.norms * sqrt(sigma.sq) + mu
    #   mean(exp(log.odds) / (1 + exp(log.odds)))
    tryCatch({
        momentsLogitnorm(mu, sqrt(sigma.sq))[1]
    },error=function(e){
#        print('there was an error using momentsLogitnorm. Using simulation to get mean of logitNorm distribution.')
#        norms = rnorm(10000, mu, sigma.sq)
        SNORMS = rnorm(10000)
        norms = SNORMS * sqrt(sigma.sq) + mu
        expits = exp(norms) / (1+exp(norms))
        mean(expits)
    })
}
