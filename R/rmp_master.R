
##--
##-- CONSTRUCTORS --##
##

#'@title Create a Measurement Predictor
#'
#'@param dataset A data frame with the columns 'measurement', 'test', 'id', 'time', as well as column for each entry in 'fixed.covariate.names'
#'@param time.transformation One of 'linear', 'linear.spline'
#'
#'@export
#' @importFrom R2jags jags
create.measurement.predictor <- function(dataset,
                                         num.latent.classes=1,
                                         time.transformation,
                                         test.minima,
                                         test.maxima,
                                         time.transformation.parameters=NULL,
                                         fixed.covariate.names=NULL,
                                         initial.class.assignments=NULL,
                                         iterations.after.burn,
                                         burn,
                                         keep,
                                         chains=2,
                                         pi.alpha=1,
                                         seed=1223)
{

##-- SET UP INPUTS TO JAGS --##
    subj.to.id = as.character(unique(dataset$id))
    n = length(subj.to.id)
    id.to.subj = 1:n
    names(id.to.subj) = subj.to.id
    subj = id.to.subj[as.character(dataset$id)]

    N = dim(dataset)[1]
    y = dataset$measurement

    d.to.test = unique(dataset$test)
    D = length(d.to.test)
    test.to.d = 1:D
    names(test.to.d) = d.to.test
    dimension = test.to.d[dataset$test]

    K = num.latent.classes

    if (is.null(fixed.covariate.names))
        n.fixed = 1
    else
        n.fixed = 1+length(fixed.covariate.names)

    fixed.covariates = cbind(rep(1,N),dataset[,fixed.covariate.names])

    test.minima = test.minima[d.to.test]
    test.maxima = test.maxima[d.to.test]
    if (any(is.na(test.minima)) || any(is.na(test.maxima)) || any(is.infinite(test.minima)) || any(is.infinite(test.maxima)))
        stop("A finite minimum and maximum must be specified for each test in the 'test' column of the dataset argument")

    if (time.transformation=='linear')
    {
        P = 2
        tt.func = function(t){c(1,t)}
        slope.t.func = function(t){c(0,1)}

        FIXED.BETA.MEANS = t(sapply(1:D, function(d){
            c((test.maxima[d]-test.minima[d])/2, 0)
        }))

        FIXED.BETA.SDS = t(sapply(1:D, function(d){
            rep((test.maxima[d]-test.minima[d])/2, 2)
        }))
    }
    else if (time.transformation=='linear.spline')
    {
        if (is.null(time.transformation.parameters) || length(time.transformation.parameters) < 1 ||
            (class(time.transformation.parameters)!='numeric' && class(time.transformation.parameters) != 'integer'))
            stop("For a time transformation of 'linear spline', the time.transformation.parameters must be passed an integer or numeric with the spline time")

        P = 2+length(time.transformation.parameters)
        tt.func = function(t){c(1,t,pmax(rep(0,length(time.transformation.parameters)),t-time.transformation.parameters))}
        slope.t.func = function(t){c(0,1,as.numeric(t>=time.transformation.parameters))}

        FIXED.BETA.MEANS = t(sapply(1:D, function(d){
            c((test.maxima[d]-test.minima[d])/2, rep(0,P-1))
        }))
        FIXED.BETA.SDS = t(sapply(1:D, function(d){
            rep((test.maxima[d]-test.minima[d])/2, P)
        }))
    }
    else if (time.transformation=='cubic.spline')
    {
        P = 4
        tt.func = function(t){c(1,t,t^2,t^3)}
        slope.t.func = function(t){c(0,1,2*t,3*t^2)}

        FIXED.BETA.MEANS = t(sapply(1:D, function(d){
            c((test.maxima[d]-test.minima[d])/2, rep(0,P-1))
        }))
        FIXED.BETA.SDS = t(sapply(1:D, function(d){
            rep((test.maxima[d]-test.minima[d])/2, P)
        }))
    }
    else
        stop("At this point, only 'linear' and 'linear.spline' are supported as time transformations")

    M = t(sapply(dataset$time, tt.func))

    IDENTITY = diag(rep(1,D*P))
    ZERO.VECTOR = rep(0,D*P)
    PI.ALPHA = if (length(pi.alpha)<K) rep(pi.alpha,K) else pi.alpha

    jags.data = list('N', 'n', 'D', 'P',
                     'M', 'y', 'subj', 'dimension',
                     'n.fixed', 'fixed.covariates',
                     'FIXED.BETA.MEANS', 'FIXED.BETA.SDS', 'IDENTITY', 'ZERO.VECTOR')

    if (K > 1)
        jags.data = c(jags.data, 'K', 'PI.ALPHA')

##-----------##
##-- INITS --##
##-----------##

    if (K>1 && !is.null(initial.class.assignments))
    {
        generators = lapply(1:K, function(k){
            dataset.for.k = dataset[sapply(dataset$id, function(id){initial.class.assignments[as.character(id)]==k}),]
            random.fixed.beta.and.sd.generator(dataset.for.k, d.to.test, P, tt.func)
        })
    }

    inits = function()
    {
        rv = list()

        if (K>1 && !is.null(initial.class.assignments))
        {
#            rv$fixed.beta = sapply(1:K, function(k){
 #               generators[[k]]$r.fixed.betas()
  #          })

            rv$eta = initial.class.assignments[subj.to.id]

            rv$individual.omega = array(0, dim=c(D*P,D*P,K))
            for (k in 1:K)
                rv$individual.omega[,,k] = solve(generators[[k]]$r.sigma())
        }

        rv
    }

##------------------------##
##-- PARAMETERS TO SAVE --##
##------------------------##

    parameters.to.save = c('fixed.beta', 'individual.sigma', 'individual.beta', 'sigma.epsilon')
    if (K > 1)
        parameters.to.save = c(parameters.to.save, 'pi', 'eta')

##-----------------##
##-- WRITE MODEL --##
##-----------------##

    model.file = tempfile()
    if (K==1)
        write.one.class.model(model.file)
    else
        write.k.class.model(model.file)

##--------------##
##-- RUN JAGS --##
##--------------##

    if (!is.na(seed))
        set.seed(seed)

    print(paste0("Running JAGS..."))

    jagsfit <- jags(data = jags.data,
                    inits = inits,
                    parameters.to.save = parameters.to.save,
                    n.iter = iterations.after.burn + burn,
                    n.burnin = burn,
                    n.thin = chains * iterations.after.burn / keep,
                    n.chains= chains,
                    model.file = model.file,
                    jags.seed=seed)
    print("...DONE running JAGS!")

    rv = jagsfit

##-------------------##
##-- PACKAGE IT UP --##
##-------------------##

    rv$D = D
    rv$P = P
    rv$tt.func = tt.func
    rv$slope.t.func = slope.t.func

    rv$time.transformation = time.transformation
    rv$dataset = dataset

    rv$dimensions = as.character(d.to.test)
    rv$subj.to.id = subj.to.id
    rv$id.to.subj = id.to.subj

    rv$fixed.covariate.names = fixed.covariate.names
    rv$n.fixed = n.fixed

    rv$K = K

    rv
}




##-- PREDICTORS --##

predict.symptoms <- function(rmp,
                             measurements,
                             predict.times)
{

}

predict.event <- function(rmp, event,
                          measurements)
{

}


#HELPERS

#rv is a list with two elements:
# $r.fixed.betas(n) - a function that generates n random sets of beta1, beta2, beta3 (intercept, slope, change from slope)
#               returns an array of dimensions [3*D, 1]
# $r.sd() - a function that generates one random set of sds - a vector of length T (where the last T-1 elements are all the same)
#

#'@title exported for testing
#'@export
random.fixed.beta.and.sd.generator <- function(dataset, d.to.test, P, tt.func)
{
    obs = list() #indexed[[id]][[test]] - with elements $vals $times
    D = length(d.to.test)

    for (j in 1:dim(dataset)[1])
    {
        id = as.character(dataset[j,'id'])
        if (is.null(obs[[id]]))
        {
            obs[[id]] = list()
            for (d in 1:D)
            {
                obs[[id]][[d.to.test[d]]] = list(vals=numeric(), times=numeric())
            }
        }

        test = dataset[j,'test']

        obs[[id]][[test]]$vals = c(obs[[id]][[test]]$vals, dataset[j, 'measurement'])
        obs[[id]][[test]]$times = c(obs[[id]][[test]]$times, dataset[j, 'time'])
    }

    n = length(obs)

    betas = sapply(1:D, function(d){
        sapply(1:n, function(i){
            X=matrix(t(sapply(obs[[i]][[d]]$times, tt.func)), ncol=P)
            Y=matrix(obs[[i]][[d]]$vals, ncol=1)
            tryCatch(as.vector(solve(t(X)%*%X)%*%t(X)%*%Y), error=function(e){rep(NA,P)})
        })
    })

    dim(betas) = c(P,n,D)
    #betas is now indexed [p,i,d]

    long.betas = matrix(NaN, nrow=n, ncol=D*P)
    for (i in 1:n)
    {
        for (d in 1:D)
        {
            long.betas[i,P*(d-1) + 1:P] = betas[,i,d]
        }
    }

    rv = list()
    rv$r.fixed.betas <- function()
    {
        boots = sample(1:n, n, replace=T)
        rv = matrix(0, ncol=1, nrow=P*D)

        for (d in 1:D)
        {
            rv[P*(d-1) + 1:P, 1] = rowMeans(betas[,boots,d], na.rm=T)
        }

        rv
    }

    rv$r.sds <- function()
    {
        boots = sample(1:n, n, replace=T)

        sapply(1:D, function(d){
            sd(unlist(sapply(boots, function(i){
                obs[[i]][[d]]$vals
            })), na.rm=T)
        })
    }

    rv$r.sigma <- function()
    {
        boots = sample(1:n, n, replace=T)
        cov(long.betas[boots,], use='pairwise.complete.obs')
    }

    rv
}
