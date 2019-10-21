#' @useDynLib MCox
"_PACKAGE"


#' Projection onto a L1 Ball
#'
#' @param data 
#' @param task_index 
#' @param time_index 
#' @param censored_index 
#' @param weight_index 
#' @param features_indices 
#' @param ... 
#'
#' @return list
#' @export
#'
#' @examples
#' \dontrun{
#'    ProjectL1(1:3, 1)
#' }
MCox = function(data, task_index, time_index, censored_index, weight_index, features_indices, ...){
    call = match.call()
    #catch all optional arguments
    opt = list(...)
    for (name in names(opt) ) assign(name, opt[[name]])
    optionalParamNames = c(
        "penaltyFactor",
        "regPower",
        "regMixing",
        "iExclude",
        "df_max",
        "df_max_ever",
        "nLambda",
        "lambda",
        "lambdaFraction",
        "algorithm",
        "backtrackingFraction",
        "threshold",
        "maxIteration"
    )
    unusedParams = setdiff(names(opt),optionalParamNames)
    if(length(unusedParams)) warning('Unused parameters: ',paste(unusedParams,collapse = ', '))
    
    # TODO: Deal with the simultaneous case (copy the features?)
    
    # weights can be missing in which case they are set to 1
    if(missing(weight_index)){ 
        data$weights = 1
        weight_index = which(colnames(data) == "weights")
        names(data)[weight_index] = "weight"
    }
    
    # take all remaining columns as the features if no columns index is supplied
    if(missing(features_indices)){ 
        features_indices = seq(ncol(data))[-c(task_index, time_index, censored_index, weight_index)]
    }
    
    # check features_indices
    unusedIndices = setdiff(features_indices, seq(ncol(data)))
    if(length(unusedIndices)>0)warning('Unused feature indices (out of range): ',paste(unusedIndices,collapse = ', '))
    features_indices = subset(features_indices, !(features_indices %in% unusedIndices))
    maskedIndices = intersect(c(task_index, time_index, censored_index, weight_index), features_indices)
    if(length(maskedIndices)>0)warning('Unused feature indices (masked by other columns): ',paste(maskedIndices,collapse = ', '))
    features_indices = subset(features_indices, !(features_indices %in% maskedIndices))
    features_indices = sort(features_indices)
    nFeatures = length(features_indices)
    feature_names = names(data)[features_indices]
    
    # extract dimensions and names
    task_name = names(data)[task_index]
    censored_name = names(data)[censored_index]
    times_name = names(data)[time_index]
    weight_name = names(data)[weight_index]
    task_names = levels(factor(data[,task_index]))
    nTasks = length(task_names)
    nObs = nrow(data)
    
    # factor to integer
    data$iTask = sapply(data[,task_index], function(task) which(task_names == task))
    
    # order by task, time and censored indicator
    data = data[order(data$iTask, -data[,time_index], -data[,censored_index]),]
    
    # complete the optional parameters
    
        # penaltyFactor is set to 1 for all features if not supplied (to evetually perform adaptive Lasso)
        if(is.null(opt$penaltyFactor)) opt$penaltyFactor = rep(1, nFeatures)
        if(length(opt$penaltyFactor) != nFeatures) stop("penaltyFactor must have the same length as the features")
        if(any(opt$penaltyFactor<0)) stop("penaltyFactor must be nonegative")
    
        # regPower (0=L_infty, 2=L_2)
        if(is.null(opt$regPower)) opt$regPower = 0
        if(!(opt$regPower %in% c(0,2))) stop("regPower must be either 0 or 2")
        # regMixing (0=L_q only, 1=L_1 only)
        if(is.null(opt$regMixing)) opt$regMixing = 0
        if(opt$regMixing < 0 || opt$regMixing > 1) stop("regMixing must be between 0 and 1")
        # df_max (df_max is current number of features in model; df_max_ever is number of featuers ever to enter the model)
        if(is.null(opt$df_max)) opt$df_max = nFeatures
        if(opt$df_max > nFeatures) opt$df_max = nFeatures
        if(is.null(opt$df_max_ever)) opt$df_max_ever = nFeatures
        if(opt$df_max_ever > nFeatures) opt$df_max_ever = nFeatures
    
        # lambda sequence preparation
        if(is.null(opt$lambda)){
            # not user-supplied then nLambda is the length and lambdaFraction determines are small is the last compared to the first
            # the first is computed at the Fortran core such that it is smallest with no features in the model
            if(is.null(opt$lambdaFraction)) opt$lambdaFraction <- 1e-3
            if(is.null(opt$nLambda)) opt$nLambda <- 100
            if(opt$lambdaFraction>=1) stop("lambdaFraction should be less than 1 when lambda is not supplied")
            if(opt$lambdaFraction<1.0E-6) stop("lambdaFraction is too small")
            if(opt$nLambda<2) stop("nLambda should be at least 1")
            if(opt$nLambda>1000) stop("nLambda should be at most 1000")
            opt$lambda = rep(1, opt$nLambda) # just to pass in the foreign function call, will be overwritten
        }else{
            # user-supplied (most likely from CV)
            opt$lambdaFraction <- as.double(1) #will trigger use of user lambda in fortran core
            if(any(opt$lambda<0)) stop("lambdas should be non-negative")
            opt$lambda <- as.double(rev(sort(opt$lambda))) #lambda should be decreasing (because of warm start)
            opt$nLambda <- as.integer(length(opt$lambda))
        }
    
        # algorithm (backtracking = 2, simple = 0)
        if(is.null(opt$algorithm)) opt$algorithm = 0
        if(!(opt$algorithm %in% c(0,2))) stop("algorithm must be either 0 or 2")
    
        # backtrackingFraction (if backtracking is used, this is by how much we shrink the stepsize every iteration)
        if(is.null(opt$backtrackingFraction)) opt$backtrackingFraction = 0.7
        if(opt$backtrackingFraction < 0 || opt$backtrackingFraction > 1) stop("backtrackingFraction must be between 0 and 1")
    
        # threshold (convergence threshold)
        if(is.null(opt$threshold)) opt$threshold = 1E-6
        if(opt$threshold < 1E-16) stop("threshold must be larger than 1e-16")
    
        # maxIteration (maximum number of gradient descent cycles)
        if(is.null(opt$maxIteration)) opt$maxIteration = 1E5
        if(opt$maxIteration <1 ) stop("maxIteration must be positive")
        

    # foreign function call
    out = .Fortran(
        "mcox_solution_path", 
        PACKAGE = "MCox",
        # dimensions
        nTasks=as.integer(nTasks),
        nObs=as.integer(nObs),
        nFeatures=as.integer(nFeatures),
        # data
        time=as.double(data[,time_index]), 
        iTask=as.integer(data$iTask), 
        censored=as.integer(data[,censored_index]), 
        weight=as.double(data[,weight_index]), 
        features=apply(as.matrix.noquote(data[,features_indices]),2,as.numeric), 
        # regularization description
        penaltyFactor=as.double(opt$penaltyFactor), 
        regPower=as.integer(opt$regPower), 
        regMixing=as.double(opt$regMixing), 
        iExclude=as.integer(rep(0, nFeatures)), 
        df_max=as.integer(opt$df_max), 
        df_max_ever=as.integer(opt$df_max_ever),
        # regularization parameter sequence
        nLambda=as.integer(opt$nLambda), 
        lambda=as.double(opt$lambda), 
        lambdaFraction=as.double(opt$lambdaFraction),
        # algorithm description
        algorithm=as.integer(opt$algorithm), 
        threshold=as.double(opt$threshold), 
        backtrackingFraction=as.double(opt$backtrackingFraction), 
        maxIteration=as.integer(opt$maxIteration),
        # outputs
        failureTimes=double(nObs),
        nLambda_done=integer(1), 
        beta_path=double(opt$nLambda*nTasks*nFeatures), 
        betaNorm_path=double(opt$nLambda*nFeatures), 
        logLik_path=double(opt$nLambda*nTasks), 
        baselineHazard_path=double(opt$nLambda*nObs), 
        # kkt conditions
        kkt_condition_zero=double(opt$nLambda*nFeatures), 
        kkt_condition_nonzero=double(opt$nLambda*nFeatures), 
        # variable selection
        nBeta=integer(opt$nLambda), 
        nBeta_ever=integer(opt$nLambda), 
        iEntered=integer(nFeatures),
        # statistics
        nUpdates=integer(1), 
        nCycles=integer(1), 
        iError=integer(5)
    )
    iError = out$iError
    # iError contains 5 entries:
    # 1: 0=OK, 1=Fatal, 2=EarlyStop,
    # 2: Error type
    # 3: Feature id
    # 4: Task id
    # 5: lambda id
    if(iError[1] == 1){
        print("Fatal error code:")
        print(iError)
        # Fatal error
        switch(iError[2],
               "1" = stop("Allocation error."),
               "2" = stop(paste("Feature ", iError[3], ", task ", iError[4], " has zero variance.", sep = "")),
               "9" = stop(paste("KKT condition of feature ", iError[3], "at the ", iError[5], "th lambda value id NaN.", sep = ""))
        )
    }
    if(iError[1] == 2){
        # Warning
        switch(iError[2],
               "1" = warning(paste("Too many features included after the ", iError[5], "th lambda value. Solution for largert lambda returned.", sep = "")),
               "2" = warning(paste("Too many features ever included after the ", iError[5], "th lambda value. Solution for largert lambda returned.", sep = "")),
               "3" = warning(paste("Too many cycles executed before the ", iError[5], "th lambda value. Solution for largert lambda returned.", sep = "")),
               "4" = warning(paste("Too many cycles executed during the ", iError[5], "th lambda value. Solution for largert lambda returned.", sep = "")),
               "5" = warning(paste("Degenerate Hessian for feature ", iError[3]," at the ", iError[5], "th lambda value. Solution for largert lambda returned.", sep = "")),
               "6" = warning(paste("Strong rule cycle did not stop after 100 cycles at the ", iError[5], "th lambda value. Solution for largert lambda returned.", sep = ""))
        )
    }
    # otherwise, no problem uncovered
    
    # Store output in array format
    nLam = out$nLambda_done
    beta = array(out$beta_path, dim=c(nFeatures,nTasks,opt$nLambda))[,,seq(nLam)]
    beta_norm = array(out$betaNorm_path, dim=c(nFeatures,opt$nLambda))[,seq(nLam)]
    kkt_condition_zero = array(out$kkt_condition_zero, dim=c(nFeatures,opt$nLambda))[,seq(nLam)]
    kkt_condition_nonzero = array(out$kkt_condition_nonzero, dim=c(nFeatures,opt$nLambda))[,seq(nLam)]
    logLik = array(out$logLik_path, dim=c(nTasks,opt$nLambda))[,seq(nLam)]
    # Baseline hazard post processing
        baselineHazard = array(out$baselineHazard_path, dim=c(nObs,opt$nLambda))[,seq(nLam)]
        # find the number of groups
        nGroups = which.min(apply(baselineHazard,1,max)) - 1
        baselineHazard = baselineHazard[seq(nGroups),]
        # find group beginning and ends
        groupsIn = which(c(1, apply(apply(baselineHazard,2,diff)>0,1,max) )== 1)
        groupsOut = c(groupsIn[-1]-1, nGroups)
        # contruct lists
        tmp = list()
        failureTimes = list()
        for(k in seq(nTasks)){
            tmp[[k]] = baselineHazard[seq(groupsIn[k],groupsOut[k]),]
            failureTimes[[k]] = out$failureTimes[seq(groupsIn[k],groupsOut[k])]
        }
        baselineHazard = tmp
    # return
    out = list(
        nLambda = nLam,
        lambda = out$lambda[seq(nLam)],
        beta = beta,
        betaNorm = beta_norm,
        kkt_condition_zero = kkt_condition_zero,
        kkt_condition_nonzero = kkt_condition_nonzero,
        logLik = logLik,
        baselineHazard = baselineHazard,
        failureTimes = failureTimes,
        nBeta = out$nBeta,
        nBeta_ever = out$nBeta_ever,
        iEntered = out$iEntered,
        nUpdates = out$nUpdates,
        nCycles = out$nCycles,
        iError = out$iError
    )
}