#' Projection onto a L1 Ball
#'
#' @param vec 
#' @param radius 
#'
#' @return projection
#' @export
#'
#' @examples
#' \dontrun{
#'    ProjectL1(1:3, 1)
#' }
ProjectL1 = function(vec, radius){
    p = as.integer(length(vec))
    out <- .Fortran(
        "ProjB1Mich", 
        PACKAGE = "MCox", 
        vector = as.double(vec), 
        dim = p, 
        radius = as.double(radius), 
        projection = double(p)
    )
    return(out$projection)
}


#' Projection onto a L1 Ball
#'
#' @param vec 
#' @param radius 
#'
#' @return projection
#' @export
#'
#' @examples
#' \dontrun{
#'    ProjectL1(1:3, 1)
#' }
Proximal = function(lambda, penalty_factor, stepsize, alpha, q, beta, gradient){
    K = as.integer(length(beta))
    out <- .Fortran(
        "Proximal", 
        PACKAGE = "MCox", 
        lam = as.double(lambda),
        pf = as.double(penalty_factor),
        sig = as.double(stepsize),
        alpah = as.double(alpha),
        reg = as.integer(q),
        ntasks = K,
        beta = as.double(beta),
        grad = as.double(gradient)
    )
    return(out$beta)
}


#' Projection onto a L1 Ball
#'
#' @param vec 
#' @param radius 
#'
#' @return projection
#' @export
#'
#' @examples
#' \dontrun{
#'    ProjectL1(1:3, 1)
#' }
Preprocessing = function(data, task_index, time_index, censored_index, weight_index, features_indeces){
    if(missing(weight_index)){
        data$weights = 1
        weight_index = which(colnames(data) == "weights")
    }
    if(missing(features_indeces)){#take all remaining columns
        features_indeces = seq(ncol(data))[-c(task_index, time_index, censored_index, weight_index)]
    }
    task_names = levels(factor(data[,task_index]))
    data$iTask = sapply(data[,task_index], function(task) which(task_names == task))
    data = data[order(data$iTask, -data[,time_index], -data[,censored_index]),]
    out <- .Fortran(
        "Preprocessing", 
        PACKAGE = "MCox", 
        time=as.double(data[,time_index]), 
        iTask=as.integer(data$iTask), 
        censored=as.integer(data[,censored_index]), 
        weight=as.double(data[,weight_index]), 
        features=apply(as.matrix.noquote(data[,features_indeces]),2,as.numeric), 
        nObs=as.integer(nrow(data)), 
        nFeatures=as.integer(length(features_indeces)), 
        nTasks=integer(1), 
        obsBoundByTaskIn=integer(100), 
        obsBoundByTaskOut=integer(100),
        groupBoundByTaskIn=integer(100), 
        groupBoundByTaskOut=integer(100),
        nGroups=integer(1), 
        obsBoundByGroupIn=integer(100), 
        obsBoundByGroupOut=integer(100),
        obsGroupId=integer(100),
        failureTimes=double(100),
        tiesTotalWeight=double(100),
        featureMean=double(length(features_indeces) * length(task_names)),
        featureStdDev=double(length(features_indeces) * length(task_names)),
        iError=integer(5)
    )
    data$groupId = out$obsGroupId[seq(nrow(data))]
    data$timeOfRiskSet = out$timeOfRiskSet[seq(nrow(data))]
    return(list(out, data))
}






#' Projection onto a L1 Ball
#'
#' @param vec 
#' @param radius 
#'
#' @return projection
#' @export
#'
#' @examples
#' \dontrun{
#'    ProjectL1(1:3, 1)
#' }
LogLikelihood = function(obj){
    out <- .Fortran(
        "LogLikelihood", 
        PACKAGE = "MCox",
        nTasks = as.integer(obj$nTasks), 
        nObs = as.integer(obj$nObs), 
        nGroups = as.integer(obj$nGroups), 
        obsBoundByGroupIn = as.integer(obj$obsBoundByGroupIn[1:obj$nGroups]), 
        obsBoundByGroupOut = as.integer(obj$obsBoundByGroupOut[1:obj$nGroups]),
        groupBoundByTaskIn = as.integer(obj$groupBoundByTaskIn[1:obj$nTasks]), 
        groupBoundByTaskOut = as.integer(obj$groupBoundByTaskOut[1:obj$nTasks]),
        weight = as.double(obj$weight[1:obj$nObs]), 
        censored = as.integer(obj$censored[1:obj$nObs]), 
        tiesTotalWeight = as.double(obj$tiesTotalWeight[1:obj$nGroups]), 
        linearPredictor = as.double(rep(0, obj$nObs)), 
        logLik = double(obj$nTasks)
    )
    return(out)
}



#' Projection onto a L1 Ball
#'
#' @param vec 
#' @param radius 
#'
#' @return projection
#' @export
#'
#' @examples
#' \dontrun{
#'    ProjectL1(1:3, 1)
#' }
PartialDerivatives = function(obj, j=1){
    out <- .Fortran(
        "PartialDerivatives", 
        PACKAGE = "MCox",
        nTasks = as.integer(obj$nTasks), 
        nObs = as.integer(obj$nObs), 
        nGroups = as.integer(obj$nGroups), 
        obsBoundByGroupIn = as.integer(obj$obsBoundByGroupIn[1:obj$nGroups]), 
        obsBoundByGroupOut = as.integer(obj$obsBoundByGroupOut[1:obj$nGroups]),
        groupBoundByTaskIn = as.integer(obj$groupBoundByTaskIn[1:obj$nTasks]), 
        groupBoundByTaskOut = as.integer(obj$groupBoundByTaskOut[1:obj$nTasks]),
        weight = as.double(obj$weight[1:obj$nObs]), 
        features = as.double(obj$features[1:obj$nObs, j]), 
        censored = as.integer(obj$censored[1:obj$nObs]), 
        tiesTotalWeight = as.double(obj$tiesTotalWeight[1:obj$nGroups]), 
        linearPredictor = as.double(rep(0, obj$nObs)), 
        gradient = double(obj$nTasks), 
        hessian = double(obj$nTasks)
    )
    return(out)
}




#' Projection onto a L1 Ball
#'
#' @param vec 
#' @param radius 
#'
#' @return projection
#' @export
#'
#' @examples
#' \dontrun{
#'    ProjectL1(1:3, 1)
#' }
GPG_Cycle = function(obj, lambda=1.0){
    out <- .Fortran(
        "gpg_cycle", 
        PACKAGE = "MCox",
        nTasks = as.integer(obj$nTasks), 
        nFeatures = as.integer(obj$nFeatures), 
        nObs = as.integer(obj$nObs), 
        nGroups = as.integer(obj$nGroups), 
        obsBoundByTaskIn = as.integer(obj$obsBoundByTaskIn[1:obj$nTasks]), 
        obsBoundByTaskOut = as.integer(obj$obsBoundByTaskOut[1:obj$nTasks]),
        groupBoundByTaskIn = as.integer(obj$groupBoundByTaskIn[1:obj$nTasks]), 
        groupBoundByTaskOut = as.integer(obj$groupBoundByTaskOut[1:obj$nTasks]),
        obsBoundByGroupIn = as.integer(obj$obsBoundByGroupIn[1:obj$nGroups]), 
        obsBoundByGroupOut = as.integer(obj$obsBoundByGroupOut[1:obj$nGroups]),
        tiesTotalWeight = as.double(obj$tiesTotalWeight[1:obj$nGroups]), 
        weight = as.double(obj$weight[1:obj$nObs]), 
        features = as.double(obj$features[1:obj$nObs, ]), 
        censored = as.integer(obj$censored[1:obj$nObs]),
        iActive = as.integer(rep(1,obj$nFeatures)),
        lambda = as.double(lambda),
        penaltyFactor = as.double(rep(1,obj$nFeatures)),
        regPower = as.integer(2),
        regMixing = as.double(0.0),
        beta = as.double(rep(0, obj$nFeatures * obj$nTasks)),
        linearPredictor = as.double(rep(0, obj$nObs)), 
        logLik = as.double(rep(0, obj$nTasks)),
        stepsize = as.double(rep(0, obj$nFeatures)), 
        nUpdates = as.integer(100), 
        gradient = double(obj$nFeatures*obj$nTasks),
        hessian = double(obj$nFeatures*obj$nTasks), 
        iError = integer(5)
    )
    return(out)
}




#' Projection onto a L1 Ball
#'
#' @param vec 
#' @param radius 
#'
#' @return projection
#' @export
#'
#' @examples
#' \dontrun{
#'    ProjectL1(1:3, 1)
#' }
GPG_Descent = function(obj, lambda=1.0){
    out <- .Fortran(
        "gpg_descent", 
        PACKAGE = "MCox",
        nTasks = as.integer(obj$nTasks), 
        nFeatures = as.integer(obj$nFeatures), 
        nObs = as.integer(obj$nObs), 
        nGroups = as.integer(obj$nGroups), 
        obsBoundByTaskIn = as.integer(obj$obsBoundByTaskIn[1:obj$nTasks]), 
        obsBoundByTaskOut = as.integer(obj$obsBoundByTaskOut[1:obj$nTasks]),
        groupBoundByTaskIn = as.integer(obj$groupBoundByTaskIn[1:obj$nTasks]), 
        groupBoundByTaskOut = as.integer(obj$groupBoundByTaskOut[1:obj$nTasks]),
        obsBoundByGroupIn = as.integer(obj$obsBoundByGroupIn[1:obj$nGroups]), 
        obsBoundByGroupOut = as.integer(obj$obsBoundByGroupOut[1:obj$nGroups]),
        tiesTotalWeight = as.double(obj$tiesTotalWeight[1:obj$nGroups]), 
        weight = as.double(obj$weight[1:obj$nObs]), 
        features = as.double(obj$features[1:obj$nObs, ]), 
        censored = as.integer(obj$censored[1:obj$nObs]),
        lambda = as.double(lambda),
        penaltyFactor = as.double(rep(1,obj$nFeatures)),
        regPower = as.integer(2),
        regMixing = as.double(0.0),
        algorithm = as.integer(1),
        threshold = as.double(1e-3),
        backtrackingFactor = as.double(0.001),
        maxIteration = as.integer(1e5),
        beta = as.double(rep(0, obj$nFeatures * obj$nTasks)),
        linearPredictor = as.double(rep(0, obj$nObs)), 
        logLik = as.double(rep(0, obj$nTasks)), 
        nCycles = integer(1), 
        nUpdates = integer(1),
        iError = integer(5)
    )
    return(out)
}