#' @useDynLib MCox
"_PACKAGE"


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
MCox = function(data, task_index, time_index, censored_index, weight_index, features_indeces){
    if(missing(weight_index)){
        data$weights = 1
        weight_index = which(colnames(data) == "weights")
    }
    if(missing(features_indeces)){#take all remaining columns
        features_indeces = seq(ncol(data))[-c(task_index, time_index, censored_index, weight_index)]
    }
    task_names = levels(factor(data[,task_index]))
    nTasks = length(task_names)
    data$iTask = sapply(data[,task_index], function(task) which(task_names == task))
    data = data[order(data$iTask, -data[,time_index], -data[,censored_index]),]
    nObs = nrow(data)
    nFeatures = length(features_indeces)
    nLambda = 10
    out <- .Fortran(
        "mcox_solution_path", 
        PACKAGE = "MCox",
        nTasks=as.integer(nTasks),
        nObs=as.integer(nObs),
        nFeatures=as.integer(nFeatures),
        time=as.double(data[,time_index]), 
        iTask=as.integer(data$iTask), 
        censored=as.integer(data[,censored_index]), 
        weight=as.double(data[,weight_index]), 
        features=apply(as.matrix.noquote(data[,features_indeces]),2,as.numeric), 
        penaltyFactor=as.double(rep(1, nFeatures)), 
        regPower=as.double(2), 
        regMixing=as.double(0), 
        iExclude=as.integer(rep(0, nFeatures)), 
        df_max=as.integer(nFeatures), 
        df_max_ever=as.integer(nFeatures),
        nLambda=as.integer(nLambda), 
        lambda=as.double(seq(nLambda,1,-1)/nLambda), 
        lambdaFraction=as.double(1e-3),
        algorithm=as.integer(0), 
        threshold=as.double(1e-5), 
        backtrackingFraction=as.double(0.8), 
        maxIteration=as.integer(1e5),
        nLambda_done=integer(1), 
        beta_path=double(nLambda*nTasks*nFeatures), 
        betaNorm_path=double(nLambda*nFeatures), 
        logLik_path=double(nLambda*nTasks), 
        baselineHazard_path=double(nLambda*nObs), 
        nBeta=integer(nLambda), 
        nBeta_ever=integer(nLambda), 
        iEntered=integer(nFeatures),
        nUpdates=integer(1), 
        nCycles=integer(1), 
        iError=integer(5)
    )
    return(out)
}