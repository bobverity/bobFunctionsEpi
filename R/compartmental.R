
# -----------------------------------
#' plotModel
#'
#' Default plot of model output from functions in this package.
#'
#' @param x model output
#'
#' @export

plotModel <- function(x, vars=NULL, ...) {
    
    # get input arguments
    args <- list(...)
    argNames <- names(args)
    
    # subset to chosen vars
    tvec <- x$time
    if (is.null(vars)) {
        vars <- setdiff(names(x), "time")
    }
    vars <- intersect(vars, names(x))
    if (length(vars)==0) {
        stop("none of chosen variables in x")
    }
    x <- subset(x, select=vars)
    
    # get basic data properties
    nvar <- length(vars)
    maxVal <- max(x,na.rm=TRUE)
    maxt <- max(tvec,na.rm=TRUE)
    
    # set defaults on undefined arguments
    if (! "xlim" %in% argNames) {
        args$xlim <- c(0, maxt)
    }
    if (! "ylim" %in% argNames) {
        args$ylim <- c(0, 1.1*maxVal)
    }
    if (! "xlab" %in% argNames) {
        args$xlab <- "time"
    }
    if (! "ylab" %in% argNames) {
        args$ylab <- ""
    }
    
    # fixed arguments, or arguments that have special meaning
    args$type <- NULL
    
    # extract colours or define defaults
    if ("col" %in% argNames) {
        colVec <- rep_len(args$col,nvar)
        args$col <- NULL
    } else {
        colVec <- 1:nvar
    }
    
    # plot with finalised list of parameters
    do.call(plot, c(list(x=0, type="n"), args))    # create empty plot
    for (i in 1:nvar) {
        do.call(lines, c(list(x=tvec, y=x[,i], col=colVec[i]), args))
    }
}

# -----------------------------------
#' simQuantiles
#'
#' Runs a given stochastic simulation function many times, computing the mean and quantiles over replicates. Note that this method will only work with simulations that have a fixed time step, i.e. synchronous or hybrid simulations, and not with asynchronous simulations. In the hybrid case the maxIterations limit cannot be reached in any simulation.
#'
#' @param FUN the stochastic simulation function to use.
#' @param args a list of arguments to the function.
#' @param reps number of times to repeat the stochastic simulation.
#' @param quantiles which quantiles to compute over replicates.
#'
#' @export

simQuantiles <- function(FUN="SIS_stochastic_hybrid", args=list(), reps=1e2, quantiles=c(0.05,0.5,0.95)) {
    
    # run function once to get dimensions and variable names
    testOutput <- do.call(FUN, args)
    varNames <- setdiff(names(testOutput),"time")
    
    # repeat simulation many times and store in array
    simArray <- array(0, dim=c(nrow(testOutput), length(varNames), reps))
    simArray[,,1] <- as.matrix(testOutput[,varNames])
    if (reps>1) {
        for (i in 2:reps) {
            simArray[,,i] <- as.matrix(do.call(FUN, args)[,varNames])
        }
    }
    
    # compute mean and quantiles over replicates and store in data frame
    df <- data.frame(time=testOutput$time)
    for (i in 1:length(varNames)) {
        m <- rowMeans(simArray[,i,,drop=FALSE], na.rm=TRUE)
        q <- apply(simArray[,i,,drop=FALSE], 1, quantile, prob=quantiles, na.rm=TRUE)
        df_new <- as.data.frame(t(rbind(m,q)))
        names(df_new) <- paste(varNames[i], c("mean", paste("Q", quantiles, sep="")), sep="_")
        df <- cbind(df, df_new)
    }
    
    # return summary data frame
    return(df)
}

# -----------------------------------
#' SIS_analytical
#'
#' Returns analytical solution to deterministic SIS model. At equilibrium the number of infectives is given by \eqn{I* = N(1 - r/beta)}.
#'
#' @param beta contact rate.
#' @param r recovery rate.
#' @param I_init initial number of infectious individuals.
#' @param N total number of individuals in population.
#' @param times vector of times at which output should be returned.
#'
#' @export

SIS_analytical <- function(beta=1, r=0.25, I_init=10, N=1e3, times=0:100) {
    
    # analytical solution to I
    I <- (beta-r)/(beta/N + (beta-r-beta*I_init/N)/I_init*exp(-(beta-r)*times))
    
    # create output object
    ret <- data.frame(time=times, S=N-I, I=I)
    
    return(ret)
}

# -----------------------------------
#' SIS_deterministic
#'
#' Returns solution to deterministic SIS model using the \code{odin} package.
#'
#' @param beta contact rate.
#' @param r recovery rate.
#' @param I_init initial number of infectious individuals.
#' @param N total number of individuals in population.
#' @param times vector of times at which output should be returned.
#'
#' @export

SIS_deterministic <- function(beta=1, r=0.25, I_init=10, N=1e3, times=0:100) {
    
    stop("function commented out due to issues with odin")
    
    # solve ode
    #mod <- SIS_deterministic_odin2(beta=beta, r=r, I_init=I_init, N=N)
    #output <- as.data.frame(mod$run(times))
    #names(output)[1] <- 'time'
    
    #return(output)
}

# -----------------------------------
#' SIS_stochastic_async
#'
#' Draw from asynchronous stochastic SIS model. Return state of the system at all time points at which any event occurs. Stop when maxIterations is reached.
#'
#' @param beta contact rate.
#' @param r recovery rate.
#' @param I_init initial number of infectious individuals.
#' @param N total number of individuals in population.
#' @param maxIterations exit if this number of iterations is reached.
#'
#' @export

SIS_stochastic_async <- function(beta=1, r=0.25, I_init=100, N=1e3, maxIterations=1e4) {
    
    # run model
    args <- list(beta=beta, r=r, I_start=I_init, N=N, maxIterations=maxIterations)
    rawOutput <- SIS_stochastic_async_cpp(args)
    
    # format output object
    ret <- as.data.frame(rawOutput)
    ret <- subset(ret, I>=0)
    ret$S <- N - ret$I
    
    return(ret)
}

# -----------------------------------
#' SIS_stochastic_hybrid
#'
#' Draw from stochastic SIS model using a compromise between a synchronous and an asynchronous algorithm. The basic algorithm is asynchronous (based on Gillespie's algorithm), but values are only stored and returned at discrete time points. The function still exits automatically at a defined maximum number of iterations.
#'
#' @param beta contact rate.
#' @param r recovery rate.
#' @param I_init initial number of infectious individuals.
#' @param N total number of individuals in population.
#' @param times vector of times at which output should be returned.
#' @param maxIterations exit if this number of iterations is reached.
#'
#' @export

SIS_stochastic_hybrid <- function(beta=1, r=0.25, I_init=100, N=1e3, times=0:100, maxIterations=1e4) {
    
    # run model
    args <- list(beta=beta, r=r, I_start=I_init, N=N, t_vec=times, maxIterations=maxIterations)
    rawOutput <- SIS_stochastic_hybrid_cpp(args)
    
    # format output object
    ret <- as.data.frame(rawOutput)
    ret[ret<0] <- NA
    ret$S <- N - ret$I
    
    return(ret)
}

# -----------------------------------
#' SIS_stochastic_sync
#'
#' Draw from synchronous stochastic SIS model. Return infectives at known time points. Note that results of the synchronous method only match up with the asynchronous method when the time step is small relative to the rates that drive the system.
#'
#' @param beta contact rate.
#' @param r recovery rate.
#' @param I_init initial number of infectious individuals.
#' @param N total number of individuals in population.
#' @param times vector of times at which output should be returned.
#'
#' @export

SIS_stochastic_sync <- function(beta=1, r=0.25, I_init=100, N=1e3, times=0:100) {
    
    # run model
    args <- list(beta=beta, r=r, I_start=I_init, N=N, t_vec=times)
    rawOutput <- SIS_stochastic_sync_cpp(args)
    
    # format output object
    ret <- as.data.frame(rawOutput)
    ret[ret<0] <- NA
    ret$S <- N - ret$I
    
    return(ret)
}

# -----------------------------------
#' SIR_deterministic
#'
#' Returns solution to deterministic SIR model using the \code{odin} package.
#'
#' @param beta contact rate.
#' @param r recovery rate.
#' @param mu natural death rate (same in all compartments). Also rate of new births into susceptible compartment.
#' @param I_init initial number of infectious individuals.
#' @param R_init initial number of recovered (immune) individuals.
#' @param N total number of individuals in population.
#' @param times vector of times at which output should be returned.
#'
#' @export

SIR_deterministic <- function(beta=1, r=0.25, mu=0.01, I_init=100, R_init=0, N=1e3, times=0:100) {
    
    stop("function commented out due to issues with odin")
    
    # solve ode
    #mod <- SIR_deterministic_odin(beta=beta, r=r, mu=mu, I_init=I_init, R_init=R_init, N=N)
    #output <- as.data.frame(mod$run(times))
    #names(output)[1] <- 'time'
    
    #return(output)
}

# -----------------------------------
#' SIR_stochastic_async
#'
#' Draw from asynchronous stochastic SIR model with optional birth/death. Return state of the system at all time points at which any event occurs. Stop when maxIterations is reached. Note that natural deaths are exactly matched by births into the susceptible state in this formulation, meaning the population size stays constant at N.
#'
#' @param beta contact rate.
#' @param r recovery rate.
#' @param mu natural death rate (same in all compartments). Also rate of new births into susceptible compartment.
#' @param I_init initial number of infectious individuals.
#' @param R_init initial number of recovered (immune) individuals.
#' @param N total number of individuals in population.
#' @param maxIterations exit if this number of iterations is reached.
#'
#' @export

SIR_stochastic_async <- function(beta=1, r=0.25, mu=0.01, I_init=100, R_init=0, N=1e3, maxIterations=1e4) {
    
    # run model
    args <- list(beta=beta, r=r, mu=mu, I_init=I_init, R_init=R_init, N=N, maxIterations=maxIterations)
    rawOutput <- SIR_stochastic_async_cpp(args)
    
    # format output object
    ret <- as.data.frame(rawOutput)
    ret <- subset(ret, I>=0)
    
    return(ret)
}

# -----------------------------------
#' SIR_stochastic_hybrid
#'
#' Draw from stochastic SIR model with optional birth/death using a compromise between a synchronous and an asynchronous algorithm. The basic algorithm is asynchronous (based on Gillespie's algorithm), but values are only stored and returned at discrete time points. The function still exits automatically at a defined maximum number of iterations. Note that natural deaths are exactly matched by births into the susceptible state in this formulation, meaning the population size stays constant at N.
#'
#' @param beta contact rate.
#' @param r recovery rate.
#' @param mu natural death rate (same in all compartments). Also rate of new births into susceptible compartment.
#' @param I_init initial number of infectious individuals.
#' @param R_init initial number of recovered (immune) individuals.
#' @param N total number of individuals in population.
#' @param times vector of times at which output should be returned.
#' @param maxIterations exit if this number of iterations is reached.
#'
#' @export

SIR_stochastic_hybrid <- function(beta=1, r=0.25, mu=0.01, I_init=100, R_init=0, N=1e3, times=0:100, maxIterations=1e4) {
    
    # run model
    args <- list(beta=beta, r=r, mu=mu, I_init=I_init, R_init=R_init, N=N, t_vec=times, maxIterations=maxIterations)
    rawOutput <- SIR_stochastic_hybrid_cpp(args)
    
    # format output object
    ret <- as.data.frame(rawOutput)
    ret[ret<0] <- NA
    
    return(ret)
}

# -----------------------------------
#' SIR_stochastic_sync
#'
#' Draw from synchronous stochastic SIR model with optional birth/death. Return state of the system at known time points. Note that natural deaths are exactly matched by births into the susceptible state in this formulation, meaning the population size stays constant at N. Results of the synchronous method only match up with the asynchronous method when the time step is small relative to the rates that drive the system.
#'
#' @param beta contact rate.
#' @param r recovery rate.
#' @param mu natural death rate (same in all compartments). Also rate of new births into susceptible compartment.
#' @param I_init initial number of infectious individuals.
#' @param R_init initial number of recovered (immune) individuals.
#' @param N total number of individuals in population.
#' @param times vector of times at which output should be returned.
#'
#' @export

SIR_stochastic_sync <- function(beta=1, r=0.25, mu=0.01, I_init=100, R_init=0, N=1e3, times=0:100) {
    
    # run model
    args <- list(beta=beta, r=r, mu=mu, I_init=I_init, R_init=R_init, N=N, t_vec=times)
    rawOutput <- SIR_stochastic_sync_cpp(args)
    
    # format output object
    ret <- as.data.frame(rawOutput)
    ret[ret<0] <- NA
    
    return(ret)
}

# -----------------------------------
#' SLIR_deterministic
#'
#' Returns solution to deterministic SLIR model, where L is an incubation (lag) stage of defined length. Solves delay differential equation using the \code{odin} package.
#'
#' @param beta contact rate.
#' @param dur_lag length of time in incubation state.
#' @param r recovery rate.
#' @param I_init initial number of infectious individuals.
#' @param R_init initial number of recovered (immune) individuals.
#' @param N total number of individuals in population.
#' @param times vector of times at which output should be returned.
#'
#' @export

SLIR_deterministic <- function(beta=0.5, dur_lag=1, r=0.25, I_init=10, R_init=0, N=1e3, times=0:100) {
    
    stop("function commented out due to issues with odin")
    
    # solve ode
    #mod <- SLIR_deterministic_odin(beta=beta, dur_lag=dur_lag, r=0.25, I_init=I_init, R_init=R_init, N=N)
    #output <- as.data.frame(mod$run(times))
    #names(output)[1] <- 'time'
    
    #return(output)
}

# -----------------------------------
#' SLIR_stochastic_async
#'
#' Draw from asynchronous stochastic SLIR model, where L is an incubation (lag) stage of defined length. Return state of the system at all time points at which any event occurs. Stop when maxIterations is reached.
#'
#' @param beta contact rate.
#' @param dur_lag length of time in incubation state.
#' @param r recovery rate.
#' @param I_init initial number of infectious individuals.
#' @param R_init initial number of recovered (immune) individuals.
#' @param N total number of individuals in population.
#' @param maxIterations exit if this number of iterations is reached.
#'
#' @export

SLIR_stochastic_async <- function(beta=1, dur_lag=1, r=0.25, I_init=100, R_init=0, N=1e3, maxIterations=1e4) {
    
    # run model
    args <- list(beta=beta, dur_lag=dur_lag, r=r, I_init=I_init, R_init=R_init, N=N, maxIterations=maxIterations)
    rawOutput <- SLIR_stochastic_async_cpp(args)
    
    # format output object
    ret <- as.data.frame(rawOutput)
    ret <- subset(ret, I>=0)
    
    return(ret)
}

# -----------------------------------
#' SLIR_stochastic_hybrid
#'
#' Draw from stochastic SLIR model, where L is an incubation (lag) stage of defined length, using a compromise between a synchronous and an asynchronous algorithm. The basic algorithm is asynchronous (based on Gillespie's algorithm), but values are only stored and returned at discrete time points. The function still exits automatically at a defined maximum number of iterations.
#'
#' @param beta contact rate.
#' @param dur_lag length of time in incubation state.
#' @param r recovery rate.
#' @param I_init initial number of infectious individuals.
#' @param R_init initial number of recovered (immune) individuals.
#' @param N total number of individuals in population.
#' @param times vector of times at which output should be returned.
#' @param maxIterations exit if this number of iterations is reached.
#'
#' @export

SLIR_stochastic_hybrid <- function(beta=1, dur_lag=1, r=0.25, I_init=100, R_init=0, N=1e3, times=0:100, maxIterations=1e4) {
    
    # run model
    args <- list(beta=beta, dur_lag=dur_lag, r=r, I_init=I_init, R_init=R_init, N=N, t_vec=times, maxIterations=maxIterations)
    rawOutput <- SLIR_stochastic_hybrid_cpp(args)
    
    # format output object
    ret <- as.data.frame(rawOutput)
    ret[ret<0] <- NA
    
    return(ret)
}


