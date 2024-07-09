# Convert output from move_seq() to moveData
as.moveData <- function(data, ID = "foo"){
  out <- data %>% 
    rename(
      step = r, 
      angle = theta_rel, 
      x = x1, 
      y = y1
    ) %>% 
    mutate(
      ID = ID
    ) %>% 
    dplyr::select(ID, step, angle, x, y)
  
  class(out) <- c("moveData", "data.frame")
  return(out)
}

# Wrapper for moveHMM::viterbi()
viterbi <- function(m, newdata = NULL){
  if (!inherits(m, "moveHMM")) 
    stop("'m' must be a moveHMM object (as output by fitHMM)")
  if (!is.null(newdata) && !inherits(newdata, "moveData")) 
    stop("'newdata' must be a moveData object (as output by prepData) or NULL")
  data <- if (!is.null(newdata)) 
    newdata
  else m$data
  nbStates <- ncol(m$mle$stepPar)
  beta <- m$mle$beta
  delta <- m$mle$delta
  stepDist <- m$conditions$stepDist
  angleDist <- m$conditions$angleDist
  stepPar <- m$mle$stepPar
  anglePar <- m$mle$anglePar
  zeroInflation <- m$conditions$zeroInflation
  knownStates <- m$knownStates
  if (nbStates == 1) 
    stop("No states to decode (nbStates=1)")
  covsCol <- which(names(data) != "ID" & names(data) != "x" & 
                     names(data) != "y" & names(data) != "step" & names(data) != 
                     "angle")
  nbCovs <- length(covsCol) - 1
  covs <- data[, covsCol]
  probs <- moveHMM:::allProbs(data, nbStates, stepDist, angleDist, stepPar, 
                    anglePar, zeroInflation, knownStates)
  trMat <- moveHMM:::trMatrix_rcpp(nbStates, beta, as.matrix(covs))
  nbAnimals <- length(unique(data$ID))
  aInd <- c(1, which(data$ID[-1] != data$ID[-nrow(data)]) + 
              1)
  allStates <- NULL
  for (zoo in 1:nbAnimals) {
    nbObs <- length(which(data$ID == unique(data$ID)[zoo]))
    obsInd <- which(!is.na(data$step) & !is.na(data$angle))
    if (zoo != nbAnimals) {
      p <- probs[aInd[zoo]:(aInd[zoo + 1] - 1), ]
      tm <- trMat[, , aInd[zoo]:(aInd[zoo + 1] - 1)]
    }
    else {
      p <- probs[aInd[zoo]:nrow(probs), ]
      tm <- trMat[, , aInd[zoo]:nrow(probs)]
    }
    xi <- matrix(NA, nbObs, nbStates)
    foo <- delta * p[1, ]
    xi[1, ] <- foo/sum(foo)
    for (i in 2:nbObs) {
      foo <- apply(xi[i - 1, ] * tm[, , i], 2, max) * p[i, ]
      xi[i, ] <- foo/sum(foo)
    }
    stSeq <- rep(NA, nbObs)
    stSeq[nbObs] <- which.max(xi[nbObs, ])
    for (i in (nbObs - 1):1) {
      stSeq[i] <- which.max(tm[, stSeq[i + 1], i + 1] * xi[i, ])
    }
    allStates <- c(allStates, stSeq)
  }
  return(allStates)
}

# Wrapper for moveHMM::stateProbs()
stateProbs <- function(model){
  moveHMM::stateProbs(model)
}

# Wrapper for moveHMM::fitHMM()
fit_HMM <- function(data, 
                    stepPar0 = c(2.5, 1.5, 1.5, 0.8), 
                    anglePar0 = c(0, pi, 0.1, 0.1), 
                    stepDist = "lnorm", angleDist = "wrpcauchy", 
                    formula = ~1, 
                    ...){
  moveHMM::fitHMM(data = data, 
                  nbStates = 2,
                  stepPar0 = stepPar0, 
                  stepDist = stepDist, 
                  angleDist = angleDist,
                  anglePar0 = anglePar0, 
                  knownStates = ifelse(data$step == max(data$step, na.rm = TRUE), 1, NA), # label the longest distance movement as state 1
                  formula = formula, 
                  ...)
}

trans_mat <- function(x){
  x$mle$gamma
}


