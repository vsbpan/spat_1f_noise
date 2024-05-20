# Wrapper for moveHMM::plotStates()
plot_states <- function(x, animals = NULL){
  moveHMM::plotStates(x, animals, ask = FALSE)
}

# S3 plotting method for moveHMM
plot.moveHMM <- function(fit, type = c("track", "sl", "ta", "residuals", "states"), ...){
  type <- match.arg(type)
  
  f <- switch(type,
              "track" = plot_track_states,
              "sl" = plot_sl_states, 
              "ta" = plot_ta_states, 
              "residuals" = moveHMM::plotPR,
              "states" = plot_states)
  f(fit, ...)
}

registerS3method("plot", class = "moveHMM", method = plot.moveHMM)

# Same code as moveHMM:::getPalette()
.getPalette <- function(nbStates){
  if (nbStates < 8) {
    pal <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
             "#0072B2", "#D55E00", "#CC79A7")
    col <- pal[1:nbStates]
  }
  else {
    hues <- seq(15, 375, length = nbStates + 1)
    col <- hcl(h = hues, l = 65, c = 100)[1:nbStates]
  }
  return(col)
}

# Modified moveHMM::getPlotData() to have a nicer log grid
.getPlotDatalog_hmm <- function(m, stepGrid){
  nbStates <- ncol(m$mle$stepPar)
  out <- list()
  w <- 1
  if (nbStates > 1) {
    states <- viterbi(m)
    w <- sapply(1:nbStates, function(s) length(which(states == 
                                                       s))/length(states))
  }
  if (m$conditions$zeroInflation) {
    zeromass <- m$mle$stepPar[nrow(m$mle$stepPar), ]
    stepPar <- as.matrix(m$mle$stepPar[-nrow(m$mle$stepPar), 
    ])
  }
  else {
    stepPar <- m$mle$stepPar
  }
  stepFun <- paste0("d", m$conditions$stepDist)
  stepDensities <- data.frame(step = stepGrid)
  for (state in 1:nbStates) {
    stepArgs <- list(stepGrid)
    for (j in 1:nrow(stepPar)) {
      stepArgs[[j + 1]] <- stepPar[j, state]
    }
    if (m$conditions$stepDist == "gamma") {
      shape <- stepArgs[[2]]^2/stepArgs[[3]]^2
      rate <- stepArgs[[2]]/stepArgs[[3]]^2
      stepArgs[[2]] <- shape
      stepArgs[[3]] <- rate
    }
    stepDensities[[paste0("state", state)]] <- w[state] * 
      do.call(stepFun, stepArgs)
    if (m$conditions$zeroInflation) {
      stepDensities[[paste0("state", state)]] <- (1 - 
                                                    zeromass[state]) * stepDensities[[paste0("state", 
                                                                                             state)]]
    }
  }
  stepDensities$total <- rowSums(stepDensities[, -1])
  out$step <- stepDensities
  return(out)
}

# A simpler version of moveHMM:::plot.moveHMM() for step lengths
plot_sl_states <- function(fit, states = NULL, col = NULL, nclass = 50, log = TRUE){
  mfrow_og <- par()$mfrow
  mar_og <- par()$mar
  
  nbStates <- ncol(fit$mle$stepPar)
  
  
  if (is.null(col) | (!is.null(col) & length(col) != nbStates)) {
    col <- .getPalette(nbStates = nbStates)
  }
  
  if(is.null(states)){
    if (nbStates > 1) {
      states <- viterbi(fit)
    } else {
      states <- rep(1, nrow(fit$data))
    }
  }
  
  
  par(mar = c(5, 4, 4, 2) - c(0, 0, 2, 1))
  legText <- c(paste("state", 1:nbStates), "total")
  lty <- c(rep(1, nbStates), 2)
  lwd <- c(rep(1, nbStates), 2)
  lineCol <- c(col, "black")
  
  if(log){
    h <- hist(log(fit$data$step), plot = FALSE, nclass = nclass)
    
    distData <- .getPlotDatalog_hmm(fit, exp(h$mids))
    
    ymax <- 1.3 * max(h$density)
    maxdens <- max(distData$step$total)
    if (maxdens > ymax & maxdens < 1.5 * ymax) {
      ymax <- maxdens
    }
    hist(log(fit$data$step), ylim = c(0, ymax), prob = TRUE, main = "", 
         xlab = "log(step length)", col = "lightgrey", border = "white", 
         nclass = nclass)
  } else {
    distData <- moveHMM::getPlotData(fit, type = "dist")
    
    h <- hist(fit$data$step, plot = FALSE, nclass = nclass)
    ymax <- 1.3 * max(h$density)
    maxdens <- max(distData$step$total)
    if (maxdens > ymax & maxdens < 1.5 * ymax) {
      ymax <- maxdens
    }
    hist(fit$data$step, ylim = c(0, ymax), prob = TRUE, main = "", 
         xlab = "step length", col = "lightgrey", border = "white", 
         nclass = nclass)
  }
  
  if(log){
    for (i in 1:(nbStates + 1)) {
      lines(log(distData$step$step), distData$step[, i + 1] * distData$step$step, 
            col = lineCol[i], 
            lty = lty[i], lwd = lwd[i])
    }
  } else {
    for (i in 1:(nbStates + 1)) {
      lines(distData$step$step, distData$step[, i + 1], 
            col = lineCol[i], 
            lty = lty[i], lwd = lwd[i])
    }
  }
  
  legend("top", legText, lwd = lwd, col = lineCol, lty = lty, 
         bty = "n")
  
  
  par(mfrow = mfrow_og)
  par(mar = mar_og)
}

# A simpler version of moveHMM:::plot.moveHMM() for turn angles
plot_ta_states <- function(fit, states = NULL, col = NULL, nclass = 50){
  mfrow_og <- par()$mfrow
  mar_og <- par()$mar
  
  
  nbStates <- ncol(fit$mle$stepPar)
  
  if (is.null(col) | (!is.null(col) & length(col) != nbStates)) {
    col <- .getPalette(nbStates = nbStates)
  }
  
  if(is.null(states)){
    if (nbStates > 1) {
      states <- viterbi(fit)
    } else {
      states <- rep(1, nrow(fit$data))
    }
  }
  
  par(mar = c(5, 4, 4, 2) - c(0, 0, 2, 1))
  distData <- moveHMM::getPlotData(fit, type = "dist")
  legText <- c(paste("state", 1:nbStates), "total")
  lty <- c(rep(1, nbStates), 2)
  lwd <- c(rep(1, nbStates), 2)
  lineCol <- c(col, "black")
  
  h1 <- hist(fit$data$angle, plot = FALSE, nclass = nclass)
  breaks <- seq(-pi, pi, length = length(h1$breaks))
  h2 <- hist(fit$data$angle, plot = FALSE, nclass = nclass)
  ymax <- 1.3 * max(h2$density)
  hist(fit$data$angle, ylim = c(0, ymax), prob = TRUE, main = "",
       xlab = "turning angle", col = "lightgrey", border = "white",
       breaks = breaks, xaxt = "n")
  axis(1, at = c(-pi, -pi/2, 0, pi/2, pi), labels = expression(-pi,
                                                               -pi/2, 0, pi/2, pi))
  for (i in 1:(nbStates + 1)) {
    lines(distData$angle$angle, distData$angle[, i + 1],
          col = lineCol[i], lty = lty[i], lwd = lwd[i])
  }
  legend("top", legText, lwd = lwd, col = lineCol, lty = lty,
         bty = "n")
  
  par(mfrow = mfrow_og)
  par(mar = mar_og)
}


# A simpler version of moveHMM:::plot.moveHMM() for tracks
plot_track_states <- function(fit, animals, states = NULL, col = NULL){
  mfrow_og <- par()$mfrow
  mar_og <- par()$mar
  
  nbStates <- ncol(fit$mle$stepPar)
  
  if (is.null(col) | (!is.null(col) & length(col) != nbStates)) {
    col <- .getPalette(nbStates = nbStates)
  }
  
  if(is.null(states)){
    if (nbStates > 1) {
      states <- viterbi(fit)
    } else {
      states <- rep(1, nrow(fit$data))
    }
  }
  
  nbAnimals <- length(unique(fit$data$ID))
  if (is.character(animals)) {
    if (any(!animals %in% unique(fit$data$ID))) {
      stop("Check 'animals' argument, ID not found")
    }
    animalsInd <- which(unique(fit$data$ID) %in% animals)
  }
  else if (is.numeric(animals)) {
    if (min(animals) < 1 | max(animals) > nbAnimals) {
      stop("Check 'animals' argument, index out of bounds")
    }
    animalsInd <- animals
  }
  else {
    animalsInd <- 1:nbAnimals
  }
  nbAnimals <- length(animalsInd)
  ID <- unique(fit$data$ID)[animalsInd]
  if (nbStates > 1) {
    par(mfrow = c(1, 1))
    par(mar = c(5, 4, 4, 2) - c(0, 0, 2, 1))
    for (zoo in 1:nbAnimals) {
      ind <- which(fit$data$ID == ID[zoo])
      s <- states[ind]
      x <- fit$data$x[ind]
      y <- fit$data$y[ind]
      if (!all(y == 0)) {
        plot(x, y, pch = 16, col = col[s], cex = 0.5, 
             asp = 1, xlab = "x", ylab = "y")
        segments(x0 = x[-length(x)], y0 = y[-length(y)], 
                 x1 = x[-1], y1 = y[-1], col = col[s[-length(s)]], 
                 lwd = 1.3)
      }
      else {
        plot(x, xlab = "time", ylab = "x", pch = 16, 
             cex = 0.5, col = col[s])
        segments(x0 = 1:(length(x) - 1), y0 = x[-length(x)], 
                 x1 = 2:length(x), y1 = x[-1], col = col[s[-length(x)]], 
                 lwd = 1.3)
      }
      mtext(paste("Animal ID:", ID[zoo]), side = 3, outer = TRUE, 
            padj = 2)
    }
  }
  par(mfrow = mfrow_og)
  par(mar = mar_og)
}

