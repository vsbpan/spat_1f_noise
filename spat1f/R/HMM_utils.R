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
viterbi <- function(model, newdata = NULL){
  moveHMM::viterbi(model, newdata = newdata)
}

# Wrapper for moveHMM::stateProbs()
stateProbs <- function(model){
  moveHMM::stateProbs(model)
}

# Wrapper for moveHMM::fitHMM()
fit_HMM <- function(data, 
                    stepPar0 = c(2.5, 1.5, 1.5, 0.8), anglePar0 = c(0, pi, 0.1, 0.1), 
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



