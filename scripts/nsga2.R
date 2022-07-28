variation_de <- function(X, P, phi = 0.5, ...) {
  phi <- 0.5
  new.solution <- X
  dimX <- dim(X)[1]
  for (i in 1:dim(X)[1]) {
    idx <- sample.int(dimX, 3,
                      replace = TRUE)#,
    # prob    = P[, i])
    new.solution[i, ] <-
      X[idx[1], ] + phi * (X[idx[2], ] - X[idx[3], ])
  }
  return (new.solution)
}

nsgaps <-
  function(X,
           problem,
           varNo,
           objDim,
           lowerBounds = rep(-Inf, varNo),
           upperBounds = rep(Inf, varNo),
           popSize = 100,
           tourSize = 2,
           maxevals = 10000,
           cprob = 0.7,
           XoverDistIdx = 5,
           mprob = 0.2,
           MuDistIdx = 10,
           saving.dir = NULL,
           ...) {
    # initializing the population
    
    Y <- evaluate_population(X       = X,
                             problem = problem,
                             nfe     = 0)$Y
    # Y   <- YV$Y
    
    parent <- cbind(X, Y)
    
    # ranking the initial population
    ranking <-
      fastNonDominatedSorting(parent[, (varNo + 1):(varNo + objDim)])
    
    # Rank index for each chromosome
    rnkIndex <- integer(popSize)
    
    i <- 1
    
    while (i <= length(ranking)) {
      rnkIndex[ranking[[i]]] <- i
      
      i <- i + 1
      
    }
    parent <- cbind(parent, rnkIndex)
    
    # crowding distance calculation
    objRange <-
      apply(parent[, (varNo + 1):(varNo + objDim)], 2, max) - apply(parent[, (varNo +
                                                                                1):(varNo + objDim)], 2, min)
    
    cd <- crowdingDist4frnt(parent, ranking, objRange)
    
    parent <- cbind(parent, apply(cd, 1, sum))
    
    nfe <- dim(parent)[1]
    iter <- 0
    if(!is.null(saving.dir)){
      write.table(data.frame(X = X, Y = Y, iter = iter, nfe = nfe, run = run), paste0(saving.dir, "/all_solutions.csv"), append = T, col.names = F, row.names = F, sep =",")
    }
    # while (iter < maxiter) {
    while (nfe < maxevals) {
      
      # main start: RA
      indexes <-
        1:dim(parent)[1] # right now all solutions are selected
      if (resource.allocation$name != "none") {
        epsilon <- 1e-50
        ra <- runif(indexes)
        indexes <-
          sample(
            x = 1:length(ra),
            size = resource.allocation$n,
            prob = ra + epsilon
          )
        indexes <- sort(indexes)
      }
      popSize <- length(indexes)
      
      temp.parent <- parent
      parent <- parent[indexes,]
      # main end: RA
      
      # tournament selection
      matingPool <- tournamentSelection(parent, popSize, tourSize)
      
      # crossover operator
      # childAfterX <-
      #   boundedSBXover(matingPool[, 1:varNo], lowerBounds, upperBounds, cprob, XoverDistIdx)
      # Only design parameters are input as the first argument
      childAfterX <- variation_de(matingPool[,1:varNo], matrix(1, nrow(matingPool[,1:varNo]), ncol(matingPool[,1:varNo])))
      childAfterX <- t(childAfterX)
      childAfterX <- t(matrix(pmax(0, pmin(childAfterX, 1)),
                              nrow  = nrow(childAfterX),
                              byrow = FALSE))
      # mutation operator
      childAfterM <-
        boundedPolyMutation(childAfterX, lowerBounds, upperBounds, mprob, MuDistIdx)
      
      # evaluate the objective fns of childAfterM
      # childAfterM <- cbind(childAfterM, t(apply(childAfterM, 1, fn)))
      Y <- evaluate_population(X       = childAfterM,
                               problem = problem,
                               nfe     = 0)$Y
      childAfterM <- cbind(childAfterM, Y)
      
      
      nfe <- nfe + dim(childAfterM)[1]
      
      
      # RA pop adjustments - start
      parent <- temp.parent
      temp.parent[indexes, 1:(varNo + objDim)] <- childAfterM
      childAfterM <- temp.parent[, 1:(varNo + objDim)]
      popSize <- dim(childAfterM)[1]
      # RA pop adjustments - end
      
      # "Rt = Pt + Qt"
      # Combine the parent with the childAfterM (No need to retain the rnkIndex and cd of parent)
      parentNext <- rbind(parent[, 1:(varNo + objDim)], childAfterM)
      # ranking again
      ranking <-
        fastNonDominatedSorting(parentNext[, (varNo + 1):(varNo + objDim)])
      
      i <- 1
      
      while (i <= length(ranking)) {
        rnkIndex[ranking[[i]]] <- i
        
        i <- i + 1
        
      }
      parentNext <- cbind(parentNext, rnkIndex)
      
      # crowded comparison again
      objRange <-
        apply(parentNext[, (varNo + 1):(varNo + objDim)], 2, max) - apply(parentNext[, (varNo +
                                                                                          1):(varNo + objDim)], 2, min)
      
      cd <- crowdingDist4frnt(parentNext, ranking, objRange)
      
      parentNext <- cbind(parentNext, apply(cd, 1, sum))
      
      parentNext.sort <-
        parentNext[order(parentNext[, varNo + objDim + 1], -parentNext[, varNo +
                                                                         objDim + 2]), ]
      
      # choose the first 'popSize' rows for next generation
      parent <- parentNext.sort[1:popSize, ]
      # cat("Gen: ",iter)
      iter <- iter + 1
      if(!is.null(saving.dir)){
        write.table(data.frame(X = X, Y = Y, iter = iter, nfe = nfe, run = run), paste0(saving.dir, "/all_solutions.csv"), append = T, col.names = F, row.names = F, sep =",")
      }
    }
    # report on nsga2 settings and results
    result = list(
      # functions = fn,
      # parameterDim = varNo,
      # objectiveDim = objDim,
      # lowerBounds = lowerBounds,
      # upperBounds = upperBounds,
      # popSize = popSize,
      # tournamentSize = tourSize,
      # generations = iter,
      # XoverProb = cprob,
      # XoverDistIndex = XoverDistIdx,
      # mutationProb = mprob,
      # mutationDistIndex = MuDistIdx,
      iter = iter,
      nfe = nfe,
      parameters = parent[, 1:varNo],
      objectives = parent[, (varNo + 1):(varNo + objDim)]#,
      # paretoFrontRank = parent[, varNo + objDim + 1],
      # crowdingDistance = parent[, varNo + objDim + 2]
    )
    
    class(result) = "nsga2R"
    
    return(result)
  }

