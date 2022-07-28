boundedSBXover = function (parent_chromosome, lowerBounds, upperBounds, cprob, 
                           mu) 
{# fix this using the pymoo code: https://github.com/anyoptimization/pymoo/blob/7b719c330ff22c10980a21a272b3a047419279c8/pymoo/operators/crossover/sbx.py#L7
  popSize = nrow(parent_chromosome)
  varNo = ncol(parent_chromosome)
  child <- parent_chromosome
  p <- 1
  
  y2 <- min(child)
  y1 <- max(child)
  
  for (i in 1:(popSize/2)) {
    if (runif(1) < cprob) {
      for (j in 1:varNo) {
        
        child1 =  child[p, j]
        child2 =  child[p + 1, j]
        
        yl <- lowerBounds[j]
        yu <- upperBounds[j]
        
        rnd = runif(1)
        if (rnd <= 0.5) {
          #### calculating for child1
          delta = (y2 - y1)
          
          if (delta < 1.0e-10){
            delta = 1.0e-10
          }
          beta = 1.0 + (2.0 * (y1 - yl) / delta)
          alpha = 2 - (beta^(-(1 + mu)))
          
          rnd = runif(1)
          if (rnd <= 1/alpha) {
            alpha = alpha * rnd
            betaq = alpha^(1/(1 + mu))
          }
          else {
            alpha = 1/(2 - alpha * rnd)
            betaq = alpha^(1/(1 + mu))
          }
          child1 = 0.5 * ((y1 + y2) - betaq * (y2 - y1))
          #### calculating for child2
          beta = 1.0 + (2.0 * (yu - y2) / delta)
          alpha = 2 - (beta^(-(1 + mu)))
          rnd = runif(1)
          if (rnd <= 1/alpha) {
            alpha = alpha * rnd
            betaq = alpha^(1/(1 + mu))
          }
          else { 
            alpha = 1/(2 - alpha * rnd)
            betaq = alpha^(1/(1 + mu))
          }
          child2 = 0.5 * ((y1 + y2) + betaq * (y2 -  y1))
          #### end child2
          if (child1 > yu) {
            child1 = yu
          }
          else if (child1 < yl) {
            child1 = yl
          }
          if (child2 > yu) {
            child2 = yu
          }
          else if (child2 < yl) {
            child2 = yl
          }
        }
        child[p, j] <- child1
        child[p + 1, j] <- child2
      }
    }
    p <- p + 2
  }
  return(child)
}


function (inputData) 
{
  popSize = nrow(inputData)
  idxDominators = vector("list", popSize)
  idxDominatees = vector("list", popSize)
  for (i in 1:(popSize - 1)) {
    for (j in i:popSize) {
      if (i != j) {
        xi = inputData[i, ]
        xj = inputData[j, ]
        if (all(xi <= xj) && any(xi < xj)) {
          idxDominators[[j]] = c(idxDominators[[j]], 
                                 i)
          idxDominatees[[i]] = c(idxDominatees[[i]], 
                                 j)
        }
        else if (all(xj <= xi) && any(xj < xi)) {
          idxDominators[[i]] = c(idxDominators[[i]], 
                                 j)
          idxDominatees[[j]] = c(idxDominatees[[j]], 
                                 i)
        }
      }
    }
  }
  noDominators <- lapply(idxDominators, length)
  rnkList <- list()
  rnkList <- c(rnkList, list(which(noDominators == 0)))
  solAssigned <- c()
  solAssigned <- c(solAssigned, length(which(noDominators == 
                                               0)))
  while (sum(solAssigned) < popSize) {
    Q <- c()
    noSolInCurrFrnt <- solAssigned[length(solAssigned)]
    for (i in 1:noSolInCurrFrnt) {
      solIdx <- rnkList[[length(rnkList)]][i]
      hisDominatees <- idxDominatees[[solIdx]]
      for (i in hisDominatees) {
        noDominators[[i]] <- noDominators[[i]] - 1
        if (noDominators[[i]] == 0) {
          Q <- c(Q, i)
        }
      }
    }
    rnkList <- c(rnkList, list(sort(Q)))
    solAssigned <- c(solAssigned, length(Q))
  }
  return(rnkList)
}

nsga2 <-
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
    while (nfe < maxevals) {
      # tournament selection
      matingPool <- tournamentSelection(parent, popSize, tourSize)
      
      # crossover operator
      childAfterX <-
        boundedSBXover(matingPool[, 1:varNo], lowerBounds, upperBounds, cprob, XoverDistIdx)
      # Only design parameters are input as the first argument
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
                               nfe     = nfe)$Y
      childAfterM <- cbind(childAfterM, Y)
      
      
      nfe <- nfe + dim(childAfterM)[1]
      
    
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
        X = parent[, 1:varNo]
        nd = find_nondominated_points(Y)
        Y = matrix(Y[nd,], ncol = objDim)
        data = data.frame(X = X[nd,], Y = Y, iter = iter, nfe = nfe, run = run)
        write.table(data, paste0(saving.dir, "/all_solutions.csv"), append = T, col.names = F, row.names = F, sep =",")
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

