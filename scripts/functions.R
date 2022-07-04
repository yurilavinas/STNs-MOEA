library(assertthat)


# concatenate the decision space values of each solutions
concat.solutions <- function(all.data) {
  assert_that(is.data.frame(all.data))
  origin.solutions <-
    apply(
      X = all.data,
      MARGIN = 1,
      FUN = paste0,
      collapse = ""
    )
  origin.solutions
}

get.objectives.values <- function(all.data) {
  assert_that(is.data.frame(all.data))
  Y = select(all.data, starts_with('Y'))
  Y
}

get.decision.values <- function(all.data) {
  assert_that(is.data.frame(all.data))
  X = select(all.data, starts_with('X'))
  X
}

## function from MOEADr
### used in uniform.decomposition
is_coprime <- function(x, y) {
  a <- x
  b <- y
  while (b != 0) {
    if (a == 1 || b == 1)
      return(TRUE)
    t <- b
    b <- a %% b
    a <- t
  }
  return(FALSE)
}

## function from MOEADr
uniform.decomposition <- function (N, n.obj) {
  assert_that(is.numeric(N))
  assert_that(is.numeric(n.obj))
  assert_that(0 < n.obj)
  assert_that(0 < N)
  
  div <- which(sapply(seq_len(N), is_coprime, N))
  H <- t(utils::combn(x = div, m = n.obj - 1))
  construct_un <- function(h, N) {
    U <- t(sapply(seq_len(N), function(x) {
      (h * x) %% N
    }))
    U <- U + N * (1 - sign(U))
  }
  cd2 <- function(U) {
    magic <- (13 / 12) ^ ncol(U)
    S1 <- (2 / nrow(U)) * sum(apply((1 + (
      abs(U - 0.5) - abs(U -
                           0.5) ^ 2
    ) / 2), 1, prod))
    S2 <- (1 / nrow(U) ^ 2) * sum(sapply(1:nrow(U), function(y) {
      tU <- matrix(U[y,],
                   nrow = nrow(U),
                   ncol = ncol(U),
                   byrow = TRUE)
      apply(1 + (abs(U - 0.5) + abs(tU - 0.5)) / 2 - abs(U -
                                                           tU) / 2,
            MARGIN = 1,
            FUN = prod)
    }))
    return(magic - S1 + S2)
  }
  min_h <- H[which.min(apply(
    H,
    MARGIN = 1,
    FUN = function(x) {
      cd2(construct_un(x, N))
    }
  )),]
  Un <- (construct_un(min_h, N) - 0.5) / N
  if (nrow(Un) == 1) {
    Un <- t(Un)
  }
  U_pow <- t(t(Un) ^ sapply(seq_len(n.obj - 1), function(x) {
    (n.obj - x) ^ -1
  }))
  pow_prod <- t(apply(
    U_pow,
    MARGIN = 1,
    FUN = function(x) {
      sapply(
        seq_len(length(x)),
        FUN = function(y) {
          prod(x[seq_len(y - 1)])
        }
      )
    }
  ))
  if (nrow(pow_prod) == 1) {
    pow_prod <- t(pow_prod)
  }
  W <- (1 - U_pow) * pow_prod
  colnames(W) <- paste("Var", 1:ncol(W), sep = "")
  return(cbind(W, VarLast = apply(
    U_pow, MARGIN = 1, FUN = prod
  )))
}

# Get estimate of 'ideal' point (minP)
getminP <- function(X){
  apply(X, MARGIN = 2, FUN = min, na.rm = TRUE)
}



# Get estimate of 'nadir' point (maxP)
getmaxP <- function(X){
  apply(X, MARGIN = 2, FUN = max, na.rm = TRUE)
}

scale_vector <- function(Y, minP=NULL, maxP=NULL) {
  
  if(is.null(minP) && is.null(maxP)){
    minP <- getminP(Y)  
    maxP <- getmaxP(Y)
  }
  
  
  
  MinP <- matrix(rep(minP, times = nrow(Y)),
                 nrow  = nrow(Y),
                 byrow = TRUE)
  MaxP <- matrix(rep(maxP, times = nrow(Y)),
                 nrow  = nrow(Y),
                 byrow = TRUE)
  
  # Perform linear scaling
  Y    <- (Y - MinP) / (MaxP - MinP + 1e-16)
  return (Y)
  
  return(Y)
}

find.representative <- function(filtered.data, W, minP, maxP) {
  assert_that(is.data.frame(filtered.data))
  assert_that(is.matrix(W))
  assert_that(5 == dim(W)[1])
  assert_that(is.vector(minP))
  # getting the objective values from the whole data
  Y = get.objectives.values(filtered.data)
  
  #scaling? confirm this
  Y = scale_vector(as.data.frame(Y), minP, maxP)
  
  selected.Y = data.frame()
  
  # this is a reference value need for the scaling of the Weight Tchebycheff
  big.minP = matrix(minP,
                    nrow  = nrow(Y),
                    ncol  = ncol(Y),
                    byrow = TRUE)
  
  # calculating the scalar values using the Weight Tchebycheff
  for (i in 1:dim(W)[1]) {
    scalar_values <- apply(t(W[i, ] * t(Y - big.minP + 1e-16)),
                           MARGIN = 1,
                           FUN    = max)
    
    # finding the solution that has the minimal scalar value
    idx = which.min(scalar_values)
    
    idx = idx[sample.int(length(idx))[1]]
    
    # combine into one data frame, with the minimal scalar value solution
    # one for each of the weight vectors
    selected.Y <- rbind(selected.Y, filtered.data[idx,])
  }
  return(selected.Y)
}


# finding the representative solutions for each of the vectors
generate.vector.data.combinatorial <- function(all.data, minP) {
  assert_that(is.data.frame(all.data))
  assert_that(is.vector(minP))
  max.iter = max(all.data$iter) - 1
  
  final.output = data.frame()
  for (j in 1:max(all.data$run)) {
    algorithm.data  = filter(all.data, run == j)
    for (i in 0:max.iter) {
      #################################
      #### getting starting points ####
      #################################
      
      #filtering data of the main iteration - where the search 'come from'
      filtered.data1 = filter(algorithm.data, iter == i)
      
      # finding the representatives solutions
      representatives1 <-
        find.representative(filtered.data1, W, minP)
      
      # getting the decision values from the representative solutions
      X1 = get.decision.values(representatives1)
      
      # concatenating the dimensions of the solutions from the whole data
      solution1 = concat.solutions(X1)
      
      # output of the iteration/generation - 1st step
      output = round(1.0 - get.objectives.values(representatives1), 6)
      output$Solution1 = solution1
      
      ###############################
      #### getting ending points ####
      ###############################
      
      #filtering data of the next iteration - where the search 'go to'
      filtered.data2 = filter(algorithm.data, iter == i + 1)
      
      # finding the representatives solutions
      representatives2 <-
        find.representative(filtered.data2, W, minP)
      
      # getting the decision values from the representative solutions
      X2 = get.decision.values(representatives2)
      
      # concatenating the dimensions of the solutions from the whole data
      solution2 = concat.solutions(X2)
      
      # output of the iteration/generation - 2nd step
      output$Solution2 = solution2
      
      #################################################
      #### adding iteration and vector information ####
      #################################################
      
      
      # adding information about the run of the algorithm
      output$run = representatives1$Run
      
      # adding generation (iteration) data
      output$Gen = i
      # adding vector column given the number of Weight vectors
      output$Vector = paste0("V", 1:dim(W)[1])
      
      
      # combining the data generated within the code above
      final.output = rbind(final.output, output)
    }
  }
  final.output
}

# load data files of one algorithm on UF/RE* problems
## add filename
read.data.continuous <- function(inpath, algorithm, problem) {
  assert_that(is.character(algorithm))
  assert_that(is.character(problem))
  all.data = data.frame()
  
  filename = paste0(inpath,
                    algorithm,
                    '_',
                    problem,
                    '/all_solutions.csv')
  
  all.data = read.table(
    file = filename,
    header = T,
    sep = ','
  )
  
  all.data
}

# finding the representative solutions for each of the vectors
## this is different because we need to partition the continuous search space
generate.vector.data.continuous <- function(all.data, minP, maxP, d = 3, width = 8) {
  assert_that(is.data.frame(all.data))
  assert_that(is.vector(minP))
  max.iter = max(all.data$iter) - 1 
  
  final.output = data.frame()
  for (j in 1:max(all.data$run)) {
    algorithm.data  = filter(all.data, run == j)
    for (i in 0:max.iter) {
      #################################
      #### getting starting points ####
      #################################
      
      #filtering data of the main iteration - where the search 'come from'
      filtered.data1 = filter(algorithm.data, iter == i)
      
      # finding the representatives solutions
      representatives1 <-
        find.representative(filtered.data1, W, minP, maxP)
      
      
      # getting the decision values from the representative solutions
      X1 = get.decision.values(representatives1)
      
      # concatenating the dimensions of the solutions from the whole data
      # creating the hyperspace partition
      solution1 = format(as.hexmode(apply(round(X1, d) * 10^d, 2, as.integer)), width=width)
      solution1 = apply(solution1, 1, paste0, collapse = "")

      
      output = get.objectives.values(representatives1)
      output$Solution1 = solution1
      
      ###############################
      #### getting ending points ####
      ###############################
      
      #filtering data of the next iteration - where the search 'go to'
      filtered.data2 = filter(algorithm.data, iter == i + 1)
      
      # finding the representatives solutions
      representatives2 <-
        find.representative(filtered.data2, W, minP, maxP)
      
      
      # getting the decision values from the representative solutions
      X2 = get.decision.values(representatives2)
      
      # concatenating the dimensions of the solutions from the whole data
      # creating the hyperspace partition
      solution2 = format(as.hexmode(apply(round(X2, d) * 10^d, 2 ,as.integer)), width=width)
      solution2 = apply(solution2, 1, paste0, collapse = "")
      
      # output of the iteration/generation - 2nd step
      output$Solution2 = solution2
      
      #################################################
      #### adding iteration and vector information ####
      #################################################
      # adding information about the run of the algorithm
      output$Run = representatives1$run
      
      # adding generation (iteration) data
      output$Gen = i
      # adding vector column given the number of Weight vectors
      output$Vector = paste0("V", 1:dim(W)[1])
      
      
      # combining the data generated within the code above
      final.output = rbind(final.output, output)
    }
  }
  final.output
}