get_mazda <- function() {
  n_var <- 222 # number of variables
  n_const = 498
  n.obj = 2
  # Creating Variable Bounds
  discrete = as.matrix(read.table(
    "~/Documents/estudos/MOEADr/inst/extdata/DiscreteValues.txt",
    col.names = paste0("V", seq_len(18)),
    sep = ",",
    fill = TRUE
  ))
  maximum = c(0, nrow = n_var)
  minimum = c(0, nrow = n_var)
  for (i in 1:n_var) {
    maximum[i] = max(discrete[i,], na.rm = TRUE)
    minimum[i] = min(discrete[i,], na.rm = TRUE)
  }
  
  # Function responsible to convert the continuous values to the nearest discrete values
  Discretize <- function(X) {
    nearest = apply(abs(discrete - X), 1, FUN = which.min)
    k = 0
    discreteValues = c(0)
    for (i in nearest) {
      k = k + 1
      discreteValues[k] = discrete[k, i]
    }
    return(discreteValues)
  }
  
  Evaluate <- function() {
    objectives <-
      scan(paste("../evaluate_mazda/pop_objs_eval.txt", sep = "/"),
           quiet = TRUE)
    objectives <- matrix(objectives, ncol = 5, byrow = TRUE)
    return(objectives[, 1:2])
  }
  
  # Definition of the problem
  evaluate_mazda <- function(X) {
    write(
      X,
      file = "../evaluate_mazda/pop_vars_eval.txt",
      ncolumns = n_var,
      sep = "\t"
    )
    system("./mazda_mop ../evaluate_mazda")
    
    objectives <-
      scan(paste("../evaluate_mazda/pop_objs_eval.txt", sep = "/"),
           quiet = TRUE)
    objectives <- matrix(objectives, ncol = 5, byrow = TRUE)
    objectives = objectives[, 1:2]
    constraints <-
      scan(paste("../evaluate_mazda/pop_cons_eval.txt", sep = "/"),
           quiet = TRUE)
    constraints <- matrix(constraints, ncol = 54, byrow = TRUE)
    
    
    objectives[, 1] = objectives[, 1] / 74
    objectives[, 2] = objectives[, 2] - 2
    
    objectives[which(rowSums(constraints)>0),] = c(1e10,1e10)
    
    return(objectives)
  }
  
  constraints <- function(X) {
    # Prepare output matrix of constraint function values
    Cmatrix <- matrix(numeric(),
                      nrow = nrow(X), # TODO: check this
                      ncol = n_const) # 296 box constraints and 36 inequality constraints
    
    # Set informative column names (be nice to your users!)
    colnames(Cmatrix) <- c(paste0("x",
                                  rep(1:n_var, times = 2),
                                  rep(c("min", "max"), each = n_var)),
                           rep(c("g1"), each = 54))
    
    # Box limits of the feasible space
    Xmin <-
      matrix(minimum,
             ncol = n_var,
             nrow = nrow(X),
             byrow = TRUE)
    Xmax <-
      matrix(maximum,
             ncol = n_var,
             nrow = nrow(X),
             byrow = TRUE)
    
    # Calculate "x_i >= 0" and "x_i <= 1" constraints
    Cmatrix[, 1:n_var]              <- Xmin - X
    Cmatrix[, (n_var + 1):(2 * n_var)] <- X - Xmax
    
    # g1 and h1 functions
    g1 <- function(X) {
      constraints <-
        scan(paste("../evaluate_mazda/pop_cons_eval.txt", sep = "/"),
             quiet = TRUE)
      constraints <- matrix(constraints, ncol = 54, byrow = TRUE)
      return(constraints)
    }
    
    # Calculate g1(x)
    Cmatrix[, (2 * n_var + 1):(2 * n_var + 54)] <- -g1(X)
    
    # Assemble matrix of *violations*
    Vmatrix <- Cmatrix
    Vmatrix[, 1:(2 * n_var + 54)] <-
      pmax(Vmatrix[, 1:(2 * n_var + 54)], 0)
    
    v = rowSums(Vmatrix)
    
    # Scaling the Penalties
    if (is.null(parent.frame(2)$iter)) {
      #First generation (No incumbent solutions)
      v[which(v != 0)] = (v[which(v != 0)] - min(v)) / (max(v) - min(v) + 1e-16) + 0.000001
    }
    else{
      # Getting the incumbent solutions
      e = parent.frame(2)
      Vtmatrix = e$Vt$Vmatrix
      vt = rowSums(Vtmatrix)
      
      # Extract max and min for scaling
      max = max(v, vt)
      min = min(v, vt)
      
      # Updating the new scaled penalties of the incumbent solutions
      e$Vt$v[which(vt != 0)] = (vt[which(vt != 0)] - min) / (max - min + 1e-16) + 0.000001
      
      # Scaling the new solution's penalties
      v[which(v != 0)] = (v[which(v != 0)] - min) / (max - min + 1e-16) + 0.000001
    }
    
    # Return necessary variables
    return(list(
      Cmatrix = Cmatrix,
      Vmatrix = Vmatrix,
      v       = v
    ))
  }
  
  out1 = list(
    evaluate = evaluate_mazda,
    constraints = constraints,
    n.var = n_var,
    n.obj = n.obj,
    n.const = n_const,
    xl = minimum,
    xu = maximum
  )
  return(out1)
}
