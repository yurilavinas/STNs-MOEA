get_moon <- function() {
  # Creating Variable Bounds
  maximum = c(1, 1)
  minimum = c(0, 0)
  n.var = 2
  n.obj = 3
  n.const = 2
  
  
  evaluate_moon <- function(X) {
    write(
      X,
      file = "../evaluate_moon/pop_vars_eval.txt",
      ncolumns = n.var,
      sep = "\t"
    )
    system("./moon_mop ../evaluate_moon/")
    
    objectives <-
      scan(paste0("../evaluate_moon/pop_objs_eval.txt"),
           quiet = TRUE)
    objectives <- matrix(objectives, ncol = n.obj, byrow = TRUE)
    
    constraints <-
      scan(paste("../evaluate_moon/pop_cons_eval.txt", sep = "/"),
           quiet = TRUE)
    constraints <- matrix(constraints, ncol = n.const, byrow = TRUE)
    
    objectives[which(rowSums(constraints)>0),] = c(1e10,1e10,1e10)
    
    return(objectives)
  }
  
  constraints <- function(X) {
    nv <- n.var # number of variables
    # Prepare output matrix of constraint function values
    Cmatrix <- matrix(numeric(),
                      nrow = nrow(X), # TODO: check this
                      ncol = 2 * nv + 2)
    
    # Set informative column names (be nice to your users!)
    colnames(Cmatrix) <- c(paste0("x",
                                  rep(1:nv, times = nv),
                                  rep(c("min", "max"), each = nv)),
                           rep(c("g1"), each = nv))
    
    # Box limits of the feasible space
    Xmin <- matrix(minimum,
                   ncol = nv,
                   nrow = nrow(X),
                   byrow = TRUE)
    Xmax <- matrix(maximum,
                   ncol = nv,
                   nrow = nrow(X),
                   byrow = TRUE)
    
    # Calculate "x_i >= 0" and "x_i <= 1" constraints
    Cmatrix[, 1:nv]              <- Xmin - X
    Cmatrix[, (nv + 1):(2 * nv)] <- X - Xmax
    
    # g1 and h1 functions
    g1 <- function(X) {
      constraints <-
        scan(paste("../evaluate_moon/pop_cons_eval.txt", sep = "/"),
             quiet = TRUE)
      constraints <- matrix(constraints, ncol = n.const, byrow = TRUE)
      return(constraints)
    }
    
    # Calculate g1(x)
    Cmatrix[, (2 * nv + 1):(2 * nv + 2)] <- -g1(X)
    
    # Assemble matrix of *violations*
    Vmatrix <- Cmatrix
    Vmatrix[, 1:(2 * nv + 2)] <-
      pmax(Vmatrix[, 1:(2 * nv + 2)], 0)        # inequality constraints
    
    v = rowSums(Vmatrix)
    # Return necessary variables
    return(list(
      Cmatrix = Cmatrix,
      Vmatrix = Vmatrix,
      v       = v
    ))
  }
  
  out1 = list(
    evaluate = evaluate_moon,
    constraints = constraints,
    n.var = n.var,
    n.obj = n.obj,
    n.const = n.const,
    xl = minimum,
    xu = maximum
  )
  
  return(out1)
}
