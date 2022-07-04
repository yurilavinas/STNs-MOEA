rm(list = ls(all = TRUE))

library(nsga2R)
library(MOEADr)
library(MOEADps)
library(reticulate)

inpath = "../"

source(paste0(inpath, 'algorithms/nsga2.R'))
source(paste0(inpath, 'algorithms/moead.R'))
source(paste0(inpath, 'scripts/utils.R'))
source(paste0(inpath, 'scripts/functions.R'))



source(paste0("MOON.R"))
source(paste0("MAZDA.R"))

source_python(paste0(inpath, 'Python/reproblem.py'))

maxeval <- 30000
popSize <- 100

repetition <- 10


fun_rw <-
  c('RE21','RE22','RE23','RE24','RE25','RE31','RE32','RE33','RE34','RE35','RE36','RE37')


for (f in fun_rw) {
  print(f)
  
  benchmark = get_reproblem(f)
  benchmark$name = f
  benchmark$n.obj = benchmark$n_objectives
  benchmark$xl = benchmark$lbound
  benchmark$xu = benchmark$ubound
  problem.1 <- function(X) {
    out = c()
    for (i in 1:dim(X)[1]) {
      out = rbind(out, unlist(benchmark$evaluate(X[i, ]))[1:benchmark$n.obj])
    }
    return(out)
  }
  
  problem <-
    list(
      name       = problem.1,
      # function that executes the MOP
      xmin       = benchmark$xl,
      # minimum parameter value for each dimension
      xmax       = benchmark$xu,
      # maximum parameter value for each dimension
      m          = benchmark$n.obj
    )

  
  n.obj = benchmark$n.obj
  
  if (n.obj == 2) {
    decomp2    <-
      list(name       = "sld", H = popSize - 1) # <-- H = 99 in the original
    
    W2  <- generate_weights(decomp = decomp2,
                            m      = n.obj)
    
    X  <- create_population(N       = nrow(W2),
                            problem = problem)
  }
  else{
    decomp2    <-
      list(name       = "sld", H = 21) # <-- H = 99 in the original
    
    W2  <- generate_weights(decomp = decomp2,
                            m      = n.obj)
    
    X  <- create_population(N       = nrow(W2),
                            problem = problem)
    popSize = nrow(W2)
  }
  
  
  scaling <- preset_moead("moead.de")$scaling
  scaling$name <- "simple"
  
  for (run in 1:repetition) {
    # print("MOEA/D")

    # dir.name <- paste0(inpath, "algorithm_data/moead_", f, "/")
    # if (!dir.exists(dir.name)) {
    #   dir.create(dir.name)
    #   if (!file.exists(paste0(dir.name, "/all_solutions.csv"))) {
    #     colname = c(paste0("X", 1:length(problem$xmin)),
    #                 paste0('Y', 1:problem$m),
    #                 'iter',
    #                 'nfe',
    #                 'run')
    #     colname = t(as.data.frame(colname))
    #     write.table(
    #       colname,
    #       paste0(dir.name, "/all_solutions.csv"),
    #       sep = ',',
    #       col.names = F,
    #       row.names = F
    #     )
    #   }
    # }
    # 
    # 
    # out1 <- moead(
    #   X = X,
    #   preset   = preset_moead("moead.de"),
    #   problem  = problem,
    #   decomp = decomp2,
    #   scaling = scaling,
    #   saving.dir = dir.name,
    #   showpars = list(show.iters = "dots", showevery = 1000),
    #   stopcrit = list(list(name    = "maxeval",
    #                        maxeval = maxeval))
    # )

    dir.name <- paste0(inpath, "algorithm_data/nsga2_", f, "/")
    if (!dir.exists(dir.name)) {
      dir.create(dir.name)
      if (!file.exists(paste0(dir.name, "/all_solutions.csv"))) {
        colname = c(paste0("X", 1:length(problem$xmin)),
                    paste0('Y', 1:problem$m),
                    'iter',
                    'nfe',
                    'run')
        colname = t(as.data.frame(colname))
        write.table(
          colname,
          paste0(dir.name, "/all_solutions.csv"),
          sep = ',',
          col.names = F,
          row.names = F
        )
      }
    }
    
    print("NSGA-II")
    out2 <-
      nsga2(
        problem = problem,
        varNo = length(problem$xmin),
        objDim = n.obj,
        X = X,
        lowerBounds = problem$xmin,
        upperBounds = problem$xmax,
        popSize = popSize,
        tourSize = 2,
        maxeval = maxeval,
        cprob = 0.9,
        XoverDistIdx = 20,
        mprob = 0.1,
        MuDistIdx = 3,
        scaling = T,
        saving.dir = dir.name
      )
    
    
  }
}
