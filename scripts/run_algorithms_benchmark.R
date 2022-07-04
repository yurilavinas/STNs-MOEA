rm(list = ls(all = TRUE))

library(nsga2R)
library(MOEADr)
library(MOEADps)
library(smoof)
library(reticulate)

inpath = "../"

source(paste0(inpath, 'algorithms/nsga2.R'))
source(paste0(inpath, 'algorithms/moead.R'))
source(paste0(inpath, 'scripts/utils.R'))
source(paste0(inpath, 'scripts/functions.R'))



source_python(paste0(inpath, 'Python/reproblem.py'))

maxeval <- 30000
popSize <- 100
dimensions <- 10
repetition <- 10

fun_benchmarks <- paste0("UF", 1:10)

for (f in fun_benchmarks) {
  print(f)
  
  num = as.integer(strsplit(f, "UF")[[1]][2])
  if (num <= 7) {
    n.obj = 2
  }
  else{
    n.obj = 3
  }
  
  problem.smoof.UF <-
    makeUFFunction(dimensions = dimensions, id = as.numeric(strsplit(f, split = "[A-Z]")[[1]][3]))
  problem.sr <- function(X) {
    t(apply(X, MARGIN = 1,
            FUN = problem.smoof.UF))
  }
  
  par.set = ParamHelpers::getParamSet(problem.smoof.UF)
  
  ## Set the input parameters for the moead() routine
  ## This reproduces the Original MOEA/D of Zhang and Li (2007)
  ## (with a few changes in the computational budget, to make it run faster)
  problem   <- list(
    name       = "problem.sr",
    xmin       = as.numeric(ParamHelpers::getLower(par.set)),
    xmax       = as.numeric(ParamHelpers::getUpper(par.set)),
    m          = n.obj
  )
  
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
      list(name       = "sld", H = 21) 
    
    W2  <- generate_weights(decomp = decomp2,
                            m      = n.obj)
    
    X  <- create_population(N       = nrow(W2),
                            problem = problem)
    popSize = nrow(W2)
  }
  
  scaling <- preset_moead("moead.de")$scaling
  scaling$name <- "simple"
  
  for (run in 1:repetition) {
    print("MOEA/D")
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
    #   saving.dir = dir.name,
    #   decomp = decomp2,
    #   scaling = scaling,
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
    out <-
      nsga2(
        X = X,
        problem = problem,
        varNo = dimensions,
        objDim = n.obj,
        maxeval = maxeval,
        lowerBounds = as.numeric(ParamHelpers::getLower(par.set)),
        upperBounds = as.numeric(ParamHelpers::getUpper(par.set)),
        popSize = popSize,
        tourSize = 2,
        cprob = 0.9,
        XoverDistIdx = 20,
        mprob = 0.1,
        MuDistIdx = 3,
        scaling = T,
        saving.dir = dir.name
      )
    
  }
}
