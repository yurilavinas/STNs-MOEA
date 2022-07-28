rm(list = ls(all = TRUE))

options(scipen=999)

library(dplyr)
library(feather)

source("functions.R")

#parameters
inpath = '../algorithm_data/'
outpath = '../STN_data/'

od = 3
width = 12
# problem
fun_UF <-
  paste0("UF", 1:10)
fun_rw <- c('RE21','RE22','RE23','RE24','RE25','RE31','RE32','RE33','RE34','RE35','RE36','RE37')

problems = c(fun_UF, fun_rw)

for (prob in problems){
  cat(prob, ': ')
  if (prob == 'UF10' || prob == 'UF8' || prob == 'UF9' || 
      prob == 'RE31'|| prob == 'RE32'|| prob == 'RE33'|| prob == 'RE34'|| prob == 'RE35'|| prob == 'RE36'|| prob == 'RE37'|| prob == 'RE38'
      ){
    n.obj = 3
  }
  else{
    n.obj = 2
  }
  
  # getting all data
  ## 1st for moead
  moead.data = read.data.continuous(
    inpath = inpath,
    algorithm = 'moead',
    problem = prob
  )
  
  # then for nsga2
  nsga2.data = read.data.continuous(
    inpath = inpath,
    algorithm = 'nsga2',
    problem = prob
  )

  # generating the 5 weight vectors for a 2, 3 or more obj problem
  # might not work for all number of objectives
  W = uniform.decomposition(N = 5, n.obj = n.obj)
  
  # getting the objective values from the representative solutions
  obj.values.moead = get.objectives.values(moead.data)
  obj.values.nsga2 = get.objectives.values(nsga2.data)
  
  # This is a reference point for the calculation of the
  ## Weighted Tchebycheff method.
  ### should be combined w/ data of all algorithms!
  minP.moead =
    getminP(obj.values.moead)
  minP.nsga2 =
    getminP(obj.values.nsga2)
  combined.minP = data.frame(rbind(minP.moead, minP.nsga2))
  
  minP = getminP(combined.minP)
  
  maxP.moead =
    getmaxP(obj.values.moead)
  maxP.nsga2 =
    getmaxP(obj.values.nsga2)
  combined.maxP = data.frame(rbind(maxP.moead, maxP.nsga2))
  
  maxP = getminP(combined.maxP)
  
  # finding the representative solutions for each of the vectors
  cat("MOEA/D... ")
  moead.data.vectors = generate.vector.data.continuous(moead.data, minP, maxP, od, width)
  cat("NSGA-II... ")
  nsga2.data.vectors  = generate.vector.data.continuous(nsga2.data, minP, maxP, od, width)
  
  write.table(
    moead.data.vectors,
    file = paste0(
      outpath,
      'moead_',
      prob,
      '.txt'
    ),
    quote = FALSE,
    row.names = FALSE
  )
  
  write.table(
    nsga2.data.vectors,
    file = paste0(
      outpath,'nsga2_',
      prob,
      '.txt'
    ),
    quote = FALSE,
    row.names = FALSE
  )
  
  print("done!")
  
}
