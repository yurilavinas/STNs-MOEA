#########################################################################
# SNCS 
# Yuri Lavinas, Gabriela Ochoa and Claus Aranha, June 2022
# STN construction - Version 2
# Input:  Text file with trajectory data of several runs
# Output: STN graph objectseaf
#########################################################################
rm(list = ls(all = TRUE))

options(scipen=999)

library(dplyr)
library(ggplot2)
library(eaf)
library(MOEADr)


source("functions.R")

#parameters
number.repetitions = 10
inpath = '../algorithm_data/'
outpath = '../overall_results/'

precision_factor = 3
# problem
# fun_UF <-
#   paste0("UF", 1:10)
# fun_rw <- c('RE21','RE22','RE23','RE24','RE25','RE31','RE32','RE33','RE34','RE35','RE36','RE38')
fun_rw <- c('RE23')
# problems = c(fun_UF, fun_rw)
problems = c(fun_rw)

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
  
  # getting the theoretical pf data
  pfname <- paste0("../pf/", prob, "_pf.txt")
  pf <- read.table(pfname, colClasses=c("double", "double"))  # Read the Pareto Front file
  colnames(pf) = paste0("Y",1:n.obj)
  
  # getting all data
  ## 1st for moead
  moead.data = read.data.continuous(
    inpath = inpath,
    algorithm = 'moead',
    problem = prob
  )
  
  ## then for nsga2
  nsga2.data = read.data.continuous(
    inpath = inpath,
    algorithm = 'nsga2',
    problem = prob)
  
  moead.data = moead.data %>% filter(iter < 150) 
  nsga2.data = nsga2.data %>% filter(iter < 150) 
  
  ## generate hv convergence
  cat('hv convergence on going...')
  maxs.moead = getmaxP(get.objectives.values(moead.data))
  maxs.nsga2 = getmaxP(get.objectives.values(nsga2.data))
  combined = rbind(maxs.moead, maxs.nsga2)
  ref.point = getmaxP(combined)
  hv.moead = rep(0, max(moead.data$iter+1))
  for(i in 0:max(moead.data$iter)){
    Y.moead.iter = moead.data %>% filter(iter <= i)
    hv.moead[i+1] = hypervolume(get.objectives.values(Y.moead.iter), ref = ref.point)
  }
  
  hv.nsga2 = rep(0, max(nsga2.data$iter+1))
  for(i in 0:max(nsga2.data$iter)){
    Y.nsga2.iter = nsga2.data %>% filter(iter <= i)
    hv.nsga2[i+1] = hypervolume(get.objectives.values(Y.nsga2.iter), ref = ref.point)
  }
  hvs = rbind(data.frame(hv = hv.moead, name = 'MOEA/D', iter = unique(moead.data$iter)), data.frame(hv = hv.nsga2, name = 'NSGA-II', iter = unique(nsga2.data$iter)))
  ggplot() + geom_point(data = hvs ,aes(iter, hv, color = name)) + theme_minimal()
  
  filename = paste0(outpath, prob, '_hvs.pdf')
  ggsave(filename, device = "pdf")
  print("done")  
  
  ## generate igd convergence
  cat('igd convergence on going...')
  igd.moead = rep(0, max(moead.data$iter+1))
  for(i in 0:max(moead.data$iter)){
    Y.moead.iter = moead.data %>% filter(iter <= i)
    igd.moead[i+1] = igd(get.objectives.values(Y.moead.iter), reference = pf)
  }
  
  igd.nsga2 = rep(0, max(nsga2.data$iter+1))
  for(i in 0:max(nsga2.data$iter)){
    Y.nsga2.iter = nsga2.data %>% filter(iter <= i)
    igd.nsga2[i+1] = igd(get.objectives.values(Y.nsga2.iter), reference = pf)
  }
  igds = rbind(data.frame(igd = igd.moead, name = 'MOEA/D', iter = unique(moead.data$iter)), data.frame(igd = igd.nsga2, name = 'NSGA-II', iter = unique(nsga2.data$iter)))
  ggplot() + geom_point(data = igds ,aes(iter, igd, color = name)) + theme_minimal()
  filename = paste0(outpath, prob, '_igds.pdf')
  ggsave(filename, device = "pdf")
  print("done")  
  
  ## generate pareto front w/ theoretical points
  cat('approximation of pareto front (all solutions) w/ theoretical PF on going...')
  Y.moead = get.objectives.values(moead.data)
  Y.nsga2 = get.objectives.values(nsga2.data)
  Y.moead$name = "MOEA/D"
  Y.nsga2$name = "NSGA-II"
  pf$name = "PF"
  Ys = rbind(Y.moead, Y.nsga2, pf)
  
  ggplot() + geom_point(data = Ys ,aes(Y1, Y2,color = name)) + theme_minimal()
  filename = paste0(outpath, prob, '_paretofront.pdf')
  ggsave(filename, device = "pdf")
  
  print("done")  

    ## generate pareto front (nondom only) w/ theoretical points
  cat('approximation of pareto front (nondom only) w/ theoretical PF on going...')
  nd = find_nondominated_points(Y.moead)
  Y.moead = Y.moead[nd,]
  nd = find_nondominated_points(Y.nsga2)
  Y.nsga2 = Y.nsga2[nd,]
  Ys = rbind(Y.moead, Y.nsga2, pf)
  
  ggplot() + geom_point(data = Ys ,aes(Y1, Y2,color = name)) + theme_minimal()
  filename = paste0(outpath, prob, '_nondom_paretofront.pdf')
  ggsave(filename, device = "pdf")
  
  print("done")  
  
  ## generate eaf comparison
  cat('eaf comparision on going...')
  filename = paste0("../overall_results/", prob, '_eafdiffplot.eps')
  setEPS()
  postscript(file = filename, width = 512, height = 256)
  eaf.moead = moead.data %>% select(run, starts_with("Y"))
  eaf.nsga2 = nsga2.data %>% select(run, starts_with("Y"))
  eafdiffplot(data.left = eaf.moead, data.right  = eaf.nsga2, full.eaf = F, col = c("white", "orange", "red"),title.left = "MOEA/D",
              title.right = "NSGA-II")
  dev.off()
  print("done")
}