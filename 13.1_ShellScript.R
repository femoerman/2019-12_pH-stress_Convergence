#clear workspace
rm(list=ls())

#import packages needed for script
library(parallel)
library(foreach)
library(doParallel)

#create list with files to run
 files <- c("/media/felix/DataDrive2/Documenten/PhD/05_student_project_2017/3_analysis/13a_cascading_Bayesian_Models_CU427.4.R",
"/media/felix/DataDrive2/Documenten/PhD/05_student_project_2017/3_analysis/13b_cascading_Bayesian_Models_SB3539.R",
"/media/felix/DataDrive2/Documenten/PhD/05_student_project_2017/3_analysis/13c_cascading_Bayesian_Models_B2086.2.R",
"/media/felix/DataDrive2/Documenten/PhD/05_student_project_2017/3_analysis/13d_cascading_Bayesian_Models_CU428.2.R",
"/media/felix/DataDrive2/Documenten/PhD/05_student_project_2017/3_analysis/15a_cascading_Bayesian_Models_K_CU427.4.R",
"/media/felix/DataDrive2/Documenten/PhD/05_student_project_2017/3_analysis/15c_cascading_Bayesian_Models_K_SB3539.R",
"/media/felix/DataDrive2/Documenten/PhD/05_student_project_2017/3_analysis/15d_cascading_Bayesian_Models_K_B2086.2.R",
"/media/felix/DataDrive2/Documenten/PhD/05_student_project_2017/3_analysis/15b_cascading_Bayesian_Models_K_CU428.2.R",
"/media/felix/DataDrive2/Documenten/PhD/05_student_project_2017/3_analysis/17a_cascading_Bayesian_Models_alfa_CU427.4.R",
"/media/felix/DataDrive2/Documenten/PhD/05_student_project_2017/3_analysis/17c_cascading_Bayesian_Models_alfa_SB3539.R",
"/media/felix/DataDrive2/Documenten/PhD/05_student_project_2017/3_analysis/17d_cascading_Bayesian_Models_alfa_B2086.2.R",
"/media/felix/DataDrive2/Documenten/PhD/05_student_project_2017/3_analysis/17b_cascading_Bayesian_Models_alfa_CU428.2.R")

#start parallelization --> detect number of cores
no_cores <- detectCores() - 1

# Initiate cluster for parallelization
cl <- makeCluster(no_cores)
registerDoParallel(cl)

#start the parallelized process --> Run all bayesian analyses in parallel
#x <- foreach(row=1:nrow(scenario), .combine=c) %dopar% (wrapper(row, scenario, run.program, input))
mclapply(files, source, mc.cores=no_cores)

#teminate the cluster
stopCluster(cl)