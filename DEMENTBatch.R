args <- commandArgs(trailingOnly=T)
job.time <- args[1]
task.ID <- as.numeric(args[2])

source("DEMENT.0.7.2.R")

out<- TraitModel(job.time,task.ID)
filename <- paste("outputs7/", out[[1]]$timestamp, ".RData", sep = "")
save.image(file=filename)
