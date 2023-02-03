library('sva')

args <- commandArgs(T)
output <- args[1]

data <- read.csv(paste0(output,'/data.csv'),header = F)
batch <- read.csv(paste0(output,'/batch.csv'),header = F)$V1
combat_data <- ComBat(data, batch=batch)
write.csv(combat_data,file=paste0(output,'/combat_data.csv'),row.names = FALSE)
