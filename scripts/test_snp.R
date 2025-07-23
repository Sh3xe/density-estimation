rm(list=ls())
source("read.mesh.R")
library(fdaPDE)

# Read and convert the data
data_csv <- read.csv("./data/snp500/original_data.csv")
df <- data.frame(data_csv)
df <- df[c("Date", "SP500")]
df$Date <- as.Date(df$Date)
df <- df[df$Date > as.Date("2000-01-01"),]
size <- length(df$Date)
space_data <- df[2:size,]$SP500 - df[1:size-1,]$SP500
time_data <- df[2:size,]$Date - df$Date[1]

# Plot
plot(x = df$Date, y=df$SP500)
hist(space_data, breaks=30, freq=FALSE)

# Convert data points to csv
write.csv(space_data, "data_space.csv")
write.csv(time_data, "data_time.csv")

# Creation of the mesh
create.mesh.1.5D(
	
)