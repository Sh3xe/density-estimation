# Set the working directory to where your files are located
# setwd("your_directory_path")

# List all files starting with 'cpp_' and ending with '.csv'
csv_files <- list.files(path = "./outputs", pattern = "^cpp_.*\\.csv$")

# Apply the transformation to each file
for (file in csv_files) {
  csv = read.csv( paste("./outputs/", file, sep="") )
	write.csv(csv, paste("./fix_out/", file, sep=""))
}
