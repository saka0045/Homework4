# The script needs to be called with:
# Rscript plot_variability.R /path/to/variability.csv
args <- commandArgs(trailingOnly=TRUE)

# Read in the variability csv file
cat("Reading in:\n")
print(args[1])
variability_table <- read.csv(args[1], header = TRUE)

# Create PDF and plot
pdf(file = "variability.pdf", width = 16, height = 8)
plot(variability_table$Base, variability_table$Smoothed.Variability, type = "l", 
     xlab = "Position (bp)", ylab = "Sequence Identity (%)", 
     main = "Sequence Identity (%) vs Position (bp)", xaxt = "n")
axis(1, at = seq(0, 1500, by = 100), las = 2)
abline(h = 75, col = "red")

cat("\nScript is done running\n")
cat("Variability graph is saved as variability.pdf\n")