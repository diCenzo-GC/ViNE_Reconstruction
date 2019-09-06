#!/usr/bin/env Rscript

# Init system
library('agricolae')

# Load plant data
data.p <- read.table('rnaSeqModel_plant.txt')
rnames<- data.p[,1]
data.p <- data.frame(data.matrix(data.p[,2:ncol(data.p)]))
row.names(data.p) <- rnames

# Load bacterial data
data.b <- read.table('rnaSeqModel_bacterium.txt')
rnames<- data.b[,1]
data.b <- data.frame(data.matrix(data.b[,2:ncol(data.b)]))
row.names(data.b) <- rnames

# Make output variables
out.p <- data.frame(matrix(, nrow = nrow(data.p), ncol = 5))
out.b <- data.frame(matrix(, nrow = nrow(data.b), ncol = 5))

# Prepare intermediate variables for plants
test.p <- data.frame(matrix(, nrow = 15, ncol = 2))
test.p[1,1] <- 'A'
test.p[2,1] <- 'A'
test.p[3,1] <- 'A'
test.p[4,1] <- 'B'
test.p[5,1] <- 'B'
test.p[6,1] <- 'B'
test.p[7,1] <- 'C'
test.p[8,1] <- 'C'
test.p[9,1] <- 'C'
test.p[10,1] <- 'D'
test.p[11,1] <- 'D'
test.p[12,1] <- 'D'
test.p[13,1] <- 'E'
test.p[14,1] <- 'E'
test.p[15,1] <- 'E'
temp.p <- data.frame(matrix(, nrow = 5, ncol = 2))

# Prepare intermediate variables for bacteria
test.b <- data.frame(matrix(, nrow = 15, ncol = 2))
test.b[1,1] <- 'A'
test.b[2,1] <- 'A'
test.b[3,1] <- 'A'
test.b[4,1] <- 'B'
test.b[5,1] <- 'B'
test.b[6,1] <- 'B'
test.b[7,1] <- 'C'
test.b[8,1] <- 'C'
test.b[9,1] <- 'C'
test.b[10,1] <- 'D'
test.b[11,1] <- 'D'
test.b[12,1] <- 'D'
test.b[13,1] <- 'E'
test.b[14,1] <- 'E'
test.b[15,1] <- 'E'
temp.b <- data.frame(matrix(, nrow = 5, ncol = 2))

# Perform the analysis for the plant data
for (n in 1:nrow(data.p)) {
  test.p[,2] <- t(data.p[n,])
  stats.out <- kruskal(test.p$X2, test.p$X1, console = FALSE)
  temp.p[,1] <- rownames(stats.out$groups)
  temp.p[,2] <- as.character(stats.out$groups$groups)
  temp.p <- temp.p[order(temp.p[,1]),]
  out.p[n,] <- t(temp.p[,2])
}
row.names(out.p) <- rownames(data.p)
  
# Perform the analysis for the bacterium data
for (n in 1:nrow(data.b)) {
  test.b[,2] <- t(data.b[n,])
  stats.out <- kruskal(test.b$X2, test.b$X1, console = FALSE, alpha = 0.2)
  temp.b[,1] <- rownames(stats.out$groups)
  temp.b[,2] <- as.character(stats.out$groups$groups)
  temp.b <- temp.b[order(temp.b[,1]),]
  out.b[n,] <- t(temp.b[,2])
}
row.names(out.b) <- rownames(data.b)

# Save and export
write.table(out.p, file = "rnaSeqModel_plant_stats.txt")
write.table(out.b, file = "rnaSeqModel_bacterium_stats.txt")
save.image(file = "rnaSeqStatsWorkspace.RData")
quit("no")
  
  