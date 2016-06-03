setwd("~/Desktop/fastDate/")
source("simFunctions.R")

# Run prefix
prefix <- format(Sys.time(), "%b%d-%H%M")
#### Set tree parameters
ntips <- 100
trHeight <- 30 
rho <- 1
lambda <- 2
mu <- 1
#### Generate an ultrametric tree
myTree <- generateFDtree(lambda, mu, ntips, trHeight)
write.tree(myTree, paste(prefix, "-origTree.tre", sep = ""))
plot(myTree)
#### Scale branch lengths
gmean <- 0.0168
gsd <- 0.0018
myTree.scaled <- scaleFDtree(gmean, gsd, ntips, myTree)
write.tree(myTree.scaled, paste(prefix, ".tre", sep = ""))
#### Pull "true" data
myNP <- getNodePriors(myTree, ntips, format = "norm", 1, nnodes = 50)
#### Make fastdate command file
path2FD <- "~/bin/speed-dating-master/src/fastdate"
ngrid <- 3000
cmd <- paste(path2FD, "--method_nodeprior --grid", ngrid, "--prior_file", myNP, "--bd_mu", mu, "--bd_lambda", lambda, "--max_age", trHeight*3, " --bd_rho ", rho, "--tree_file", paste(prefix, ".tre", sep = ""),"--rate_mean", gmean, "--rate_variance", gsd, "--out_file", paste(prefix, ".out", sep = ""))
write(cmd, paste(prefix, ".sh", sep = ""))
# RUN FASTDATE

## LOOK AT RESULTS
### How well do they match set node priors?
list.results <- compareNodes2Tree("Jun03-1240.out", "Jun03-1240-50nodes.txt")
# Make plot to compare nodepriors to estimated
compare <- list.results[[1]]
plot(compare[,1], pch = 16, ylab = "Node Depth", bty = "n", ylim = c(0, 60))
points(compare[,2], col = "red", pch = 16)
for (i in 1:dim(compare)[1]) {
  lines(c(compare[i,1], compare[i,2])~ c(i,i), lty = 3)
}
legend("topright", legend = c("simulated", "estimated"), col = c("black", "red"), bty = "n", pch = 16)

## compare 
node.num <- (ntips + 1) : (2 * ntips - 1)
node.noPrior <- setdiff(node.num, results[[2]])
NoPr <- compareNoPrior("Jun03-1240.out", node.noPrior, "Jun03-1240-origTree.tre")
# Make plot to compare nodes without prior to estimated
plot(NoPr[,1], pch = 16, ylab = "Node Depth", bty = "n", ylim = c(0, 60))
points(NoPr[,2], col = "red", pch = 16)
for (i in 1:dim(NoPr)[1]) {
  lines(c(NoPr[i,1], NoPr[i,2])~ c(i,i), lty = 3)
}
legend("topright", legend = c("simulated", "estimated"), col = c("black", "red"), bty = "n", pch = 16)

