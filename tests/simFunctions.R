library(phytools); library(phyclust); library(gtools); library(phangorn); library(phylobase)

#### FUNCTIONS
generateFDtree <- function(blam, dmu, ntaxa, simDepth) {
  for (i in 1:10000) {
    print(i)
    origTree <- pbtree(b = blam, d = dmu, n = ntaxa, scale = simDepth, nsim = 1, method = "continuous", extant.only = TRUE)
    if (!is.null(origTree)) {
      rootedge <- which(origTree$edge[,1] == (ntaxa + 1))
      singOut <- rootedge[which(origTree$edge[rootedge,2] < (ntaxa + 1))]
      if (length(singOut) == 1) {
        return(origTree)
        break()
      }
    }
  }
}

scaleFDtree <- function(gamma.mean, gamma.sd, ntaxa, origTree) {
  simgam <- rgamma(n = 10000, shape = gamma.mean^2/gamma.sd, rate = gamma.mean/gamma.sd)
  bl.scalars <- sample(size = (2*ntaxa - 2), x = simgam)
  scaleTree <- origTree
  scaleTree$edge.length <- scaleTree$edge.length * bl.scalars
  return(scaleTree)
}


getNodePriors <- function(tree, ntaxa, format = c("exp", "norm"), choose.sd, offset = 0.75, pre = prefix, nnodes) {
  temp <- phylo4(tree)
  treeheight <- max(nodeHeights(tree))
  nodePriors <- nPexp <- NULL
  for (j in (ntaxa+1):(2*ntaxa - 1)) {
    tipNames <- names(descendants(temp, j, type = "tips"))
    tipNames <- c(tipNames[1], tail(tipNames, 1))
    if (getMRCA(tree, c(tipNames)) == j) {
      heights <- treeheight - nodeheight(tree, j)
      if (format == "exp")
        nodePrior <- paste(paste(tipNames, collapse = ","), " exp (", round(heights, 4), "," , round(offset*heights, 4), ")", sep = "")
      if (format == "norm") 
        nodePrior <- paste(paste(tipNames, collapse = ","), " norm (", round(heights, 4), "," , choose.sd, ",", round(offset*heights, 4), ")", sep = "")
      nodePriors <- c(nodePriors, nodePrior)
    } 
  }
  which.nodes <- sample(1:length(nodePriors), nnodes)
  write(nodePriors[which.nodes], paste(pre, "-", nnodes, "nodes.txt", sep = ""))
  return(paste(pre, "-", nnodes, "nodes.txt", sep = ""))
}

compareNodes2Tree <- function(path2tr, path2node) {
  a <- read.tree(path2tr)
  treeHeight <- max(nodeHeights(a))
  np <- read.table(path2node)
  comptab <- mrcas <- NULL
  for (i in 1:dim(np)[1]) {
    nodes <- unlist(strsplit(as.vector(np[i,1]), ","))
    mean <- gsub("\\(", "", strsplit(as.vector(np[i,3]), ",")[[1]][1])
    mrca <- getMRCA(a, nodes)
    mrcas <- c(mrcas, mrca)
    est <- treeHeight - nodeheight(a, mrca)
    comptab <- rbind(comptab, as.numeric(c(mean, est)))
  }
  return(list(comptab, mrcas))
}


compareNoPrior <- function(outTr, mrcaNoPr, origTree) {
  a <- read.tree(outTr)
  b <- read.tree(origTree)
  treeHeight2 <- max(nodeHeights(a))
  treeHeight1 <- max(nodeHeights(b))
  comptab <- NULL
  for (i in 1:length(mrcaNoPr)) {
    sim <- treeHeight1 - nodeheight(b, mrcaNoPr[i])
    est <- treeHeight2 - nodeheight(a, mrcaNoPr[i])
    comptab <- rbind(comptab, c(sim, est))
  }
  return(comptab)
}
