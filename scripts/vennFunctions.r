Venn2 <- function(set1, set2, names, title)
{
  stopifnot( length(names) == 2)
 
  # Form universe as union of all three sets
  universe <- sort( unique( c(set1, set2) ) )
 
  Counts <- matrix(0, nrow=length(universe), ncol=2)
  colnames(Counts) <- names
 
  for (i in 1:length(universe))
  {
    Counts[i,1] <- universe[i] %in% set1
    Counts[i,2] <- universe[i] %in% set2
  }
 
  vennDiagram( vennCounts(Counts),cex=c(1.2,1.1,1), lwd=2, circle.col=c(4,15,3,2,6), counts.col=1, main=title )
}

Venn3 <- function(set1, set2, set3, names, title)
{
  stopifnot( length(names) == 3)
 
  # Form universe as union of all three sets
  universe <- sort( unique( c(set1, set2, set3) ) )
 
  Counts <- matrix(0, nrow=length(universe), ncol=3)
  colnames(Counts) <- names
 
  for (i in 1:length(universe))
  {
    Counts[i,1] <- universe[i] %in% set1
    Counts[i,2] <- universe[i] %in% set2
    Counts[i,3] <- universe[i] %in% set3
  }
 
  vennDiagram( vennCounts(Counts),cex=c(1.2,1.1,1), lwd=2, circle.col=c(4,15,3,2,6), counts.col=1, main=title )
}


Venn4 <- function(set1, set2, set3, set4, names, title)
{
  stopifnot( length(names) == 4)
 
  # Form universe as union of all three sets
  universe <- sort( unique( c(set1, set2, set3, set4) ) )
 
  Counts <- matrix(0, nrow=length(universe), ncol=4)
  colnames(Counts) <- names
 
  for (i in 1:length(universe))
  {
    Counts[i,1] <- universe[i] %in% set1
    Counts[i,2] <- universe[i] %in% set2
    Counts[i,3] <- universe[i] %in% set3
    Counts[i,4] <- universe[i] %in% set4

  }
 
  vennDiagram( vennCounts(Counts),cex=c(1,1.1,1), lwd=2, circle.col=c(4,15,3,2,6), counts.col=1, main=title )
}

Venn5 <- function(set1, set2, set3, set4, set5, names, title)
{
  stopifnot( length(names) == 5)
 
  # Form universe as union of all three sets
  universe <- sort( unique( c(set1, set2, set3, set4, set5) ) )
 
  Counts <- matrix(0, nrow=length(universe), ncol=5)
  colnames(Counts) <- names
 
  for (i in 1:length(universe))
  {
    Counts[i,1] <- universe[i] %in% set1
    Counts[i,2] <- universe[i] %in% set2
    Counts[i,3] <- universe[i] %in% set3
    Counts[i,4] <- universe[i] %in% set4
    Counts[i,5] <- universe[i] %in% set5
  }
 
  vennDiagram( vennCounts(Counts),cex=c(1.2,1.1,1), lwd=2, circle.col=c(4,15,3,2,6), counts.col=1, main=title )
}


