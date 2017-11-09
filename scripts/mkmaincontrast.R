mkmaincontrast <- function(dmat,fcts=c("R","G","S","T"),contrast=list(fct="R",lev=c("E","M")),constraint=list(fct=setdiff(fcts,contrast$fct),lev=c("M","C","H")))
{

# usage:
# mkmaincontrast(dmat,fct,fcts,contrast,constraint)
#
# arguments
#
# dmat design matrix of group model the group labels must come from design table. Read it into workspace before call this function
# fcts=c("R","G","S","T")  #for "regime","gender","status","tissue"   the names of the main factor and the order must the same as the 
#  code for group in design table. for example the code "E_F_C_H" is one instance (group) of "regime_gender_status_tissue". For psudooscura 
#  data, using dafault value for argument:  fcts,   
# contrast=list(fct="R",lev=c("E","M"))  a list to specify the contrast. contrast$fct define which main factoe is concerned 
#  contrast$lev specifies which two levels are contrast, any the first level is take as reference
#  the default contrast set the M/E contrast in regime factor
# constraint=list(fct=setdiff(fcts,contrast$fct),lev=c("M","C","H"))  a list to specify the constraint put on the contrast
#  contain two character vectors, constrant$fac and constraint$lev. The former give the factors which constraint may be applied
#  the constraint%lev gives the range of the constaint, which must be either "all" or a string for the one level of the corresponding
#  factor. The two vectors must have same length and correct correspondence.
#
# The default argument values are for your data, and contrast two regime level M/E for Male Head samples only
# 
# output value:
# A single column matrix for specifying the argument "contrast" in glmLRT function (edgeR package)
#


  pra <- colnames(dmat)
  pra <- substr(pra,nchar("group")+1,100L)
  x <- unlist(strsplit(pra[1],split="_"))
  stopifnot(length(x) == length(fcts))          # stop if the number of factor shown by fcts and parameter names are different
  x <- t(matrix(unlist(strsplit(pra,split="_")),nrow=length(fcts)))
  colnames(x) <- fcts
  #
  # check the argument contrast
  if(is.null(contrast$fct) | any(is.null(contrast$lev))) stop("NULL component exist in contrast list") 
  if(!any(fcts %in% contrast$fct) | length(contrast$fct) > 1) stop("check main contrast factor, wrong input")
  xc <- which(fcts %in% contrast$fct)
  if(length(contrast$lev) != 2 | sum(contrast$lev %in% x[,xc]) !=2)  stop("check levels of contrast factor, wrong input")
  #
  # check the argument constraint
  if(length(constraint$lev) != length(constraint$fct))  stop("check levels of constraint, wrong input")
  for(k in 1:length(constraint$fct)) {
    xc <- which(fcts %in% constraint$fct[k])
    if(sum(constraint$lev[k] %in% x[,xc]) !=1 & tolower(constraint$lev[k]) != "all") stop("check levels of contrast factor, wrong input")
  }
  ## 
  ## make contrast vector
  cvect <- rep(0,nrow(x))
  xc <- which(fcts %in% contrast$fct)
  cvect[x[,xc] %in% contrast$lev[1]] <- -1
  cvect[x[,xc] %in% contrast$lev[2]] <- 1
  ## get the constrained set (row number of dmat) 
  constr <- rep(1,nrow(x))
  for(k in 1:length(constraint$fct)) {
    xc <- which(fcts %in% constraint$fct[k]) 
    if(sum(constraint$lev[k] %in% x[,xc]) ==1) constr <- constr * (x[,xc] %in% constraint$lev[k])
  }
  ## update the contrast vector
  cvect <- cvect * constr
  ## get the scale right
  cvect[cvect > 0] <- cvect[cvect>0]/sum(cvect>0)
  cvect[cvect < 0] <- cvect[cvect<0]/sum(cvect<0)
  names(cvect) <- colnames(dmat)
  print(contrast)	
  print(constraint)
  print(cvect)
  cmat <- matrix(cvect,nrow=length(cvect),,dimnames=list(colnames(dmat),"specified"))
} 

## end of the function
#
##### application example
#
# contrast V/C by constrained in E fomale head samples
# contrast <- list(fct="S",lev=c("C","V"))
# constraint <- list(fct=c("R","G","T"),lev=c("E","F","H"))
# fcts <- c("R","G","S","T")
##  load or read "dmat"
## path <- "X:/codefun"     #path for to the folder where function mkmaincontrast.R was saved
## source(paste0(path,"/mamaincontrast.R"))   
# cvector <- mkanycontrast(dmat,fcts=fcts,contrast=contrast,constraint=constraint) 
# contres <-glmLRT(fitres,contrast=cvector)
# 
##### extract the results 
# test_results_table <- contres$table
# contrast_made <- contres$comparison
##### add FDR value
# test_results_table <- cbind(test_results_table, FDR=p.adjust(test_results_table$PVal,method="fdr"))
#


## Email
# Hi Paris,
# 
# I coded a function to generate the contrast matrix for you mentioned contrast. It also can be used for any similar contrasts. By changing the argument "contrast" and "constraint" properly, the function give the contrast matrix for your specified contrast.
# 
# I have detailed comment lines with the function code, also an application example for contrasting V/C two levels within E Male Head samples. 
# 
# Actually, you can contrast two levels in any a factor, while the range of the other factors may be constrained. By changing the values for the factors you can specify any constraint you wanted. In addition, if you do not want a factor to be contained, just simply set its constraint value as "all". 
# 
# You said you want to contrast wo levels M/E in regime factor, with the constraint on other three factors is Male Head. the contrast argument of function mkmaincontrast can be set as contrast=list(fct="R", lev = c("E","V"))   where 
# fct="R" specifies the factor of concerning is "Regime"
# lev=c("E","M") specifies the two levels "E" "M" are compared using "E" (the first element of the vector) as reference   
# The three constrained factors are "G" for gender  "S" for status (of courtship) and "T" for tissue, if as you have mentioned for Male Head tissue only, the constraint argument of the function can be specified as:
# constraint=list(fct=c("G","S","T"),lev=("M","all","H").
# The argument value means "gender" factor is constrained in "M" (Male) only, the "S" status factor value is "all" which says not be constrained, and "T" (tissue) factors has value "H" (head) so it is constrained in Head samples.
# The above constraint and contrast values are the default value of the function. So, if the design matrix "dmat" is  available in R workspace, you can simply the obtain the contrast matrix for the contrast by 
# c1mat <- mkmaincontrast(dmat)
# the c1mat can be used to specify the argument "contrast" in glmLRT, that is
# glmLRT(glmfitobj,contrast=c1mat)
# 
# I put the function "mkmaincontrast.R" in  http://cgr.liv.ac.uk/illum/ID0095_7c86f14de6d6f17d/code/
# 
# Let me know if you encounter any problem in using this function.
# 
# 
# 
# 
