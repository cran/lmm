\name{example.lmm}
\alias{example.lmm}
\title{
lmm library example command file
}

\description{
These data, taken from an article by Weil, A.T., Zinberg, N.E. and
Nelson, J.M. (1968; Clinical and psychological effects of marihuana
in man; Science, 162, 1234-1242), come from a pilot study to
investigate the clinical and psychological effects of marijuana use
in human subjects. Nine subjects subjects each received three
treatments---low-dose, high-dose, and placebo. Under each treatment,
changes in heart rate (beats per minute) were measured 15 and 90
minutes after administration. NA denotes a missing value.
}

\seealso{
\code{\link[lmm]{ecmeml.lmm}}, \code{\link[lmm]{ecmerml.lmm}},
\code{\link[lmm]{fastml.lmm}}, \code{\link[lmm]{fastrml.lmm}},
\code{\link[lmm]{fastmcmc.lmm}}, \code{\link[lmm]{fastmode.lmm}}, 
\code{\link[lmm]{mgibbs.lmm}}
}

\examples{
#
#   -----------------------------------------------------------------
#                         15 minutes                 90 minutes
#                   ----------------------     ----------------------
#                   Placebo   Low   High      Placebo   Low   High
#   -----------------------------------------------------------------
#   Subject 1          16     20     16          20     -6     -4       
#           2          12     24     12          -6      4     -8
#           3           8      8     26          -4      4      8
#           4          20      8     NA          NA     20     -4
#           5           8      4     -8          NA     22     -8
#           6          10     20     28         -20     -4     -4
#           7           4     28     24          12      8     18
#           8          -8     20     24          -3      8    -24
#           9          NA     20     24           8     12     NA
#   -----------------------------------------------------------------
#
########################################################################
# Below we show how to fit a traditional compound symmetry model
# with a fixed effect for each column (occasion) and a random 
# intercept for each subject. First we enter the data.
#
y <- c(16,20,16,20,-6,-4,
    12,24,12,-6,4,-8,
    8,8,26,-4,4,8,
    20,8,20,-4,
    8,4,-8,22,-8,
    10,20,28,-20,-4,-4,
    4,28,24,12,8,18,
    -8,20,24,-3,8,-24,
    20,24,8,12)
occ <- c(1,2,3,4,5,6,
      1,2,3,4,5,6,
      1,2,3,4,5,6,
      1,2,5,6,
      1,2,3,5,6,
      1,2,3,4,5,6,
      1,2,3,4,5,6,
      1,2,3,4,5,6,
      2,3,4,5)
subj <- c(1,1,1,1,1,1,
       2,2,2,2,2,2,
       3,3,3,3,3,3,
       4,4,4,4,
       5,5,5,5,5,
       6,6,6,6,6,6,
       7,7,7,7,7,7,
       8,8,8,8,8,8,
       9,9,9,9)
########################################################################
# Now we must specify the model. 
# If the six measurements per subject were ordered in time, we might
# consider using a model with time of measurement entered with linear
# (or perhaps higher-order polynomial) effects. But because the
# six measurements are not clearly ordered, let's use a model that has
# an intercept and five dummy codes to allow the population means for
# the six occasions to be estimated freely. We will also allow the
# intercept to randomly vary by subject. For a subject i with no
# missing values, the covariate matrices will be
#
#                   1 1 0 0 0 0              1
#                   1 0 1 0 0 0              1
#           Xi =    1 0 0 1 0 0       Zi =   1
#                   1 0 0 0 1 0              1
#                   1 0 0 0 0 1              1
#                   1 0 0 0 0 0              1
#
# The Xi's and Zi's are combined into a single matrix called
# pred. The pred matrix has length(y) rows. Each column of Xi and Zi
# must be represented in pred. Because Zi is merely the first column
# of Xi, we do not need to enter that column twice. So pred is simply
# the matrices Xi (i=1,...,9), stacked upon each other.
#
pred <- cbind(int=rep(1,49),dummy1=1*(occ==1),dummy2=1*(occ==2),
    dummy3=1*(occ==3),dummy4=1*(occ==4),dummy5=1*(occ==5))
xcol <- 1:6
zcol <- 1
########################################################################
# Now find ML estimates using the ECME procedure and the faster
# scoring algorithm
#
ecmeml.result <- ecmeml.lmm(y,subj,pred,xcol,zcol)
fastml.result <- fastml.lmm(y,subj,pred,xcol,zcol)
#
# In this example, ECME converged in 212 cycles, but the fast
# algorithm took only 8. The results can be viewed by printing the
# various components of "ecmeml.result" and "fastml.result". 
# For example, extract the ML estimate of the fixed effects beta.
#
beta.hat <- fastml.result$beta
#
# Because of the dummy codes used in the Xi's, the first element of
# beta (the intercept) estimates the mean for the last occasion,
# and the other elements of beta estimate the differences in means
# between the first five occasions and the last one. So we can find
# the estimated means for the six occasions like this:
#
muhat <- c(beta.hat[2]+beta.hat[1], beta.hat[3]+beta.hat[1],
   beta.hat[4]+beta.hat[1], beta.hat[5]+beta.hat[1],
   beta.hat[6]+beta.hat[1], beta.hat[1]) 
#
# The functions for RML estimation work exactly the same way:
#
ecmerml.result <- ecmerml.lmm(y,subj,pred,xcol,zcol)
fastrml.result <- fastrml.lmm(y,subj,pred,xcol,zcol)
#
#######################################################################
# The function "fastrml.lmm" calculates the improved variance
# estimates for random effects described in Section 4 of Schafer
# (1998). The code below reproduces Table 2, which compares 
# 95\% interval estimates under the new method to conventional
# empirical Bayes intervals.
#
b.hat <- as.vector(fastrml.result$b.hat)
se.new <- sqrt(as.vector(fastrml.result$cov.b.new))
se.old <- sqrt(as.vector(fastrml.result$cov.b))
table2 <- cbind(round(b.hat,3),
   round(cbind(b.hat-2*se.old,b.hat+2*se.old,
      b.hat-2*se.new,b.hat+2*se.new),2),
   round(100*(se.new-se.old)/se.old))
dimnames(table2) <- list(paste("Subject",format(1:9)),
   c("Est.","Lower.old","Upper.old","Lower.new","Upper.new",
   "Increase (\%)"))
print(table2)
#
#######################################################################
# The functions "mgibbs.lmm" and "fastmcmc.lmm" perform the MCMC
# procedures described in Section 5. The code below runs each
# algorithm for 5,000 cycles, and then displays autocorrelation
# plots like those of Figure 1.
#
prior <- list(a=3*100,b=3,c=3,Dinv=3*5)
gibbs.result <- mgibbs.lmm(y,subj,pred,xcol,zcol,prior=prior,
   seed=1234,iter=5000)
fmcmc.result <- fastmcmc.lmm(y,subj,pred,xcol,zcol,prior=prior,
   seed=2345,iter=5000)
# 
# Before doing this, make sure that a graphics device is active:
library(ts)
par(mfrow=c(2,1))
#
# there were problems with debian Linux; so add ylim
acf(log(gibbs.result$psi.series[1,1,]),lag.max=10, ylim=0:1)
acf(log(fmcmc.result$psi.series[1,1,]),lag.max=10, ylim=0:1)

#######################################################################
}

\keyword{models}
