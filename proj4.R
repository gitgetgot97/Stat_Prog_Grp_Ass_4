# s1704264 Benedict Troy
# s2451950 Scott Jenkins
# s2453724 Emily Segers

# Github Repo:
# https://github.com/gitgetgot97/Stat_Prog_Grp_Ass_4

# Contribution:
# None of us contributed significantly more, or less, than the others. 
# Hence we deserve equal marks on this assignment. 

## Newton Optimizer

## This programme contains an R function, newt, implementing Newton's method for
## minimization of functions.

###############################################################################

## Test Function rb

rb <- function(th,k=2) {
k*(th[2]-th[1]^2)^2 + (1-th[1])^2
}

th <- c(2,1)

rb(th)

## gb returns the gradient of rb
gb <- function(th,k=2) {
c(-2*(1-th[1])-k*4*th[1]*(th[2]-th[1]^2),k*2*(th[2]-th[1]^2))
}

gb(th)

## hb returns the hessian of hb
hb <- function(th,k=2) {
h <- matrix(0,2,2)
h[1,1] <- 2-k*2*(2*(th[2]-th[1]^2) - 4*th[1]^2)
h[2,2] <- 2*k
h[1,2] <- h[2,1] <- -4*k*th[1]
h
}

hb(th)



## SJ: Trying to create grad and hess functions
## For reference: there is package for this, but can't use it, because not in base R
## https://cran.r-project.org/web/packages/EKMCMC/index.html

rb <- function(th,k=2) {
  k*(th[2]-th[1]^2)^2 + (1-th[1])^2
}

## Check the function for different theta vectors
th <- c(2,1)
th_1 <- c(4,1)
rb(th); rb(th_1)

## Want to differentiate rb wrt th[1] and th[2], grad is a vector of these. hess is the second partial derivatives

deriv(rb(th), "th[1]")           ## these aren't working as I would expect. says grad = 0.
deriv(rb(th), "th[2]")
deriv(rb(th_1), "th[1]")
eval(deriv(rb(th_1), "th[1]"))

gb(th); gb(th_1)
hb(th); hb(th_1)


# A more simple example in one variable
f = expression(x^2 + 5*x + 4)
deriv(f, 'x')                                    ## QUESSTION: WE WANT .grad[,'x'] from here - how do we access it??
x <- 1
eval(f)   ## eval function to find the value at a given x
eval(deriv(f,'x'))  
deriv(deriv(f,'x'),'x')

## Have a good weekend both!
## Scott



newt<- function(theta,func,grad,hess=NULL,...,tol=1e-8,fscale=1,maxit=100,max.half=20,eps=1e-6){
  
}## End of newt function

## theta is a vector of initial values for the optimization parameters.
## func is the objective function to minimize. Its first argument is the vector of optimization parameters. 
## Remaining arguments will be passed from newt using '...'.
## grad is the gradient function. It has the same arguments as func 
## but returns the gradient vector of the objective w.r.t. the elements of parameter vector.
## hess is the Hessian matrix function. 
## It has the same arguments as func but returns the Hessian matrix of the objective w.r.t. the elements of parameter vector. 
## If not supplied then newt should obtain an approximation to the Hessian by finite differencing of the gradient vector (hint: (t(A)+A)/2 is exactly symmetric for any matrix A).
## ... any arguments of func, grad and hess after the first (the parameter vector) are passed using this.
## tol the convergence tolerance.
## fscale a rough estimate of the magnitude of func near the optimum - used in convergence testing.
## maxit the maximum number of Newton iterations to try before giving up.
## max.half the maximum number of times a step should be halved before concluding that the step has failed to improve the objective.
## eps the finite difference intervals to use when a Hessian function is not provided.


## newt should return a list containing:
## f the value of the objective function at the minimum.
## theta the value of the parameters at the minimum.
## iter the number of iterations taken to reach the minimum.
## g the gradient vector at the minimum (so the user can judge closeness to numerical zero).
## Hi the inverse of the Hessian matrix at the minimum (useful if the objective is a negative log likelihood).



