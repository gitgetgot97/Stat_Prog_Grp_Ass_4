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

th <- c(100,200)

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



rb <- function(th,k=2) {
  k*(th[2]-th[1]^2)^2 + (1-th[1])^2
}

## Check the function for different theta vectors
th <- c(2,1)
th_1 <- c(4,1)
rb(th); rb(th_1)


gb(th); gb(th_1)
hb(th); hb(th_1)




newt<- function(theta,func,grad,hess=NULL,...,tol=1e-8,fscale=1,maxit=100,max.half=20,eps=1e-6){
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
  
  iter <- 0 ## initialise count of number of iterations
  g <- grad(theta, ...) ## calculate the gradient vector for initial guess
  hess0 <- hess(theta, ...) ## calculate hessian matrix for initial guess
  I <- diag(length(theta)) ## create appropriately sized identity to add to hessian if it is not positive definite
  
  
  while(iter <= maxit 
        & abs(max(grad(theta, ...))) > (tol*(abs(func(theta,...))+fscale))){
    ## execute an iteration of Newton's method if we have not reached the maximum number of iterations
    ## and if the gradient of the objective is above the tolerance level
    
    mult = 1e-8 ## initial multiple of identity to add to non-positive definite hessian
    while(inherits(try(chol(hess0),silent = TRUE), "try-error") & mult < 10){
      ## check if the hessian is positive definite by attempting a Cholesky decomposition
      ## if it is not positive definite add identity matrices to it until it is
      hess0 <- hess0 + mult*I
      mult <- 10*mult
    }
    
    R <- chol(hess0) ## Cholesky decomposition of the hessian
    hess_inv<-chol2inv(R) ## use the decomposition to calculate the inverse
    
    delta <- -(hess_inv)%*%g ## find the delta that minimises the quadratic approximation
    temp_theta <- theta + delta ## store current guess for theta + delta
    half_iter <- 1 ## initialise count for the number of half steps used
    
    while(func(theta, ...) < func(temp_theta, ...) & half_iter <= max.half){
      ## If our initial guess is far from the objective, we may overshoot
      ## if the updated theta gives a worse objective value, add delta/2
      ## divide delta by 2 repeatedly until a better objective is foun
      temp_theta <- theta + delta*(1/2 ^ half_iter)
      if(abs(func(temp_theta,...)) == Inf){
        ## if the new half-step returns a non-finite value, reset temp_theta to its previous value
        ## this will simply move to the next half-step
        temp_theta <- theta + delta*(1/2 ^ (half_iter - 1))
      }
      half_iter <- half_iter + 1
    }
    
    if(half_iter == max.half){
      ## alert user that the maximum number of half step were used
      warning("Max half steps reached")
    }
    
    theta <- temp_theta
    f <- func(theta, ...)
    g <- grad(theta, ...)
    hess0 <- hess(theta, ...)
    iter <- iter + 1
  }
  
  R <- chol(hess0)
  Hi <- chol2inv(R)
  if(inherits(try(chol(hess0),silent = TRUE), "try-error")){
    warning("Hessian at minimum is not positive definite")
  }
  output <- list(theta = theta, f=f, iter=iter, g=g, Hi=Hi)
  return(output) 
}## End of newt function



trial <- newt(th, rb, gb, hess=hb)





