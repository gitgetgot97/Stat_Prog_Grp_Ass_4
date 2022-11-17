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

##############################################################################

finite_diff <- function(theta, grad, eps, ...){
  ## comment what this function does...
  ## Produces the hessian if user does not provide
  grad_0 <- grad(theta, ...)
  Hfd = matrix(0, length(theta), length(theta))
  
  for(i in 1:length(theta)){
    
    th1 <- theta; th1[i] <- th1[i] + eps
    grad_1 <- grad(th1, ...)
    Hfd[i,] <- (grad_1 - grad_0)/eps
  }
  return (t(Hfd) + Hfd)/2
} ## end of finite_diff function


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
  f <- func(theta, ...)
  
  if(!is.finite(func(theta, ...)) | sum(is.finite(grad(theta,...))) < length(theta)){
   stop("Objective or derivative in not finite at initial theta, please enter a different theta.") 
  }
  
  
  if(!is.null(hess)){
    hess0 <- hess(theta, ...) ## calculate hessian matrix for initial guess
  }else{
    hess0 <- finite_diff(theta, grad,eps,...)
  }
  I <- diag(length(theta)) ## create appropriately sized identity to add to hessian if it is not positive definite
  
  
  while(iter < maxit & abs(max(g)) > (tol*(abs(f+fscale)))){
    ## execute an iteration of Newton's method if we have not reached the maximum number of iterations
    ## and if the gradient of the objective is above the tolerance level
    
    mult = 1e-8 ## initial multiple of identity to add to non-positive definite hessian
    while(inherits(try(chol(hess0),silent = TRUE), "try-error")){
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
      if(is.finite(func(temp_theta,...))){
        ## if the new half-step returns a non-finite value, reset temp_theta to its previous value
        ## this will simply move to the next half-step
        temp_theta <- theta + delta*(1/2 ^ (half_iter - 1))
      }
      half_iter <- half_iter + 1
    }
    
    if(half_iter == max.half){
      ## alert user that the maximum number of half step were used
      warning("Max half steps reached without improving objective. Return most recent output.")
      ## break out of the larger while loop - objective won't be improved further
      break
    }
    
    
    theta <- temp_theta
    f <- func(theta, ...)
    g <- grad(theta, ...)
    
    if(!is.null(hess)){
      hess0 <- hess(theta, ...)
    }else{
      hess0 <- finite_diff(theta, grad, eps, ...)
    }
    iter <- iter + 1
  } ## End of outer while loop
  
  
  if(iter == maxit){
    ## alert user that the maximum number of iterations has been reached without convergence
    warning("Max iterations reached without convergence. Return most recent output.")
  }
  
  
  R <- try(chol(hess0), silent=TRUE)
  Hi <- try(chol2inv(R), silent=TRUE)
  if(inherits(try(chol(hess0),silent = TRUE), "try-error")){
    warning("Hessian at minimum is not positive definite")
    Hi = NULL
  }
  output <- list(theta = theta, f=f, iter=iter, g=g, Hi=Hi)
  return(output) 
}## End of newt function


######################################
## Testing

rb <- function(th,k=2) {
  k*(th[2]-th[1]^2)^2 + (1-th[1])^2
}

## gb returns the gradient of rb
gb <- function(th,k=2) {
  c(-2*(1-th[1])-k*4*th[1]*(th[2]-th[1]^2),k*2*(th[2]-th[1]^2))
}

## hb returns the hessian of hb
hb <- function(th,k=2) {
  h <- matrix(0,2,2)
  h[1,1] <- 2-k*2*(2*(th[2]-th[1]^2) - 4*th[1]^2)
  h[2,2] <- 2*k
  h[1,2] <- h[2,1] <- -4*k*th[1]
  h
}

th <- c(1,1)

trial <- newt(th, rb, gb) ## we need to know what to do when the minimum is entered as a first guess

#######################################################
## Another test

f3 <- function(x0) {
  x0^3 - 3*x0^2 + 1 + 100*sin(x0+1)
}

f3_prime <- function(x0) {
  3*x0^2 - 6*x0 + 100*cos(x0+1)
}

trial <- newt(10, f3, f3_prime)   ## 10 is our initial x0 guess

##################
## Another test

doogle <- function(x, alpha1=3, alpha2=1, alpha3=0.1, alpha4=0.3, alpha5=10){
  return(alpha1*x[1]^2 + alpha2*x[2]^2 + alpha3*x[3]^2 + alpha4*x[4]^2 + alpha5*x[5]^2)
}

doogle_diff <- function(x, alpha1=3, alpha2=1, alpha3=0.1, alpha4=0.3, alpha5=10){
  c(2*alpha1*x[1], 2*alpha2*x[2], 2*alpha3*x[3], 2*alpha4*x[4], 2*alpha5*x[5]^2)
}

xguess = c(10,10,12,3,12) 

doogle_test <- newt(xguess, doogle, doogle_diff)

###############
##Another Test

saddle <- function(y){
  return(y[1]^2 - y[2]^2)
}

saddle_diff <- function(y){
  c(2*y[1], 2*y[2])
}

y_guess <- c(0,0)

saddle_test <- newt(y_guess, saddle, saddle_diff) ### this lands us at a saddle point

## 1) Commenting
## 2) try test functions for fun. Aim to try a func which warns on half steps and on non finite grad.