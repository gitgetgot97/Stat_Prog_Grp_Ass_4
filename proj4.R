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
## minimization of functions. Detailed description within the function

##############################################################################

finite_diff <- function(theta, grad, eps, ...){
  ## When the hessian is not provided, this function provides an estimate by
  ## finite differencing
  grad_0 <- grad(theta, ...)     ## store the initial gradient as grad_0
  Hfd = matrix(0, length(theta), length(theta)) ## initialise matrix of correct size
  
  for(i in 1:length(theta)){
    ## creates the ith row of the hessian
    th1 <- theta; th1[i] <- th1[i] + eps  ## shift theta by small value
    grad_1 <- grad(th1, ...)              ## compute gradient at shifted theta
    Hfd[i,] <- (grad_1 - grad_0)/eps      ## insert values into hessian matrix
  }
  return (t(Hfd) + Hfd)/2                 ## ensures that output is symmetric
} ## end of finite_diff function


newt<- function(theta,func,grad,hess=NULL,...,tol=1e-8,fscale=1,maxit=100,max.half=20,eps=1e-6){
  ## theta is a vector of initial values for the optimization parameters
  ## func is the objective function to minimize. Its first argument is theta 
  ## Remaining arguments will be passed from newt using '...'.
  ## grad is the gradient function. It has the same arguments as func...
  ## but returns the gradient vector of the objective w.r.t. the elements of parameter vector.
  ## hess is the Hessian matrix function. It has the same arguments as func...
  ## but returns the Hessian matrix of the objective w.r.t. the elements of parameter vector. 
  ## If not supplied then newt finds an approximation to the Hessian by using our finite_diff function
  ## ... any arguments of func, grad and hess after the first (the parameter vector) are passed using this.
  ## tol is the convergence tolerance.
  ## fscale is a rough estimate of the magnitude of func near the optimum - used in convergence testing.
  ## maxit is the maximum number of Newton iterations to try before giving up.
  ## max.half is the maximum number of times a step is halved before...
  ## concluding that the step has failed to improve the objective.
  ## eps is the finite difference intervals to use when a Hessian function is not provided.
  
  
  ## newt returns a list containing:
  ## f: the value of the objective function at the minimum.
  ## theta: the value of the parameters at the minimum.
  ## iter: the number of iterations taken to reach the minimum.
  ## g: the gradient vector at the minimum (so the user can judge closeness to numerical zero).
  ## Hi: the inverse of the Hessian matrix at the minimum.
  
  iter <- 0 ## initialise count of number of iterations
  g <- grad(theta, ...) ## calculate the gradient vector for initial guess
  f <- func(theta, ...) ## evaluates objective function at initial guess
  
  ## Check that objective function and gradient are finite at initial guess
  if(!is.finite(func(theta, ...)) | sum(is.finite(grad(theta,...))) < length(theta)){
   stop("Objective or derivative in not finite at initial theta, please enter a different theta.") 
  }
  
  
  if(!is.null(hess)){
    hess0 <- hess(theta, ...) ## calculate hessian matrix for initial guess
  }else{
    ## Use the finite_diff function if hessian is not provided
    hess0 <- finite_diff(theta, grad,eps,...)
  }
  
  I <- diag(length(theta)) ## create appropriately sized identity to add to hessian if it is not positive definite
  
  while(iter < maxit & abs(max(g)) > (tol*(abs(f+fscale)))){
    ## execute an iteration of Newton's method if we have not reached the maximum number of iterations
    ## and if the gradient of the objective is above the tolerance level
    mult = 1e-8 ## initial multiple of identity to add to non-positive definite hessian
    mat_norm = norm(hess0, "F")  ## use the Frobenius norm to improve efficiency of perturbations
    while(inherits(try(chol(hess0),silent = TRUE), "try-error")){
      ## check if the hessian is positive definite by attempting a Cholesky decomposition
      ## while it is not positive definite, add increasing multiples of the identity matrix until it is
      hess0 <- hess0 + mult*mat_norm*I
      mult <- 10*mult      ## increase multiplier before next iteration
    }
    
    R <- chol(hess0) ## Cholesky decomposition of the hessian
    hess_inv<-chol2inv(R) ## use the decomposition to calculate the inverse
    
    delta <- -(hess_inv)%*%g  ## find the delta which minimises the quadratic approximation
    temp_theta <- theta + delta ## store new guess
    half_iter <- 1 ## initialise count for the number of half steps used
    
    while(func(theta, ...) < func(temp_theta, ...) & half_iter <= max.half){
      ## If our initial guess is far from the objective, we may overshoot
      ## if the updated theta gives a worse objective value, add delta/2 instead
      ## reduce step size by half on each iteration until a better objective is found
      temp_theta <- theta + delta*(1/2 ^ half_iter)
      if(!is.finite(func(temp_theta,...))){
        ## if the new half-step returns a non-finite value, reset temp_theta to theta
        ## this will simply move to the next half-step
        temp_theta <- theta ## re-enter the while loop to try next (smaller) half step
      }
      half_iter <- half_iter + 1 ## increase half_iter counter
    }## end while loop
    
    if(half_iter == max.half){
      ## alert user that the maximum number of half step were used without improvement
      warning("Max half steps reached without improving objective. Return most recent output.")
      ## break out of the larger while loop - objective won't be improved further
      break
    }
    
    # updating arguments for next iteration. Final values will be output
    theta <- temp_theta     
    f <- func(theta, ...)
    g <- grad(theta, ...)
    
    if(!is.null(hess)){
      hess0 <- hess(theta, ...) ## calculate hessian matrix for new guess
    }else{
      ## Use the finite_diff function if hessian is not provided
      hess0 <- finite_diff(theta, grad,eps,...)
    }

    iter <- iter + 1  ## update iteration counter
  } ## End of outer while loop
  
  if(iter == maxit){
    ## alert user that the maximum number of iterations has been reached without convergence
    warning("Max iterations reached without convergence. Return most recent output.")
  }
  # try to calculate hessian matrix inverse by cholesky decomposition
  R <- try(chol(hess0), silent=TRUE)
  Hi <- try(chol2inv(R), silent=TRUE)
  if(inherits(try(chol(hess0),silent = TRUE), "try-error")){
    # if we can not do this, alert user that the hessian is not +ve definite 
    warning("Hessian at minimum is not positive definite")
    Hi = NULL  ## In this case return NULL
  }
  output <- list(theta = theta, f=f, iter=iter, g=g, Hi=Hi) ## create output list
  return(output) 
}## End of newt function
