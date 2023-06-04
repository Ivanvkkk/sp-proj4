## name ; university user name ; 
## Shiyu Cai ; s2307982
## Yifan Wu ; s2316499
## Yifan Cai ; s2301729

## github repo: https://github.com/eavancai/proj4.git

## Yifan Wu: the code work from begin to calculate step for theta and debugging
## Shiyu Cai: the code work after calculating step for theta
## Yifan Cai: test the function, write outlines and comments and debugging

## outline:(The process of Newton's method for minimizing function)

## Input the function to minimize, parameters and etc.
## Check if the gradient is finite:
## If the gradient is not finite, we should stop and reset the inputs.
## If the gradient is finite, we should check if the hessian function provided:
## If not, we use find_hessian function to generate hessian matrix.

## We set two limitations for calculating the minimization:
## The iterations should lower than or equal to maximun number of iterations
## We need to judge the converge 

## We should check whether the hessian matrix is positive value:
## if not, we need convert it to positive value

## Use the step formula generated form Taylor's theorem 
## We can get the next theta value after knowing the step 

## To ensure that we are on the way to find the minimum:
## Check if the function value at next theta is less than the one 
## at the previous theta value:
## If the function value becomes larger, we should to change step by halving.
## We should also set limitation to the halving times.
## We will repeat the process of calculate the step value and next few steps 
## until much close to the minimum of function.
## To get much close to the minimum, we also need to check few things and 
## set limitations.

# test function
rb <- function(th, k = 2) {
  k * (th[2]-th[1]^2)^2 + (1 - th[1])^2
}

#test gradient
gb <- function(th, k = 2) {
  c(-2*(1-th[1])-k*4*th[1]*(th[2]-th[1]^2), k*2*(th[2]-th[1]^2))
}

## test hessian
hb <- function(th,k=2) {
  h <- matrix(0,2,2)
  h[1,1] <- 2-k*2*(2*(th[2]-th[1]^2) - 4*th[1]^2)
  h[2,2] <- 2*k
  h[1,2] <- h[2,1] <- -4*k*th[1]
  h 
}

find_hessian <- function(theta, grad, eps = 1e-6){
  '
  find_hessian(theta, grad, eps)
    This function is to calculate the hessian matrix for theta with 
    gradient funcion(grad).
  
  input:
    theta(vector): a vector of values for the parameters
    grad: the gradient function
    eps(float): the finite intervals 
    
  output: 
    Hfd(matrix): the generated hessian matrix  
  '
  theta_len <- length(theta) ## find the length of theta
  gll0 <- grad(theta) ## find the gradient of the initial theta
  Hfd <- matrix(0, theta_len, theta_len) ## initialize the matrix
  
  for (i in 1:length(theta)) { ## loop over parameters
    th1 <- theta; th1[i] <- th1[i] + eps ## increase theta by eps
    gll1 <- grad(th1) ## compute resulting gradient
    Hfd[i,] <- (gll1 - gll0)/eps ## approximate second derivative
  }
  Hfd <- (t(Hfd)+Hfd)/2 ## ensure the matrix is symmetric
  Hfd
}

newt <- function(theta, func, grad, hess=NULL,..., tol = 1e-8, fscale=1, 
                 maxit=100, max.half=20, eps=1e-6){
  '
  newt(theta, func, grad, hess=NULL,..., tol, fscale, maxit, max.half, eps)
    This function is to implement Newton method for minimizing of functions.
  
  input:
    theta(vector): a vector of values for the parameters
    func: the function to be minimised
    grad: the gradient function
    hess: the hessian function
    ...: other arguments of func, grad and hess 
    tol: the convegence tolerance
    fscale: a rough estimate of the magnitude of func near the optimum 
    maxit(int): the maximum number of Newton iterations to try before giving up
    max.half(int): the maximum number of times a step should be halved before 
    concluding that the step has failed to
    eps(float): the finite intervals 
    
  output: 
    a list of :
      f: the value of the objective function at the minimum
      theta: the value of the parameters at the minimum
      iter: the number of iterations taken to reach the minimum
      g: the gradient vector at the minimum 
      Hi: the inverse of the Hessian matrix at the minimum 
  '
  
  ## if the objective function are not finite then stop the function
  ifelse(is.infinite(func(theta, ...)), 
         stop("the objective are not finite at the initial theta"),'')
  
  ## if the derivative function are not finite then stop the function
  ifelse(is.infinite(grad(theta)), 
         stop("the derivatives are not finite at the initial theta"),'')
  
  iter = 0 ## initialize the iteration 
  
  ## Jump out of the loop in both cases:
  ## 1.the number of iteration is greater than the maximum number of iterations
  ## 2.the convergence is reached
  while (iter <= maxit ){
    if(abs(norm(grad(theta, ...),type='2')) <= 
       (abs(norm(func(theta, ...),type='2'))+fscale)*tol){
      break
    }
    ## determine if the hessian matrix is empty
    if (is.null(hess) == TRUE){
    ## when hessian is not supplied, use find_hessian function to calculate it
      hess_value = find_hessian(theta, grad, eps) 
    }else{
      hess_value = hess(theta, ...)
      hess_value = (t(hess_value) + hess_value)/2
    }
    
    ## try cholesky decomposition to judge if the matrix is positive definite
    ## turn the non-positive definite matrices into positive definite
    k = 1e-8
    I = diag(length(theta))
    while (TRUE){
      if (typeof(try(chol(hess_value + k*I),TRUE)) == "double"){
        break
      }
      k = k * 10
    }
    ## use the new positive definite matrix
    hess_value <- hess_value + k*I 
    
    ## compute step according to Taylor’s theorem
    step_theta <- -1 * chol2inv(chol(hess_value)) %*% grad(theta, ...)
    ## update theta
    theta_new <- theta + step_theta
    
    ## set initial times of halving step_theta
    time_of_half <- 0
    ##half the step if below 2 conditions occur within max.half times:
    ##1.func(theta_new) > func(theta)
    ##2.func(theta_new) is not finite
    while (time_of_half < max.half) {
      if (func(theta_new, ...) > func(theta, ...) | func(theta_new, ...) 
          == Inf | func(theta_new, ...) == -Inf){
        step_theta <- step_theta / 2
        theta_new <- theta + step_theta
        time_of_half <- time_of_half + 1
      }
      else break
    }
    ## If the step fails to reduce the objective despite 
    ## trying max.half step halvings
    if (func(theta_new, ...) > func(theta, ...) | func(theta_new, ...) 
        == Inf | func(theta_new, ...) == -Inf){
      warning('the step fails to reduce the objective despite 
              trying max.half step halvings')
    }
    ##update theta
    theta <- theta_new
    ##update times of iteration
    iter <- iter + 1
  }
  ifelse(iter == maxit & abs(grad(theta, ...)) > 
           ((func(theta, ...) + fscale)*tol), 
         warning('maxit is reached without convergence'), '')
  
  if (is.null(hess) == TRUE){
    hess_value = find_hessian(theta, grad)
  }else{
    hess_value = hess(theta)
    hess_value = (t(hess_value) + hess_value)/2
  }
  
  ## try cholesky decomposition to judge if the matrix is positive definite
  if (typeof(try(chol(hess_value),TRUE)) != "double"){
    stop(' Hessian is not positive definite')
  }
  
  list(f = func(theta, ...), theta = theta, iter = iter, 
       g = grad(theta, ...), Hi = chol2inv(chol(hess_value)))
  
}