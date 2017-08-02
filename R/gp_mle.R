# =========================== gp_mle ===========================

gp_mle <- function(gp_data) {
  # Maximum likelihood estimation for the generalized Pareto distribution
  #
  # Performs maximum likelihood estimation for the generalized Pareto
  # distribution.  Uses the function \code{gpdmle} associated with
  # Grimshaw (1993), which returns MLEs of sigma and k = - \xi.
  #
  # Grimshaw, S. D. (1993) Computing Maximum Likelihood Estimates
  #   for the Generalized Pareto Distribution.  Technometrics, 35(2), 185-191.
  #   and Computing (1991) 1, 129-133. http://dx.doi.org/10.1007/BF01889987.
  #
  # Args:
  #   gp_data : A numeric vector containing positive values, assumed to be a
  #             random sample from a generalized Pareto distribution.
  #
  # Returns:
  #   A list with components
  #     mle  : A numeric vector.  MLEs of GP parameters sigma and xi.
  #     nllh : A numeric scalar.  The negated log-likelihood at the MLE.
  #
  # Call Grimshaw (1993) function, note: k is -xi, a is sigma
  pjn <- grimshaw_gp_mle(gp_data)
  temp <- list()
  temp$mle <- c(pjn$a, -pjn$k)  # mle for (sigma,xi)
  sc <- rep(temp$mle[1], length(gp_data))
  xi <- temp$mle[2]
  temp$nllh <- sum(log(sc)) + sum(log(1 + xi * gp_data / sc) * (1 / xi + 1))
  gp_info <- gp_obs_info(temp$mle, gp_data)
  temp$cov <- solve(gp_info)
  temp$se <- sqrt(diag(temp$cov))
  return(temp)
}

# =========================== gp_obs_info ===========================

gp_obs_info <- function(gp_pars, y) {
  # Observed information for the generalized Pareto distribution
  #
  # Calculates the observed information matrix for a random sample \code{y}
  # from the generalized Pareto distribution, i.e. the negated Hessian matrix
  # of the generalized Pareto log-likelihood, evaluated at \code{gp_pars}.
  #
  # Args:
  #   gp_pars : A numeric vector. Parameters sigma and xi of the
  #   generalized Pareto distribution.
  #   y       : A numeric vector. A sample of positive data values.
  #
  # Returns:
  #   A 2 by 2 numeric matrix.  The observed information matrix.
  #
  s <- gp_pars[1]
  x <- gp_pars[2]
  i <- matrix(NA, 2, 2)
  i[1,1] <- -sum((1 - (1 + x) * y * (2 * s + x * y) / (s + x * y) ^ 2) / s ^ 2)
  i[1,2] <- i[2,1] <- -sum(y * (1 - y / s) / (1 + x * y / s) ^ 2 / s ^ 2)
  i[2,2] <- sum(2 * log(1 + x * y / s) / x ^ 3 - 2 * y / (s + x * y) / x ^ 2 -
                  (1 + 1 / x) * y ^ 2 / (s + x * y) ^ 2)
  return(i)
}

# =========================== grimshaw_gp_mle ===========================

grimshaw_gp_mle <- function(x) {
  #  Argument for function:
  #
  #  x     the sample values from the GPD
  #
  #
  #
  #  Returned from the function:
  #
  #  k       the mle of k
  #  a       the mle of a
  #
  n<-length(x)
  xq<-sort(x)
  xbar<-mean(x)
  sumx2<-sum(x^2)/n
  x1<-xq[1]
  xn<-xq[n]
  #
  #  Find the local maxima/minima of the likelihood.
  #
  #
  #  Initialize epsilon as the accuracy criterion
  #
  epsilon<-10^(-6)/xbar
  #
  #  The local maxima/minima must be found numerically by
  #  finding the zero(s) of h().
  #
  #  Algorithm for finding the zero(s) of h().
  #
  #
  #  Any roots that exist must be within the interval (lobnd,hibnd).
  #
  lobnd<-2*(x1-xbar)/x1^2
  if(lobnd>=0){
    lobnd<- -epsilon
  }
  hibnd<-(1/xn)-epsilon
  if(hibnd<=0){
    hibnd<-epsilon
  }
  #
  #  If h''(0)>0, look for one negative and one positive zero of h().
  #  If h''(0)<0, look for two negative and two positive zeros of h().
  #
  secderiv<-sumx2-2*xbar^2  #{ Evaluate h''(0). }
  if(secderiv>0){
    #
    #
    #  Look for one negative and one positive zero of h().
    #
    #
    thzeros<-cbind(c(0,0),c(0,0))
    nzeros<-2
    #
    #  Begin with the initial value at lobnd.
    #
    hlo<-(1+sum(log(1-lobnd*x))/n)*(sum(1/(1-lobnd*x))/n)-1
    if(hlo<0){
      thlo<-lobnd       #{  Orient the search so h(thlo)<0  }
      thhi<- -epsilon
    }
    else{
      thlo<- -epsilon
      thhi<-lobnd
    }
    thzero<-lobnd    #{  Initial value for modified Newton-Raphson is lobnd. }
    dxold<-abs(thhi-thlo)
    dx<-dxold
    temp1<-sum(log(1-thzero*x))/n
    temp2<-sum(1/(1-thzero*x))/n
    temp3<-sum(1/(1-thzero*x)^2)/n
    h<-(1+temp1)*(temp2)-1
    hprime<-(temp3-temp2^2-temp1*(temp2-temp3))/thzero

    #
    #  Newton-Raphson Algorithm to find the zero of the function h()
    #  for a given initial starting point.
    j<-1
    maxiter<-100  #{Maximum number of mod. Newton-Raphson iterations}
    while(j<=maxiter){
      #
      #  Determine whether it is better to use Bisection (if N-R is
      #  out of range or not decreasing fast enough) or Newton-Raphson.
      #
      c1<-(((thzero-thhi)*hprime-h)*((thzero-thlo)*hprime-h)>=0)
      c2<-(abs(2*h)>abs(dxold*hprime))
      if(c1+c2>=1){
        dxold<-dx
        dx<-(thhi-thlo)/2
        thzero<-thlo+dx
        if(thlo==thzero){       #{Change in root is negligible}
          j<-1000
        }
      }
      else{
        dxold<-dx
        dx<-h/hprime
        temp<-thzero
        thzero<-thzero-dx
        if(temp==thzero){       #{Change in root is negligible}
          j<-1001
        }
      }
      #
      #  Determine if convergence criterion is met.
      #
      if(abs(dx)<epsilon*abs(thlo+thhi)/2){
        j<-999
      }
      temp1<-sum(log(1-thzero*x))/n
      temp2<-sum(1/(1-thzero*x))/n
      temp3<-sum(1/(1-thzero*x)^2)/n
      h<-(1+temp1)*(temp2)-1
      hprime<-(temp3-temp2^2-temp1*(temp2-temp3))/thzero
      if(h<0){            #{Maintain the bracket on the root}
        thlo<-thzero
      }
      else{
        thhi<-thzero
      }
      j<-j+1
    }

    if(j>maxiter+1){
      thzeros[1,]<-cbind(thzero,j)
    }

    #
    #  Begin with the initial value at hibnd.
    #
    hlo<-(1+sum(log(1-epsilon*x))/n)*(sum(1/(1-epsilon*x))/n)-1
    if(hlo<0){
      thlo<-epsilon       #{  Orient the search so h(thlo)<0  }
      thhi<-hibnd
    }
    else{
      thlo<-hibnd
      thhi<-epsilon
    }
    thzero<-hibnd    #{  Initial value for modified Newton-Raphson is hibnd. }
    dxold<-abs(thhi-thlo)
    dx<-dxold
    temp1<-sum(log(1-thzero*x))/n
    temp2<-sum(1/(1-thzero*x))/n
    temp3<-sum(1/(1-thzero*x)^2)/n
    h<-(1+temp1)*(temp2)-1
    hprime<-(temp3-temp2^2-temp1*(temp2-temp3))/thzero

    #
    #  Newton-Raphson Algorithm to find the zero of the function h()
    #  for a given initial starting point.
    j<-1
    maxiter<-100  #{Maximum number of mod. Newton-Raphson iterations}
    while(j<=maxiter){
      #
      #  Determine whether it is better to use Bisection (if N-R is
      #  out of range or not decreasing fast enough) or Newton-Raphson.
      #
      c1<-(((thzero-thhi)*hprime-h)*((thzero-thlo)*hprime-h)>=0)
      c2<-(abs(2*h)>abs(dxold*hprime))
      if(c1+c2>=1){
        dxold<-dx
        dx<-(thhi-thlo)/2
        thzero<-thlo+dx
        if(thlo==thzero){       #{Change in root is negligible}
          j<-1000
        }
      }
      else{
        dxold<-dx
        dx<-h/hprime
        temp<-thzero
        thzero<-thzero-dx
        if(temp==thzero){       #{Change in root is negligible}
          j<-1001
        }
      }
      #
      #  Determine if convergence criterion is met.
      #
      if(abs(dx)<epsilon*abs(thlo+thhi)/2){
        j<-999
      }
      temp1<-sum(log(1-thzero*x))/n
      temp2<-sum(1/(1-thzero*x))/n
      temp3<-sum(1/(1-thzero*x)^2)/n
      h<-(1+temp1)*(temp2)-1
      hprime<-(temp3-temp2^2-temp1*(temp2-temp3))/thzero
      if(h<0){            #{Maintain the bracket on the root}
        thlo<-thzero
      }
      else{
        thhi<-thzero
      }
      j<-j+1
    }


    if(j>maxiter+1){
      thzeros[2,]=cbind(thzero,j)
    }
  }

  else{
    #
    #
    #  Look for two negative and two positive zeros of h().
    #
    #
    thzeros<-matrix(rep(0,8),ncol=2)
    nzeros<-4
    #
    #  Begin with the initial value at lobnd.
    #
    hlo<-(1+sum(log(1-lobnd*x))/n)*(sum(1/(1-lobnd*x))/n)-1
    if(hlo<0){
      thlo<-lobnd       #{  Orient the search so h(thlo)<0  }
      thhi<- -epsilon
    }
    else{
      thlo<- -epsilon
      thhi<-lobnd
    }
    thzero<-lobnd    #{  Initial value for modified Newton-Raphson is lobnd. }
    dxold<-abs(thhi-thlo)
    dx<-dxold
    temp1<-sum(log(1-thzero*x))/n
    temp2<-sum(1/(1-thzero*x))/n
    temp3<-sum(1/(1-thzero*x)^2)/n
    h<-(1+temp1)*(temp2)-1
    hprime<-(temp3-temp2^2-temp1*(temp2-temp3))/thzero

    #
    #  Newton-Raphson Algorithm to find the zero of the function h()
    #  for a given initial starting point.
    j<-1
    maxiter<-100  #{Maximum number of mod. Newton-Raphson iterations}
    while(j<=maxiter){
      #
      #  Determine whether it is better to use Bisection (if N-R is
      #  out of range or not decreasing fast enough) or Newton-Raphson.
      #
      c1<-(((thzero-thhi)*hprime-h)*((thzero-thlo)*hprime-h)>=0)
      c2<-(abs(2*h)>abs(dxold*hprime))
      if(c1+c2>=1){
        dxold<-dx
        dx<-(thhi-thlo)/2
        thzero<-thlo+dx
        if(thlo==thzero){       #{Change in root is negligible}
          j<-1000
        }
      }
      else{
        dxold<-dx
        dx<-h/hprime
        temp<-thzero
        thzero<-thzero-dx
        if(temp==thzero){       #{Change in root is negligible}
          j<-1001
        }
      }
      #
      #  Determine if convergence criterion is met.
      #
      if(abs(dx)<epsilon*abs(thlo+thhi)/2){
        j<-999
      }
      temp1<-sum(log(1-thzero*x))/n
      temp2<-sum(1/(1-thzero*x))/n
      temp3<-sum(1/(1-thzero*x)^2)/n
      h<-(1+temp1)*(temp2)-1
      hprime<-(temp3-temp2^2-temp1*(temp2-temp3))/thzero
      if(h<0){            #{Maintain the bracket on the root}
        thlo<-thzero
      }
      else{
        thhi<-thzero
      }
      j<-j+1
    }

    if(j>maxiter+1){
      thzeros[1,]<-cbind(thzero,j)
    }
    #
    #  Look at the derivative to determine where the second root lies.
    #   If h'(0)>0, second root lies between thzero and -epsilon.
    #   If h'(0)<0, second root lies between lobnd and thzero.
    #
    temp1<-sum(log(1-thzero*x))/n
    temp2<-sum(1/(1-thzero*x))/n
    temp3<-sum(1/(1-thzero*x)^2)/n
    hprime<-(temp3-temp2^2-temp1*(temp2-temp3))/thzero
    if(hprime>0){
      #
      #  h'(0)>0, so the second zero lies between thzero and -epsilon.
      #
      #
      #  Establish Initial Values.
      #
      thlo<-thzero
      thhi<- -epsilon
      thzero<-thhi
      dx<-thlo-thhi

      j<-1
      maxiter<-100  #{Maximum number of bisection iterations}
      while(j<=maxiter){
        dx<-.5*dx
        thmid<-thzero+dx
        hmid<-(1+sum(log(1-thmid*x))/n)*(sum(1/(1-thmid*x))/n)-1
        if(hmid<0){
          thzero<-thmid
        }
        if(hmid==0){   #{Zero of h() has been found}
          j<-999
        }
        if(abs(dx)<epsilon*abs(thlo+thhi)/2){
          j<-999
        }
        j<-j+1
      }

      if(j>maxiter+1){
        thzeros[2,]<-cbind(thzero,j)
      }
    }
    else{
      #
      #  h'(0)<0, so the second zero lies between lobnd and thzero.
      #
      #
      #  Establish Initial Values.
      #
      thlo<-lobnd
      thhi<-thzero
      thzero<-thlo
      dx<-thhi-thlo

      j<-1
      maxiter<-100  #{Maximum number of bisection iterations}
      while(j<=maxiter){
        dx<-.5*dx
        thmid<-thzero+dx
        hmid<-(1+sum(log(1-thmid*x))/n)*(sum(1/(1-thmid*x))/n)-1
        if(hmid<0){
          thzero<-thmid
        }
        if(hmid==0){   #{Zero of h() has been found}
          j<-999
        }
        if(abs(dx)<epsilon*abs(thlo+thhi)/2){
          j<-999
        }
        j<-j+1
      }

      if(j>maxiter+1){
        thzeros[2,]<-cbind(thzero,j)
      }
    }
    #
    #  Begin with the initial value at hibnd.
    #
    hlo<-(1+sum(log(1-epsilon*x))/n)*(sum(1/(1-epsilon*x))/n)-1
    if(hlo<0){
      thlo<-epsilon       #{  Orient the search so h(thlo)<0  }
      thhi<-hibnd
    }
    else{
      thlo<-hibnd
      thhi<-epsilon
    }
    thzero<-hibnd    #{  Initial value for modified Newton-Raphson is hibnd. }
    dxold<-abs(thhi-thlo)
    dx<-dxold
    temp1<-sum(log(1-thzero*x))/n
    temp2<-sum(1/(1-thzero*x))/n
    temp3<-sum(1/(1-thzero*x)^2)/n
    h<-(1+temp1)*(temp2)-1
    hprime<-(temp3-temp2^2-temp1*(temp2-temp3))/thzero

    #
    #  Newton-Raphson Algorithm to find the zero of the function h()
    #  for a given initial starting point.
    j<-1
    maxiter<-100  #{Maximum number of mod. Newton-Raphson iterations}
    while(j<=maxiter){
      #
      #  Determine whether it is better to use Bisection (if N-R is
      #  out of range or not decreasing fast enough) or Newton-Raphson.
      #
      c1<-(((thzero-thhi)*hprime-h)*((thzero-thlo)*hprime-h)>=0)
      c2<-(abs(2*h)>abs(dxold*hprime))
      if(c1+c2>=1){
        dxold<-dx
        dx<-(thhi-thlo)/2
        thzero<-thlo+dx
        if(thlo==thzero){       #{Change in root is negligible}
          j<-1000
        }
      }
      else{
        dxold<-dx
        dx<-h/hprime
        temp<-thzero
        thzero<-thzero-dx
        if(temp==thzero){       #{Change in root is negligible}
          j<-1001
        }
      }
      #
      #  Determine if convergence criterion is met.
      #
      if(abs(dx)<epsilon*abs(thlo+thhi)/2){
        j<-999
      }
      temp1<-sum(log(1-thzero*x))/n
      temp2<-sum(1/(1-thzero*x))/n
      temp3<-sum(1/(1-thzero*x)^2)/n
      h<-(1+temp1)*(temp2)-1
      hprime<-(temp3-temp2^2-temp1*(temp2-temp3))/thzero
      if(h<0){            #{Maintain the bracket on the root}
        thlo<-thzero
      }
      else{
        thhi<-thzero
      }
      j<-j+1
    }

    if(j>maxiter+1){
      thzeros[3,]<-cbind(thzero,j)
    }
    #
    #  Look at the derivative to determine where the second root lies.
    #   If h'(0)>0, second root lies between thzero and hibnd.
    #   If h'(0)<0, second root lies between epsilon and thzero.
    #
    temp1<-sum(log(1-thzero*x))/n
    temp2<-sum(1/(1-thzero*x))/n
    temp3<-sum(1/(1-thzero*x)^2)/n
    hprime<-(temp3-temp2^2-temp1*(temp2-temp3))/thzero
    if(hprime>0){
      #
      #  h'(0)>0, so the second zero lies between thzero and hibnd.
      #
      #
      #  Establish Initial Values.
      #
      thlo<-thzero
      thhi<-hibnd
      thzero<-thhi
      dx<-thlo-thhi

      j<-1
      maxiter<-100  #{Maximum number of bisection iterations}
      while(j<=maxiter){
        dx<-.5*dx
        thmid<-thzero+dx
        hmid<-(1+sum(log(1-thmid*x))/n)*(sum(1/(1-thmid*x))/n)-1
        if(hmid<0){
          thzero<-thmid
        }
        if(hmid==0){   #{Zero of h() has been found}
          j<-999
        }
        if(abs(dx)<epsilon*abs(thlo+thhi)/2){
          j<-999
        }
        j<-j+1
      }

      if(j>maxiter+1){
        thzeros[4,]<-cbind(thzero,j)
      }
    }
    else{
      #
      #  h'(0)<0, so the second zero lies between epsilon and thzero.
      #
      #
      #  Establish Initial Values.
      #
      thlo<-epsilon
      thhi<-thzero
      thzero<-thlo
      dx<-thhi-thlo

      j<-1
      maxiter<-100  #{Maximum number of bisection iterations}
      while(j<=maxiter){
        dx<-.5*dx
        thmid<-thzero+dx
        hmid<-(1+sum(log(1-thmid*x))/n)*(sum(1/(1-thmid*x))/n)-1
        if(hmid<0){
          thzero<-thmid
        }
        if(hmid==0){   #{Zero of h() has been found}
          j<-999
        }
        if(abs(dx)<epsilon*abs(thlo+thhi)/2){
          j<-999
        }
        j<-j+1
      }

      if(j>maxiter+1){
        thzeros[4,]<-cbind(thzero,j)
      }
    }
  }
  #
  #  Of the candidate zero(s) of h(), determine whether they correspond
  #  to a local maximum or minimum of the log-likelihood.
  #
  #  Eliminate any non-convergent roots}
  thetas<-thzeros[thzeros[,2]>maxiter+1,]
  nzeros<-nrow(thetas)
  proll<-rep(0,nzeros)
  mles<-matrix(rep(0,4*nzeros),ncol=4)
  i<-1
  while(i<=nzeros){
    temp1<-sum(log(1-thetas[i,1]*x))
    mles[i,1]<- -temp1/n
    mles[i,2]<-mles[i,1]/thetas[i,1]
    mles[i,3]<- -n*log(mles[i,2])+(1/mles[i,1]-1)*temp1
    mles[i,4]<-999
    i<-i+1
  }
  ind<-1:length(mles[,4])
  ind<-ind[mles[,4]==999]
  if(sum(ind)==0){   #{ Check to see if there are any local maxima. }
    nomle<-0          #{ If not, there is no mle. }
  }
  else{
    nomle<-1
  }
  if(nomle!=0){
    mles<-mles[ind,]
    nmles<-nrow(mles)
    #
    #  Add the boundary value where k=1 to the candidates for the
    #  maximum of the log-likelihood.
    #
    mles<-rbind(mles,c(1,xn,-n*log(xn),999))
    nmles<-nmles+1
    #
    #  Choose of the candidate mles whichever has the largest
    #  log-likelihood.
    #
    maxlogl<-max(mles[,3])
    ind<-order(mles[,3])
    ind<-ind[nmles]
    k<-mles[ind,1]
    #  label(k,'GPD mle of k')
    a<-mles[ind,2]
    #  label(a,'GPD mle of a')
  }
  else{
    #
    #  No Maximum Likelihood Estimators were found.
    #
    k<-NA
    a<-NA
  }
  return(list(k=k,a=a))
}
