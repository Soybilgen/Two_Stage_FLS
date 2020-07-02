TVPss <- function(x,pars,T,nvar) {
  # Tt: State transition matrix
  cx=dim(x)[2]
  # if (cx==1) Tt<- pars[2*nvar+2*cx+1] else Tt<-matrix(diag(pars[(2*nvar+2*cx+1):(2*nvar+3*cx)]),ncol=cx)
  if (cx==1) Tt<- 1 else Tt<-matrix(diag(rep(1,cx)),ncol=cx)
  # Z_t: Measurement observation vector y_t=Z_r\times state
  xones=matrix(rep((x),nvar),ncol=nvar)
  Zt <- array(data=t(xones),dim=c(nvar,cx,T))
  # ct: measurement constant (this corresponds to intercept term in regression)
  ct <- matrix(pars[1:nvar])
  # dt: transition intercept terms
  dt <- matrix(rep(0,cx),ncol=1)
  # GGt: variance of measurement equation noise
  if (nvar==1) GGt<- matrix(pars[(nvar+1):(nvar+nvar)]^2) else GGt <- matrix(diag(pars[(nvar+1):(nvar+nvar)]^2),ncol=nvar)
  # HHt: variance of state noise
  if (cx==1) HHt=matrix(pars[(2*nvar+1):(2*nvar+cx)]^2) else HHt=diag(pars[(2*nvar+1):(2*nvar+cx)]^2)
  #  initial condtions
  a0 <-(pars[(2*nvar+cx+1):(2*nvar+2*cx)])
  if (cx==1) P0<- 1 else P0 <- diag(rep(1,cx))
  return(list(a0 = a0, P0 = P0, ct = ct, dt = dt, Zt = Zt, Tt = Tt, GGt = GGt,
              HHt = HHt))
}

objective <- function(theta,y,x,T,nvar) {
  yt=t(y)
  sp <- TVPss(x,theta,T,nvar)
  ans <- fkf(a0 = sp$a0, P0 = sp$P0, dt = sp$dt, ct = sp$ct, Tt = sp$Tt,
             Zt = sp$Zt, HHt = sp$HHt, GGt = sp$GGt, yt = yt)
  return(-ans$logLik)
}

lagmatrix<-function(data,lags){
  k=dim(data)[2]
  T=dim(data)[1]
  nlag=length(lags)
  dataf=matrix(rep(NaN,(k*nlag*T)),nrow=T)
  ii=0
  for (i in lags){
    ii=ii+1
    dataf[(i+1):T,((ii-1)*k+1):(ii*k)]=data[1:(T-i),]
  }
  return(dataf)
}

flskf <-  function (A, b, mu=1, ncap=length(b))
{
  m <- nrow (A)
  n <- ncol (A)
  X <- array (0,c(n,ncap))
  R <- array (0,c(n,n))
  diag (R) <- 1/mu
  a <- t(A[1,,drop=FALSE])
  phi <- c((1/mu) * crossprod(a) + 1)
  w <- (1/phi) * R %*% a
  rho <- b[1]
  X[,1] <- rho * w
  for (j in 2:ncap) {
    R = R - phi * tcrossprod(w)
    diag (R) <- diag (R) + 1/mu
    a <- t(A[j,,drop=FALSE])
    phi <- c(t(a) %*% R %*% a)
    phi <- phi + 1
    rho <- c(b[j] - crossprod (a, X[,j-1]))
    w <- (1/phi) * R %*% a
    X[,j] <- X[,j-1] + rho*w
  }
  X
}

fls <-  function (A, b, mu=1, ncap=length(b), smoothed=TRUE){
  m <- nrow (A)
  n <- ncol (A)
  M <- array (0,c(n,n,ncap))
  E <- array (0,c(n,ncap))
  X <- array (0,c(n,ncap))
  R <- matrix(0,n,n)
  diag(R) <- diag(R) + mu
  for (j in 1:ncap) {
    Z <- solve(qr(R + tcrossprod(A[j,]),LAPACK=TRUE),diag(1.0,n));
    M[,,j] <- mu*Z             # (5.7b)
    v <- b[j]*A[j,]
    if(j==1) p <- rep(0,n)
    else p <- mu*E[,j-1]
    w <- p + v
    E[,j] <- Z %*% w           # (5.7c)
    R <- -mu*mu*Z
    diag(R) <- diag(R) + 2*mu
  }
  # Calculate eqn (5.15) FLS estimate at ncap
  Q <- -mu*M[,,ncap-1]
  diag(Q) <- diag(Q) + mu
  Ancap <- A[ncap,,drop=FALSE]
  C <- Q + t(Ancap) %*% Ancap
  d <- mu*E[,ncap-1,drop=FALSE] + b[ncap]*t(Ancap)
  X[,ncap] <- C %*% d
  X[,ncap] <- solve(qr(C,LAPACK=TRUE),d)
  if (smoothed) {
    # Use eqn (5.16) to obtain smoothed FLS estimates for 
    # X[,1], X[,2], ..., X[,ncap-1]
    for (j in 1:(ncap-1)) {
      l <- ncap - j
      X[,l] <- E[,l] + M[,,l] %*% X[,l+1]
    }
  }
  else {
    X <- X[,ncap]
  }
  X
}

bbfls<-function(X,y,mu){
  Xt=t(X)
  k=ncol(X)
  N=nrow(X)
  G=matrix(rep(0,(k*N*N)),ncol=N)
  A=matrix(rep(0,(k*N*N*k)),ncol=(N*k))
  for (i in 1:nrow(X)){
    
    G[(1+((i-1)*k)):(i*k),i]=Xt[,i]
    if (i==1){
      As=Xt[,i]%*%t(Xt[,i])+mu*diag(k)
      A[(1+((i-1)*k)):(i*k),(1+((i)*k)):(2*i*k)]=-mu*diag(k)
      A[(1+((i-1)*k)):(i*k),(1+((i-1)*k)):(i*k)]=As
    } else if (i==N){
      As=Xt[,i]%*%t(Xt[,i])+mu*diag(k)
      A[(1+((i-1)*k)):(i*k),(1+((i-1)*k)):(i*k)]=As
      A[(1+((i-1)*k)):(i*k),(1+((i-2)*k)):((i-1)*k)]=-mu*diag(k)
    }else {
      As=Xt[,i]%*%t(Xt[,i])+2*mu*diag(k)
      A[(1+((i-1)*k)):(i*k),(1+((i-1)*k)):(i*k)]=As
      A[(1+((i-1)*k)):(i*k),(1+((i-2)*k)):((i-1)*k)]=-mu*diag(k)
      A[(1+((i-1)*k)):(i*k),(1+((i)*k)):((i+1)*k)]=-mu*diag(k)
    }
    A[(1+((i-1)*k)):(i*k),]
  }
  bfls0=solve(A)%*%G%*%y
  
  bfls=matrix(rep(0,N*k),ncol=k)
  res=y
  for (j in 1:k){
    bfls[,j]=bfls0[(seq((1+(j-1)),N*k,k))]
    res=res-X[,j]*bfls[,j]
  }
  
  sigma2hatbfls0=sum(res^2)/(N-k)
  
  stdbfls0=sqrt(diag(sigma2hatbfls0*(solve(A)%*%G%*%t(G)%*%t(solve(A)))))
  stdbfls=matrix(rep(0,N*k),ncol=k)
  for (j in 1:k){
    stdbfls[,j]=stdbfls0[(seq((1+(j-1)),N*k,k))]
  }
  
  confintbfls.low=bfls-qt(0.975,N-k)*stdbfls
  confintbfls.up=bfls+qt(0.975,N-k)*stdbfls
  
  list(bfls,stdbfls,confintbfls.low,confintbfls.up)
}