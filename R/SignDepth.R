#### KSign ####

#' \code{K}-sign depth for univariate data
#'
#' C++ and R implementations of the algorithm for
#' the computation of the \code{K}-sign depth for
#' one-dimensional data (residuals). The C++ implementation
#' is the fastest one that should be used in practice; its
#' complexity is O(n). For comparison,
#' also the (very slow) naive R version of complexity O(n^K)
#' is provided.
#'
#' @param x A vector of -1 and +1 of length \code{n}. These
#' are typically the signs of the residuals (-1 for negative sign,
#' +1 for positive sign).
#'
#' @param K Positive integer, the parameter of the depth.
#'
#' @param naive Logical indicator whether the naive implementation
#' in R should be used. By default set to \code{FALSE}. Use \code{FALSE}
#' in all situations, the naive regime is only for checking and comparison.
#'
#' @return If \code{naive=FALSE}:
#' A vector of length \code{2K} of number of occurences
#' of alternating patterns in the sequence \code{x}. The vector
#' is named \code{+j} or \code{-j}, where the initial sign
#' stands for the starting sign of the alternating sequence
#' (-1 for - and +1 for +), and \code{j} is the length of the sequence.
#' \code{+3} thus means the number of triples of the form \code{+-+}
#' in the sequence \code{x} etc. The final \code{K}-sign depth is
#' the sum of the final two entries of this vector (\code{+K} and \code{-K})
#' divided by \code{choose(n,K)}, where \code{n} is the length of \code{x}.
#' The function also provides all \code{j}-sign depths for \code{j<K};
#' one just needs to sum the two corresponding entries in the vector.
#'
#' If \code{naive=TRUE}: Total number of \code{K}-sign alternating
#' sequences. This corresponds to the sum of the two elements in the
#' last row of the output matrix if \code{naive=FALSE}.
#'
#' @examples
#' n = 500
#' x = sign(rnorm(n))
#' K = 3
#' sum(KSign(x,K)[K,])/choose(n,K) # K-sign depth
#'
#' n = 50
#' x = sign(rnorm(n))
#' K = 3
#' sum(KSign(x,K)[K,])==KSign(x,K,naive=TRUE)

KSign = function(x,K,naive = FALSE){
  if(any(x^2!=1)) stop("x must be a vector of +1 and -1")
  if(naive){
    # find all K-element subsets of x and check condition one-by-one
    n = length(x)
    C = combn(n,K)
    m = 0
    for(i in 1:ncol(C)){
      if(all(x[C[-1,i]]!=x[C[-K,i]])) m = m+1
    }
    return(m)
  } else {
    res = altern_Craw2(x, K)
    res = res[-1,]
    colnames(res) = c("last-","last+")
    rownames(res) = 1:K
    return(res)
  }
}

#### component Ldepth ####

#' Component \code{L}-depths for multivariate data
#'
#' C++ implementations of the algorithm for
#' the computation of (both simplified and full) component
#' \code{L}-depth for multi-dimensional data (residuals).
#' The complexity of the algorithms is O(n).

#' @param X The bivariate (or multivariate) dataset represented
#' as a numerical matrix with \code{p} columns and \code{n} rows.
#' The rows of \code{X} are typically the residuals.
#'
#' @param L Positive integer, parameter of the depth.
#'
#' @param type Type of the depth to compute. Can take values \code{simplified}
#' or \code{full}. By default, \code{simplified} is computed.
#'
#' @return The depth value; a numerical value between 0 and 1.
#'
#' @seealso \link{KSign} for the function called in the case of a full
#' version of the depth.
#'
#' @examples
#' n = 500
#' p = 3
#' X = matrix(rnorm(p*n),ncol=p)
#' L = 2
#'
#' componentLDepth(X,L,"simplified")*(n-L)
#' # the same evaluated component-wise
#' for(j in 1:p) print(componentLDepth(X[,j],L,"simplified")*(n-L))
#'
#' componentLDepth(X,L,"full")*choose(n,L+1)
#' # the same evaluated component-wise
#' for(j in 1:p) print(componentLDepth(X[,j],L,"full")*choose(n,L+1))

componentLDepth = function(X, L=1, type=c("simplified", "full")){
  K = L
  type = match.arg(type)
  if(is.vector(X)) X = matrix(X,ncol=1)
  p = ncol(X)
  n = nrow(X)
  X = sign(X)
  if(type=="simplified"){
    res = Inf
    for(j in 1:p){
      res = min(res,compKdepth(X[,j],L))
    }
  return(res/(n-L))
  }
  if(type=="full"){
    res = Inf
    for(j in 1:p) res = min(res,sum(KSign(X[,j],L+1)[L+1,]))
  return(res/choose(n,L+1))
  }
}

#### simplLSimplex ####

#' Simplified \code{L}-simplex depth for bivariate data
#'
#' C++ and R implementations of the algorithm for
#' the computation of the simplified \code{L}-simplex depth for
#' two-dimensional data (residuals). The C++ implementation
#' is the fastest one that should be used in practice; the complexity
#' of both algorithms is O(n). For comparison,
#' also two versions of (slower) naive R versions are provided.
#'
#' @param X The bivariate dataset represented as a numerical matrix
#' with two columns and \code{n} rows. The rows of \code{X} are
#' typically the residuals.
#'
#' @param L Positive integer 1 or 2, parameter of the depth. If
#' \code{L} is 1, only the simplified 1-simplex depth is computed. If \code{L}
#' is 2, both simplified 1-simplex depth and simplified 2-simplex depth are
#' computed.
#'
#' @param method String with the method used to compute the depth. Possible
#' values are \code{"Cpp"} for the fast C++ implementation, \code{"naive"}
#' for a slower R implementation of the same algorithm,
#' or \code{"naive0"} for the simplest
#' implementation of a loop via \code{depth.simplicial} from package
#' \code{ddalpha}. By default set to \code{"Cpp"}. Use \code{"Cpp"}
#' in all situations, the naive regime is only for checking and comparison.
#'
#' @param echo Logical indicator whether the combinations of triples
#' of consecutive indices that result in triangles containing the origin
#' should be given. By default set to \code{FALSE}.
#'
#' @return A vector of length \code{2}. In the first element we have
#' the relative number of consecutive triangles (\code{L=1}) formed from
#' the rows of  \code{X} that contain the origin, and in the second
#' element we have the relative number of consecutive pairs of triangles
#' (\code{L=2}) that both contain the origin. Both numbers lie between
#' 0 and 1.
#'
#' @seealso \link{fullLSimplex} for a C++ implementation of the full
#' version of the 2-simplex depth, and \link{LSimplex} for a wrapper
#' function that can call both the full version of the depth
#' (\link{fullLSimplex}), and the simplified version of the depth
#' (\link{simplLSimplex}).
#'
#' @examples
#' n = 500
#' X = matrix(rnorm(2*n),ncol=2)
#' L = 2
#' simplLSimplex(X,L)

simplLSimplex = function(X,L=2,method=c("Cpp","naive0","naive1"),
                    echo=FALSE){
  K = L
  method = match.arg(method)
  if(method=="naive0"){
    n = nrow(X)
    D1 = 0
    D2 = 0
    for(i in 1:(n-2)){
      if(ddalpha::depth.simplicial(x=c(0,0),data=X[i:(i+2),], exact = T)){
        D1 = D1 + 1
        if(echo) print(c(i,i+1,i+2))
        if((K==2)&(i<=n-3)){
          D2 = D2 + ddalpha::depth.simplicial(x=c(0,0),data=X[(i+1):(i+3),], exact = T)
        }
      }
    }
    return(c(D1,D2)/c(n-2,n-3))
  }
  XR = X[,2] / X[,1]
  if(any(is.nan(XR))) warning("Some elements of X are numerically equal to 0,
  computation might be incorrect")
  sgn = sign(X[,1])
  n = length(sgn)
  if(method=="naive1"){
    rnk = rank(XR) # rank here
    alt = c(1,-1,1)
    D1 = 0;
    if(K==2) D2 = 0 else D2 = NA
    for(i in 1:(n-2)){
      chck = (alt == sgn[i:(i+2)][order(rnk[i:(i+2)])])
      if(all(chck) | (all(!chck))){
        D1 = D1 + 1
        if(echo) print(c(i,i+1,i+2))
        if((K==2)&(i<=n-3)){
          chck = (alt == sgn[(i+1):(i+3)][order(rnk[(i+1):(i+3)])])
          if(all(chck) | (all(!chck))) D2 = D2 + 1
        }
      }
    }
    return(c(D1/(n-2),D2/(n-3)))
  }
  if(method=="Cpp"){
    return(kSD_Craw2(XR, sgn, K, echo)/c(n-2,n-3))
  }
}

#### fullLSimplex ####

#' Full-\code{1}-simplex and full-\code{2}-simplex depth for bivariate data
#'
#' C++ implementations of the algorithm for
#' the computation of the full-\code{2}-simplex depth for
#' two-dimensional data (residuals). The C++ implementation
#' is the fastest one that should be used in practice; the complexity
#' of both algorithms is however only O(n^4).
#'
#' @param X The bivariate dataset represented as a numerical matrix
#' with two columns and \code{n} rows. The rows of \code{X} are
#' typically the residuals.
#'
#' @param naive Logical indicator whether the naive implementation
#' in C++ should be used. By default set to \code{FALSE} if \code{n} is
#' greater than 125, otherwise \code{TRUE} is used. It was observed that
#' for \code{n} small, the naive implementation tends to be faster. As \code{n}
#' grows, the non-naive C++ is superior.
#'
#' @return A vector of length \code{2}, both elements are numeric values
#' between 0 and 1. First element is the full-\code{1}-simplex depth,
#' that is the same as the two-dimensional simplicial depth of the origin.
#' The second element it the full-\code{2}-simplex depth, that is the
#' relative number of non-consecutive pairs of triangles that share two
#' vertices and both contain the origin.
#'
#' @seealso \link{simplLSimplex} for an implementation of the simplified
#' version of the 1-simplex and 2-simplex depth.
#'
#' @examples
#' n = 100
#' X = matrix(rnorm(2*n),ncol=2)
#' fullLSimplex(X)
#' fullLSimplex(X, naive=TRUE)
#' fullLSimplex(X, naive=FALSE)
#'
#' n = 20
#' X = matrix(rnorm(2*n),ncol=2)
#' fullLSimplex(X)

fullLSimplex = function(X, naive = NULL){
  n = nrow(X)
  XR = X[,2] / X[,1]
  if(any(is.nan(XR))) warning("Some elements of X are numerically equal to 0,
  computation might be incorrect")
  rnk = rank(XR) # rank here
  sgn = sign(X[,1])
  if(is.null(naive)) naive = (n<150)
  #
  if(naive){
    # full-2-simplex
    res2 = kSD_Cfull2(sgn,rnk)
    # full-1-simplex
    S = sign(X[order(XR),1])
    res = sum(SignDepth::KSign(S, 3)[3,])
    res1 = res
  } else {
    # S<- sign(X[,2])
    # R<-rank(acos(sign(X[,2])*X[,1]/sqrt(X[,1]^2+X[,2]^2)))
    res = kSD_Cfull4(sgn,rnk)
    res1 = res[1]
    res2 = res[2]
  }
  res = c(res1/choose(n,3), res2/choose(n,4))
  names(res) = c("full1Simplex", "full2Simplex")
  return(res)
}

#### LSimplex ####

#' Simplified and full \code{L}-simplex depth for bivariate data
#'
#' Wrapper function for the fastest C++ implementations of the algorithms for
#' the computation of the simplified \code{L}-simplex depths, and the
#' full \code{2}-simplex depth for two-dimensional data (residuals).
#'
#' @param X The bivariate dataset represented as a numerical matrix
#' with two columns and \code{n} rows. The rows of \code{X} are
#' typically the residuals.
#'
#' @param L Positive integer 1 or 2, parameter of the depth. If
#' \code{L} is 1, only the (simplified or full) 1-simplex depth is computed.
#' If \code{L} is 2, both (simplified or full) 1-simplex depth and (simplified
#' or full) 2-simplex depth are computed.
#'
#' @param type Type of the depth to compute. Can take values \code{simplified}
#' or \code{full}.
#'
#' @return If \code{L=1}, a single numerical value with the (simplified or full)
#' \code{1}-simplex depth. If \code{L=2}, a vector of length \code{2}, both 
#' elements are numeric values between 0 and 1. The first element is the 
#' (simplified or full) \code{1}-simplex depth. The second element it the 
#' (simplified or full) \code{2}-simplex depth. 
#'
#' @seealso \link{simplLSimplex} and \link{fullLSimplex} for an implementation
#' of the simplified and full versions of the 1-simplex and 2-simplex depth,
#' respectively.
#'
#' @examples
#' n = 100
#' X = matrix(rnorm(2*n),ncol=2)
#' LSimplex(X, type="simpl")
#' LSimplex(X, type="full")

LSimplex = function(X, L=2, type=c("simplified", "full")){
  K = L
  type = match.arg(type)
  if(type=="simplified") return(simplLSimplex(X,K)[1:L])
  if(type=="full") if(K==2) return(fullLSimplex(X)[1:L]) else {
    x = X[order((X[,2]) / (X[,1])),1]
    res = sum(KSign(sign(x), 3)[3,])/choose(length(x),3)
    return(res[1:L])
    # return(ddalpha::depth.simplicial(c(0,0),data=X,exact=TRUE))
  }
}

#### LSimplex_MC ####

#' A Monte Carlo approximation of the null distribution of
#' \code{L}-simplex depths for bivariate data
#'
#' C++ implementation of a Monte Carlo approximation of the distribution
#' of the simplified \code{1}-simplex depth, simplified \code{2}-simplex depth,
#' and the full \code{1}-simplex depth and full \code{2}-simplex depth
#' for two-dimensional data (residuals).
#'
#' @param n Sample size, positive integer.
#'
#' @param B Number of independent Monte Carlo replicates to consider.
#' By default is set to 100000.
#'
#' @param full Logical variable that says whether also the
#' full \code{L}-simplex depths should be considered. The procedure is
#' much faster without the full \code{L}-simplex depths.
#'
#' @return A large matrix of dimensions \code{B}-times-\code{j}, where
#' \code{j} is 2 if only the simplified \code{L}-simplex depths are considered
#' (\code{full=FALSE}), 4 if also the full \code{L}-simplex depths are
#' considered (\code{full=TRUE}): first column for the simplified
#' \code{1}-simplex depth,
#' second column for the simplified \code{2}-simplex depth, third column to the
#' full \code{1}-simplex depth (that is, the usual 2-dimensional simplicial
#' depth of the origin), and fourth column for the full \code{2}-simplex depth.
#'
#' @seealso \link{LSimplex}, \link{simplLSimplex}, and \link{fullLSimplex} for
#' the implementations of the functions used to generate the Monte Carlo
#' approximation.
#'
#' @examples
#' n = 100
#' B = 50
#' (LMC = LSimplex_MC(n,B))
#'
#' # standardizing the columns to (approximately) pivotal quantities
#'
#' mns = c(1/4, 1/12, 1/4, 1/12)
#' sds = c(1/4*sqrt(11/3)/sqrt(n-2),
#'         1/12*sqrt(169/10)/sqrt(n-3),
#'         1/n,
#'         1/n)
#'
#' (LMC_stand = t((t(LMC)-mns)/sds)) # standardized kMC

LSimplex_MC = function(n,B=1e5,full=TRUE,naive=NULL){
  rnks = replicate(B,sample(n,n))
  rnks = lapply(seq_len(ncol(rnks)), function(i) rnks[,i])
  sgns = matrix((rbinom(n*B,1,prob=1/2)-1/2)*2,nrow=n)
  sgns = lapply(seq_len(ncol(sgns)), function(i) sgns[,i])
  if(is.null(naive)) naive = (n<125)
  if(full){
    if(naive){
      res = matrix(
        unlist(
          purrr::map2(rnks,sgns,function(x,y)
            c(kSD_Craw2(x, y, 2, FALSE)/c(n-2,n-3),
              sum(KSign(y, 3)[3,])/choose(n,3),
              kSD_Cfull2(y,x)/choose(n,4)))),
        ncol=4, byrow=TRUE)
    } else {
      res = matrix(
        unlist(
          purrr::map2(rnks,sgns,function(x,y)
            c(kSD_Craw2(x, y, 2, FALSE)/c(n-2,n-3),
              kSD_Cfull4(y,x)/c(choose(n,3),choose(n,4))))),
        ncol=4, byrow=TRUE)
    }
    colnames(res) = c("1Simplex", "2Simplex",
                      "full1Simplex", "full2Simplex")
    return(res)
  } else {
    res = matrix(
      unlist(
        purrr::map2(rnks,sgns,function(x,y)
          c(kSD_Craw2(x, y, 2, FALSE)/c(n-2,n-3)))),
      ncol=2, byrow=TRUE)
    colnames(res) = c("1Simplex", "2Simplex")
    return(res)
  }
}
