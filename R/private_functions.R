sigdig <- function(x,digits=3) {
  m <- digits - ceiling(log10(max(abs(x),na.rm=T)))
  if (m > 0) {
    round(x,m)
  } else {
    signif(round(x),digits)
  }
}

pvar <- function(mu=0,V=NULL,weights=NULL) {
  if (!is.null(V)) {
    n <- nrow(V)
  } else {
    n <- length(mu)
  }
  if (is.null(weights)) {
    weights <- rep(1,n)
  }
  weights <- weights/sum(weights)
  if (!is.null(V)) {
    x <- sum(diag(V)*weights) - matrix(weights,nrow=1) %*% V %*% matrix(weights,ncol=1) + sum(mu^2*weights) - sum(mu*weights)^2
  } else {
    x <- sum(mu^2*weights) - sum(mu*weights)^2
  }
  as.numeric(x)
}

kron <- function(eigen.A, B) {
  #returns Q such that QtQ is solution
  #eigen.B <- eigen(B)
  tmp <- svd(B)
  eigen.B <- list(values=tmp$d,vectors=tmp$u)
  V1 <- kronecker(Diagonal(x=sqrt(eigen.A$values)),Diagonal(x=sqrt(eigen.B$values)))
  V1.inv <- kronecker(Diagonal(x=1/sqrt(eigen.A$values)),Diagonal(x=1/sqrt(eigen.B$values)))
  V2 <- as(as(Matrix(eigen.B$vectors,dimnames=list(rownames(B),rownames(B))),"generalMatrix"),"unpackedMatrix")
  V3 <- kronecker(eigen.A$vectors,V2,make.dimnames=T)
  return(list(mat=tcrossprod(V1,V3), inv=tcrossprod(V1.inv,V3)))
}

Keff <- function(r2,alpha) {
  m <- nrow(r2)
  if (m > 1) {
    Q <- sqrt(r2)
    Q[upper.tri(Q,diag=T)] <- NA
    rmax <- apply(Q[-1,],1,max,na.rm=T)
    kappa <- sqrt(1-rmax^(-1.31*log10(alpha)))
    return(1+sum(kappa))
  } else {
    return(1)
  }
}

get_x <- function(map) {
  #takes a map with chrom and position and returns x axis values for plotting multiple chromosomes
  a <- tapply(map[,2],map[,1],max)
  n <- length(a)
  m <- tapply(map[,2],map[,1],length)
  b <- c(0,apply(array(1:(n-1)),1,function(k){sum(a[1:k])}))
  x <- map[,2] + rep(b,times=m)
  return(x)
}

coerce_dpo <- function(x) {
  x2 <- try(as(x,"dpoMatrix"),silent=TRUE)
  if (class(x2)=="dpoMatrix") {
    tmp <- try(Cholesky(x2), silent = TRUE)
    if (class(tmp)=="Cholesky") {
      return(x2)
    } else {
      tmp <- try(chol(x2),silent=TRUE)
      if (class(tmp)=="Cholesky") {
        return(x2)
      }
    }
  }

  d <- Diagonal(x=1/sqrt(diag(x)))
  x2 <- crossprod(d,x%*%d)
  eg <- eigen(x2,symmetric=TRUE)
  thresh <- .Machine$double.eps*10
  repeat {
    lambda <- ifelse(eg$values < thresh,thresh,eg$values)
    K <- Matrix(Diagonal(x=sqrt(diag(x))) %*% eg$vectors %*% Diagonal(x=sqrt(lambda)))
    x3 <- tcrossprod(K)
    dimnames(x3) <- dimnames(x)
    tmp <- chol(x3)
    if (class(tmp)=="Cholesky") {
      return(x3)
    }
    thresh <- thresh*10
  }
}

f.cov.trait <- function(vc,traits,us) {
  n.trait <- length(traits)
  if (us)  {
    cov.mat <- matrix(0,nrow=n.trait,n.trait)
    dimnames(cov.mat) <- list(traits,traits)
    tmp <- expand.grid(traits,traits)
    tmp <- apply(tmp,2,as.character)[,c(2,1)]
    tmp2 <- apply(tmp,1,paste,collapse=":")
    ix <- sapply(as.list(tmp2),grep,x=rownames(vc),fixed=T)
    iv <- which(sapply(ix,length)>0)
    ix <- unlist(ix[iv])
    tmp <- tmp[iv,]
    cov.mat[cbind(tmp[,1],tmp[,2])] <- vc[ix,1]
    cov.mat[upper.tri(cov.mat,diag=F)] <- cov.mat[lower.tri(cov.mat,diag=F)]
  } else {
    iu <- apply(array(traits),1,grep,x=rownames(vc),fixed=T)
    cov.mat <- diag(vc[iu,1])
    dimnames(cov.mat) <- list(traits,traits)
    if (vc[1,4]!="B") {
      cov.mat[1,2] <- cov.mat[2,1] <- vc[1,1]*sqrt(vc[iu[1],1]*vc[iu[2],1])
    }
  }
  return(coerce_dpo(cov.mat))
}

f.cov.loc <- function(vc,locs) {
  n.loc <- length(locs)
  vcnames <- rownames(vc)
  if (n.loc==2) {
    iu <- apply(array(locs),1,grep,x=vcnames,fixed=T)
    cov.mat <- diag(vc[iu,1])
    dimnames(cov.mat) <- list(locs,locs)
    cov.mat[1,2] <- cov.mat[2,1] <- vc[1,1]*sqrt(vc[iu[1],1]*vc[iu[2],1])
    fa.mat <- t(chol(cov.mat))
  } else {
    psi <- diag(vc[apply(array(paste0(locs,"!var")),1,grep,x=vcnames,fixed=T),1])
    fa.mat <- matrix(0,nrow=n.loc,ncol=2)
    rownames(fa.mat) <- locs
    ix <- apply(array(paste0(locs,"!fa1")),1,grep,x=vcnames,fixed=T)
    fa.mat[,1] <- vc[ix,1]
    ix <- apply(array(paste0(locs,"!fa2")),1,grep,x=vcnames,fixed=T)
    fa.mat[,2] <- vc[ix,1]
    cov.mat <- tcrossprod(fa.mat) + psi
  }

  #rotate
  tmp <- svd(fa.mat)
  fa.mat <- tmp$u %*% diag(tmp$d)
  #scale
  D <- diag(1/sqrt(diag(cov.mat)))
  dimnames(D) <- list(locs,locs)

  return(list(cov.mat=coerce_dpo(cov.mat),
              loadings=D%*%fa.mat))
}

#
# f.cov.har <- function(vc, hars) {
#   n.loc <- length(hars)
#   vcnames <- rownames(vc)
#
#
#
#
#   if (n.loc==2) {
#     iu <- apply(array(hars),1,grep,x=vcnames,fixed=T)
#     cov.mat <- diag(vc[iu,1])
#     dimnames(cov.mat) <- list(hars,hars)
#     cov.mat[1,2] <- cov.mat[2,1] <- vc[1,1]*sqrt(vc[iu[1],1]*vc[iu[2],1])
#     fa.mat <- t(chol(cov.mat))
#   } else {
#     psi <- diag(vc[apply(array(paste0(hars,"!var")),1,grep,x=vcnames,fixed=T),1])
#     fa.mat <- matrix(0,nrow=n.loc,ncol=2)
#     rownames(fa.mat) <- hars
#     ix <- apply(array(paste0(hars,"!fa1")),1,grep,x=vcnames,fixed=T)
#     fa.mat[,1] <- vc[ix,1]
#     ix <- apply(array(paste0(hars,"!fa2")),1,grep,x=vcnames,fixed=T)
#     fa.mat[,2] <- vc[ix,1]
#     cov.mat <- tcrossprod(fa.mat) + psi
#   }
#
#   #rotate
#   tmp <- svd(fa.mat)
#   fa.mat <- tmp$u %*% diag(tmp$d)
#   #scale
#   D <- diag(1/sqrt(diag(cov.mat)))
#   dimnames(D) <- list(locs,locs)
#
#   return(list(cov.mat=coerce_dpo(cov.mat),
#               loadings=D%*%fa.mat))
# }
#



cov_to_cor <- function(x) {
  d <- Diagonal(x=1/sqrt(diag(x)))
  x2 <- crossprod(d,x%*%d)
  dimnames(x2) <- dimnames(x)
  return(as.matrix(x2))
}

direct_sum <- function(x) {
  n <- length(x)
  m <- sapply(x,nrow)
  m.cumulative <- apply(array(1:n),1,function(k){sum(m[1:k])})

  z <- expand.grid(col=1:m[1],row=1:m[1])
  out <- data.frame(row=z$row,col=z$col,value=as.vector(x[[1]]))
  if (n > 1) {
    for (i in 2:n) {
      z <- expand.grid(col=(1:m[i])+m.cumulative[i-1],row=(1:m[i])+m.cumulative[i-1])
      tmp <- data.frame(row=z$row,col=z$col,value=as.vector(x[[i]]))
      out <- rbind(out,tmp)
    }
  }
  out <- out[out$col <= out$row,]
  dimO <- max(m.cumulative)
  dname <- as.character(1:dimO)
  return(new("dsTMatrix",uplo="L",i=out[,1]-1L,j=out[,2]-1L,x=out[,3],
             Dim=c(dimO,dimO),Dimnames=list(dname,dname),factors=list()))
}




G_fa <- function (var, term, G.param, cond.fac = "") {
  est <- G.param[[term]][[var]]$initial
  k <- strsplit(G.param[[term]][[var]]$facnam, split = "k = ",
                fixed = TRUE)[[1]][2]
  k <- as.numeric(substr(k, start = 1, stop = nchar(k) - 1))
  nlevs <- length(G.param[[term]][[var]]$levels) - k
  if (nlevs * (k + 1) != length(est)) {
    stop(paste("Number of levels of ", var, " and the number of variance parameters do not agree",
               sep = ""))
  }
  specvar <- diag(est[grepl("!var", names(est))])
  loadings <- matrix(est[grepl("!fa", names(est))], ncol = k)
  G <- tcrossprod(loadings) + specvar
  return(G)
}


mat_Gvar <- function (var, term, G.param, kspecial, cond.fac = "",  ...) {
  if (cond.fac != "") cond.fac <- paste(cond.fac, "_", sep = "")
  G <- switch(kspecial$cortype,
              id = asremlPlus:::G.id(var = var, term = term, G.param = G.param, cond.fac = cond.fac),
              ar1 = asremlPlus:::G.ar1(var = var, term = term, G.param = G.param, cond.fac = cond.fac),
              ar2 = asremlPlus:::G.ar2(var = var, term = term, G.param = G.param,
                                       cond.fac = cond.fac),
              ar3 = asremlPlus:::G.ar3(var = var, term = term, G.param = G.param, cond.fac = cond.fac),
              sar = asremlPlus:::G.sar(var = var, term = term, G.param = G.param,
                                       cond.fac = cond.fac),
              sar2 = asremlPlus:::G.sar2(var = var, term = term, G.param = G.param, cond.fac = cond.fac),
              ma1 = asremlPlus:::G.ma1(var = var, term = term, G.param = G.param,
                                       cond.fac = cond.fac),
              ma2 = asremlPlus:::G.ma2(var = var,term = term, G.param = G.param, cond.fac = cond.fac),
              arma = asremlPlus:::G.arma(var = var, term = term, G.param = G.param,
                                         cond.fac = cond.fac),
              exp = asremlPlus:::G.exp(var = var,term = term, G.param = G.param, cond.fac = cond.fac),
              gau = asremlPlus:::G.gau(var = var, term = term, G.param = G.param,
                                       cond.fac = cond.fac),
              cor = asremlPlus:::G.cor(var = var,term = term, G.param = G.param, cond.fac = cond.fac),
              corb = asremlPlus:::G.corb(var = var, term = term, G.param = G.param,
                                         cond.fac = cond.fac),
              corg = asremlPlus:::G.corg(var = var,term = term, G.param = G.param, cond.fac = cond.fac),
              us = asremlPlus:::G.us(var = var, term = term, G.param = G.param,
                                     cond.fac = cond.fac), spl = G.spl(var = var,
                                                                       term = term, G.param = G.param, cond.fac = cond.fac),
              fa = G_fa(var = var, term = term, G.param = G.param, cond.fac = cond.fac),
              rr = asremlPlus:::G.rr(var = var, term = term,  G.param = G.param, cond.fac = cond.fac),
              ante = asremlPlus:::G.us(var = var, term = term, G.param = G.param, cond.fac = cond.fac),
              G.unknown(kspecial))
  if (kspecial$final == "v") {
    if (kspecial$cortype == "id")
      vpname <- paste(cond.fac, term, "!", var,
                      sep = "")
    else vpname <- paste(cond.fac, term, "!", var,
                         "!var", sep = "")
    if (!is.na(G.param[[term]][[var]]$initial[vpname]))
      G <- G.param[[term]][[var]]$initial[vpname] *
        G
  }
  if (kspecial$final == "h") {
    if (grepl("+", term, fixed = TRUE))
      vpnames <- paste(cond.fac, term, "!", G.param[[term]][[var]]$model,
                       "(", G.param[[term]][[var]]$facnam, ")_",
                       G.param[[term]][[var]]$levels, sep = "")
    else vpnames <- paste(cond.fac, term, "!",
                          var, "_", G.param[[term]][[var]]$levels,
                          sep = "")
    D <- G.param[[term]][[var]]$initial[vpnames]
    D <- sqrt(diag(D, nrow = length(G.param[[term]][[var]]$levels)))
    G <- D %*% G %*% D
  }
  return(G)
}




# A function to extract the residual variance covariance matrix at a grouping level
getVarCov.asreml <- function(asreml.obj, which.matrix = "V", extra.matrix = NULL,  ignore.terms = NULL, ignore.vars = NULL,
                             fixed.spline.terms = NULL, bound.exclusions = c("F","B", "S", "C")){
  require(asremlPlus)
  asr4 <- asremlPlus:::isASRemlVersionLoaded(4, notloaded.fault = TRUE)
  if (!asr4) stop("This function requires asreml4 or later.")
  validasr <- asremlPlus:::validAsreml(asreml.obj)
  if (is.character(validasr)) stop(validasr)
  mat.options <- c("V", "G", "R")
  mat.opt <- mat.options[asremlPlus:::check.arg.values(which.matrix, mat.options)]
  ran.specials <- c("at", "spl", "dev", "grp",
                    "fa", "rr", "ante")
  res.specials <- c("dsum")
  common.specials <- c("id", "diag", "us",
                       "ar1", "ar2", "ar3", "sar", "sar2",
                       "ma1", "ma2", "arma", "exp",
                       "gau", "cor", "corb", "corg")
  other.specials <- c("sph", "chol", "ante",
                      "sfa", "facv")
  design.specials <- c("spl", "dev", "grp",
                       "fa", "rr")
  if (!all(bound.exclusions %in% c("F", "B", "S", "C")))  stop("At least one bound.type is not one of those allowed with ASReml-R version 3")
  call <- asreml.obj$call
  if (!("data" %in% names(call))) stop("estimateV.asreml assumes that data has been set in call to asreml")
  dat <- eval(call$data)
  dat <- asreml.obj$mf
  n <- nrow(dat)
  V <- matrix(0, nrow = n, ncol = n)
  G.param <- asreml.obj$G.param
  ranterms <- names(G.param)
  R.param <- asreml.obj$R.param
  resterms <- names(R.param)
  if (!(is.null(ignore.terms)))
    if (any(is.na(match(ignore.terms, c(ranterms, resterms)))))
      stop(paste("The following terms are not amongst the variance parameters: ",
                 paste(ignore.terms[is.na(match(ignore.terms,
                                                c(ranterms, resterms)))], collapse = ", "),
                 sep = ""))
  if (!is.null(extra.matrix) & (nrow(V) != n | ncol(V) != n))
    stop("V must be square and of order equal to the number of rows in data")

  if (mat.opt %in% c("V", "G")) {
    if (!is.null(fixed.spline.terms)) {
      if (any(is.na(match(fixed.spline.terms, ranterms))))
        stop(paste("The following spline terms are not amongst the random variance parameters: ",
                   paste(fixed.spline.terms[is.na(match(fixed.spline.terms,
                                                        ranterms))], collapse = ", "), sep = ""))
      ranterms <- ranterms[-match(fixed.spline.terms, ranterms)]
      G.param <- G.param[ranterms]
    }
    if (!is.null(ignore.terms)) {
      kignore <- na.omit(match(ignore.terms, ranterms))
      if (length(kignore) > 0) {
        ranterms <- ranterms[-na.omit(match(ignore.terms,
                                            ranterms))]
        G.param <- G.param[ranterms]
      }
    }
    for (term in ranterms) {
      bound <- FALSE
      if (!is.null(bound.exclusions) & term %in% names(asreml.obj$vparameters.con)) {
        bound <- (asreml.obj$vparameters.con[[term]] %in% bound.exclusions)
        if (bound)
          warning(paste(term, "not included in V because its bound is",
                        asreml.obj$vparameters.con[[term]], sep = " "))
      }
      if (!bound) {
        if (term == "units" | (G.param[[term]]$variance$model ==
                               "idv" & (all(unlist(lapply(G.param[[term]][2:length(G.param[[term]])],
                                                          function(x) {
                                                            is.id <- (x$model == "id")
                                                            return(is.id)
                                                          })))))) {
          if (term == "units") {
            V <- V + diag(G.param$units$variance$initial,
                          nrow = n)
          } else {
            Z <- model.matrix(as.formula(paste("~ - 1 + ", term)), data = dat)
            V <- V + G.param[[term]]$variance$initial * Z %*% t(Z)
          }
        } else {
          termvars <- names(G.param[[term]])[-1]
          ignore_vars_term <- ignore.vars[[term]]
          termvars <- setdiff(termvars, ignore_vars_term)
          cond.fac <- ""
          if (G.param[[term]]$variance$model == "idv") {
            G <- G.param[[term]]$variance$initial
          } else {
            G <- 1
          }
          has.fa <- ((any(unlist(lapply(G.param[[term]][2:length(G.param[[term]])],
                                        function(x) {
                                          is.fa <- (x$model %in% c("fa", "rr"))
                                          return(is.fa)
                                        })))))
          if (has.fa) full.levs <- req.levs <- term
          for (var in termvars) {
            kspecial <- asremlPlus:::checkSpecial(var = var, term = term,
                                                  G.param = G.param, specials = c(ran.specials, common.specials), residual = FALSE)
            G <- kronecker(G, mat_Gvar(var = var, term = term,
                                       G.param = G.param, kspecial = kspecial,
                                       cond.fac = cond.fac, residual = FALSE))
            if (has.fa) {
              var.levs <- G.param[[term]][[var]]$levels
              if (any(!is.na(G.param[[term]][[var]]$initial))) {
                vp.levs <- as.data.frame(do.call(rbind,
                                                 strsplit(names(G.param[[term]][[var]]$initial),
                                                          split = "!", fixed = TRUE)),
                                         stringsAsFactors = FALSE)
                vp.levs <- vp.levs[vp.levs$V3 == "var",
                                   2]
              } else {
                vp.levs <- var.levs
              }
              full.levs <- as.vector(outer(var.levs,
                                           full.levs, paste, sep = "!"))
              req.levs <- as.vector(outer(vp.levs, req.levs,
                                          paste, sep = "!"))
            }
          }
          if ((any(unlist(lapply(G.param[[term]][2:length(G.param[[term]])],
                                 function(x) {
                                   need.design <- (x$model %in% design.specials)
                                   return(need.design)
                                 }))))) {
            if (is.null(asreml.obj$design)) {
              asreml::asreml.options(design = TRUE)
              asreml.obj <- eval(call)
            }
            if (has.fa) {
              cols <- c(1:length(full.levs))[match(req.levs, full.levs)]
              cols <- paste(var, cols, sep = "")
              Z <- asreml.obj$design[, match(cols, colnames(asreml.obj$design))]
            }
            else {
              Z <- as.matrix(getTermDesignMatrix(term,
                                                 asreml.obj))
            }
            Z <- as.matrix(Z)
            # V <- V + Z %*% G %*% t(as.matrix(Z))
          } else {
            if (grepl("+", term, fixed = TRUE)) {
              str.terms <- strsplit(term, split = "+",
                                    fixed = TRUE)[[1]]
              cols <- NULL
              for (str.term in str.terms) {
                cols <- c(cols, paste(str.term, 1:asreml.obj$noeff[str.term],
                                      sep = ""))
              }
              Z <- asreml.obj$design[, match(cols, colnames(asreml.obj$design))]
              # V <- V + Z %*% G %*% t(as.matrix(Z))
            } else {
              Z <- model.matrix(as.formula(paste("~ -1 +", var)), data = droplevels(dat))
              # V <- V + Z %*% G %*% t(Z)
            }
          }
        }
      }
    }
  }
  for (term in resterms) {
    if (R.param[[term]]$variance$model == "idv")
      foundvar <- TRUE
    else foundvar <- FALSE
  }
  if (mat.opt %in% c("V", "R")) {
    if (!is.null(ignore.terms)) {
      kignore <- na.omit(match(ignore.terms, resterms))
      if (length(kignore) > 0) {
        resterms <- resterms[-na.omit(match(ignore.terms,
                                            resterms))]
        R.param <- R.param[resterms]
      }
    }
    nosections <- length(R.param)
    Rlist <- vector(mode = "list", length = nosections)
    names(Rlist) <- resterms
    if (nosections > 1) {
      cond.fac <- strsplit(labels(asreml.obj$formulae$residual)[1],
                           split = "|", fixed = TRUE)[[1]][2]
      if (grepl(",", cond.fac, fixed = TRUE))
        cond.fac <- strsplit(cond.fac, ",", fixed = TRUE)[[1]][1]
      if (grepl(")", cond.fac, fixed = TRUE))
        cond.fac <- strsplit(cond.fac, ")", fixed = TRUE)[[1]][1]
      if (grepl("~", cond.fac, fixed = TRUE))
        cond.fac <- strsplit(cond.fac, "~", fixed = TRUE)[[1]][2]
      cond.fac <- trimws(cond.fac)
    } else {
      cond.fac <- ""
    }
    for (term in resterms) {
      if ((R.param[[term]]$variance$model == "idv" |
           R.param[[term]]$variance$model == "id") &
          (all(unlist(lapply(R.param[[term]][2:length(R.param[[term]])],
                             function(x) {
                               is.id <- (x$model == "id")
                               return(is.id)
                             }))))) {
        if (R.param[[term]]$variance$model == "idv")
          Rlist[[term]] <- R.param[[term]]$variance$initial *
            diag(1, nrow = R.param[[term]]$variance$size)
        else Rlist[[term]] <- diag(1, nrow = R.param[[term]]$variance$size)
      } else {
        termvars <- names(R.param[[term]])[-1]
        ignore_vars_term <- ignore.vars[[term]]
        termvars <- setdiff(termvars, ignore_vars_term)
        if (R.param[[term]]$variance$model == "idv")
          R <- R.param[[term]]$variance$initial
        else R <- 1
        for (var in termvars) {
          kspecial <- asremlPlus:::checkSpecial(var = var, term = term,
                                                G.param = R.param, specials = c(res.specials, common.specials, other.specials), residual = TRUE)
          if (kspecial$final == "v" | kspecial$final == "h") foundvar <- TRUE
          R <- kronecker(R, mat_Gvar(var = var, term = term,
                                     G.param = R.param, kspecial = kspecial, cond.fac = cond.fac))
        }
        Rlist[[term]] <- as.matrix(R)
      }
    }
    if (which.matrix == "V") {
      if (length(Rlist) == 1) {
        V <- V + Rlist[[1]]
      } else {
        V <- V + dae:::mat.dirsum(Rlist)
      }
    } else {
      R <- dae::mat.dirsum(Rlist)
    }
  }
  if (!foundvar)
    V <- asreml.obj$sigma2 * V
  if (!is.null(extra.matrix))
    V <- V + extra.matrix

  if (which.matrix == "G") {
    return(G)
  } else if (which.matrix == "R") {
    return(R)
  } else {
    return(V)
  }
}


