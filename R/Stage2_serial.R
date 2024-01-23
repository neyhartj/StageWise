#' Stage 2 analysis of multi-environment trials
#'
#' Stage 2 analysis of multi-environment trials
#'
#' Stage 2 of the two-stage approach described by Damesa et al. 2017, using ASReml-R for variance component estimation. The variable \code{data} has three mandatory column: id, env, BLUE. Optionally, \code{data} can have a column labeled "loc", which changes the main effect for genotype into a separable genotype-within-location effect, using a FA2 covariance model for the locations. Optionally, \code{data} can have a column labeled "trait", which uses an unstructured covariance model. The multi-location and multi-trait analyses cannot be combined. Missing data are allowed in the multi-trait but not the single-trait analysis. The argument \code{geno} is used to partition genetic values into additive and non-additive components. Any individuals in \code{data} that are not present in \code{geno} are discarded.
#'
#' The argument \code{vcov} is used to partition the macro- and micro-environmental variation, which are called GxE and residual in the output. \code{vcov} is a named list of variance-covariance matrices for the BLUEs within each environment, with id for rownames (single trait) or id:trait. The order in \code{vcov} and \code{data} should match. Both \code{data} and \code{vcov} can be created using the function \code{\link{Stage1}}.
#'
#' Because ASReml-R can only use relationship matrices defined in the global environment, this function creates and then removes global variables when either \code{vcov} or \code{geno} is used. By default, the workspace memory for ASReml-R is set at 500mb. If you get an error about insufficient memory, try increasing it. ASReml-R version 4.1.0.148 or later is required.
#'
#' @references Damesa et al. 2017. Agronomy Journal 109: 845-857. doi:10.2134/agronj2016.07.0395
#'
#' @param data data frame of BLUEs from Stage 1 (see Details)
#' @param vcov named list of variance-covariance matrices for the BLUEs
#' @param har.cor.str the covariance structure for the random genotype x harvest year effect
#' @param geno output from \code{\link{read_geno}}
#' @param fix.eff.marker markers in \code{geno} to include as additive fixed effect covariates
#' @param silent TRUE/FALSE, whether to suppress ASReml-R output
#' @param workspace Memory limit for ASRreml-R variance estimation
#' @param non.add one of the following: "none","g.resid","dom"
#' @param max.iter maximum number of iterations for asreml
#'
#' @return List containing
#' \describe{
#' \item{aic}{AIC}
#' \item{vars}{variance components for \code{\link{blup_prep}}, as variable of class \code{\link{class_var}}}
#' \item{params}{Estimates and SE for fixed effects and variance components}
#' \item{random}{Random effect predictions}
#' \item{loadings}{scaled loadings for the FA2 multi-loc model}
#' }
#'
#' @importFrom stats model.matrix var
#' @importFrom methods new
#' @importFrom rlang .data
#' @import Matrix
#' @import ggplot2
#' @import ggrepel
#'
#' @export

SerialStage2 <- function(data, vcov = NULL, har.cor.str = c("idt", "corv", "corh", "ar1", "ar1h", "us"),
                         geno = NULL, fix.eff.marker = NULL, silent = TRUE, rescale = FALSE,
                         workspace="500mb", non.add = c("g.resid", "none", "dom"), max.iter=20) {

  stopifnot(inherits(data,"data.frame"))
  stopifnot(requireNamespace("asreml"))
  non.add <- match.arg(non.add)
  har.cor.str <- match.arg(har.cor.str)
  if (non.add=="dom")
    stopifnot(class(geno)=="class_genoD")
  library(asreml)

  # Error handling
  required_cols <- c("id", "har", "trl", "subtrl", "BLUE")
  stopifnot(required_cols %in% colnames(data))
  data$id <- as.character(data$id)
  data$har <- as.character(data$har)
  data$trl <- as.character(data$trl)
  data$id.trl <- apply(data[,c("id","trl")],1,paste,collapse=":")
  data$id.subtrl <- apply(data[,c("id","subtrl")],1,paste,collapse=":")
  data$id.har <- apply(data[,c("id","har")],1,paste,collapse=":")

  missing <- which(is.na(data$BLUE))
  if (length(missing) > 0) {
    data <- data[!(data$id.har %in% data$id.har[missing]),]
  }

  diagG <- diagD <- numeric(0)
  dom <- NULL
  if (!is.null(geno)) {
    stopifnot(inherits(geno,"class_geno"))
    id <- sort(intersect(data$id,rownames(geno@G)))
    n <- length(id)
    data <- data[data$id %in% id,]
    id.weights <- table(factor(data$id,levels=id))
    .GlobalEnv$asremlG <- geno@G[id,id]
    meanG <- as.numeric(pvar(V=.GlobalEnv$asremlG,weights=id.weights))
    diagG <- mean(diag(.GlobalEnv$asremlG))

    if (non.add=="dom") {
      .GlobalEnv$asremlD <- geno@D[id,id]
      meanD <- as.numeric(pvar(V=.GlobalEnv$asremlD,weights=id.weights))

      dom.covariate <- as.numeric(geno@coeff.D[id,] %*% matrix(1,nrow=ncol(geno@coeff.D),ncol=1))
      dom.covariate <- dom.covariate/(geno@scale*(geno@ploidy-1))
      names(dom.covariate) <- id
      data$dom <- dom.covariate[data$id]

      diagD <- mean(diag(.GlobalEnv$asremlD))
      dom <- "dom"
    }
  } else {
    id <- sort(unique(data$id))
    n <- length(id)
  }

  if (!is.null(fix.eff.marker)) {
    stopifnot(!is.null(geno))
    n.mark <- length(fix.eff.marker)
    stopifnot(fix.eff.marker %in% colnames(geno@coeff))
    dat2 <- data.frame(id=rownames(geno@coeff),as.matrix(geno@coeff[,fix.eff.marker]))
    colnames(dat2) <- c("id",fix.eff.marker)
    data <- merge(data,dat2,by="id")
  } else {
    n.mark <- 0
  }


  # harvest years
  hars <- sort(unique(data$har))
  n.har <- length(hars)
  if (n.har == 1 & (non.add == "g.resid")) {
    stop("Need more than one harvest year for g.resid model")
  }
  # Trials
  trls <- unique(data$trl)
  n.trl <- length(trls)
  subtrls <- unique(data$subtrl)
  n.subtrl <- length(subtrls)

  # Make sure to fill in the sequence of harvest years as factors
  if (grepl("ar1", har.cor.str)) {
    data$har <- fct_expand_seq(data$har)
  } else {
    data$har <- factor(data$har, levels = hars)
  }


  data$id <- factor(data$id,levels=id)
  # Observed id x har combinations
  id.har <- sort(unique(data$id.har))
  tmp <- expand.grid(har = levels(data$har), id = id)
  id.har.all <- apply(tmp[,c("id","har")],1,paste,collapse=":")
  # Set levels of id.har to all factorial levels
  data$id.har <- factor(data$id.har, levels = id.har.all)
  id.trl <- unique(data$id.trl)
  data$id.trl <- factor(data$id.trl, levels = id.trl)
  id.subtrl <- data$id.subtrl
  data$id.subtrl <- factor(data$id.subtrl, levels = id.subtrl)
  data$trl <- factor(data$trl, levels = trls)


  # Sort the data by trial then harvest year;
  # expand the data if the cor.str is ar1 or ar1h
  if ("loc" %in% names(data)) {
    data <- data[order(data$trl, data$loc, data$har),]
    data_tmp <- cbind(tmp, id.har = id.har.all)
  } else {
    if (grepl("ar1", har.cor.str)) {
      data_tmp <- cbind(tmp, id.har = id.har.all)
      data_tmp$id.har <- factor(data_tmp$id.har, levels = id.har.all)
      data <- merge(x = data_tmp, y = data, all.x = TRUE)
      hars <- as.character(sort(unique(data$har)))
    }
    data <- data[order(data$trl, data$id, data$har),]

  }

  # non-missing data per harvest year
  har.miss <- tapply(X = data$BLUE, INDEX = data$har, FUN = function(x) sum(!is.na(x)))
  har.miss[is.na(har.miss)] <- 0

  # Calculate g x har weights
  gH.weights <- table(na.omit(data)$id.har)

  out <- vector("list", 4)
  names(out) <- c("aic","vars","fixed","random")
  n.loc <- 1
  n.trait <- 1

  if ("trait" %in% colnames(data)) {
    data$trait <- as.character(data$trait)
    traits <- sort(unique(data$trait))
    n.trait <- length(traits)
    stopifnot(n.trait > 1)

    data$env.id.trait <- apply(data[,c("env.id","trait")],1,paste,collapse=":")
    env.id.trait <- unique(data$env.id.trait)
    data$env.id.trait <- factor(data$env.id.trait,levels=env.id.trait)

    if (!is.null(vcov)) {
      dname <- strsplit(rownames(vcov[[1]])[1:n.trait],split=":",fixed=T)
      traits <- sapply(dname,"[[",2)
    }
    data$Trait <- factor(data$trait,levels=traits)
    data <- data[order(data$env.id,data$Trait),]

  } else {

    traits <- ""  #NULL
    if ("loc" %in% colnames(data)) {
      data$loc <- factor(as.character(data$loc))
      data <- data[order(data$loc),]
      locations <- levels(data$loc)
      n.loc <- length(locations)
      stopifnot(n.loc > 1)
      loc.weights <- table(data$loc)
      data$gL <- paste(as.character(data$id),as.character(data$loc),sep=":")
      tmp <- expand.grid(loc=locations,id=id)
      gL <- apply(tmp[,c("id","loc")],1,paste,collapse=":")
      data$gL <- factor(data$gL,levels=gL)
      gL.weights <- table(data$gL)
    }
  }

  if (is.null(geno)) {
    tmp <- c("har","fixed.marker","additive","add x loc","dominance","heterosis","genotype","g x loc","g x har", "g x har corr", "Stage1.error","residual")
    # tmp <- c("env","fixed.marker","additive","add x loc","dominance","heterosis","genotype","g x loc","g x env","Stage1.error","residual")
  } else {
    tmp <- c("har","fixed.marker","additive","add x loc", "add x har","dominance","heterosis","g.resid","g x loc","g x har", "g x har corr", "Stage1.error","residual")
    # tmp <- c("env","fixed.marker","additive","add x loc","dominance","heterosis","g.resid","g x loc","g x env","Stage1.error","residual")
  }
  vars <- array(data=numeric(0),dim=c(n.trait,n.trait,length(tmp)),dimnames=list(traits,traits,tmp))

  if (n.trait==1) {
    # Heterogeneous residual per location
    if (n.loc==1) {
      model <- "asreml(data=data,fixed=BLUE~FIX-1,random=~RANDOM,residual=~idv(units)"
    } else {
      model <- "asreml(data=data,fixed=BLUE~FIX-1,random=~RANDOM,residual=~dsum(~idv(units)|loc)"
    }
    if (n.mark > 0 | (non.add=="dom")) {
      if (n.loc==1) {
        model <- sub("FIX",paste(c("har",fix.eff.marker,dom),collapse="+"),model,fixed=T)
      } else {
        model <- sub("FIX",paste(c("har",paste0(c(fix.eff.marker,dom),":loc")),collapse="+"),model,fixed=T)
      }
    } else {
      model <- sub("FIX","har",model,fixed=T)
    }
  } else {
    model <- "asreml(data=data,fixed=BLUE~FIX-1,random=~RANDOM,residual=~id(env.id):us(Trait)"
    if (n.mark > 0 | (non.add=="dom")) {
      model <- sub("FIX",paste(paste0(c("har",fix.eff.marker,dom),":Trait"),collapse="+"),model,fixed=T)
    } else {
      model <- sub("FIX","har:Trait",model,fixed=T)
    }
  }

  # Possibilities for formulating the random effects formula:
  # 1. genotype data provided (GRM) or not
  # 2. harvest year covariance
  # 3. single location or multiple location
  # 4. single trial or multiple trials
  # 5. single trait or multiple traits
  #

  # 1. How to code genotypes
  if (is.null(geno)) {
    random.effects.id <- "id"
  } else {
    random.effects.id <- "vm(id,source=asremlG,singG='PSD')"
  }

  # 2. How to code harvest year covariance
  if (har.cor.str == "idt") {
    random.effects.har <- NULL
    random.effects.idhar <- random.effects.id
  } else {
    if (n.har == 2) {
      random.effects.har <- "corh(har)"
    } else {
      random.effects.har <- paste0(har.cor.str, "(har)")
    }
    # Combine id and har effects
    random.effects.idhar <- paste0(random.effects.id, ":", random.effects.har)
  }


  # 3. How to code single locations
  if (n.loc > 1) {
    random.effects.loc.har <- paste0("loc:", random.effects.har)
    if (n.loc == 2) {
      random.effects.loc.gen.har <- paste0(random.effects.id, ":corh(loc):", random.effects.har)
    } else {
      random.effects.loc.gen.har <- paste0(random.effects.id, ":fa(loc, 2):", random.effects.har)
    }
  } else {
    random.effects.loc.har <- random.effects.loc.gen.har <- NULL
  }

  # 4. How to code single or multiple trials
  if (n.trl > 1) {
    random.effects.loc.trl.har <- paste0("trl:", random.effects.har)
  } else {
    random.effects.loc.trl.har <- NULL
  }

  # 5. How to code single or multiple traits
  if (n.trait > 1) {
    stop ("Multiple traits are not yet supported")
  }

  # Combine random effects strings
  random.effects <- paste(c(random.effects.idhar, random.effects.loc.har, random.effects.loc.gen.har, random.effects.loc.trl.har),
                          collapse = " + ")


  n.gE <- sum(!is.na(data$BLUE))/n.trait

  if (!is.null(vcov)) {
    stopifnot(is.list(vcov))
    stopifnot(subtrls %in% names(vcov))
    vcov <- vcov[subtrls]
    vcov <- mapply(FUN=function(x, y){
      # Remove missing
      # xmiss <-
      tmp <- paste(rownames(x), y, sep=":")
      dimnames(x) <- list(tmp, tmp)
      return(x)},x=vcov,y=as.list(names(vcov)))

    if (n.trait==1) {
      ix <- lapply(vcov,function(y){which(rownames(y) %in% id.subtrl)})
      random.effects <- paste0(random.effects,"+vm(id.subtrl,source=asremlOmega)")
    } else {
      ix <- lapply(vcov,function(y){which(rownames(y) %in% id.subtrl.trait)})
      random.effects <- paste0(random.effects,"+vm(id.subtrl.trait,source=asremlOmega)")
    }

    omega.list <- vector("list", n.subtrl)
    for (j in 1:n.subtrl) {
      omega.list[[j]] <- as(vcov[[j]][ix[[j]],ix[[j]]],"dpoMatrix")
    }
    .GlobalEnv$asremlOmega <- direct_sum(lapply(omega.list,solve))
    dname <- lapply(omega.list,rownames)
    dimnames(.GlobalEnv$asremlOmega) <- list(unlist(dname),unlist(dname))
    attr(.GlobalEnv$asremlOmega,"INVERSE") <- TRUE

    Omega <- bdiag(omega.list)
    ix <- split(1:(n.gE * n.trait), rep(1:n.trait, times = n.gE))
    for (i in 1:n.trait)
      for (j in i:n.trait) {
        Omega_ij <- Omega[ix[[i]],ix[[j]]]
        vars[i,j,"Stage1.error"] <- mean(diag(Omega_ij)) - mean(Omega_ij)
      }
  }

  asreml::asreml.options(workspace=workspace,maxit=max.iter,trace=!silent)
  model <- sub(pattern = "RANDOM", replacement = random.effects, model, fixed = TRUE)

  model_initial <- eval(parse(text=paste0(model,",start.values = TRUE, na.action = na.method(x = 'ignore'))")))
  start.table <- model_initial$vparameters.table

  if (!is.null(vcov)) {
    k <- grep("Omega",start.table$Component,fixed=T)
    start.table$Value[k] <- 1
    start.table$Constraint[k] <- "F"

    # Hold residual variance to 1
    k <- grep("units",start.table$Component,fixed=T)
    start.table$Value[k] <- 1
    start.table$Constraint[k] <- "F"

    # model <- paste0(model, ", R.param = start.table")

  }

  if (!is.null(geno) & (n.loc > 1)) {
    k <- grep(":loc!loc!cor",start.table$Component,fixed=T)
    k2 <- grep("asremlG",start.table$Component,fixed=T)
    k <- setdiff(k,k2)
    start.table$Value[k] <- 0.9999
    start.table$Constraint[k] <- "F"
  }

  # # Expand grid
  # data1 <- expand.grid(trl = trls, har = levels(data$har), id = levels(data$id), stringsAsFactors = TRUE)
  # data1$id.har <- as.factor(id.har.all)
  # data1 <- merge(x = data1, y = data[,c("trl", "har", "id", "BLUE")], all.x = TRUE, all.y = FALSE)

  ans <- eval(parse(text = paste0(model,", G.param = start.table, na.action = na.method(x = 'ignore'))")))
  # ans <- eval(parse(text = paste0(model,", family = asr_gaussian(dispersion = 1.0))")))
  if (!ans$converge) {
    stop("ASReml-R did not converge. Try increasing max.iter")
  }

  ## Try adjusting the R matrix instead

#
#   model_initial <- asreml(data=data1, fixed=BLUE~har-1, random=~id:ar1(har),
#                           residual= ~ trl.id.har, start.values = TRUE, na.action = na.method(x = 'ignore'))
#   start.table <- model_initial$vparameters.table
#   R_param <- model_initial$R.param
#
#
#   if (!is.null(vcov)) {
#     k <- grep("Omega",start.table$Component,fixed=T)
#     start.table$Value[k] <- 1
#     start.table$Constraint[k] <- "F"
#
#     # Hold residual variance to 1
#     k <- grep("units",start.table$Component,fixed=T)
#     start.table$Value[k] <- 1
#     start.table$Constraint[k] <- "F"
#
#     model <- paste0(model, ", R.param = start.table")
#
#   }
#
#   if (!is.null(geno) & (n.loc > 1)) {
#     k <- grep(":loc!loc!cor",start.table$Component,fixed=T)
#     k2 <- grep("asremlG",start.table$Component,fixed=T)
#     k <- setdiff(k,k2)
#     start.table$Value[k] <- 0.9999
#     start.table$Constraint[k] <- "F"
#   }
#


  # while (!ans$converge) {
  #   cat("ASReml-R failed to converge. Do you wish to continue running? y/n \n")
  #   input <- readLines(n=1)
  #   if (input=="y") {
  #     ans <- asreml::update.asreml(ans)
  #   } else {
  #     return()
  #   }
  # }

  sans <- summary(ans,coef=TRUE)
  out$aic <- as.numeric(sans$aic)

  B <- matrix(0,nrow=n.trait,ncol=n.trait)
  dimnames(B) <- list(traits,traits)

  #fixed effects
  if (n.trait==1) {
    beta <- sans$coef.fixed
    beta.names <- rownames(beta)
    beta.names <- gsub("har_","",beta.names)
    ix <- match(levels(data$har),beta.names)
    out$params$har <- data.frame(har=levels(data$har),
                                 estimate=as.numeric(beta[ix,1]),
                                 SE=as.numeric(beta[ix,2]))
    har_weights <- tapply(X = data$BLUE, INDEX = data$har, FUN = function(x) sum(!is.na(x)))
    vars[1,1,"har"] <- pvar(mu=out$params$har$estimate,weights=as.numeric(har_weights))

    if (non.add=="dom") {
      ix <- grep("dom",beta.names,fixed=T)
      if (n.loc==1) {
        vars[1,1,"heterosis"] <- pvar(mu=dom.covariate*beta[ix,1],weights=id.weights)
        out$params$heterosis <- data.frame(estimate=as.numeric(beta[ix,1]),
                                           SE=as.numeric(beta[ix,2]))
      } else {
        beta.names[ix] <- gsub("loc_","",beta.names[ix])
        beta.names[ix] <- gsub("dom","",beta.names[ix])
        beta.names[ix] <- gsub(":","",beta.names[ix])
        ix <- ix[match(locations,beta.names[ix])]
        mu <- as.numeric(kronecker(matrix(dom.covariate,ncol=1),beta[ix,1]))
        vars[1,1,"heterosis"] <- pvar(mu=mu,weights=gL.weights)
        out$params$heterosis <- data.frame(loc=locations,
                                           estimate=as.numeric(beta[ix,1]),
                                           SE=as.numeric(beta[ix,2]))
      }
    }
    if (n.mark > 0) {
      if (n.loc==1) {
        ix <- match(fix.eff.marker,beta.names)
        out$params$marker <- data.frame(marker=fix.eff.marker,
                                        estimate=as.numeric(beta[ix,1]),
                                        SE=as.numeric(beta[ix,2]))
        mu <- as.numeric(as.matrix(geno@coeff[id,fix.eff.marker])%*%beta[ix,1])
        vars[1,1,"fixed.marker"] <- pvar(mu=mu,weights=id.weights)
      } else {
        ix <- grep("loc_",beta.names)
        tmp <- strsplit(beta.names[ix],split=":",fixed=T)
        beta.names[ix] <- sapply(tmp,function(z){
          tmp2 <- apply(array(z),1,grep,pattern="loc_")
          paste(z[order(sapply(tmp2,length))],collapse=":")
        })
        beta.names[ix] <- gsub("loc_","",beta.names[ix])
        loc.marker <- expand.grid(loc=locations,marker=fix.eff.marker,
                                  stringsAsFactors = FALSE)
        loc.marker <- loc.marker[,c(2,1)]
        ix <- match(apply(loc.marker,1,paste,collapse=":"),beta.names)
        x <- out$params$marker <- data.frame(marker=loc.marker$marker,
                                       loc=loc.marker$loc,
                                       estimate=as.numeric(beta[ix,1]),
                                       SE=as.numeric(beta[ix,2]))
        mu_gL <- as.numeric(kronecker(as.matrix(geno@coeff[id,fix.eff.marker]),diag(n.loc))%*%beta[ix,1])
        vars[1,1,"fixed.marker"] <- pvar(mu=mu_gL,weights=gL.weights)
      }
    }
  } else {
    #multi-trait

    beta <- sans$coef.fixed
    ix <- grep("env_",rownames(beta),fixed=T)
    rownames(beta) <- gsub("env_","",rownames(beta))
    rownames(beta) <- gsub("Trait_","",rownames(beta))
    tmp <- strsplit(rownames(beta)[ix],split=":",fixed=T)
    out$params$env <- data.frame(env=sapply(tmp,"[[",1),
                                 trait=sapply(tmp,"[[",2),
                                 estimate=as.numeric(beta[ix,1]),
                                 SE=as.numeric(beta[ix,2]))

    weights <- as.numeric(table(data$env[data$trait==traits[1]]))
    weights <- weights/sum(weights)
    for (i in 1:n.trait) {
       for (j in i:n.trait) {
        mu1 <- out$params$env$estimate[out$params$env$trait==traits[i]]
        mu2 <- out$params$env$estimate[out$params$env$trait==traits[j]]
        vars[i,j,"env"] <- sum(mu1*mu2*weights) - sum(mu1*weights)*sum(mu2*weights)
      }
    }

    if (non.add=="dom") {
      gamma <- (geno@ploidy/2 - 1)/(geno@ploidy - 1)
      ix <- grep("dom",rownames(beta),fixed=T)
      out$params$heterosis <- data.frame(trait=traits,
                                         estimate=as.numeric(beta[ix,1]),
                                         SE=as.numeric(beta[ix,2]))

      weights <- id.weights/sum(id.weights)
      for (i in 1:n.trait)
        for (j in i:n.trait) {
          mu1d <- beta[ix[i],1]*dom.covariate
          mu2d <- beta[ix[j],1]*dom.covariate
          vars[i,j,"heterosis"] <- sum(mu1d*mu2d*weights) - sum(mu1d*weights)*sum(mu2d*weights)
          B[j,i] <- B[i,j] <- gamma^2*(mean(mu1d*mu2d) - mean(mu1d)*mean(mu2d))
        }
    } else {
      gamma <- 0
    }

    if (n.mark > 0) {
      trait.marker <- expand.grid(trait=traits,marker=fix.eff.marker,stringsAsFactors = FALSE)
      ix <- match(apply(trait.marker,1,paste,collapse=":"),rownames(beta))
      x <- out$params$marker <- data.frame(marker=trait.marker$marker,
                                       trait=trait.marker$trait,
                                       estimate=as.numeric(beta[ix,1]),
                                       SE=as.numeric(beta[ix,2]))
      weights <- id.weights/sum(id.weights)
      for (i in 1:n.trait) {
        for (j in i:n.trait) {
          b1 <- matrix(x[x$trait==traits[i],"estimate"],ncol=1)
          mu1 <- as.numeric(as.matrix(geno@coeff[id,fix.eff.marker])%*%b1)
          b2 <- matrix(x[x$trait==traits[j],"estimate"],ncol=1)
          mu2 <- as.numeric(as.matrix(geno@coeff[id,fix.eff.marker])%*%b2)
          vars[i,j,"fixed.marker"] <- sum(mu1*mu2*weights) - sum(mu1*weights)*sum(mu2*weights)
          if (non.add=="dom") {
            mu1 <- mu1+gamma*mu1d
            mu2 <- mu2+gamma*mu2d
          }
          B[j,i] <- B[i,j] <- mean(mu1*mu2) - mean(mu1)*mean(mu2)
        }
      }
    }
  }

  #variances
  vc <- sans$varcomp
  # Rescale if called
  if (rescale) {
    vc$component <- ans$vparameters
  }

  if (!is.null(vcov)) {
    vc <- vc[-which(vc$bound=="F" & !grepl("har!cor", row.names(vc))),]
  } else {
    vc <- vc[-which(row.names(vc) == "units!units"),]
  }

  vc.names <- rownames(vc)
  name2 <- vc.names
  name2 <- gsub("vm(id, source = asremlG, singG = \"PSD\")","additive",name2,fixed=T)
  name2 <- gsub("vm(id, source = asremlD)","dominance",name2,fixed=T)
  name2 <- gsub("units!units","residual",name2,fixed=T)
  name2 <- gsub("units","residual",name2,fixed=T)
  name2 <- gsub("har.id","residual",name2,fixed=T)
  name2 <- gsub("Trait!Trait!","",name2,fixed=T)
  name2 <- gsub("Trait!Trait_","",name2,fixed=T)
  out$params$vc <- data.frame(name=name2,
                              estimate=vc$component, SE=vc$std.error)

  Imat <- Diagonal(n=n,x=1)
  dimnames(Imat) <- list(id,id)

  if (n.trait==1) {
    if (!is.null(vcov)) {
      resid.vc <- Matrix(vars[1,1,"Stage1.error"], nrow = 1, ncol = 1)
    } else {
      resid.vc <- Matrix(vc[grep("units",vc.names,fixed=T)[1],1],ncol=1)
    }

    if (n.loc > 1) {
      rownames(resid.vc) <- locations
      if (!is.null(vcov)) {
        vars[1,1,"g x env"] <- pvar(V=diag(as.numeric(resid.vc)),weights=loc.weights)
      } else {
        vars[1,1,"residual"] <- pvar(V=diag(as.numeric(resid.vc)),weights=loc.weights)
      }

      if (is.null(geno)) {
        if (n.loc > 2) {
          iz <- grep("id:fa",vc.names,fixed=T)
        } else {
          iz <- grep("id:loc",vc.names,fixed=T)
        }
        cov.ans <- f.cov.loc(vc=vc[iz,], locations)
        geno1.vc <- cov.ans$cov.mat
        geno2.vc <- Matrix(NA,nrow=0,ncol=0)
        model <- 0L
        vars[1,1,"genotype"] <- mean(geno1.vc[upper.tri(geno1.vc,diag=FALSE)])*(1 - 1/n)

        K <- kronecker(Imat,geno1.vc,make.dimnames = T)
        vars[1,1,"g x har"] <- pvar(V = K, weights=gH.weights) - vars[1,1,"genotype"]

      } else {
        cov.ans <- f.cov.loc(vc[grep("source = asremlG",vc.names,fixed=T),],locations)
        geno1.vc <- cov.ans$cov.mat
        geno2.vc <- Matrix(NA,nrow=0,ncol=0)
        vars[1,1,"additive"] <- mean(geno1.vc[upper.tri(geno1.vc,diag=FALSE)]) * meanG
        K <- kronecker(.GlobalEnv$asremlG,geno1.vc,make.dimnames = T)
        vars[1,1,"add x loc"] <- pvar(V=K,weights=gL.weights) - vars[1,1,"additive"]
        model <- 1L

        if (non.add=="g.resid") {
          tmp <- Matrix(vc[grep("id:loc",vc.names,fixed=T),1],ncol=1)
          rownames(tmp) <- locations
          geno2.vc <- coerce_dpo(tcrossprod(sqrt(tmp)))
          K <- kronecker(Imat,geno2.vc,make.dimnames = T)
          vars[1,1,"g.resid"] <- pvar(V=K,weights=gL.weights)
          model <- 2L
        }
        if (non.add=="dom") {
          tmp <- Matrix(vc[grep("source = asremlD):loc",vc.names,fixed=T),1],ncol=1)
          rownames(tmp) <- locations
          geno2.vc <- coerce_dpo(tcrossprod(sqrt(tmp)))
          K <- kronecker(.GlobalEnv$asremlD,geno2.vc,make.dimnames = T)
          vars[1,1,"dominance"] <- pvar(V=K,weights=gL.weights)
          model <- 3L
        }

      }
      out <- c(out,loadings=list(cov.ans$loadings))

    } else {

      #one location
      if (is.null(vcov)) {
        vars[1,1,"residual"] <- as.numeric(resid.vc)*(1 - 1/n.gE)
      }

      if (is.null(geno)) {
        # Get the g x har VCOV matrix
        model <- 0L
        # Create a matrix of the cor structure is none
        if (har.cor.str == "idt") {
          cov.ans <- diag(ans$vparameters["id"], nrow = n.har, ncol = n.har)
        } else {
          cov.ans <- getVarCov.asreml(asreml.obj = ans, which.matrix = "G",
                                      ignore.terms = setdiff(names(ans$G.param), "id:har"),
                                      ignore.vars = list("id:har" = c("id")))
        }
        geno1.vc <- cov.ans * ans$sigma2
        dimnames(geno1.vc) <- list(hars, hars)
        geno1.vc <- coerce_dpo(geno1.vc)
        geno2.vc <- Matrix(NA,nrow=0,ncol=0)
        if (har.cor.str == "idt") {
          vars[1,1,"genotype"] <- mean(diag(geno1.vc))
        } else {
          vars[1,1,"genotype"] <- mean(geno1.vc[upper.tri(geno1.vc,diag=FALSE)])*(1 - 1/n)
          K <- kronecker(Imat,geno1.vc,make.dimnames = T)
          vars[1,1,"g x har"] <- pvar(V = K, weights = gH.weights) - vars[1,1,"genotype"]
        }

        vars[1,1,"g x har corr"] <- vc["id:har!har!cor",1]

      } else {
        # Create a matrix of the cor structure is none
        if (har.cor.str == "idt") {
          cov.ans <- diag(ans$vparameters["id"], nrow = n.har, ncol = n.har)
        } else {
          terms <- names(ans$G.param)
          terms_ignore <- terms[!grepl("id", terms)]
          vars_ignore <- lapply(ans$G.param[setdiff(terms, terms_ignore)],
                                function(x) names(x)[grepl("id", names(x))])
          cov.ans <- getVarCov.asreml(asreml.obj = ans, which.matrix = "G",
                                      ignore.terms = terms_ignore,
                                      ignore.vars = vars_ignore)
          dimnames(cov.ans) <- list(hars,hars)

        }

        geno1.vc <- cov.ans
        geno2.vc <- Matrix(NA,nrow=0,ncol=0)
        vars[1,1,"additive"] <- mean(geno1.vc[upper.tri(geno1.vc,diag=FALSE)]) * meanG
        K <- kronecker(.GlobalEnv$asremlG,geno1.vc,make.dimnames = T)
        vars[1,1,"add x har"] <- pvar(V=K,weights=gH.weights) - vars[1,1,"additive"]
        model <- 1L

        if (non.add=="g.resid") {
          tmp <- Matrix(vc[grep("id:.*har",vc.names,fixed=T),1],ncol=1)
          rownames(tmp) <- hars
          geno2.vc <- coerce_dpo(tcrossprod(sqrt(tmp)))
          K <- kronecker(Imat,geno2.vc,make.dimnames = T)
          vars[1,1,"g.resid"] <- pvar(V=K,weights=gL.weights)
          model <- 2L
        }
        if (non.add=="dom") {
          tmp <- Matrix(vc[grep("source = asremlD):loc",vc.names,fixed=T),1],ncol=1)
          rownames(tmp) <- locations
          geno2.vc <- coerce_dpo(tcrossprod(sqrt(tmp)))
          K <- kronecker(.GlobalEnv$asremlD,geno2.vc,make.dimnames = T)
          vars[1,1,"dominance"] <- pvar(V=K,weights=gL.weights)
          model <- 3L
        }



        model <- 1L
        geno1.vc <- Matrix(vc[grep("source = asremlG",vc.names,fixed=T),1],ncol=1)
        geno2.vc <- Matrix(NA,nrow=0,ncol=0)
        vars[1,1,"additive"] <- as.numeric(geno1.vc)*meanG
        if (non.add=="g.resid") {
          geno2.vc <- Matrix(vc["id",1],ncol=1)
          model <- 2L
          vars[1,1,"g.resid"] <- as.numeric(geno2.vc)*(1 - 1/n)
        }
        if (non.add=="dom") {
          geno2.vc <- Matrix(vc[grep("source = asremlD",vc.names,fixed=T),1],ncol=1)
          model <- 3L
          vars[1,1,"dominance"] <- as.numeric(geno2.vc)*meanD
        }
      }
    }
  } else {
    #multi-trait

    iv <- grep("env.id:Trait!",vc.names,fixed=T)
    resid.vc <- f.cov.trait(vc[iv,],traits,us=TRUE)
    if (!is.null(vcov)) {
      tmp <- "g x env"
    } else {
      tmp <- "residual"
    }

    for (i in 1:n.trait)
      for (j in i:n.trait)
        vars[i,j,tmp] <- resid.vc[i,j]*(1 - 1/n.gE)

    vc <- vc[-iv,]
    if (is.null(geno)) {
      model <- 0L
      geno2.vc <- Matrix(NA,nrow=0,ncol=0)
      geno1.vc <- f.cov.trait(vc[grep("id:Trait!",rownames(vc),fixed=T),],
                              traits,us=(n.trait>2))
      for (i in 1:n.trait)
        for (j in i:n.trait)
          vars[i,j,"genotype"] <- geno1.vc[i,j]*(1 - 1/n)

    } else {
      model <- 1L
      geno2.vc <- Matrix(NA,nrow=0,ncol=0)
      geno1.vc <- f.cov.trait(vc[grep("source = asremlG",rownames(vc),fixed=T),],
                            traits,us=(n.trait>2))
      for (i in 1:n.trait)
        for (j in i:n.trait)
          vars[i,j,"additive"] <- geno1.vc[i,j]*meanG

      if (non.add=="g.resid") {
        geno2.vc <- f.cov.trait(vc[grep("id:Trait!",rownames(vc),fixed=T),],
                              traits,us=(n.trait>2))
        model <- 2L
        for (i in 1:n.trait)
          for (j in i:n.trait)
            vars[i,j,"g.resid"] <- geno2.vc[i,j]*(1 - 1/n)

      }
      if (non.add=="dom") {
        geno2.vc <- f.cov.trait(vc[grep("source = asremlD",rownames(vc),fixed=T),],
                            traits,us=(n.trait>2))
        model <- 3L
        for (i in 1:n.trait)
          for (j in i:n.trait)
            vars[i,j,"dominance"] <- geno2.vc[i,j]*meanD

      }
    }
  }

  if (n.mark==0) {
    fix.eff.marker <- character(0)
  }

  if (n.trait > 1) {
    for (i in seq_len(dim(vars)[3]))
      vars[,,i][lower.tri(vars[,,i])] <- vars[,,i][upper.tri(vars[,,i])]
  }

  out$vars <- new(Class="class_var",geno1=geno1.vc,geno2=geno2.vc,model=model,
                  resid=resid.vc,diagG=diagG,diagD=diagD,
                  vars=vars,B=B,fix.eff.marker=fix.eff.marker)

  #random effects
  if (n.trait==1) {
    # id har names
    id.har.names <- expand.grid(id = id, har = hars, stringsAsFactors = FALSE)
    u <- sans$coef.random
    if (n.loc > 1) {
      id.loc.names <- expand.grid(loc=locations,id=id,stringsAsFactors = F)[,c(2,1)]
      if (!is.null(geno)) {
        if (n.loc==2) {
          id.loc <- expand.grid(paste0("loc_",locations),
                                paste0("vm(id, source = asremlG, singG = \"PSD\")_",id),
                                stringsAsFactors = F)[,c(2,1)]
        } else {
          id.loc <- expand.grid(paste0("fa(loc, 2)_",locations),
                                paste0("vm(id, source = asremlG, singG = \"PSD\")_",id),
                                stringsAsFactors = F)[,c(2,1)]
        }
        ix1 <- match(apply(id.loc,1,paste,collapse=":"),rownames(u))
        out$random <- data.frame(id.loc.names,add=as.numeric(u[ix1,1]))

        if (non.add=="g.resid") {
          id.loc <- expand.grid(paste0("loc_",locations),
                                paste0("id_",id),stringsAsFactors = F)[,c(2,1)]
          ix2 <- match(apply(id.loc,1,paste,collapse=":"),rownames(u))
          out$random$g.iid=as.numeric(u[ix2,1])
        }
        if (non.add=="dom") {
          id.loc <- expand.grid(paste0("loc_",locations),
                                paste0("vm(id, source = asremlD)_",id),stringsAsFactors = F)[,c(2,1)]
          ix2 <- match(apply(id.loc,1,paste,collapse=":"),rownames(u))
          out$random$dom=as.numeric(u[ix2,1])
        }
      } else {
        if (n.loc==2) {
          id.loc <- expand.grid(paste0("loc_",locations),paste0("id_",id),stringsAsFactors = F)[,c(2,1)]
        } else {
          id.loc <- expand.grid(paste0("fa(loc, 2)_",locations),paste0("id_",id),stringsAsFactors = F)[,c(2,1)]
        }
        ix <- match(apply(id.loc,1,paste,collapse=":"),rownames(u))
        out$random <- data.frame(id.loc.names,g.iid=as.numeric(u[ix,1]))
      }
    } else {
      if (!is.null(geno)) {
        ix1 <- match(paste("vm(id, source = asremlG, singG = \"PSD\")",id,sep="_"),rownames(u))
        out$random <- data.frame(id=id,add=as.numeric(u[ix1,1]))
        if (non.add=="g.resid") {
          ix2 <- match(paste("id",id,sep="_"),rownames(u))
          out$random$g.iid=as.numeric(u[ix2,1])
        }
        if (non.add=="dom") {
          ix2 <- match(paste("vm(id, source = asremlD)",id,sep="_"),rownames(u))
          out$random$dom=as.numeric(u[ix2,1])
        }
      } else {
        # Get id - har BLUPs
        id.har.use <- expand.grid(paste0("id_",id), paste0("har_", hars), stringsAsFactors = F)
        ix <- match(apply(id.har.use,1, paste, collapse=":"), rownames(u))
        out$random <- data.frame(id.har.names, g.iid = as.numeric(u[ix,1]))
      }
    }
  } else {
    #multi-trait
    u <- sans$coef.random
    id.trait.names <- expand.grid(trait=traits,id=id,stringsAsFactors = F)[,c(2,1)]

    if (is.null(geno)) {
      id.trait <- expand.grid(paste0("Trait_",traits),paste0("id_",id),stringsAsFactors = F)[,c(2,1)]
      ix <- match(apply(id.trait,1,paste,collapse=":"),rownames(u))
      out$random <- data.frame(id.trait.names,g.iid=as.numeric(u[ix,1]))
    } else {
      id.trait <- expand.grid(paste0("Trait_",traits),paste0("vm(id, source = asremlG, singG = \"PSD\")_",id),
                              stringsAsFactors = F)[,c(2,1)]
      ix1 <- match(apply(id.trait,1,paste,collapse=":"),rownames(u))
      out$random <- data.frame(id.trait.names,add=as.numeric(u[ix1,1]))
      if (non.add=="g.resid") {
        id.trait <- expand.grid(paste0("Trait_",traits),paste0("id_",id),stringsAsFactors = F)[,c(2,1)]
        ix2 <- match(apply(id.trait,1,paste,collapse=":"),rownames(u))
        out$random$g.iid=as.numeric(u[ix2,1])
      }
      if (non.add=="dom") {
        id.trait <- expand.grid(paste0("Trait_",traits),paste0("vm(id, source = asremlD)_",id),
                                stringsAsFactors = F)[,c(2,1)]
        ix2 <- match(apply(id.trait,1,paste,collapse=":"),rownames(u))
        out$random$dom=as.numeric(u[ix2,1])
      }
    }
  }

  if (!is.null(geno))
    rm("asremlG",envir = .GlobalEnv)
  if (!is.null(vcov))
    rm("asremlOmega",envir = .GlobalEnv)
  if (non.add=="dom")
    rm("asremlD",envir = .GlobalEnv)

  return(out)
}
