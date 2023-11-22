#' Stage 1 analysis of multi-year trials in perennial crops
#'
#' Computes genotype BLUEs for each experiment
#'
#' The input file must have one column labeled "id" for the individuals, one labeled "trl" for each trial, and one labeled "har" for the harvest year or other repeated factor. The data for each trial are analyzed independently with a linear mixed model. Although not used in Stage1, to include a genotype x location effect in \code{\link{Stage2}}, a column labeled "loc" should be present in the input file.
#'
#' Argument \code{effects} is used to specify other i.i.d. effects besides genotype and has three columns: name, fixed, factor. The "name" column is a string that must match a column in the input file. The fixed column is a logical variable to indicate whether the effect is fixed (TRUE) or random (FALSE). The factor column is a logical variable to indicate whether the effect is a factor (TRUE) or numeric (FALSE).
#'
#' Argument \code{solver} specifies which software to use for REML. Current options are "asreml" and "spats". For "spats", the argument \code{spline} must be a vector of length two, with the names of the x and y variables (respectively) for the 2D spline.
#'
#' The heritability and residuals in the output are based on a random effects model for id.
#'
#' Missing response values are omitted for single-trait analysis but retained for multi-trait analysis (unless both traits are missing), to allow for prediction in Stage 2.
#'
#' Argument \code{workspace} is a vector of length two containing the workspace and pworkspace limits for ASReml-R, with default values of 500mb. If you get an error about insufficient memory, try increasing the appropriate value (workspace for variance estimation and pworkspace for BLUE computation).
#'
#' For multiple traits, only "asreml" is supported, and only the BLUE model is run, so the returned object does not contain H2.
#'
#' If the input file has a column "expt", this indicates multiple experiments within environment, which may be needed when using spatial analyses. Each experiment is first analyzed separately, and then the BLUEs from all experiments in one env are jointly analyzed to compute a single BLUE per env. The estimation errors from each experiment are propagated into the multi-expt model.
#'
#' @param filename Name of CSV file. If missing, assumes that \code{data} is passed.
#' @param data A \code{data.frame} with the columns described above for \code{filename}. If missing, assumes \code{filename} is passed.
#' @param traits trait names (see Details)
#' @param effects data frame specifying other effects in the model (see Details)
#' @param cor.effects data frame specifying serially correlated effects in the model (see Details)
#' @param solver one of the following: "asreml","spats"
#' @param spline vector of variable names for 2D spline with SpATS
#' @param residual Optional custom residual formula (for asreml only)
#' @param silent TRUE/FALSE, whether to suppress REML output
#' @param workspace memory limits for ASRreml-R
#'
#' @return List containing
#' \describe{
#' \item{blues}{data frame of BLUEs}
#' \item{vcov}{list of variance-covariance matrices for the BLUEs, one per experiment (env)}
#' \item{fit}{data frame with broad-sense H2 (plot basis) and/or AIC}
#' \item{resid}{For single trait, list of diagnostic plots and data frame of residuals. For multi-trait, list of resid var-cov matrices.}
#' }
#'
#' @importFrom utils combn read.csv
#' @import ggplot2
#' @import SpATS
#' @importFrom ggpubr ggarrange
#' @importFrom rlang .data
#' @importFrom stats resid
#' @importFrom spam as.dgCMatrix.spam
#' @import Matrix
#'
#' @export
#'
SerialStage1 <- function(filename, data, traits, effects = NULL, cor.effects = NULL,
                         solver=c("asreml"), spline=NULL, residual = NULL, silent = TRUE,

                         pred.terms, workspace=c("500mb","500mb")) {

  if (missing(data)) {
    data <- read.csv(file=filename,check.names=F)
  } else if (!missing(data)) {
    stopifnot(is.data.frame(data))
  } else {
    stop("One of 'filename' or 'data' must be passed.")
  }

  solver <- match.arg(solver)
  solver <- toupper(solver)
  stopifnot(traits %in% colnames(data))

  if (!is.null(residual)) stopifnot(inherits(residual, "formula"))
  n.trait <- length(traits)

  stopifnot(requireNamespace("asreml"))
  library(asreml)
  asreml::asreml.options(maxit=30,workspace=workspace[1],pworkspace=workspace[2],trace=FALSE)

  if (solver=="SPATS") {
    if (n.trait > 1) {
      stop("Use asreml for multiple traits")
    }
    stopifnot(requireNamespace("SpATS"))
    stopifnot(spline %in% colnames(data))
  }

  if (!is.null(effects)) {
    # Check effect names in data
    for (nm in effects$name) {
      # If interaction; split
      nm_split <- strsplit(x = nm, split = ":")[[1]]
      check_nm <- sapply(nm_split, function(nm) nm %in% colnames(data))
      if (any(!check_nm)) stop(cat("The following 'effects' terms were not found in 'data':", names(check_nm)[!check_nm]))
    }
  }
  if (!is.null(cor.effects)) {
    # Check effect names in data
    for (nm in cor_effects$name) {
      # If interaction; split
      nm_split <- strsplit(x = nm, split = ":")[[1]]
      check_nm <- sapply(nm_split, function(nm) nm %in% colnames(data))
      if (any(!check_nm)) stop(cat("The following 'effects' terms were not found in 'data':", names(check_nm)[!check_nm]))
    }
  }
  stopifnot(c("id","trl", "har") %in% colnames(data))

  if (!is.element("expt",colnames(data))) {
    data$expt <- data$trl
    no.expt <- TRUE
  } else {
    no.expt <- FALSE
    tmp <- split(data$trl, data$expt)
    tmp2 <- sapply(tmp,function(z){length(unique(z))})
    if (any(tmp2 > 1)) {
      stop("expt names must be unique within env")
    }
  }

  data$trl <- as.character(data$trl)
  data$expt <- as.character(data$expt)
  data$id <- as.character(data$id)
  expt.og <- unique(data$expt)

#   iz <- apply(as.matrix(data[,traits]),1,function(z){!all(is.na(z))})
#   data <- data[iz,]

  tmp <- split(data[,c("id", "har")], data$expt)
  replicated <- sapply(tmp,function(x){
    y <- table(table(x))
    any(as.integer(names(y)) > 1)
  })
  expts <- names(which(replicated))
  data <- data[data$expt %in% expts,]

  if (nrow(data)==0) {
    stop("No experiments with replication")
  }
  trls <- unique(data$trl)
  expt.in.trl <- lapply(split(data$expt,factor(data$trl,levels=trls)),unique)
  expt.missing <- setdiff(expt.og,expts)

  if (length(expt.missing)>1) {
    cat("Some experiments removed due to missing data or lack of replication:\n")
    cat(paste0(paste(expt.missing,collapse="\n"),"\n"))
  }
  n.expt <- length(expts)
  n.trl <- length(trls)

  effect.table <- matrix("",nrow=2,ncol=2)
  rownames(effect.table) <- c("blue","blup")
  colnames(effect.table) <- c("fixed","random")

  if (is.null(residual)) {
    if (n.trait == 1) {
      residual <- "id(units)"
    } else {
      residual <- "id(units):us(trait)"
    }
  } else {
    residual <- as.character(residual)[2]
  }

  # Number of harvests per experiment
  har_exp_tab <- table(data$expt, data$har)
  har_exp_tab <- ifelse(har_exp_tab > 1, 1, 0)

  # COunt the number of hars
  mult.har <- any(rowSums(har_exp_tab) > 1)

  if (n.trait == 1) {
    if (solver=="ASREML") {
      model <- sub("traits",traits,
                 "asreml(data=data1,na.action=na.method(x='omit'),fixed=traits~FIX,random=~RANDOM,residual=~RESIDUAL)",fixed=T)
      model <- sub(pattern = "RESIDUAL", replacement = residual, x = model)

      if (mult.har) {
        effect.table[1,1] <- "id + har + id:har"
        effect.table[2,1] <- "har"
        effect.table[2,2] <- "+ id + id:har"
      } else {
        effect.table[1,1] <- "id"
        effect.table[2,2] <- "id"
      }
    }
    if (solver=="SPATS")
      model <- sub("traits",traits,
                   "SpATS(data=data1,response='traits',genotype='id',fixed=~FIX,random=~RANDOM,spatial=~SAP(spline.x,spline.y),genotype.as.random=GARgar",fixed=T)
  } else {
    model <- sub("response",paste(traits,collapse=","),
                 "asreml(data=data1,na.action=na.method(x='omit'),fixed=cbind(response)~FIX,random=~RANDOM,residual=~RESIDUAL)",fixed=T)
    model <- sub(pattern = "RESIDUAL", replacement = residual, x = model)
    effect.table[1,1] <- "id:trait"
  }

  factor.vars <- "id"

  if (!is.null(effects)) {
    if (n.trait > 1) {
      effects$name2 <- apply(array(effects$name),1,paste,"trait",sep=":")
    } else {
      effects$name2 <- effects$name
    }
    for (i in 1:nrow(effects)) {
      term_i <- effects$name2[i]
      if (effects$fixed[i]) {
        # If the effect is fixed but contains "id", move it to random for the BLUP model
        if (grepl(pattern = factor.vars, x = term_i)) {
          # Keep it as fixed in the blue model
          effect.table[1, 1] <- paste(effect.table[1,1], term_i, sep="+")
          # Do not make it fixed in the blup model; make it random
          effect.table[2, 2] <- paste(effect.table[2,2], term_i, sep="+")
        } else {
          effect.table[,1] <- paste(effect.table[,1],term_i,sep="+")
        }
      } else {
        effect.table[,2] <- paste(effect.table[,2],term_i,sep="+")
      }
    }



    factor.vars <- c(factor.vars,effects$name[which(effects$factor)])
    numeric.vars <- effects$name[which(!effects$factor)]
  } else {
    numeric.vars <- character(0)
  }
  n.numeric <- length(numeric.vars)

  # Parse the cor.effects table
  if (!is.null(cor.effects)) {
    if (n.trait > 1) {
      cor.effects$name2 <- apply(array(cor.effects$name),1,paste,"trait",sep=":")
    } else {
      cor.effects$name2 <- cor.effects$name
    }
    for (i in 1:nrow(cor.effects)) {
      cor_term_i <- paste0(cor.effects$name2[i], ":", cor.effects$str[i], "(har)")
      effect.table[,2] <- paste(effect.table[,2], cor_term_i, sep="+")
    }
    factor.vars <- c(factor.vars,cor.effects$name[which(cor.effects$factor)])
  }

  #eliminate leading "+"
  effect.table <- apply(effect.table,c(1,2),function(z){if(substr(z,1,1)=="+"){substr(z,2,nchar(z))}else{z}})

  #BLUE model
  if (effect.table[1,2]=="") {
    blue.model <- sub("random=~RANDOM,","",model,fixed=T)
  } else {
    blue.model <- sub("RANDOM",effect.table[1,2],model,fixed=T)
  }
  if (effect.table[1,1]=="") {
    blue.model <- sub("fixed=~FIX,","",blue.model,fixed=T)
  } else {
    blue.model <- sub("FIX",effect.table[1,1],blue.model,fixed=T)
  }

  #BLUP model
  if (solver=="ASREML" & effect.table[2,1]=="") {
    effect.table[2,1] <- "1"
  }
  if (effect.table[2,2]=="") {
    blup.model <- sub("random=~RANDOM,","",model,fixed=T)
  } else {
    blup.model <- sub("RANDOM",effect.table[2,2],model,fixed=T)
  }
  if (effect.table[2,1]=="") {
    blup.model <- sub("fixed=~FIX,","",blup.model,fixed=T)
  } else {
    blup.model <- sub("FIX",effect.table[2,1],blup.model,fixed=T)
  }
  if (solver=="SPATS") {
    blup.model <- sub("GARgar","TRUE",blup.model,fixed=T)
    blue.model <- sub("GARgar","FALSE",blue.model,fixed=T)
    blup.model <- sub("spline.x",spline[1],blup.model,fixed=T)
    blup.model <- sub("spline.y",spline[2],blup.model,fixed=T)
    blue.model <- sub("spline.x",spline[1],blue.model,fixed=T)
    blue.model <- sub("spline.y",spline[2],blue.model,fixed=T)
    if (silent) {
      blup.model <- paste0(blup.model,",control=list(monitoring=0))")
      blue.model <- paste0(blue.model,",control=list(monitoring=0))")
    } else {
      blup.model <- paste0(blup.model,")")
      blue.model <- paste0(blue.model,")")
    }
  }

  resid.blup <- NULL
  vcov <- vector("list",n.expt)
  names(vcov) <- expts

  if ("loc" %in% colnames(data)) {
    fit <- data[!duplicated(data$expt),c("loc","trl","expt")]
  } else {
    fit <- data[!duplicated(data$expt),c("trl","expt")]
  }

  fit <- fit[match(expts,fit$expt),]
  if (n.trait==1) {
    fit$H2 <- as.numeric(NA)
    if (solver=="ASREML")
      fit$AIC <- as.numeric(NA)
  } else {
    fit$AIC <- as.numeric(NA)
  }
  if (no.expt) fit <- fit[,-match("expt",colnames(fit))]

  blue.out <- NULL
  blup.resid <- NULL
  resid.vc <- vector("list",n.expt)
  names(resid.vc) <- expts
  if (solver=="SPATS") {
    spatial.plot <- vector("list",n.expt)
    names(spatial.plot) <- expts
  }
  if (!silent)  cat(sub("X",paste(traits,collapse=" "),"Traits: X\n"))
  for (j in 1:n.expt) {
    if (!silent) cat(sub("X",expts[j],"expt: X\n"))

    ix <- which(data$expt==expts[j])
    data1 <- as.data.frame(data[ix,])
    hars <- unique(data1$har)
    n.har <- length(hars)

    # Split and unlist factor.vars
    factor.vars <- unique(unlist(strsplit(factor.vars, ":")))

    for (q in 1:length(factor.vars)) {
      eval(parse(text="data1[,factor.vars[q]] <- factor(as.character(data1[,factor.vars[q]]))"))
    }
    if (n.numeric > 0) {
      for (q in 1:n.numeric) {
        eval(parse(text="data1[,numeric.vars[q]] <- as.numeric(data1[,numeric.vars[q]])"))
      }
    }

    #BLUP model
    if (n.trait==1) {
      ans <- try(eval(parse(text=blup.model)),silent=TRUE)
      # If class error; print the error
      if (class(ans) == "try-error") {
        cat("Error in BLUP model: ", ans)
        next
      } else if (solver == "ASREML" && (!ans$converge)) {
        cat("BLUP model failed to converge.\n")
        next
      }
      residuals <- resid(ans)
      blup.resid <- rbind(blup.resid,
                          data.frame(id=as.character(data1$id),expt=expts[j],resid=residuals))
      if (solver=="ASREML") {
        # vc <- summary(ans)$varcomp
        # Vg <- vc[match("id",rownames(vc)),1]
        # Vge <- vc[match("id:har",rownames(vc)),1]
        # Ve <- vc[match("units!units",rownames(vc)),1]
        #         fit$H2[j] <- round(Vg/(Vg+Ve+Vge),2)

        # Calculate heritability using the variance of genotype differences according to Cullis 2006
        vc <- summary(ans)$varcomp
        Vg <- vc[match("id",rownames(vc)),1]
        vdBLUP.mat <- predict(ans, classify="id", only="id", sed=TRUE)$sed^2
        vdBLUP.avg <- mean(vdBLUP.mat[upper.tri(vdBLUP.mat, diag=FALSE)]) # take mean of upper triangle
        H2Cullis <- 1 - (vdBLUP.avg / 2 / Vg)
        fit$H2[j] <- round(H2Cullis, 2)
      }
      if (solver=="SPATS") {
        fit$H2[j] <- round(as.numeric(getHeritability(ans)),2)
        x.coord <- ans$data[,ans$terms$spatial$terms.formula$x.coord]
        y.coord <- ans$data[,ans$terms$spatial$terms.formula$y.coord]
        fit.spatial.trend <- obtain.spatialtrend(ans)
        p1.data <- data.frame(x=x.coord,y=y.coord,z=residuals)
        p1 <- ggplot(p1.data,aes(x=.data$x,y=.data$y,fill=.data$z)) + geom_tile() + scale_fill_viridis_c(name="") + xlab(spline[1]) + ylab(spline[2]) + ggtitle("Residuals")
        p2.data <-data.frame(expand.grid(x=fit.spatial.trend$col.p, y=fit.spatial.trend$row.p), z=as.numeric(t(fit.spatial.trend$fit)))
        p2 <- ggplot(p2.data,aes(x=.data$x,y=.data$y,fill=.data$z)) + geom_tile() + scale_fill_viridis_c(name="") + xlab(spline[1]) + theme(axis.text.y = element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank()) + ggtitle("Spatial Trend")
        spatial.plot[[j]] <- ggarrange(p1,p2,common.legend = TRUE,legend = "right")
      }
    }

    #BLUE model
    asreml::asreml.options(trace=!silent)
    ans <- try(eval(parse(text=blue.model)),silent=TRUE)
    if ((class(ans)=="try-error") || ((solver=="ASREML")&&(!ans$converge))) {
      cat("BLUE model failed to converge.\n")
    } else {
      asreml::asreml.options(trace=FALSE)
      vc <- summary(ans)$varcomp

      if (n.trait==1) {
        if (solver=="ASREML") {
          predans <- asreml::predict.asreml(ans, classify=pred.terms,vcov = TRUE)
          tmp <- predans$pvals[,c(strsplit(pred.terms, ":")[[1]], "predicted.value")]
          colnames(tmp) <- c(strsplit(pred.terms, ":")[[1]], "BLUE")
          vcov[[j]] <- predans$vcov
          fit$AIC[j] <- round(as.numeric(summary(ans)$aic),1)
        }
        if (solver=="SPATS") {
          predans <- predict.SpATS(object=ans,which="id",predFixed="marginal",
                                   return.vcov.matrix = TRUE)
          tmp <- predans[,c("id","predicted.values")]
          colnames(tmp) <- c("id","BLUE")
          tmp2 <- Matrix(spam::as.dgCMatrix.spam(attr(predans,"vcov")))
          vcov[[j]] <- as(forceSymmetric(tmp2),"packedMatrix")
        }
        for (col in strsplit(pred.terms, ":")[[1]]) {
          tmp[[col]] <- as.character(tmp[[col]])
        }
        tmp_dimnames <- tmp[strsplit(pred.terms, ":")[[1]]]
        vcov_dimnames <- apply(X = tmp_dimnames, MARGIN = 1, FUN = paste0, collapse = ":")

        dimnames(vcov[[j]]) <- list(vcov_dimnames, vcov_dimnames)
        blue.out <- rbind(blue.out,data.frame(expt=expts[j],tmp))

      } else {
        vc <- vc[-which(vc$bound=="F" & round(vc$component)==1L),c("component","std.error")]
        vc.names <- rownames(vc)
        iv <- grep("units:trait!",vc.names,fixed=T)
        resid.vc[[j]] <- f.cov.trait(vc[iv,],traits,us=TRUE)

        fit$AIC[j] <- round(as.numeric(summary(ans)$aic),1)
        predans <- asreml::predict.asreml(ans,classify="id:trait",vcov = TRUE)
        tmp <- predans$pvals[,c("id","trait","predicted.value")]
        colnames(tmp) <- c("id","trait","BLUE")
        vcov[[j]] <- predans$vcov
        tmp$id <- as.character(tmp$id)
        id.trait <- apply(tmp[,1:2],1,paste,collapse=":")
        dimnames(vcov[[j]]) <- list(id.trait,id.trait)
        blue.out <- rbind(blue.out,data.frame(expt=expts[j],tmp))
      }
    }
  }

  ik <- which(!sapply(vcov,is.null))
  vcov <- vcov[ik]
  expt.in.trl <- expt.in.trl[ik]
  blue.out$trl <- data$trl[match(blue.out$expt,data$expt)]

  nee <- sapply(expt.in.trl,length)
  iu <- which(nee==1)
  if (length(iu) > 0) {
    tmp <- names(vcov)[iu]
    names(vcov)[iu] <- blue.out$trl[match(tmp,blue.out$expt)]
  }

  for (j in which(nee > 1)) {

    omega.list <- vcov[expt.in.trl[[j]]]
    vcov2 <- mapply(FUN=function(x,y){
      tmp <- paste(y,rownames(x),sep=":")
      dimnames(x) <- list(tmp,tmp)
      return(x)},x=omega.list,y=as.list(expt.in.trl[[j]]))

    .GlobalEnv$asremlOmega <- direct_sum(lapply(vcov2,solve))
    dname <- lapply(vcov2,rownames)
    dimnames(.GlobalEnv$asremlOmega) <- list(unlist(dname),unlist(dname))
    attr(.GlobalEnv$asremlOmega,"INVERSE") <- TRUE

    ix <- which(blue.out$trl==trls[j])
    data2 <- blue.out[ix,]
    data2$expt.id <- factor(paste(data2$expt,data2$id,sep=":"))
    data2$id <- factor(data2$id)
    data2$expt <- factor(data2$expt)
    # Convert expt to numeric
    data2$expt_n <- as.numeric(data2$expt)
    # Weights vector
    wts <- do.call("c", map(vcov2, diag))
    data2$smith.w <- 1 / wts

    ## Prepare the model but do not fit it
    initial <- asreml(fixed = BLUE ~ expt+ id,
                      random = ~ vm(expt.id, source = asremlOmega),
                      residual = ~ idv(units),
                      data = data2,
                      na.action = na.method(y = "omit", x = "omit"),
                      # family = asr_gaussian(dispersion = 1.0),
                      # weights = smith.w,
                      start.values = TRUE)

    # Edit the components
    start.table <- initial$vparameters.table
    k <- grep("Omega",start.table$Component,fixed=T)
    start.table$Value[k] <- 1
    start.table$Constraint[k] <- "F"

    # Fit
    ans3 <- asreml(fixed = BLUE ~ expt + id,
                   random = ~ vm(expt.id, source = asremlOmega),
                   residual = ~ idv(units),
                   data = data2,
                   na.action = na.method(y = "omit", x = "omit"),
                   # family = asr_gaussian(dispersion = 1.0),
                   # weights = smith.w,
                   G.param = start.table)

    # Predict overall id means
    predans3 <- asreml::predict.asreml(ans3, classify="id", vcov = TRUE)
    blue3 <- data.frame(expt=NA,predans3$pvals[,c("id","predicted.value")],env=envs[j])

    colnames(blue3) <- c("expt","id","BLUE","env")
    vcov3 <- predans3$vcov
    dimnames(vcov3) <- list(blue3$id,blue3$id)
    vcov2 <- vcov[-match(expt.in.trl[[j]], names(vcov))]
    vcov <- c(vcov2, vcov3)
    names(vcov) <- c(names(vcov2),envs[j])
    blue.out <- rbind(blue.out[-ix,],blue3)
    rm("asremlOmega",envir = .GlobalEnv)
  }

  if (n.trait==1) {
    blue.out <- blue.out[,c("trl", strsplit(pred.terms, ":")[[1]], "BLUE")]
  } else {
    blue.out <- blue.out[,c("trl",strsplit(pred.terms, ":")[[1]], "trait","BLUE")]
  }
  if ("har" %in% colnames(blue.out)) {
    blue.out <- blue.out[order(blue.out$trl, blue.out$id, blue.out$har),]
  } else {
    blue.out <- blue.out[order(blue.out$trl, blue.out$id),]
  }
  if ("loc" %in% colnames(data)) {
    blue.out$loc <- as.character(data$loc[match(blue.out$trl, data$trl)])
  }

  if (n.trait==1) {
    p1 <- ggplot(data=blup.resid,aes(y=.data$resid,x=.data$expt)) + ylab("Residual") + xlab("") +
      stat_boxplot(outlier.color="red") + theme_bw() + theme(axis.text.x=element_text(angle=90,vjust=0.5))
    p2 <- ggplot(data=blup.resid,aes(sample=.data$resid)) + stat_qq() + stat_qq_line() + facet_wrap(~expt) + theme_bw() + xlab("Expected") + ylab("Observed")
    if (solver=="ASREML")
      return(list(blues=blue.out,vcov=vcov,fit=fit, vc = vc,
                  resid=list(boxplot=p1,qqplot=p2,table=blup.resid)))
    if (solver=="SPATS")
      return(list(blues=blue.out,vcov=vcov,fit=fit,
                  resid=list(boxplot=p1,qqplot=p2,spatial=spatial.plot,table=blup.resid)))
  } else {
    #Multi-trait
    return(list(blues=blue.out,vcov=vcov,fit=fit,resid=resid.vc))
  }
}
