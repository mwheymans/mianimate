
rm(list = ls())

library(haven)
library(mice)


dataset <- read_sav(file="Backpain50 MI missing.sav")
dataset <- dataset[, c("Pain", "Tampascale")]
data = dataset
m = 1
method = c("", "norm")
maxit = 10 
predictorMatrix = matrix(c(0,1,1,0), 2,2,byrow=T)
dimnames(predictorMatrix) <- list(c(names(dataset)), c(names(dataset)))

where = NULL 
visitSequence = NULL
blots = NULL
post = NULL 
defaultMethod = c("norm.predict")

printFlag = TRUE
seed = 1232
data.init = NULL

#### Extra

check.data <- function(data, method) {
  check.dataform(data)
  
}

check.dataform <- function(data) {
  if (!(is.matrix(data) || is.data.frame(data)))
    stop("Data should be a matrix or data frame", call. = FALSE)
  if (ncol(data) < 2)
    stop("Data should contain at least two columns", call. = FALSE)
  data <- as.data.frame(data)
  mat <- sapply(data, is.matrix)
  if (any(mat)) stop("Cannot handle columns with class matrix: ", 
                     colnames(data)[mat])
  
  dup <- duplicated(colnames(data))
  if (any(dup)) stop("Duplicate names found: ", 
                     paste(colnames(data)[dup], collapse = ", "))
  
  data
}

check.m <- function(m) {
  m <- m[1L]
  if (!is.numeric(m))
    stop("Argument m not numeric", call. = FALSE)
  m <- floor(m)
  if (m < 1L)
    stop("Number of imputations (m) lower than 1.", call. = FALSE)
  m
}

check.cluster <- function(data, predictorMatrix) {
  # stop if the cluster variable is a factor
  isclassvar <- apply(predictorMatrix == -2, 2, any)
  for (j in colnames(predictorMatrix)) {
    if (isclassvar[j] && lapply(data, is.factor)[[j]]) 
      stop("Convert cluster variable ", j, " to integer by as.integer()")
  }
  TRUE
}

make.predictorMatrix <- function(data, blocks = make.blocks(data)) {
  data <- check.dataform(data)
  predictorMatrix <- matrix(1, nrow = length(blocks), ncol = ncol(data))
  dimnames(predictorMatrix) <- list(names(blocks), colnames(data))
  for (i in row.names(predictorMatrix)) 
    predictorMatrix[i, grep(i, colnames(predictorMatrix))] <- 0
  predictorMatrix
}

check.predictorMatrix <- function(predictorMatrix, 
                                  data,
                                  blocks = NULL) {
  data <- check.dataform(data)
  
  if (!is.matrix(predictorMatrix))
    stop("predictorMatrix not a matrix", call. = FALSE)
  if (any(dim(predictorMatrix) == 0L))
    stop("predictorMatrix has no rows or columns", call. = FALSE)
  
  # if we have no blocks, restrict to square predictorMatrix
  if (is.null(blocks)) {
    if (nrow(predictorMatrix) != ncol(predictorMatrix))
      stop(paste("If no blocks are specified, predictorMatrix must", 
                 "have same number of rows and columns"), 
           call. = FALSE)
    if (is.null(dimnames(predictorMatrix))) {
      if (ncol(predictorMatrix) == ncol(data)) 
        dimnames(predictorMatrix) <- list(colnames(data), colnames(data))
      else
        stop("Missing row/column names in predictorMatrix", call. = FALSE)
    }
    for (i in row.names(predictorMatrix))
      predictorMatrix[i, grep(i, colnames(predictorMatrix))] <- 0
    return(predictorMatrix)
  }
  
  # check conforming arguments
  if (nrow(predictorMatrix) > length(blocks))
    stop(paste0("predictorMatrix has more rows (", nrow(predictorMatrix), 
                ") than blocks (", length(blocks), ")"),
         call. = FALSE)
  
  # borrow rownames from blocks if needed
  if (is.null(rownames(predictorMatrix)) && 
      nrow(predictorMatrix) == length(blocks))
    rownames(predictorMatrix) <- names(blocks)
  if (is.null(rownames(predictorMatrix)))
    stop("Unable to set row names of predictorMatrix", call. = FALSE)
  
  # borrow blocknames from predictorMatrix if needed
  if (is.null(names(blocks)) &&
      nrow(predictorMatrix) == length(blocks))
    names(blocks) <- rownames(predictorMatrix)
  if (is.null(names(blocks)))
    stop("Unable to set names of blocks", call. = FALSE)
  
  # check existence of row names in blocks
  found <- rownames(predictorMatrix) %in% names(blocks)
  if (!all(found))
    stop("Names not found in blocks: ", 
         paste(rownames(predictorMatrix)[!found], collapse = ", "), 
         call. = FALSE)
  
  # borrow colnames from data if needed
  if (is.null(colnames(predictorMatrix)) && 
      ncol(predictorMatrix) == ncol(data)) 
    colnames(predictorMatrix) <- names(data)
  if (is.null(colnames(predictorMatrix))) 
    stop("Unable to set column names of predictorMatrix", call. = FALSE)
  
  # check existence of variable names on data
  found <- colnames(predictorMatrix) %in% names(data) 
  if (!all(found))
    stop("Names not found in data: ", 
         paste(colnames(predictorMatrix)[!found], collapse = ", "), 
         call. = FALSE)
  
  list(predictorMatrix = predictorMatrix,
       blocks = blocks)
}
make.where <- function(data, 
                       keyword = c("missing", "all", "none", "observed")) {
  keyword <- match.arg(keyword)
  
  data <- check.dataform(data)
  where <- switch(keyword,
                  missing = is.na(data),
                  all = matrix(TRUE, nrow = nrow(data), ncol = ncol(data)),
                  none = matrix(FALSE, nrow = nrow(data), ncol = ncol(data)), 
                  observed = !is.na(data))
  
  dimnames(where) <- dimnames(data)
  where
}


check.where <- function(where, data, blocks) {
  if (is.null(where)) 
    where <- make.where(data, keyword = "missing")
  
  if (!(is.matrix(where) || is.data.frame(where)))
    if (is.character(where)) return(make.where(data, keyword = where))
  else
    stop("Argument `where` not a matrix or data frame", call. = FALSE)
  if (!all(dim(data) == dim(where)))
    stop("Arguments `data` and `where` not of same size", call. = FALSE)
  
  where <- as.logical(as.matrix(where))
  if (anyNA(where))
    stop("Argument `where` contains missing values", call. = FALSE)
  
  where <- matrix(where, nrow = nrow(data), ncol = ncol(data))
  dimnames(where) <- dimnames(data)
  where[, !colnames(where) %in% unlist(blocks)] <- FALSE
  where
}

check.method <- function(method, data, where, blocks, defaultMethod) {
  
  if (is.null(method)) return(make.method(data = data, 
                                          where = where,
                                          blocks = blocks,
                                          defaultMethod = defaultMethod))
  nimp <- nimp(where, blocks)
  
  # expand user's imputation method to all visited columns
  # single string supplied by user (implicit assumption of two columns)
  if (length(method) == 1) {
    if (is.passive(method))
      stop("Cannot have a passive imputation method for every column.")
    method <- rep(method, length(blocks))
    method[nimp == 0] <- ""
  }
  
  # check the length of the argument
  if (length(method) != length(blocks))
    stop("Length of method differs from number of blocks", call. = FALSE)
  
  # add names to method
  names(method) <- names(blocks)
  
  # check whether the requested imputation methods are on the search path
  active.check <- !is.passive(method) & nimp > 0 & method != ""
  passive.check <- is.passive(method) & nimp > 0 & method != ""
  check <- all(active.check) & any(passive.check)
  if (check) {
    fullNames <- rep.int("mice.impute.passive", length(method[passive.check]))
  } else {
    fullNames <- paste("mice.impute", method[active.check], sep = ".")
    if (length(method[active.check]) == 0) fullNames <- character(0)
  }
  notFound <- !vapply(fullNames, exists, logical(1), 
                      mode = "function", inherits = TRUE)
  if (any(notFound)) {
    stop(paste("The following functions were not found:",
               paste(fullNames[notFound], collapse = ", ")))
  }
  
  # type checks on built-in imputation methods
  for (j in names(blocks)) {
    vname <- blocks[[j]]
    y <- data[, vname, drop = FALSE]
    mj <- method[j]
    mlist <- list(m1 = c("logreg", "logreg.boot", "polyreg", "lda", "polr"), 
                  m2 = c("norm", "norm.nob", "norm.predict", "norm.boot",
                         "mean", "2l.norm", "2l.pan",
                         "2lonly.pan", "quadratic", "ri"), 
                  m3 = c("norm", "norm.nob", "norm.predict", "norm.boot",
                         "mean", "2l.norm", "2l.pan", 
                         "2lonly.pan", "quadratic", "logreg", "logreg.boot"))
    cond1 <- sapply(y, is.numeric)
    cond2 <- sapply(y, is.factor) & sapply(y, nlevels) == 2 
    cond3 <- sapply(y, is.factor) & sapply(y, nlevels) > 2
    if (any(cond1) && mj %in% mlist$m1)
      warning("Type mismatch for variable(s): ", 
              paste(vname[cond1], collapse = ", "),
              "\nImputation method ", mj, " is for categorical data.",
              call. = FALSE)
    if (any(cond2) && mj %in% mlist$m2)
      warning("Type mismatch for variable(s): ", 
              paste(vname[cond2], collapse = ", "),
              "\nImputation method ", mj, " is not for factors.", 
              call. = FALSE)
    if (any(cond3) && mj %in% mlist$m3)
      warning("Type mismatch for variable(s): ", 
              paste(vname[cond3], collapse = ", "),
              "\nImputation method ", mj, " is not for factors with >2 levels.",
              call. = FALSE)
  }
  method[nimp == 0] <- ""
  unlist(method)
}
check.visitSequence <- function(visitSequence = NULL, 
                                data, where = NULL, blocks) {
  
  if (is.null(names(blocks)) || any(is.na(names(blocks))))
    stop("Missing names in `blocks`.")
  
  if (is.null(visitSequence)) return(make.visitSequence(data, blocks))
  
  if (is.null(where)) where <- is.na(data)
  nimp <- nimp(where, blocks)
  if (length(nimp) == 0) visitSequence <- nimp
  
  if (length(visitSequence) == 1 && is.character(visitSequence)) {
    code <- match.arg(visitSequence, c("roman", "arabic", "monotone",
                                       "revmonotone"))
    visitSequence <- switch(
      code, 
      roman = names(blocks)[nimp > 0],
      arabic = rev(names(blocks)[nimp > 0]),
      monotone = names(blocks)[order(nimp)],
      revmonotone = rev(names(blocks)[order(nimp)])
    )
  }
  
  # legacy handling
  if (is.numeric(visitSequence)) 
    visitSequence <- colnames(data)[visitSequence]
  
  # check against names(blocks)
  visitSequence <- visitSequence[is.element(visitSequence, names(blocks))]
  
  # remove any blocks without missing data
  visitSequence <- names((nimp > 0L)[visitSequence])
  visitSequence
}
check.post <- function(post, data) {
  
  if(is.null(post)) return(make.post(data))
  
  # check
  if (length(post) != ncol(data))
    stop("length(post) does not match ncol(data)", call. = FALSE)
  
  # change
  if (is.null(names(post))) names(post) <- colnames(data)
  
  post
}
check.blots <- function(blots, data, blocks = NULL) {
  data <- check.dataform(data)
  
  if (is.null(blots)) return(make.blots(data, blocks))
  
  blots <- as.list(blots)
  for (i in seq_along(blots)) blots[[i]] <- as.list(blots[[i]])
  
  if (length(blots) == length(blocks) && is.null(names(blots)))
    names(blots) <- names(blocks)
  blots
}
is.passive <- function(string)
{
  return("~" == substring(string, 1, 1))
}

make.method <- function(data, where = make.where(data), blocks = make.blocks(data),
                        defaultMethod = c("pmm", "logreg", "polyreg", "polr")) {
  
  assign.method <- function(y) {
    if (is.numeric(y)) return(1)
    if (nlevels(y) == 2) return(2)
    if (is.ordered(y) && nlevels(y) > 2) return(4)
    if (nlevels(y) > 2) return(3)
    if (is.logical(y)) return(2)
    return(1)
  }
  
  # assign methods based on type, 
  # use method 1 if there is no single method within the block
  method <- rep("", length(blocks))
  names(method) <- names(blocks)
  for (j in names(blocks)) {
    yvar <- blocks[[j]]
    y <- data[, yvar]
    def <- sapply(y, assign.method)
    k <- ifelse(all(diff(def) == 0), k <- def[1], 1)
    method[j] <- defaultMethod[k]
  }
  nimp <- nimp(where, blocks)
  method[nimp == 0] <- ""
  method
}

check.method <- function(method, data, where, blocks, defaultMethod) {
  
  if (is.null(method)) return(make.method(data = data, 
                                          where = where,
                                          blocks = blocks,
                                          defaultMethod = defaultMethod))
  nimp <- nimp(where, blocks)
  
  # expand user's imputation method to all visited columns
  # single string supplied by user (implicit assumption of two columns)
  if (length(method) == 1) {
    if (is.passive(method))
      stop("Cannot have a passive imputation method for every column.")
    method <- rep(method, length(blocks))
    method[nimp == 0] <- ""
  }
  
  # check the length of the argument
  if (length(method) != length(blocks))
    stop("Length of method differs from number of blocks", call. = FALSE)
  
  # add names to method
  names(method) <- names(blocks)
  
  # check whether the requested imputation methods are on the search path
  active.check <- !is.passive(method) & nimp > 0 & method != ""
  passive.check <- is.passive(method) & nimp > 0 & method != ""
  check <- all(active.check) & any(passive.check)
  if (check) {
    fullNames <- rep.int("mice.impute.passive", length(method[passive.check]))
  } else {
    fullNames <- paste("mice.impute", method[active.check], sep = ".")
    if (length(method[active.check]) == 0) fullNames <- character(0)
  }
  notFound <- !vapply(fullNames, exists, logical(1), 
                      mode = "function", inherits = TRUE)
  if (any(notFound)) {
    stop(paste("The following functions were not found:",
               paste(fullNames[notFound], collapse = ", ")))
  }
  
  # type checks on built-in imputation methods
  for (j in names(blocks)) {
    vname <- blocks[[j]]
    y <- data[, vname, drop = FALSE]
    mj <- method[j]
    mlist <- list(m1 = c("logreg", "logreg.boot", "polyreg", "lda", "polr"), 
                  m2 = c("norm", "norm.nob", "norm.predict", "norm.boot",
                         "mean", "2l.norm", "2l.pan",
                         "2lonly.pan", "quadratic", "ri"), 
                  m3 = c("norm", "norm.nob", "norm.predict", "norm.boot",
                         "mean", "2l.norm", "2l.pan", 
                         "2lonly.pan", "quadratic", "logreg", "logreg.boot"))
    cond1 <- sapply(y, is.numeric)
    cond2 <- sapply(y, is.factor) & sapply(y, nlevels) == 2 
    cond3 <- sapply(y, is.factor) & sapply(y, nlevels) > 2
    if (any(cond1) && mj %in% mlist$m1)
      warning("Type mismatch for variable(s): ", 
              paste(vname[cond1], collapse = ", "),
              "\nImputation method ", mj, " is for categorical data.",
              call. = FALSE)
    if (any(cond2) && mj %in% mlist$m2)
      warning("Type mismatch for variable(s): ", 
              paste(vname[cond2], collapse = ", "),
              "\nImputation method ", mj, " is not for factors.", 
              call. = FALSE)
    if (any(cond3) && mj %in% mlist$m3)
      warning("Type mismatch for variable(s): ", 
              paste(vname[cond3], collapse = ", "),
              "\nImputation method ", mj, " is not for factors with >2 levels.",
              call. = FALSE)
  }
  method[nimp == 0] <- ""
  unlist(method)
}
handles.format <- function(fn) {
  # determine whether function fn handles the `format` argument
  f <- get(fn)
  handles.arg(f, "format")
}
handles.arg <- function(f, a = "data") {
  # determine whether function f handles argument a
  if (!is.function(f)) return(FALSE)
  a %in% names(formals(f))
}
edit.setup <- function(data, setup,
                       allow.na = FALSE, 
                       remove.constant = TRUE,
                       remove.collinear = TRUE, 
                       remove_collinear = TRUE) {
  # legacy handling
  if (!remove_collinear) remove.collinear <- FALSE
  
  # edits the imputation model setup
  # When it detec constant or collinear variables, write in loggedEvents 
  # and continues imputation with reduced model
  
  pred <- setup$predictorMatrix
  meth <- setup$method
  vis <- setup$visitSequence
  post <- setup$post
  
  # FIXME: this function is not yet adapted to blocks
  if (ncol(pred) != nrow(pred) || length(meth) != nrow(pred) 
      || ncol(data) != nrow(pred))
    return(setup)
  
  varnames <- colnames(data)
  
  # remove constant variables but leave passive variables untouched
  for (j in seq_len(ncol(data))) {
    if (!is.passive(meth[j])) {
      d.j <- data[, j]
      v <- if (is.character(d.j)) NA else var(as.numeric(d.j), na.rm = TRUE)
      constant <- if (allow.na) {
        if (is.na(v)) FALSE else v < 1000 * .Machine$double.eps
      } else {
        is.na(v) || v < 1000 * .Machine$double.eps
      }
      didlog <- FALSE
      if (constant && any(pred[, j] != 0) && remove.constant) {
        out <- varnames[j]
        pred[, j] <- 0
        updateLog(out = out, meth = "constant")
        didlog <- TRUE
      }
      if (constant && meth[j] != "" && remove.constant) {
        out <- varnames[j]
        pred[j, ] <- 0
        if (!didlog)
          updateLog(out = out, meth = "constant")
        meth[j] <- ""
        vis <- vis[vis != j]
        post[j] <- ""
      }
    }
  }
  
  ## remove collinear variables
  ispredictor <- apply(pred != 0, 2, any)
  if (any(ispredictor)) {
    droplist <- find.collinear(data[, ispredictor, drop = FALSE], ...)
  } else {
    droplist <- NULL
  }
  if (length(droplist) > 0) {
    for (k in seq_along(droplist)) {
      j <- which(varnames %in% droplist[k])
      didlog <- FALSE
      if (any(pred[, j] != 0) && remove.collinear) {
        # remove as predictor
        out <- varnames[j]
        pred[, j] <- 0
        updateLog(out = out, meth = "collinear")
        didlog <- TRUE
      }
      if (meth[j] != "" && remove.collinear) {
        out <- varnames[j]
        pred[j, ] <- 0
        if (!didlog)
          updateLog(out = out, meth = "collinear")
        meth[j] <- ""
        vis <- vis[vis != j]
        post[j] <- ""
      }
    }
  }
  
  if (all(pred == 0L)) stop("nothing left to impute")
  
  setup$predictorMatrix <- pred
  setup$visitSequence <- vis
  setup$post <- post
  setup$method <- meth
  return(setup)
}
obtain.design <- function(data, formula = ~ .) {
  
  mf <- model.frame(formula, data = data, na.action = na.pass)
  model.matrix(formula, data = mf)
}

update.design <- function(design, data, varname = ".") {
  # Updates columns of the design matrix related to variable
  # varname in data
  
  varname <- as.character(varname[1])
  idx <- attr(design, "assign") %in% grep(varname, names(data))
  
  # variable j not found
  if (varname == "" || !any(idx)) return(design)
  
  # create model frame of variable j only
  fj <- as.formula(paste("~", varname))
  mfj <- model.frame(fj, data = data, na.action = na.pass)
  design[, idx] <- model.matrix(fj, data = mfj)[, -1, drop = FALSE]
  design
}

check.df <- function(x, y, ry) {
  # if needed, writes the df warning message to the log
  df <- sum(ry) - ncol(x) - 1
  mess <- paste("df set to 1. # observed cases:", sum(ry), " # predictors:", ncol(x) + 1)
  if (df < 1 && sum(ry) > 0)
    updateLog(out = mess, frame = 4)
}

###### End Extra

#function (data, m = 5, method = NULL, predictorMatrix, where = NULL, 
#          blocks, visitSequence = NULL, formulas, blots = NULL, post = NULL, 
#          defaultMethod = c("pmm", "logreg", "polyreg", "polr"), maxit = 5, 
#          printFlag = TRUE, seed = NA, data.init = NULL, ...) 
#{
  call <- match.call()
  if (!is.na(seed)) 
    set.seed(seed)
  data <- check.dataform(data)
  m <- check.m(m)
  mp <- missing(predictorMatrix)
  mb <- TRUE #missing(blocks)
  mf <- TRUE #missing(formulas)
  if (mp & mb & mf) {
    blocks <- make.blocks(colnames(data))
    predictorMatrix <- make.predictorMatrix(data, blocks)
    formulas <- make.formulas(data, blocks)
  }
  if (!mp & mb & mf) {
    predictorMatrix <- check.predictorMatrix(predictorMatrix, 
                                             data)
    blocks <- make.blocks(colnames(predictorMatrix), partition = "scatter")
    formulas <- make.formulas(data, blocks, predictorMatrix = predictorMatrix)
  }
  if (mp & !mb & mf) {
    blocks <- check.blocks(blocks, data)
    predictorMatrix <- make.predictorMatrix(data, blocks)
    formulas <- make.formulas(data, blocks)
  }
  if (mp & mb & !mf) {
    formulas <- check.formulas(formulas, data)
    blocks <- construct.blocks(formulas)
    predictorMatrix <- make.predictorMatrix(data, blocks)
  }
  if (!mp & !mb & mf) {
    blocks <- check.blocks(blocks, data)
    z <- check.predictorMatrix(predictorMatrix, data, blocks)
    predictorMatrix <- z$predictorMatrix
    blocks <- z$blocks
    formulas <- make.formulas(data, blocks, predictorMatrix = predictorMatrix)
  }
  if (!mp & mb & !mf) {
    formulas <- check.formulas(formulas, data)
    predictorMatrix <- check.predictorMatrix(predictorMatrix, 
                                             data)
    blocks <- construct.blocks(formulas, predictorMatrix)
  }
  if (mp & !mb & !mf) {
    blocks <- check.blocks(blocks, data, calltype = "formula")
    formulas <- check.formulas(formulas, blocks)
    predictorMatrix <- make.predictorMatrix(data, blocks)
  }
  if (!mp & !mb & !mf) {
    blocks <- check.blocks(blocks, data)
    formulas <- check.formulas(formulas, data)
    predictorMatrix <- check.predictorMatrix(predictorMatrix, 
                                             data, blocks)
  }
  chk <- check.cluster(data, predictorMatrix)
  where <- check.where(where, data, blocks)
  visitSequence <- check.visitSequence(visitSequence, data = data, 
                                       where = where, blocks = blocks)
  method <- check.method(method = method, data = data, where = where, 
                         blocks = blocks, defaultMethod = defaultMethod)
  post <- check.post(post, data)
  blots <- check.blots(blots, data, blocks)
  state <- list(it = 0, im = 0, dep = "", meth = "", log = FALSE)
  loggedEvents <- data.frame(it = 0, im = 0, dep = "", meth = "", 
                             out = "")
  setup <- list(method = method, predictorMatrix = predictorMatrix, 
                visitSequence = visitSequence, post = post)
  setup <- edit.setup(data, setup)
  method <- setup$method
  predictorMatrix <- setup$predictorMatrix
  visitSequence <- setup$visitSequence
  post <- setup$post
  nmis <- apply(is.na(data), 2, sum)
  source("initialize.imp.R")
  imp <- initialize.imp(data, m, where, blocks, visitSequence, 
                        method, nmis, data.init)
  from <- 1
  to <- from + maxit - 1
  source("sampler2.R")
  q <- sampler2(data, m, where, imp, blocks, method, visitSequence, 
               predictorMatrix, formulas, blots, post, c(from, to), 
               printFlag)
  q
  
  colMeans(q)
  
  
  if (!state$log) 
    loggedEvents <- NULL
  if (state$log) 
    row.names(loggedEvents) <- seq_len(nrow(loggedEvents))
  midsobj <- list(data = data, imp = q$imp, m = m, where = where, 
                  blocks = blocks, call = call, nmis = nmis, method = method, 
                  predictorMatrix = predictorMatrix, visitSequence = visitSequence, 
                  formulas = formulas, post = post, blots = blots, seed = seed, 
                  iteration = q$iteration, lastSeedValue = .Random.seed, 
                  chainMean = q$chainMean, chainVar = q$chainVar, loggedEvents = loggedEvents, 
                  version = packageVersion("mice"), date = Sys.Date())
  oldClass(midsobj) <- "mids"
  if (!is.null(midsobj$loggedEvents)) 
    warning("Number of logged events: ", nrow(midsobj$loggedEvents), 
            call. = FALSE)
  return(midsobj)
}
#rm(list = ls())

dataset <- read_sav(file="Backpain50 MI missing.sav")

names(dataset)

dataset <- dataset[, c("Pain", "Tampascale")]

imp.regress <- mice(dataset, m=1, maxit=10, method=c("", "norm"), seed = 1232)
imp.regress$chainMean
data.frame(imp.regress$chainMean[2,,1])

