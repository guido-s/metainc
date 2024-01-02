# Auxiliary functions
#
# Package: metainc
# Author: Guido Schwarzer <guido.schwarzer@uniklinik-freiburg.de>
# License: GPL (>= 2)
#

chklogical <- function(x, name = NULL) {
  ##
  ## Check whether argument is logical
  ##
  if (is.null(name))
    name <- deparse(substitute(x))
  ##
  if (is.numeric(x))
    x <- as.logical(x)
  ##
  if (length(x) !=  1 || !is.logical(x) || is.na(x))
    stop("Argument '", name, "' must be a logical.", call. = FALSE)
  ##
  invisible(NULL)
}

chkchar <- function(x, length = 0, name = NULL, nchar = NULL, single = FALSE,
                    NULL.ok = FALSE) {
  if (is.null(x) & NULL.ok)
    return(invisible(NULL))
  ##
  if (!missing(single) && single)
    length <- 1
  if (is.null(name))
    name <- deparse(substitute(x))
  ##
  if (length && length(x) != length) {
    if (length == 1)
      stop("Argument '", name, "' must be a character string.",
           call. = FALSE)
    else
      stop("Argument '", name, "' must be a character vector of length ",
           length, ".",
           call. = FALSE)
  }
  ##
  if (length == 1) {
    if (!is.null(nchar) && !(nchar(x) %in% nchar))
      if (length(nchar) == 1 && nchar == 1)
        stop("Argument '", name, "' must be a single character.",
             call. = FALSE)
    else
      stop("Argument '", name, "' must be a character string of length ",
           if (length(nchar) == 2)
             paste0(nchar, collapse = " or ")
           else
             paste0(nchar, collapse = ", "),
           ".",
           call. = FALSE)
  }
  ##
  if (!is.character(x) & !is.numeric(x))
    stop("Argument '", name, "' must be a character vector.")
  else {
    if (!is.null(nchar) & any(!(nchar(x) %in% nchar)))
      if (length(nchar) == 1 && nchar == 1)
        stop("Argument '", name, "' must be a vector of single characters.",
             call. = FALSE)
    else
      stop("Argument '", name, "' must be a character vector where ",
           "each element has ",
           if (length(nchar) == 2)
             paste0(nchar, collapse = " or ")
           else
             paste0(nchar, collapse = ", "),
           " characters.",
           call. = FALSE)
  }
  ##
  invisible(NULL)
}

formatCI <- function(lower, upper,
                     bracket.left = gs("CIbracket"),
                     separator = gs("CIseparator"),
                     bracket.right,
                     justify.lower = "right",
                     justify.upper = justify.lower,
                     lower.blank = gs("CIlower.blank"),
                     upper.blank = gs("CIupper.blank"),
                     ...) {
  
  ## Change layout of CIs
  ##
  chkchar(bracket.left, length = 1)
  chkchar(separator, length = 1)
  if (!missing(bracket.right))
    chkchar(bracket.right, length = 1)
  ##
  if (missing(bracket.left)) {
    bracktype <- setchar(bracket.left, c("[", "(", "{", ""))
    ##
    if (bracktype == "[") {
      bracketLeft <- "["
      bracketRight <- "]"
    }
    else if (bracktype == "(") {
      bracketLeft <- "("
      bracketRight <- ")"
    }
    else if (bracktype == "{") {
      bracketLeft <- "{"
      bracketRight <- "}"
    }
    else if (bracktype == "") {
      bracketLeft <- ""
      bracketRight <- ""
    }
    ##
    bracket.left <- bracketLeft
  }
  ##
  if (missing(bracket.right))
    bracket.right <- bracketRight
  
  format.lower <- format(lower, justify = justify.lower)
  format.upper <- format(upper, justify = justify.upper)
  ##
  if (!lower.blank)
    format.lower <- rmSpace(format.lower)
  if (!upper.blank)
    format.upper <- rmSpace(format.upper)
  ##
  if (separator == "-")
    format.upper <-
    paste0(ifelse(substring(format.upper, 1, 1) == "-", " ", ""),
           format.upper)
  ##
  res <- ifelse(lower != "NA" & upper != "NA",
                paste0(bracket.left,
                       format.lower,
                       separator,
                       format.upper,
                       bracket.right),
                "")
  ##
  res
}

formatN <- function(x, digits = 2, text.NA = "--", big.mark = "",
                    format.whole.numbers = TRUE,
                    monospaced = FALSE) {
  
  outdec <- options()$OutDec  
  
  if (!monospaced) {
    if (format.whole.numbers) {
      res <- format(ifelse(is.na(x),
                           text.NA,
                           formatC(x, decimal.mark = outdec,
                                   format = "f", digits = digits,
                                   big.mark = big.mark)
      )
      )
    }
    else {
      res <- format(ifelse(is.na(x),
                           text.NA,
                           ifelse(is_wholenumber(x),
                                  x,
                                  formatC(x, decimal.mark = outdec,
                                          format = "f", digits = digits,
                                          big.mark = big.mark)
                           )
      )
      )
    }
  }
  else {
    x <- round(x, digits)
    res <- ifelse(is.na(x),
                  text.NA,
                  format(x, decimal.mark = outdec, big.mark = big.mark))
  }
  ##
  res <- rmSpace(res, end = TRUE)
  ##
  res
}

formatPT <- function(x, lab = FALSE, labval = "p", noblanks = FALSE,
                     digits = 4, zero = TRUE, scientific = FALSE,
                     lab.NA = "--", big.mark = "",
                     JAMA = FALSE) {
  
  if (is.null(x))
    return("")
  
  outdec <- options()$OutDec
  
  n.zeros <- digits - 1
  n.zeros[n.zeros < 0] <- 0
  
  if (!scientific) {
    if (lab) {
      if (!JAMA)
        res <- format(ifelse(is.na(x) | is.nan(x),
                             paste(labval, "=", lab.NA),
                             ifelse(x == 0,
                                    paste(labval, "= 0"),
                                    ifelse(x < 1 / 10^digits,
                                           paste0(labval, " < 0", outdec,
                                                  paste(rep("0",
                                                            n.zeros), collapse = ""),
                                                  "1"),
                                           paste(paste(labval, "="),
                                                 formatC(round(x, digits),
                                                         decimal.mark = outdec,
                                                         big.mark = big.mark,
                                                         format = "f", digits = digits)
                                           )
                                    )
                             )
        )
        )
      else
        res <- format(ifelse(is.na(x) | is.nan(x),
                             paste(labval, "=", lab.NA),
                             ifelse(x < 0.001,
                                    paste0(labval, " < 0", outdec,
                                           paste(rep("0", 2), collapse = ""), "1"),
                                    ifelse(x >= 0.001 & x < 0.01,
                                           paste(paste(labval, "="),
                                                 formatC(x,
                                                         decimal.mark = outdec,
                                                         big.mark = big.mark,
                                                         format = "f", digits = 3)),
                                           ifelse(x >= 0.01 & x <= 0.99,
                                                  paste(paste(labval, "="),
                                                        formatC(x,
                                                                decimal.mark = outdec,
                                                                big.mark = big.mark,
                                                                format = "f", digits = 2)),
                                                  paste(paste(labval, ">"),
                                                        formatC(0.99,
                                                                decimal.mark = outdec,
                                                                big.mark = big.mark,
                                                                format = "f", digits = 2)))
                                    )
                             )
        )
        )
      
    }
    else {
      if (!JAMA)
        res <- format(ifelse(is.na(x) | is.nan(x),
                             lab.NA,
                             ifelse(x == 0,
                                    0,
                                    ifelse(x < 1 / 10^digits,
                                           paste0("< 0", outdec,
                                                  paste(rep("0", n.zeros), collapse = ""),
                                                  "1"),
                                           formatC(round(x, digits),
                                                   decimal.mark = outdec,
                                                   big.mark = big.mark,
                                                   format = "f", digits = digits)
                                    )
                             )
        ),
        justify = "right")
      else
        res <- format(ifelse(is.na(x) | is.nan(x),
                             lab.NA,
                             ifelse(x < 0.001,
                                    paste0("< 0", outdec,
                                           paste(rep("0", 2), collapse = ""), "1"),
                                    ifelse(x >= 0.001 & x < 0.01,
                                           formatC(x,
                                                   decimal.mark = outdec,
                                                   big.mark = big.mark,
                                                   format = "f", digits = 3),
                                           ifelse(x >= 0.01 & x <= 0.99,
                                                  formatC(x,
                                                          decimal.mark = outdec,
                                                          big.mark = big.mark,
                                                          format = "f", digits = 2),
                                                  paste(">",
                                                        formatC(0.99,
                                                                decimal.mark = outdec,
                                                                big.mark = big.mark,
                                                                format = "f", digits = 2)))
                                    )
                             )
        ),
        justify = "right")
    }
  }
  else {
    if (lab)
      res <- format(ifelse(is.na(x) | is.nan(x),
                           paste(labval, "=", lab.NA),
                           paste(labval, "=",
                                 formatC(x, decimal.mark = outdec,
                                         big.mark = big.mark,
                                         format = "e", digits = digits)
                           )
      )
      )
    else
      res <- formatC(x, decimal.mark = outdec,
                     big.mark = big.mark, format = "e", digits = digits)
  }
  ##
  if (noblanks)
    res <- gsub(" ", "", res)
  if (!zero)
    res <- gsub("0\\.", "\\.", res)
  ##
  ## Treat NaNs as NAs
  ##
  res[grep("NaN", res)] <- lab.NA
  
  res
}

rmSpace <- function(x, end = FALSE, pat = " ") {
  
  if (!end) {
    while (any(substring(x, 1, 1) == pat, na.rm = TRUE)) {
      sel <- substring(x, 1, 1) == pat
      x[sel] <- substring(x[sel], 2)
    }
  }
  else {
    last <- nchar(x)
    
    while (any(substring(x, last, last) == pat, na.rm = TRUE)) {
      sel <- substring(x, last, last) == pat
      x[sel] <- substring(x[sel], 1, last[sel] - 1)
      last <- nchar(x)
    }
  }
  
  x
}

catch <- function(argname, matchcall, data, encl)
  eval(matchcall[[match(argname, names(matchcall))]], data, enclos = encl)


chklength <- function(x, k, text = "") {
  if (length(x) != k)
    stop(text, call. = FALSE)
  #
  invisible(NULL)
}


chknumeric <- function(x, min, max, zero = FALSE, length = 0,
                       name = NULL, single = FALSE, integer = FALSE) {
  if (!missing(single) && single)
    length <- 1
  ##
  ## Check numeric variable
  ##
  if (is.null(name))
    name <- deparse(substitute(x))
  ##
  x <- x[!is.na(x)]
  if (length(x) == 0)
    return(NULL)
  ##
  if (!is.numeric(x))
    stop("Non-numeric value for argument '", name, "'.",
         call. = FALSE)
  ##
  if (length && length(x) != length)
    stop("Argument '", name, "' must be a numeric of length ", length, ".",
         call. = FALSE)
  ##
  if (!missing(min) & missing(max)) {
    if (zero & min == 0 & any(x <= min, na.rm = TRUE))
      stop("Argument '", name, "' must be positive.",
           call. = FALSE)
    else if (any(x < min, na.rm = TRUE))
      stop("Argument '", name, "' must be larger equal ",
           min, ".", call. = FALSE)
  }
  ##
  if (missing(min) & !missing(max)) {
    if (zero & max == 0 & any(x >= max, na.rm = TRUE))
      if (zero & max == 0 & any(x >= max, na.rm = TRUE))
        stop("Argument '", name, "' must be negative.",
             call. = FALSE)
    else if (any(x > max, na.rm = TRUE))
      stop("Argument '", name, "' must be smaller equal ",
           min, ".", call. = FALSE)
  }
  ##
  if ((!missing(min) & !missing(max)) &&
      (any(x < min, na.rm = TRUE) | any(x > max, na.rm = TRUE)))
    stop("Argument '", name, "' must be between ",
         min, " and ", max, ".", call. = FALSE)
  ##
  if (integer && any(!is_wholenumber(x))) {
    if (length(x) == 1)
      stop("Argument '", name, "' must be an integer.",
           call. = FALSE)
    else
      stop("Argument '", name, "' may only contain integers.",
           call. = FALSE)
  }
  ##
  invisible(NULL)
}


chknull <- function(x)
  if (is.null(x))
    stop("Argument '", deparse(substitute(x)), "' is NULL.", call. = FALSE)


chkclass <- function(x, class, name = NULL) {
  ##
  ## Check class of R object
  ##
  if (is.null(name))
    name <- deparse(substitute(x))
  ##
  n.class <- length(class)
  if (n.class == 1)
    text.class <- paste0('"', class, '"')
  else if (n.class == 2)
    text.class <- paste0('"', class, '"', collapse = " or ")
  else
    text.class <- paste0(paste0('"', class[-n.class], '"', collapse = ", "),
                         ', or ', '"', class[n.class], '"')
  ##
  if (!inherits(x, class))
    stop("Argument '", name, "' must be an object of class ", text.class, ".",
         call. = FALSE)
  ##
  invisible(NULL)
}


setchar <- function(x, val, text, list = FALSE, name = NULL,
                    stop.at.error = TRUE, addtext = "",
                    return.NULL = TRUE, nchar.equal = FALSE,
                    setNA = FALSE) {
  val <- unique(val)
  #
  if (is.null(name))
    name <- deparse(substitute(x))
  nval <- length(val)
  #
  if (is.numeric(x)) {
    numeric.x <- TRUE
    idx <- x
    idx[idx < 1] <- NA
    idx[idx >= nval + 1] <- NA
  }
  else {
    numeric.x <- FALSE
    #
    if (length(unique(tolower(x))) != length(unique(x)) |
        length(unique(tolower(val))) != length(unique(val)))
      idx <- charmatch(x, val, nomatch = NA)
    else
      idx <- charmatch(tolower(x), tolower(val), nomatch = NA)
  }
  #
  if ((anyNA(idx) || any(idx == 0)) && !setNA) {
    if (list)
      first <- "List element '"
    else
      first <- "Argument '"
    #
    if (missing(text)) {
      if (numeric.x) {
        if (nval == 1)
          vlist <- "1"
        else if (nval == 2)
          vlist <- "1 or 2"
        else
          vlist <- paste("between 1 and", nval)
      }
      else {
        if (nval == 1)
          vlist <- paste0('"', val, '"')
        else if (nval == 2)
          vlist <- paste0('"', val, '"', collapse = " or ")
        else
          vlist <- paste0(paste0('"', val[-nval], '"', collapse = ", "),
                          ', or ', '"', val[nval], '"')
      }
      #
      if (stop.at.error)
        stop(first, name, "' must be ", vlist, addtext, ".", call. = FALSE)
      else {
        if (return.NULL)
          return(NULL)
        else
          return(x)
      }
    }
    else {
      if (stop.at.error)
        stop(first, name, "' ", text, ".", call. = FALSE)
      else {
        if (return.NULL)
          return(NULL)
        else
          return(x)
      }
    }
  }
  #
  if (is.null(x))
    return(NULL)
  else
    res <- val[idx]
  #
  if (nchar.equal && nchar(res) != nchar(x))
    res <- x
  #
  res
}

is_relative_effect <- function(x)
  x %in% c("HR", "OR", "RR", "IRR", "ROM", "DOR")

is_wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
  if (is.numeric(x))
    res <- abs(x - round(x)) < tol
  else
    res <- NA
  ##
  res
}
