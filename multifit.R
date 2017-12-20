multifit <- function(mod, multief, data = NULL, formula = NULL, args = NULL, criterion = "AIC",
                     signif = TRUE, alpha = 0.05, print_sum = FALSE, plot_est = FALSE,
                     xlab = "Radio [m]", labels = NULL, type = "b", pch = c(1, 16)){
  
  # Arguments checkings
  if(!is.character(mod) || length(mod) != 1) stop("Argument mod must be a character of length 1")
  if(!is.character(multief)) stop("Argument multief must be a character vector")
  if(!is.null(args)){ if(!is.character(args)) stop("Argument args must be NULL or a character vector") }
  if(is.character(criterion)){
    if(length(criterion) == 1){
      if(!criterion %in% c("AIC", "BIC", "R2")){
        stop("Argument criterion must be one of the following: 'AIC', 'BIC' or 'R2', or a user-defined function")
      }
    } else {
      if(!is.function(tryCatch(eval(parse(text = criterion[1])), error = function(e) FALSE))){
        stop("Argument criterion does not refers to an existing function of the global environment")
      } else {
        if(!criterion[2] %in% c("max", "min")){
          stop("For user-defined functions, the second element of argument criterion must be 'max' or 'min'")
        }
      }
    }
  } else {
    stop("Argument criterion must be a character vector (one of 'AIC', 'BIC' and 'R2', or a user-defined function)")
  }
  if(!is.logical(signif) || length(signif) != 1) stop("Argument signif must be logical")
  if(!is.numeric(alpha) || !(alpha > 0 && alpha <= 1) || length(alpha) != 1) stop("Argument alpha must be a number between 0 and 1")
  if(!is.logical(print_sum) || length(print_sum) != 1) stop("Argument print_sum must be logical")
  if(!is.logical(plot_est) || length(plot_est) != 1) stop("Argument plot_est must be logical")
  if(!is.null(labels)){ if(!is.character(as.character(labels))) stop("Argument labels must be NULL or a character vector") }
  if(!is.numeric(pch)) stop("Argument pch must be numeric")

  # Check if the function to be applied exists in the environment
  if(!exists(mod)) stop(paste("Could not find function '", mod, "'. Make sure that the necessary package is loaded", sep = ""))
  
  # Check if the provided names for the effects exist as columns of the provided data.frame
  if(!any(multief %in% colnames(data)))
    stop("Could not find any of the multi-effects as columns of the provided data. Are you sure is correctly written?")
  if(!all(multief %in% colnames(data)))
    warning("Could not find some of the multi-effects as columns of the provided data")
  
  # Objects definition
  multief      <- factor(multief, levels = multief)
  initial.args <- paste(args, collapse = ",")
  fits         <- vector("list", length(multief))
  fits.GoF     <- rep(0, length(multief))
  if(!is.null(formula)) 
    init.formula <- deparse(substitute(formula))
  sum.list     <- vector("list", length(multief))
  p.values     <- rep(0, length(multief))
  e.values     <- rep(0, length(multief))
  mod_warn     <- vector("list", length(multief))
  effect_warns <- NULL
  mod_mess     <- vector("list", length(multief))
  effect_mess  <- NULL
  mod_errors   <- vector()
  ele          <- FALSE

  # Check if 'multief' was a defined as a predictor variable
  if(!grepl("multief", deparse(substitute(formula)))){
    if(!grepl("multief", initial.args))
      stop("The formula, or the fixed effects defined in 'args', must include the expression 'multief' as a predictor variable")
  }
  
  # Function. Summary table of the predictor variable at each spatial scale
  lands_table <- function(multief, data){
    if(!is.null(data)){
      table <- data.frame(spatial_scale = as.character(multief), 
                          n             = rep(0, length(multief)),
                          min           = rep(0, length(multief)),
                          max           = rep(0, length(multief)),
                          range         = rep(0, length(multief)),
                          mean          = rep(0, length(multief)),
                          median        = rep(0, length(multief)))
      for(i in 1:length(multief)){
        unique_vec <- unique(data[, as.character(multief[i])])
        table[table$spatial_scale == multief[i], "n"]      <- length(na.omit(unique_vec))
        table[table$spatial_scale == multief[i], "min"]    <- min(na.omit(unique_vec))
        table[table$spatial_scale == multief[i], "max"]    <- max(na.omit(unique_vec))
        table[table$spatial_scale == multief[i], "range"]  <- range(na.omit(unique_vec))[2] - range(na.omit(unique_vec))[1]
        table[table$spatial_scale == multief[i], "mean"]   <- mean(na.omit(unique_vec))
        table[table$spatial_scale == multief[i], "median"] <- median(na.omit(unique_vec))
      }
    } else {
      table <- NA
      message("lands_summary could not be created as no data was specified")
    }
    return(table)
  }

  # Function. Extract the 'Good of Fitness' (AIC, BIC or R2) value for the models
  extractGoF <- function(fitted_model, criterion){
    aic   <- function(x){ tail(extractAIC(x), n = 1) }
    value <- NULL
    
    if(length(criterion) == 1){
      if(criterion == "AIC") 
        value <- AIC(fitted_model)
      if(criterion == "BIC")
        value <- BIC(fitted_model)
      if(criterion == "R2")
        value <- summary(fitted_model)$r.squared
    } else {
      value <- do.call(criterion[1], list(fitted_model))
    }

    if(is.null(value) || is.na(value)){
      if(length(criterion) == 1){
        if(grepl("quasipoisson", args) && criterion %in% c("AIC", "BIC")){
          stop("Quasipoisson models do not display AIC/BIC values.\n 
               We recommend running the models with a poisson family to get the AICs/BICs, proceed to the multi-scale analyis, 
               and then reran the selected model aside with a quasipoisson family to see the output.")
        } else {
          stop(paste("Could not get the", criterion, "value for the specified models. Try another criterion?"), call. = FALSE)
        }
      } else {
        stop(paste("Could not get the '", criterion[1], "' value for the specified models. Try another criterion?", sep = ""), call. = FALSE)
      }
    } else {
      return(value)
    }
  }

  # Function. Extract coefficients of the models
  extract_coeff <- function(summary, i){
    coeff   <- summary$coefficients
    # If the model is a zero-inflated one (zeroingl or hurdle)
    if(is.list(coeff) && names(coeff) == c("count", "zero")){
      opt <- names(coeff)
      if(ele == FALSE){
        ele <- readline(prompt = paste("Please, choose the type of model from where the coefficients will be extracted (", 
                                       paste(opt, ": ", 1:length(opt), sep = "", collapse = " - "),"): ", sep = ""))
        if(ele %in% as.character(1:length(opt))){
          ele   <- as.numeric(ele)
          coeff <- coeff[[ele]]
          ele   <<- ele
        } else {
          coeff <- coeff[[1]]
          ele <<- 1
          message("Could not recognise the inputted value. Option one was taken")
        }
      } else {
        coeff <- coeff[[ele]]
      }
    }
    
    # Trying to extract p.values...
    p.types <- c("Pr(>|z|)", "Pr(>|t|)")
    p.label <- p.types[p.types %in% colnames(coeff)]
    if(length(p.label) > 0){
      p.value <- coeff[as.character(multief[i]), p.label]
    } else {
      coeff   <- summary$tTable
      p.types <- "p-value"
      p.label <- p.types[p.types %in% colnames(coeff)]
      if(length(p.label) > 0){
        p.value <- coeff[as.character(multief[i]), p.label]
      } else {
        p.value <- NA
      }
    }
    
    # Trying to extract e.values...
    if("Estimate" %in% colnames(coeff)){
      e.value <- coeff[as.character(multief[i]), "Estimate"]
    } else {
      if("Value" %in% colnames(coeff)){
        e.value <- coeff[as.character(multief[i]), "Value"]
      } else {
        e.value <- NA
      }
    }
    out <- c(p.value, e.value)
    return(out)
  }
  
  # Function. Manage model's errors, warnings and messages
  running <- function(expr, i) {
    warns <- mess <- NULL
    warn_handler <- function(w) {
      warns <<- c(warns, list(w))
      invokeRestart("muffleWarning")
    }
    mess_handler <- function(m) {
      mess <<- c(mess, list(m))
      NULL
    }
    val <- suppressMessages(tryCatch(withCallingHandlers(expr, warning = warn_handler, message = mess_handler), error = function(e) e))
    out <- list(value = val, warnings = warns, messages = mess)
    return(out)
  } 

  # Run all models
  for(i in 1:length(multief)){
    
    # Replace the string 'multief' with the effect of each particular spatial scale
    args       <- gsub("multief", paste(multief[i]), initial.args)
    args       <- gsub("\"", "", args)
    if(!is.null(formula)){
      formula    <- gsub("multief", paste(multief[i]), init.formula)
      formula    <- gsub("\"", "", formula)
      expression <- paste(mod, "(", formula, ", data =", deparse(substitute(data)), ",", args, ")")
    } else {
      expression <- paste(mod, "(data =", deparse(substitute(data)), ",", args, ")")
    }

    new_fit <- running(eval(parse(text = expression)), i)
    if(!any(class(new_fit$value) %in% c("simpleError", "error"))){
      fits[[i]]          <- new_fit$value
      fits.GoF[i]        <- extractGoF(fits[[i]], criterion)
      sum.list[[i]]      <- summary(fits[[i]])
      coeff              <- extract_coeff(sum.list[[i]], i)
      p.values[i]        <- coeff[1]
      e.values[i]        <- coeff[2]
    } else {
      fits[[i]]   <- new_fit$value
      fits.GoF[i] <- sum.list[[i]] <- p.values[i] <- NA
      mod_errors  <- append(mod_errors, as.character(multief[i]))
    }
    names(fits)[i] <- as.character(multief[i])
    if(!is.null(new_fit$warnings)){
      mod_warn[[i]] <- new_fit$warnings
      effect_warns  <- append(effect_warns, as.character(multief[i]))
    } else {
      mod_warn[[i]] <- NA
    }
    if(!is.null(new_fit$messages)){
      mod_mess[[i]] <- new_fit$messages
      effect_mess   <- append(effect_mess, as.character(multief[i]))
    } else {
      mod_mess[[i]] <- NA
    }
    names(mod_warn)[i] <- names(mod_mess)[i] <- as.character(multief[i])
  }
  if(length(mod_errors) > 0){
    mod_errors <- paste(mod_errors, collapse = ", ")
    message(paste("The following model/s threw an error and could not be included in the analysis: ", mod_errors, sep = ""))
  }
    
  # Function. Significance extraction of estimates
  f_sign <- function(p.values){
    if(any(!is.na(p.values))){
      p.sign <- ifelse(p.values < alpha, TRUE, FALSE)
    } else {
      p.sign <- rep(TRUE, length(p.values))
    }
    return(p.sign)
  }

  # If there are GoF values to work with...
  if(any(!is.na(fits.GoF))){
    get_mfrow <- par()$mfrow
    if(plot_est){
      if(!all(is.na(e.values))){
        par(mfrow = c(1, 2))
      } else {
        plot_est <- FALSE
        message("Could not plot the estimates of the models cause they could not be extracted")
      }
    }
    # Plot the model selection result
    if(is.null(labels)){ labels <- multief }
    if(signif){
      if(!all(is.na(p.values))){
        p.sign <- f_sign(p.values)
      } else {
        p.sign <- rep(1, length(p.values))
        message("p.values could not be extracted")
      }
    } else {
      p.sign <- rep(0, length(p.values))
    }
    temp_df <- data.frame(x = 1:length(multief), y = fits.GoF, signif = as.numeric(p.sign + 1))
    plot(y ~ x, xlab = xlab, xaxt = "n", ylab = criterion[1], 
         type = type, pch = pch[signif], data = temp_df)
    axis(1, at = 1:length(multief), labels = labels)
    title(main = "Model selection")

    # Plot estimates, if required
    if(plot_est){
      temp_df <- data.frame(x = 1:length(multief), y = e.values, signif = as.numeric(p.sign + 1))
      plot(y ~ x, xlab = xlab, xaxt = "n", ylab = "Estimate", 
           type = type, pch = pch[signif], data = temp_df)
      abline(0, 0, col = "gray")
      axis(1, at = 1:length(multief), labels = labels)
      title(main = "Estimates")
    }
    
    # Record plot and restablish original par options
    plot <- recordPlot()
    par(mfrow = get_mfrow)
    
    # Generate output data.frame
    summ           <- data.frame(multief, fits.GoF, Estimates = e.values, p.values = p.values)
    names(summ)[2] <- criterion[1]
    
    # Output message of the selected model
    if(length(criterion) == 1){
      if(criterion %in% c("AIC", "BIC")){
        decision_func <- "which.min"
      } else {
        decision_func <- "which.max"
      }
    } else {
      decision_func <- paste("which", criterion[2], sep = ".")
    }
    out <- as.character(summ$multief[do.call(decision_func, list(fits.GoF))])
    sum <- summary(fits[[do.call(decision_func, list(fits.GoF))]])
    message("The model including the effect '", out, "' was the best model according to the specified criterion (", 
            criterion[1], ", ", decision_func, ")")

  } else {
    # If all models threw errors...
    warning("Fatal errors were found in all the models. Check them out by typing $Models to the generated object", call. = FALSE)
    plot <- NULL
    sum  <- NULL
    summ <- NULL
  }
  
  # Generate summary table of landscape attribute at each spatial scale
  lands_summary <- lands_table(multief, data)

  # Print summary of the 'best' model?
  if(print_sum && !is.null(sum)) print(sum)
  
  # Warn the user of possible warnings and messages during the workflow
  if(!is.null(effect_warns))
    message(paste("Warnings were found in the models with the following effects:", paste(effect_warns, collapse = ", ")))
  if(!is.null(effect_mess))
    message(paste("Messages were found in the models with the following effects:", paste(effect_mess, collapse = ", ")))

  # Return a list with relevant information of the analysis
  out_obj <- list(lands_summary = lands_summary, summary = summ, plot = plot, models = fits, warnings = mod_warn, messages = mod_mess)
  invisible(out_obj)
}
