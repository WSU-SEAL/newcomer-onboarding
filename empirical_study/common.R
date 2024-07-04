get.automated.spearman <- function(dataset, metrics, spearman.threshold, verbose = F){
  
  .get.higher.correlation <- function(index1, index2, metrics, metric.correlations, verbose = T, count){
    metric1 <- metrics[index1]
    metric2 <- metrics[index2]
    if(verbose)
      cat(paste0('Step ', count, ' - {', metric1, ', ', metric2, '} > '))# are correlated with r = ', metric.correlations[metric1, metric2], '\n'))
    metric1.correlation <- mean(metric.correlations[metric1, !metrics %in% c(metric1, metric2)])
    metric2.correlation <- mean(metric.correlations[metric2, !metrics %in% c(metric1, metric2)])
    if(metric1.correlation <= metric2.correlation){
      return(index2)
    } else {
      return(index1)
    }
  }
  
  metric.correlations <- abs(rcorr(as.matrix(dataset[, metrics]), type = 'spearman')$r)
  
  above.threshold <- which((metric.correlations >= spearman.threshold), arr.ind = TRUE)
  row.names(above.threshold) <- NULL
  above.threshold <- as.data.frame(above.threshold, row.names = NULL)
  above.threshold <- above.threshold[above.threshold$row != above.threshold$col, ]
  above.threshold$correlation <- 100
  for(i in 1:nrow(above.threshold)){
    above.threshold$correlation[i] <- metric.correlations[above.threshold$row[i], above.threshold$col[i]]
  }
  above.threshold <- above.threshold[order(-above.threshold$correlation), ]
  
  exclude.metrics <- {}
  count <- 1
  repeat{
    if(nrow(above.threshold) == 0)
      break
    tmp <- above.threshold[1, ]
    exclude.index <- .get.higher.correlation(tmp$row, tmp$col, metrics, metric.correlations, verbose, count)
    exclude.metrics <- c(exclude.metrics, metrics[exclude.index])
    if(verbose){
      cat(paste0(metrics[exclude.index], '\n'))
      count <- count + 1
    }
    above.threshold <- above.threshold[-which((above.threshold$row == exclude.index) | (above.threshold$col == exclude.index)), ]
  }
  selected.metrics <- metrics[!metrics %in% exclude.metrics]
  return(selected.metrics)
}


remove.constant.categorical <-
  function(dataset,
           metrics) {
    # Check constant metrics
    constant <-
      apply(dataset[, metrics], 2, function(x)
        max(x) == min(x))
    constant <- names(constant[constant == TRUE])
    # Remove constant metrics
    if (length(constant) > 0) {
      print("Constant")
      print(constant)
      metrics <- metrics[!metrics %in% constant]
    }
    
    # Check categorical metrics
    category <- sapply(dataset[, metrics], class)
    category <- names(category[category == "character"])
    # Remove categorical metrics from Spearman Analysis
    if (length(category) > 0) {
      print("Category:")
      print(category)
      metrics <- metrics[!metrics %in% category]
    }
    
    return(metrics)
  }


stepwise.vif <-
  function (dataset,
            metrics,
            vif.threshold = 5,
            verbose = F)
  {
    dataset$dummy <- rnorm(nrow(dataset))
    output <- metrics
    step.count <- 1
    output.results <- list()
    repeat {
      vif.scores <- vif(lm(as.formula(paste0(
        "dummy~", paste0(output,
                         collapse = "+")
      )), data = dataset))
      na.coefficients <- Reduce('|', is.nan(vif.scores))
      if (na.coefficients) {
        stop("NA coefficient in a regression model.")
      }
      output.results[[step.count]] <-
        sort(vif.scores, decreasing = F)
      vif.scores <- vif.scores[vif.scores >= vif.threshold]
      if (length(vif.scores) == 0)
        break
      drop.var <-
        names(vif.scores[vif.scores == max(vif.scores)])[1]
      if (verbose) {
        print(paste0(
          "Step ",
          step.count,
          " - Exclude ",
          drop.var,
          " (VIF = ",
          max(vif.scores),
          ")"
        ))
      }
      step.count <- step.count + 1
      output <- output[!output %in% drop.var]
    }
    names(output.results) <- paste0("Iteration ", 1:step.count)
    names(output.results)[length(output.results)] <- "Final"
    return(output)
  }



AutoSpearman <-
  function(dataset,
           metrics,
           spearman.threshold = 0.7,
           vif.threshold = 5,
           verbose = T) {
    # Remove constant metrics and categorical metrics
    metrics <- remove.constant.categorical(dataset, metrics)
    print(metrics)
    
    
    spearman.metrics <-
      get.automated.spearman(dataset, metrics, spearman.threshold, verbose)
    AutoSpearman.metrics <-
      stepwise.vif(dataset, spearman.metrics, vif.threshold, verbose)
    
    return(AutoSpearman.metrics)
  }

get_p_code<- function(pvalue){
  if(length(pvalue) == 0)
    return ("$-$")
  
  p_code =""
  if (pvalue <0.001){
    p_code ="$^{***}$"
  }
  else if (pvalue <0.01){ 
    p_code ="$^{**}$"
  }
  else if (pvalue <0.05){
    p_code ="$^{*}$"
  }
  return (p_code)
}

get_OR_code<- function(or_value){
  if(length(or_value) == 0)
    return ("$-$")
  
  or_code =""
  if (or_value >1){
    or_code =paste("\\biasToWomen{",or_value,"}", sep="")
  }
  else or_code =paste("\\biasToMen{",or_value,"}", sep="")
  
  return (or_code)
}


get_impact_code<- function(impact_value){
  if(length(impact_value) == 0)
    return ("$-$")
  
  impact_code =""
  if (impact_value <1){
    impact_code =paste("\\biasToWomen{",impact_value,"}", sep="")
  }
  else impact_code =paste("\\biasToMen{",impact_value,"}", sep="")
  
  return (impact_code)
}

get_z_code<- function(z_value){
  if(length(z_value) == 0)
    return ("$-$")
  
  z_code =""
  if (z_value <0){
    z_code =paste("\\biasToWomen{",z_value,"}", sep="")
  }
  else z_code =paste("\\biasToMen{",z_value,"}", sep="")
  
  return (z_code)
}

get_z_code_review<- function(z_value){
  if(length(z_value) == 0)
    return ("$-$")
  
  z_code =""
  if (z_value >0){
    z_code =paste("\\biasToWomen{",z_value,"}", sep="")
  }
  else z_code =paste("\\biasToMen{",z_value,"}", sep="")
  
  return (z_code)
}

get_GN_code<- function(or_value){
  if((length(or_value) == 0) ||is.na(or_value))
    return ("$-$")
  
  
  or_code =""
  if (or_value >1){
    or_code =paste("\\biasToGendered{",or_value,"}", sep="")
  }
  else or_code =or_value
  
  return (or_code)
}
