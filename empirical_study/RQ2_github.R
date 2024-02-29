library(lme4)
library(car)
library(Hmisc)
library(fmsb)
library(rms)
library(pROC)
library(epiDisplay)
library(MASS)
library(pscl)
library(DescTools)
library(effects)
library(tidyverse)


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

library(rstudioapi)
cur_dir = dirname(getSourceEditorContext()$path)
setwd(cur_dir)
getwd()

# dynamic_vars = c('log_insertions','log_deletions','log_total_churn','num_patch',
#                  'log_author_promptness',
#                  'is_bugfix', 'is_gfi','change_entropy',
#                  'num_changed_files','doc_file_ratio','num_newfile',
#                  'title_length',
#                  'desc_length','title_readability','desc_readability')


dynamic_vars = c('log_insertions','log_deletions','log_total_churn','num_patch',
                 'promptness',
                 'is_bugfix', 'is_gfi','change_entropy',
                 'num_changed_files','doc_file_ratio','num_newfile',
                 'title_length',
                 'desc_length','title_readability','desc_readability')

result_vars = c('log_insertions','log_deletions','log_total_churn','num_patch',
                'promptness',
                'is_bugfix', 'is_gfi','change_entropy',
                'num_changed_files','doc_file_ratio','num_newfile',
                'title_length',
                'desc_length','title_readability','desc_readability')

rows <- result_vars


RQ1_latex_table_code =""

project="github"

r_sq = "\textit{R-squared}"


datafile =paste("./", project, ".csv", sep="")
print(datafile)

dataset = read.csv(datafile, header = TRUE)
print(nrow(dataset))

dataset$author_promptness <- as.numeric(as.character(dataset$author_promptness))
 
dataset <- na.omit(dataset)

dataset$author_promptness[dataset$author_promptness<0] <- 0



dataset$doc_file_ratio=  dataset$doc_file_count/(dataset$num_changed_files+1)
dataset$log_author_promptness=log(dataset$author_promptness+1)
dataset$log_insertions=log(dataset$additions+1)
dataset$log_deletions=log(dataset$deletions+1)
dataset$log_total_churn=log(dataset$total_churn+1)
dataset$promptness=10/(dataset$log_author_promptness+1)

newcomerDS =dataset[dataset$is_newcomer==1,]
otherDS =dataset[dataset$is_newcomer==0,]

  
for(k in 1:3){
  newcomer_indv =newcomerDS[newcomerDS$project_group==k,]
  
  #summary(newcomer_indv)
  
  survived_vars = AutoSpearman(dataset = newcomer_indv, metrics = dynamic_vars,verbose = T)
  print(survived_vars)
  
  formula_string ="is_merged ~ is_bugfix "
  
  
  
  for (variable in survived_vars){
      if (variable !="is_bugfix")
        formula_string =paste(formula_string, "+ ", variable,  sep=" ")
  
  }
    
    
  
  
  #print (formula_string)
  
  newcomer_glm_acceptance <- glm(as.formula(formula_string) , data=newcomer_indv,  x=T, y=T, 
                            family = binomial)
  
  #summary(newcomer_glm_acceptance)
  
  model_summary =coef(summary(newcomer_glm_acceptance))
  
  for(i in 1:length(result_vars)){
    odds_ratio =exp(coef(newcomer_glm_acceptance)[result_vars[i]])
    #print(format(round(odds_ratio, 2), nsmall = 2))
    if(is.na(odds_ratio)||is.infinite(odds_ratio)){
      rows[i] <- paste(rows[i],'--',sep=" & ")
    }
    else{
      pvalue =round(model_summary[grepl(paste(result_vars[i],'$',sep=""),row.names(model_summary)), 4],4)
      rows[i] <- paste(rows[i],format(round(odds_ratio, 2), nsmall = 2),sep=" & ")
      rows[i] <- paste(rows[i],get_p_code(pvalue),sep="")
    }
    
    #print(pvalue)
  }
  
  
  
  veal=round(PseudoR2(newcomer_glm_acceptance,which = "VeallZimmermann"),3)
  #print(project)
  print(veal)
  r_sq <- paste(r_sq,veal,sep=" & ")
  
  model_coefficients <- coef(newcomer_glm_acceptance)
  #print(model_summary)
  #print(exp(model_coefficients))
  
  non_newcomer_indv=otherDS[otherDS$project_group==k,]
  
  summary(non_newcomer_indv)
  
  non_new_survived_vars = AutoSpearman(dataset = non_newcomer_indv, metrics = dynamic_vars,verbose = T)
  print(non_new_survived_vars)
  
  non_new_formula_string ="is_merged ~ is_bugfix "
  
  
  for (variable in non_new_survived_vars){
    if (variable !="is_bugfix" && !grepl(variable, non_new_formula_string)){
      non_new_formula_string =paste(non_new_formula_string, "+ ", variable,  sep=" ")
      #print(non_new_formula_string)
    }
    
  }
  
  summary(non_newcomer_indv)
  
  
  non_newcomer_glm_acceptance <- glm(as.formula(non_new_formula_string) , data=non_newcomer_indv,  x=T, y=T, 
                                 family = binomial)
  
  summary(non_newcomer_glm_acceptance)
  veal=round(PseudoR2(non_newcomer_glm_acceptance,which = "VeallZimmermann"),3)
  print(veal)
  
  odds_ratio =exp(coef(non_newcomer_glm_acceptance))
  print(odds_ratio)
  
  for(i in 1:length(result_vars)){
    odds_ratio =exp(coef(non_newcomer_glm_acceptance)[result_vars[i]])
    print(format(round(odds_ratio, 2), nsmall = 2))
    if(is.na(odds_ratio)||is.infinite(odds_ratio)){
      rows[i] <- paste(rows[i],'--',sep=" & ")
    }
    else{
      pvalue =round(model_summary[grepl(paste(result_vars[i],'$',sep=""),row.names(model_summary)), 4],4)
      rows[i] <- paste(rows[i],format(round(odds_ratio, 2), nsmall = 2),sep=" & ")
      rows[i] <- paste(rows[i],get_p_code(pvalue),sep="")
    }
  }
  
  #print(rows)
  
  veal=round(PseudoR2(non_newcomer_glm_acceptance,which = "VeallZimmermann"),3)
  print(veal)
  
  r_sq <- paste(r_sq,veal,sep=" & ")
  
  #model_coefficients <- coef(non_newcomer_glm_acceptance)
  #print(model_summary)
  #print(exp(model_coefficients))
  
}

print(rows)
print(r_sq) 




















#for GitHub large non-newcomers




datafile = "./nonnewcomer_large.csv"
print(datafile)

dataset = read.csv(datafile, header = TRUE)
summary(dataset)


dynamic_vars = c('log_insertions','log_deletions','log_total_churn','num_patch',
                 'log_author_promptness',
                 'is_bugfix', 'is_gfi','change_entropy',
                 'num_changed_files','doc_file_ratio','num_newfile',
                 'title_length',
                 'desc_length','title_readability','desc_readability')

non_new_survived_vars = AutoSpearman(dataset = dataset, metrics = dynamic_vars,verbose = T)
print(non_new_survived_vars)

non_new_formula_string ="is_merged ~ is_bugfix "


for (variable in non_new_survived_vars){
  if (variable !="is_bugfix" && !grepl(variable, non_new_formula_string)){
    non_new_formula_string =paste(non_new_formula_string, "+ ", variable,  sep=" ")
    #print(non_new_formula_string)
  }
  
}

#summary(non_newcomer_indv)

#non_new_formula_string <- "is_merged ~ is_bugfix + log_insertions +  log_deletions +  num_patch  +  log_author_promptness +  is_gfi +  change_entropy +  doc_file_ratio" #+ " +   num_newfile# title_readability +  desc_readability ++ title_length +  desc_length  +  "

non_newcomer_glm_acceptance <- glm(as.formula(non_new_formula_string) , data=dataset,  x=T, y=T, 
                                   family = binomial)

summary(non_newcomer_glm_acceptance)

veal=round(PseudoR2(non_newcomer_glm_acceptance,which = "VeallZimmermann"),3)
print(veal)

odds_ratio =exp(coef(non_newcomer_glm_acceptance))
print(odds_ratio)





