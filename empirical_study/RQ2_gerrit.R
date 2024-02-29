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


projects=c("android" ="Android", 
           "libreoffice" ="LibreOffice",
           "ovirt" ="oVirt", 
           "typo3" ="Typo3", 
           "wikimedia"="Wikimedia"
           )

# dynamic_vars = c('log_insertions','log_deletions','log_total_churn','number_patches',
#                  'log_author_promptness', 'is_bugfix',
#                  'file_count','doc_file_ratio','num_new_files', 'change_entropy',
#                  'title_length',
#                  'description_length','title_redability','description_readability')


dynamic_vars = c('log_insertions','log_deletions','log_total_churn','number_patches',
                 'promptness', 'is_bugfix',
                 'file_count','doc_file_ratio','num_new_files', 'change_entropy',
                 'title_length',
                 'description_length','title_redability','description_readability')

result_vars = c('log_insertions','log_deletions','log_total_churn','number_patches',
                'promptness', 'is_bugfix',
                'file_count','doc_file_ratio','num_new_files', 'change_entropy',
                'title_length',
                'description_length','title_redability','description_readability')

 # rows = c('IT','DT','TC','NP',
 #          'AD', 'IB',
 #          'FC','DF','NN', 'CE',
 #          'LT',
 #          'LD','TR','DR')

rows <- result_vars

#print(rows)

RQ1_latex_table_code =""

project="ovirt"

r_sq = "\textit{R-squared}"


for (project in names(projects)) {
  #p_name =projects[project]
  datafile =paste("gerrit_", project, ".csv", sep="")
  #print(datafile)
  dataset = read.csv(datafile, header = TRUE)
  print(nrow(dataset))
  
 
  dataset$doc_file_ratio=  dataset$doc_file_count/dataset$file_count
  dataset$log_author_promptness=log(dataset$author_promptness+1)
  dataset$log_insertions=log(dataset$insertions+1)
  dataset$log_deletions=log(dataset$deletions+1)
  dataset$log_total_churn=log(dataset$total_churn+1)
  dataset$promptness=10/(dataset$log_author_promptness+1)
  
  newcomerDS =dataset[dataset$is_newcomer==1,]
  otherDS =dataset[dataset$is_newcomer==0,]
  
  #print(dynamic_vars)
  
  #summary(newcomerDS)
  
  survived_vars = AutoSpearman(dataset = newcomerDS, metrics = dynamic_vars,verbose = T)
  #print(survived_vars)
  
  formula_string ="is_accepted ~ is_bugfix "
  
  
  
  for (variable in survived_vars){
      if (variable !="is_bugfix")
        formula_string =paste(formula_string, "+ ", variable,  sep=" ")
  
  }
    
    
  
  
  #print (formula_string)
  
  newcomer_glm_acceptance <- glm(as.formula(formula_string) , data=newcomerDS,  x=T, y=T, 
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
  #print(veal)
  r_sq <- paste(r_sq,veal,sep=" & ")
  
  model_coefficients <- coef(newcomer_glm_acceptance)
  #print(model_summary)
  #print(exp(model_coefficients))
  
  otherDS =dataset[dataset$is_newcomer==0,]
  
  non_new_survived_vars = AutoSpearman(dataset = otherDS, metrics = dynamic_vars,verbose = T)
  #print(non_new_survived_vars)
  
  non_new_formula_string ="is_accepted ~ is_bugfix + "
  
  
  
  for (variable in non_new_survived_vars){
    if (variable !="is_bugfix")
    non_new_formula_string =paste(formula_string, "+ ", variable,  sep=" ")
    
    
  }
  
  #print(non_new_formula_string)
  
  
  non_newcomer_glm_acceptance <- glm(as.formula(non_new_formula_string) , data=otherDS,  x=T, y=T, 
                                 family = binomial)
  
  
  #summary(non_newcomer_glm_acceptance)
  
  #model_summary =coef(summary(non_newcomer_glm_acceptance))
  for(i in 1:length(result_vars)){
    odds_ratio =exp(coef(non_newcomer_glm_acceptance)[result_vars[i]])
    #print(format(round(odds_ratio, 2), nsmall = 2))
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
  #print(veal)
  
  r_sq <- paste(r_sq,veal,sep=" & ")
  
  #model_coefficients <- coef(non_newcomer_glm_acceptance)
  #print(model_summary)
  #print(exp(model_coefficients))
  
}



print(rows)
print(r_sq) 


























#draw hierarchical cluster
project="android"

datafile =paste("gerrit_", project, ".csv", sep="")
print(datafile)
dataset = read.csv(datafile, header = TRUE)

dataset$doc_file_ratio=  dataset$doc_file_count/dataset$file_count
dataset$log_author_promptness=log(dataset$author_promptness+1)
dataset$log_insertions=log(dataset$insertions+1)
dataset$log_deletions=log(dataset$deletions+1)
dataset$log_total_churn=log(dataset$total_churn+1)

newcomerDS =dataset[dataset$is_newcomer==1,]


summary(newcomerDS)


newcomerDS <- rename(newcomerDS, IT=log_insertions, DT=log_deletions,TC=log_total_churn, NP=number_patches,
                     AD=log_author_promptness, IB=is_bugfix,
                     FC=file_count, DF=doc_file_ratio, NN=num_new_files,
                     CE=change_entropy,  LT=title_length,
                     LD=description_length, TR=title_redability, DR=description_readability)

summary(newcomerDS)

dynamic_vars = c('IT','DT','TC','NP',
                 'AD','IB',
                 'FC','DF',
                 'NN',
                 'CE','LT',
                 'LD','TR','DR')

setEPS()                                             # Set postscript arguments
postscript("cluster.eps") 

vc <- varclus(~ ., data=newcomerDS[,dynamic_vars], trans="abs")
#Plot hierarchical clusters and the spearman's correlation threshold of 0.7
plot(vc)
threshold <- 0.7
abline(h=1-threshold, col = "red", lty = 2)
dev.off()   






