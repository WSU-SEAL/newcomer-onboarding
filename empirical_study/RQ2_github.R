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

library(lmtest)

source_here("common.R")



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


result_vars = c( 'log_total_churn', 'is_gfi', 'num_newfile', 'is_bugfix', 'num_changed_files', 'change_entropy',
                 'doc_file_ratio',   'log_insertions','log_deletions', 
                 'title_length',  'desc_length',  'title_readability','desc_readability',
                 'promptness', 'num_patch' )

rows_nc <- result_vars

rows_non_nc <- result_vars

RQ1_latex_table_code =""

project="github"

r_sq_nc = "\textit{R-squared}"
r_sq_non_nc = "\textit{R-squared}"


datafile =paste("./", project, ".csv", sep="")
print(datafile)

dataset = read.csv(datafile, header = TRUE)
print(nrow(dataset))

dataset <- na.omit(dataset)

dataset$author_promptness <- as.numeric(as.character(dataset$author_promptness))
 


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
  
  
  av = anova(newcomer_glm_acceptance, type="II", test="LR")
  #summary(newcomer_glm_acceptance)
  
  model_summary =coef(summary(newcomer_glm_acceptance))
  
  print(exp(coef(newcomer_glm_acceptance)))
  
  for(i in 1:length(result_vars)){
    odds_ratio =exp(coef(newcomer_glm_acceptance)[result_vars[i]])
    #print(format(round(odds_ratio, 2), nsmall = 2))
    if(is.na(odds_ratio)||is.infinite(odds_ratio)){
      rows_nc[i] <- paste(rows_nc[i],'--',sep=" & ")
      rows_nc[i] <- paste(rows_nc[i],'--',sep=" & ")
    }
    else{
      pvalue =round(model_summary[grepl(paste(result_vars[i],'$',sep=""),row.names(model_summary)), 4],4)
      
      if(pvalue<0.05){
        if (odds_ratio<1) {
          rows_nc[i] <- paste( rows_nc[i], "& \\", "negative{" ,sep="")
        }
        else {
          rows_nc[i] <- paste(rows_nc[i], "& \\", "positive{" ,sep="")
        }
        rows_nc[i] <- paste(rows_nc[i],format(round(odds_ratio, 2), nsmall = 2),sep="  ") 
        rows_nc[i] <- paste(rows_nc[i],get_p_code(pvalue),sep="")
        rows_nc[i] <- paste(rows_nc[i], "}" ,sep=" ")
      }
      
      else rows_nc[i] <- paste(rows_nc[i],format(round(odds_ratio, 2), nsmall = 2),sep="& ")
      
      
      rows_nc[i] <- paste(rows_nc[i], format(round(100*av[result_vars[i], 2]/av['NULL', 4], 2), nsmall = 2) ,sep=" & ")
    }
    
    #print(pvalue)
  }
  
  
  
  veal=round(PseudoR2(newcomer_glm_acceptance,which = "VeallZimmermann"),3)
  print(lrtest(newcomer_glm_acceptance))
  #print(project)
  print(veal)
  r_sq_nc <- paste(r_sq_nc,veal,sep=" & ")
  
  model_coefficients <- coef(newcomer_glm_acceptance)
  #print(model_summary)
  #print(exp(model_coefficients))
  
}

print(rows_nc)
print(r_sq_nc) 


for(k in 1:3){
  
  non_newcomer_indv=otherDS[otherDS$project_group==k,]
  
  #summary(non_newcomer_indv)
  
  non_new_survived_vars = AutoSpearman(dataset = non_newcomer_indv, metrics = dynamic_vars,verbose = T)
  print(non_new_survived_vars)
  
  non_new_formula_string ="is_merged ~ is_bugfix "
  
  
  for (variable in non_new_survived_vars){
    if (variable !="is_bugfix" && !grepl(variable, non_new_formula_string)){
      non_new_formula_string =paste(non_new_formula_string, "+ ", variable,  sep=" ")
      #print(non_new_formula_string)
    }
    
  }
  
  #summary(non_newcomer_indv)
  
  
  non_newcomer_glm_acceptance <- glm(as.formula(non_new_formula_string) , data=non_newcomer_indv,  x=T, y=T, 
                                 family = binomial)
  
  null_model <- glm(as.formula("is_merged ~ 1 ") , data=non_newcomer_indv,  x=T, y=T, 
                    family = binomial)
  
  av = anova(non_newcomer_glm_acceptance, type="II", test="LR")
  
  #summary(non_newcomer_glm_acceptance)
  veal=round(PseudoR2(non_newcomer_glm_acceptance,which = "VeallZimmermann"),3)
  print(veal)
  
  odds_ratio =exp(coef(non_newcomer_glm_acceptance))
  print(odds_ratio)
  
  for(i in 1:length(result_vars)){
    odds_ratio =exp(coef(non_newcomer_glm_acceptance)[result_vars[i]])
    print(format(round(odds_ratio, 2), nsmall = 2))
    if(is.na(odds_ratio)||is.infinite(odds_ratio)){
      rows_non_nc[i] <- paste(rows_non_nc[i],'--',sep=" & ")
      rows_non_nc[i] <- paste(rows_non_nc[i],'--',sep=" & ")
    }
    else{
      pvalue =round(model_summary[grepl(paste(result_vars[i],'$',sep=""),row.names(model_summary)), 4],4)
      
      if(pvalue<0.05){
        if (odds_ratio<1) {
          rows_non_nc[i] <- paste( rows_non_nc[i], "& \\", "negative{" ,sep="")
        }
        else {
          rows_non_nc[i] <- paste(rows_non_nc[i], "& \\", "positive{" ,sep="")
        }
        rows_non_nc[i] <- paste(rows_non_nc[i],format(round(odds_ratio, 2), nsmall = 2),sep="  ") 
        rows_non_nc[i] <- paste(rows_non_nc[i],get_p_code(pvalue),sep="")
        rows_non_nc[i] <- paste(rows_non_nc[i], "}" ,sep=" ")
      }
      
      else rows_non_nc[i] <- paste(rows_non_nc[i],format(round(odds_ratio, 2), nsmall = 2),sep="& ")
      
      rows_non_nc[i] <- paste(rows_non_nc[i], format(round(100*av[result_vars[i], 2]/av['NULL', 4], 2), nsmall = 2) ,sep=" & ")
    }
  }
  
  #print(rows)
  
  veal=round(PseudoR2(non_newcomer_glm_acceptance,which = "VeallZimmermann"),3)
  print (lrtest(non_newcomer_glm_acceptance))
  print(veal)
  
  r_sq_non_nc <- paste(r_sq_non_nc,veal,sep=" & ")
  
  #model_coefficients <- coef(non_newcomer_glm_acceptance)
  #print(model_summary)
  #print(exp(model_coefficients))
  
}

print(rows_non_nc)
print(r_sq_non_nc) 
























