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

source("common.R")

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



dynamic_vars = c('log_insertions','log_deletions','log_total_churn','number_patches',
                 'promptness', 'is_bugfix',
                 'file_count','doc_file_ratio','num_new_files', 'change_entropy',
                 'title_length',
                 'description_length','title_redability','description_readability')

result_vars = c( 'log_total_churn', 'num_new_files', 'is_bugfix', 'file_count', 'change_entropy',
                 'doc_file_ratio',   'log_insertions','log_deletions', 
                 'title_length',  'description_length',  'title_redability','description_readability',
                 'promptness', 'number_patches' )


rows_nc <- result_vars

rows_non_nc <- result_vars


#print(rows)

RQ1_latex_table_code =""

project="ovirt"

r_sq_nc = "\textit{R-squared}"

r_sq_non_nc = "\textit{R-squared}"

#Newcomer analysis table
for (project in names(projects)) {
  #p_name =projects[project]
  datafile =paste("gerrit_", project, ".csv", sep="")
  
  print(project)
  
  dataset = read.csv(datafile, header = TRUE)
  print(nrow(dataset))
  
  
  dataset$doc_file_ratio=  dataset$doc_file_count/dataset$file_count
  dataset$log_author_promptness=log(dataset$author_promptness+1)
  dataset$log_insertions=log(dataset$insertions+1)
  dataset$log_deletions=log(dataset$deletions+1)
  dataset$log_total_churn=log(dataset$total_churn+1)
  dataset$promptness=10/(dataset$log_author_promptness+1)
  
  newcomerDS =dataset[dataset$is_newcomer==1,]
  #otherDS =dataset[dataset$is_newcomer==0,]
  
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
  
  null_model <- glm(as.formula("is_accepted ~ 1 ") , data=newcomerDS,  x=T, y=T, 
                                 family = binomial)
  av = anova(newcomer_glm_acceptance, type="II", test="LR")
  #summary(newcomer_glm_acceptance)
  
  model_summary =coef(summary(newcomer_glm_acceptance))
  
  for(i in 1:length(result_vars)){
    odds_ratio =exp(coef(newcomer_glm_acceptance)[result_vars[i]])
    #print(format(round(odds_ratio, 2), nsmall = 2))
    if(is.na(odds_ratio)||is.infinite(odds_ratio)){
      rows_nc[i] <- paste(rows_nc[i],'--',sep=" & ")
      rows_nc[i] <- paste(rows_nc[i],'--',sep=" & ")
    }
    else{
      pvalue =round(model_summary[grepl(paste(result_vars[i],'$',sep=""),row.names(model_summary)), 4],4)
      
      if(length(pvalue) ==0)
        pvalue =0.5
      
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
  
  print (lrtest(newcomer_glm_acceptance, null_model)) #Log likelihood test
  #print(project)
  #print(veal)
  r_sq_nc <- paste(r_sq_nc,veal,sep=" & ")
  
  print("Newcomer odds ratios:")
  model_coefficients <- coef(newcomer_glm_acceptance)
  #print(model_summary)
  print(exp(model_coefficients))
  

  
  
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
  
  null_model <- glm(as.formula("is_accepted ~ 1 ") , data=otherDS,  x=T, y=T, 
                    family = binomial)
  
  print('Anova print')
  av = anova(non_newcomer_glm_acceptance, type="II", test="LR")
  #print(av[`Resid. Dev`])
  
  #print(100*av['is_bugfix', 2]/av['NULL', 4])
  
  #summary(non_newcomer_glm_acceptance)
  
  #model_summary =coef(summary(non_newcomer_glm_acceptance))
  for(i in 1:length(result_vars)){
    odds_ratio =exp(coef(non_newcomer_glm_acceptance)[result_vars[i]])
    #print(format(round(odds_ratio, 2), nsmall = 2))
    
    if(is.na(odds_ratio)||is.infinite(odds_ratio)){
      rows_non_nc[i] <- paste(rows_non_nc[i],'--',sep=" & ")
      rows_non_nc[i] <- paste(rows_non_nc[i],'--',sep=" & ")
    }
    else{
      pvalue =round(model_summary[grepl(paste(result_vars[i],'$',sep=""),row.names(model_summary)), 4],4)
      
       if(length(pvalue) ==0)
         pvalue =0.5
      
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
  print (lrtest(non_newcomer_glm_acceptance, null_model))
  #print(veal)
  
  r_sq_non_nc <- paste(r_sq_non_nc,veal,sep=" & ")
  
  #model_coefficients <- coef(non_newcomer_glm_acceptance)
  #print(model_summary)
  print("non-newcomer odds ratios")
  print(exp(model_coefficients))
  #break
  
}


print("newcomer")
print(r_sq_nc)
print(rows_nc)


print("non-newcomer")
print(r_sq_non_nc)
print(rows_non_nc)


























