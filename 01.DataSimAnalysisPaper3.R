
#============================================================
# Load packages
#============================================================

library(tidyverse)
library(dplyr)
library(TAM)
library(coda)
library(CDM)
library(rstan)
library(R.utils)
library(data.table)
library("parallel")
library("doParallel")

#=============================================================
# TRT Model script 
#=============================================================
 
# model 

twoTRT_code <- rstan::stan_model(model_code ="
  data{
    int <lower=0> n_student; //  number of individuals
    int <lower=0> n_item; //  number of items
    int <lower=0> n_testlet; //  number of testlets
    int <lower=0> ind_testlet[n_item]; //  testlet membership
    int<lower=0,upper=1> Y[n_item,n_student]; //array of responses
  }

   parameters {
    real<lower=0> alpha [n_item]; //item discrimination
    real beta [n_item]; //item difficulty parameter
    real mu_beta; // mean difficulty
    real<lower=0> sigma_beta; //difficulty sd
    vector[n_student] theta; // latent trait
    vector[n_student] te[n_testlet]; // testlet effect
    real<lower=0> sigma_te[n_testlet]; // sd of testlet effect
   }


   model{
    theta ~ normal(0,1);
    alpha ~ lognormal(0,1);
    beta ~ normal(mu_beta,sigma_beta);
    mu_beta ~ normal(0,5);
    sigma_beta ~ cauchy(0,5);
    
    for (i in 1:n_testlet){
      te[i] ~ normal(0, sigma_te[i]);
      sigma_te[i] ~ cauchy(0,5);
    }
    
     for (j in 1:n_item){
     Y[j] ~ bernoulli_logit(alpha[j]*(theta-beta[j]-te[ind_testlet[j]]));
  }
  }
  ")

twoTRT_code_DIFlow <- rstan::stan_model(model_code ="
  data{
    int <lower=0> n_student; //  number of individuals
    int <lower=0> n_item; //  number of items
    int <lower=0> n_testlet; //  number of testlets
    int <lower=0> ind_testlet[n_item]; //  testlet membership
    int<lower=0,upper=1> Y[n_item,n_student]; //array of responses
  }

   parameters {
    real<lower=0> alpha [n_item]; //item discrimination
    real beta [n_item]; //item difficulty parameter
    real mu_beta; // mean difficulty
    real<lower=0> sigma_beta; //difficulty sd
    vector[n_student] theta; // latent trait
    vector[n_student] te[n_testlet]; // testlet effect
    real<lower=0> sigma_te[n_testlet]; // sd of testlet effect
   }


   model{
    theta ~ normal(-1.8,1);
    alpha ~ lognormal(0,1);
    beta ~ normal(mu_beta,sigma_beta);
    mu_beta ~ normal(0,5);
    sigma_beta ~ cauchy(0,5);
    
    for (i in 1:n_testlet){
      te[i] ~ normal(0, sigma_te[i]);
      sigma_te[i] ~ cauchy(0,5);
    }
    
     for (j in 1:n_item){
     Y[j] ~ bernoulli_logit(alpha[j]*(theta-beta[j]-te[ind_testlet[j]]));
  }
  }
  ")
  

#=============================================================
## Simulation Conditions
#=============================================================
 
DIF_conditions = list(
  # Baseline
  cc01 = list(mma = 0, mmb = 0,   pp = 0, dep = 0,   type = "No DIF",  DIFcntry = "No DIF"),
  # Uniform DIF 
  cc02 = list(mma = 0, mmb = .31, pp = 1, dep = .10, type = "uniform", DIFcntry = "Low"),
  cc03 = list(mma = 0, mmb = .45, pp = 1, dep = .10, type = "uniform", DIFcntry = "Low"),
  cc04 = list(mma = 0, mmb =1.40, pp = 1, dep = .10, type = "uniform", DIFcntry = "Low"),
  cc05 = list(mma = 0, mmb = .31, pp = 2, dep = .10, type = "uniform", DIFcntry = "Low"),
  cc06 = list(mma = 0, mmb = .45, pp = 2, dep = .10, type = "uniform", DIFcntry = "Low"),
  cc07 = list(mma = 0, mmb =1.40, pp = 2, dep = .10, type = "uniform", DIFcntry = "Low"),
  cc08 = list(mma = 0, mmb = .31, pp = 3, dep = .10, type = "uniform", DIFcntry = "Low"),
  cc09 = list(mma = 0, mmb = .45, pp = 3, dep = .10, type = "uniform", DIFcntry = "Low"),
  cc10 = list(mma = 0, mmb =1.40, pp = 3, dep = .10, type = "uniform", DIFcntry = "Low"),
  
  cc11 = list(mma = 0, mmb = .31, pp = 1, dep = .30, type = "uniform", DIFcntry = "Low"),
  cc12 = list(mma = 0, mmb = .45, pp = 1, dep = .30, type = "uniform", DIFcntry = "Low"),
  cc13 = list(mma = 0, mmb =1.40, pp = 1, dep = .30, type = "uniform", DIFcntry = "Low"),
  cc14 = list(mma = 0, mmb = .31, pp = 2, dep = .30, type = "uniform", DIFcntry = "Low"),
  cc15 = list(mma = 0, mmb = .45, pp = 2, dep = .30, type = "uniform", DIFcntry = "Low"),
  cc16 = list(mma = 0, mmb =1.40, pp = 2, dep = .30, type = "uniform", DIFcntry = "Low"),
  cc17 = list(mma = 0, mmb = .31, pp = 3, dep = .30, type = "uniform", DIFcntry = "Low"),
  cc18 = list(mma = 0, mmb = .45, pp = 3, dep = .30, type = "uniform", DIFcntry = "Low"),
  cc19 = list(mma = 0, mmb =1.40, pp = 3, dep = .30, type = "uniform", DIFcntry = "Low"),
  
  cc20 = list(mma = 0, mmb = .31, pp = 1, dep = .50, type = "uniform", DIFcntry = "Low"),
  cc21 = list(mma = 0, mmb = .45, pp = 1, dep = .50, type = "uniform", DIFcntry = "Low"),
  cc22 = list(mma = 0, mmb =1.40, pp = 1, dep = .50, type = "uniform", DIFcntry = "Low"),
  cc23 = list(mma = 0, mmb = .31, pp = 2, dep = .50, type = "uniform", DIFcntry = "Low"),
  cc24 = list(mma = 0, mmb = .45, pp = 2, dep = .50, type = "uniform", DIFcntry = "Low"),
  cc25 = list(mma = 0, mmb =1.40, pp = 2, dep = .50, type = "uniform", DIFcntry = "Low"),
  cc26 = list(mma = 0, mmb = .31, pp = 3, dep = .50, type = "uniform", DIFcntry = "Low"),
  cc27 = list(mma = 0, mmb = .45, pp = 3, dep = .50, type = "uniform", DIFcntry = "Low"),
  cc28 = list(mma = 0, mmb =1.40, pp = 3, dep = .50, type = "uniform", DIFcntry = "Low"),
  
  cc29 = list(mma = 0, mmb = .31, pp = 1, dep = .10, type = "uniform", DIFcntry = "Med"),
  cc30 = list(mma = 0, mmb = .45, pp = 1, dep = .10, type = "uniform", DIFcntry = "Med"),
  cc31 = list(mma = 0, mmb =1.40, pp = 1, dep = .10, type = "uniform", DIFcntry = "Med"),
  cc32 = list(mma = 0, mmb = .31, pp = 2, dep = .10, type = "uniform", DIFcntry = "Med"),
  cc33 = list(mma = 0, mmb = .45, pp = 2, dep = .10, type = "uniform", DIFcntry = "Med"),
  cc34 = list(mma = 0, mmb =1.40, pp = 2, dep = .10, type = "uniform", DIFcntry = "Med"),
  cc35 = list(mma = 0, mmb = .31, pp = 3, dep = .10, type = "uniform", DIFcntry = "Med"),
  cc36 = list(mma = 0, mmb = .45, pp = 3, dep = .10, type = "uniform", DIFcntry = "Med"),
  cc37 = list(mma = 0, mmb =1.40, pp = 3, dep = .10, type = "uniform", DIFcntry = "Med"),

  cc38 = list(mma = 0, mmb = .31, pp = 1, dep = .30, type = "uniform", DIFcntry = "Med"),
  cc39 = list(mma = 0, mmb = .45, pp = 1, dep = .30, type = "uniform", DIFcntry = "Med"),
  cc40 = list(mma = 0, mmb =1.40, pp = 1, dep = .30, type = "uniform", DIFcntry = "Med"),
  cc41 = list(mma = 0, mmb = .31, pp = 2, dep = .30, type = "uniform", DIFcntry = "Med"),
  cc42 = list(mma = 0, mmb = .45, pp = 2, dep = .30, type = "uniform", DIFcntry = "Med"),
  cc43 = list(mma = 0, mmb =1.40, pp = 2, dep = .30, type = "uniform", DIFcntry = "Med"),
  cc44 = list(mma = 0, mmb = .31, pp = 3, dep = .30, type = "uniform", DIFcntry = "Med"),
  cc45 = list(mma = 0, mmb = .45, pp = 3, dep = .30, type = "uniform", DIFcntry = "Med"),
  cc46 = list(mma = 0, mmb =1.40, pp = 3, dep = .30, type = "uniform", DIFcntry = "Med"),

  cc47 = list(mma = 0, mmb = .31, pp = 1, dep = .50, type = "uniform", DIFcntry = "Med"),
  cc48 = list(mma = 0, mmb = .45, pp = 1, dep = .50, type = "uniform", DIFcntry = "Med"),
  cc49 = list(mma = 0, mmb =1.40, pp = 1, dep = .50, type = "uniform", DIFcntry = "Med"),
  cc50 = list(mma = 0, mmb = .31, pp = 2, dep = .50, type = "uniform", DIFcntry = "Med"),
  cc51 = list(mma = 0, mmb = .45, pp = 2, dep = .50, type = "uniform", DIFcntry = "Med"),
  cc52 = list(mma = 0, mmb =1.40, pp = 2, dep = .50, type = "uniform", DIFcntry = "Med"),
  cc53 = list(mma = 0, mmb = .31, pp = 3, dep = .50, type = "uniform", DIFcntry = "Med"),
  cc54 = list(mma = 0, mmb = .45, pp = 3, dep = .50, type = "uniform", DIFcntry = "Med"),
  cc55 = list(mma = 0, mmb =1.40, pp = 3, dep = .50, type = "uniform", DIFcntry = "Med"),
  
  # cc56 = list(mma = 0, mmb = .31, pp = 1, dep = .10, type = "uniform", DIFcntry = "High"),
  # cc57 = list(mma = 0, mmb = .45, pp = 1, dep = .10, type = "uniform", DIFcntry = "High"),
  # cc58 = list(mma = 0, mmb =1.40, pp = 1, dep = .10, type = "uniform", DIFcntry = "High"),
  # cc59 = list(mma = 0, mmb = .31, pp = 2, dep = .10, type = "uniform", DIFcntry = "High"),
  # cc60 = list(mma = 0, mmb = .45, pp = 2, dep = .10, type = "uniform", DIFcntry = "High"),
  # cc61 = list(mma = 0, mmb =1.40, pp = 2, dep = .10, type = "uniform", DIFcntry = "High"),
  # cc62 = list(mma = 0, mmb = .31, pp = 3, dep = .10, type = "uniform", DIFcntry = "High"),
  # cc63 = list(mma = 0, mmb = .45, pp = 3, dep = .10, type = "uniform", DIFcntry = "High"),
  # cc64 = list(mma = 0, mmb =1.40, pp = 3, dep = .10, type = "uniform", DIFcntry = "High"),
  # 
  # cc65 = list(mma = 0, mmb = .31, pp = 1, dep = .30, type = "uniform", DIFcntry = "High"),
  # cc66 = list(mma = 0, mmb = .45, pp = 1, dep = .30, type = "uniform", DIFcntry = "High"),
  # cc67 = list(mma = 0, mmb =1.40, pp = 1, dep = .30, type = "uniform", DIFcntry = "High"),
  # cc68 = list(mma = 0, mmb = .31, pp = 2, dep = .30, type = "uniform", DIFcntry = "High"),
  # cc69 = list(mma = 0, mmb = .45, pp = 2, dep = .30, type = "uniform", DIFcntry = "High"),
  # cc70 = list(mma = 0, mmb =1.40, pp = 2, dep = .30, type = "uniform", DIFcntry = "High"),
  # cc71 = list(mma = 0, mmb = .31, pp = 3, dep = .30, type = "uniform", DIFcntry = "High"),
  # cc72 = list(mma = 0, mmb = .45, pp = 3, dep = .30, type = "uniform", DIFcntry = "High"),
  # cc73 = list(mma = 0, mmb =1.40, pp = 3, dep = .30, type = "uniform", DIFcntry = "High"),
  # 
  # cc74 = list(mma = 0, mmb = .31, pp = 1, dep = .50, type = "uniform", DIFcntry = "High"),
  # cc75 = list(mma = 0, mmb = .45, pp = 1, dep = .50, type = "uniform", DIFcntry = "High"),
  # cc76 = list(mma = 0, mmb =1.40, pp = 1, dep = .50, type = "uniform", DIFcntry = "High"),
  # cc77 = list(mma = 0, mmb = .31, pp = 2, dep = .50, type = "uniform", DIFcntry = "High"),
  # cc78 = list(mma = 0, mmb = .45, pp = 2, dep = .50, type = "uniform", DIFcntry = "High"),
  # cc79 = list(mma = 0, mmb =1.40, pp = 2, dep = .50, type = "uniform", DIFcntry = "High"),
  # cc80 = list(mma = 0, mmb = .31, pp = 3, dep = .50, type = "uniform", DIFcntry = "High"),
  # cc81 = list(mma = 0, mmb = .45, pp = 3, dep = .50, type = "uniform", DIFcntry = "High"),
  # cc82 = list(mma = 0, mmb =1.40, pp = 3, dep = .50, type = "uniform", DIFcntry = "High"),
  # Non-uniform DIF
  cc56 = list(mma = .10, mmb = 0, pp = 1, dep = .10, type = "non-uniform", DIFcntry = "Low"),
  cc57 = list(mma = .39, mmb = 0, pp = 1, dep = .10, type = "non-uniform", DIFcntry = "Low"),
  cc58 = list(mma = .91, mmb = 0, pp = 1, dep = .10, type = "non-uniform", DIFcntry = "Low"),
  cc59 = list(mma = .10, mmb = 0, pp = 2, dep = .10, type = "non-uniform", DIFcntry = "Low"),
  cc60 = list(mma = .39, mmb = 0, pp = 2, dep = .10, type = "non-uniform", DIFcntry = "Low"),
  cc61 = list(mma = .91, mmb = 0, pp = 2, dep = .10, type = "non-uniform", DIFcntry = "Low"),
  cc62 = list(mma = .10, mmb = 0, pp = 3, dep = .10, type = "non-uniform", DIFcntry = "Low"),
  cc63 = list(mma = .39, mmb = 0, pp = 3, dep = .10, type = "non-uniform", DIFcntry = "Low"),
  cc64 = list(mma = .91, mmb = 0, pp = 3, dep = .10, type = "non-uniform", DIFcntry = "Low"),

  cc65 = list(mma = .10, mmb = 0, pp = 1, dep = .30, type = "non-uniform", DIFcntry = "Low"),
  cc66 = list(mma = .39, mmb = 0, pp = 1, dep = .30, type = "non-uniform", DIFcntry = "Low"),
  cc67 = list(mma = .91, mmb = 0, pp = 1, dep = .30, type = "non-uniform", DIFcntry = "Low"),
  cc68 = list(mma = .10, mmb = 0, pp = 2, dep = .30, type = "non-uniform", DIFcntry = "Low"),
  cc69 = list(mma = .39, mmb = 0, pp = 2, dep = .30, type = "non-uniform", DIFcntry = "Low"),
  cc70 = list(mma = .91, mmb = 0, pp = 2, dep = .30, type = "non-uniform", DIFcntry = "Low"),
  cc71 = list(mma = .10, mmb = 0, pp = 3, dep = .30, type = "non-uniform", DIFcntry = "Low"),
  cc72 = list(mma = .39, mmb = 0, pp = 3, dep = .30, type = "non-uniform", DIFcntry = "Low"),
  cc73 = list(mma = .91, mmb = 0, pp = 3, dep = .30, type = "non-uniform", DIFcntry = "Low"),

  cc74 = list(mma = .10, mmb = 0, pp = 1, dep = .50, type = "non-uniform", DIFcntry = "Low"),
  cc75 = list(mma = .39, mmb = 0, pp = 1, dep = .50, type = "non-uniform", DIFcntry = "Low"),
  cc76 = list(mma = .91, mmb = 0, pp = 1, dep = .50, type = "non-uniform", DIFcntry = "Low"),
  cc77 = list(mma = .10, mmb = 0, pp = 2, dep = .50, type = "non-uniform", DIFcntry = "Low"),
  cc78 = list(mma = .39, mmb = 0, pp = 2, dep = .50, type = "non-uniform", DIFcntry = "Low"),
  cc79 = list(mma = .91, mmb = 0, pp = 2, dep = .50, type = "non-uniform", DIFcntry = "Low"),
  cc80 = list(mma = .10, mmb = 0, pp = 3, dep = .50, type = "non-uniform", DIFcntry = "Low"),
  cc81 = list(mma = .39, mmb = 0, pp = 3, dep = .50, type = "non-uniform", DIFcntry = "Low"),
  cc82 = list(mma = .91, mmb = 0, pp = 3, dep = .50, type = "non-uniform", DIFcntry = "Low"),

  cc83 = list(mma = .10, mmb = 0, pp = 1, dep = .10, type = "non-uniform", DIFcntry = "Med"),
  cc84 = list(mma = .39, mmb = 0, pp = 1, dep = .10, type = "non-uniform", DIFcntry = "Med"),
  cc85 = list(mma = .91, mmb = 0, pp = 1, dep = .10, type = "non-uniform", DIFcntry = "Med"),
  cc86 = list(mma = .10, mmb = 0, pp = 2, dep = .10, type = "non-uniform", DIFcntry = "Med"),
  cc87 = list(mma = .39, mmb = 0, pp = 2, dep = .10, type = "non-uniform", DIFcntry = "Med"),
  cc88 = list(mma = .91, mmb = 0, pp = 2, dep = .10, type = "non-uniform", DIFcntry = "Med"),
  cc89 = list(mma = .10, mmb = 0, pp = 3, dep = .10, type = "non-uniform", DIFcntry = "Med"),
  cc90 = list(mma = .39, mmb = 0, pp = 3, dep = .10, type = "non-uniform", DIFcntry = "Med"),
  cc91 = list(mma = .91, mmb = 0, pp = 3, dep = .10, type = "non-uniform", DIFcntry = "Med"),

  cc92 = list(mma = .10, mmb = 0, pp = 1, dep = .30, type = "non-uniform", DIFcntry = "Med"),
  cc93 = list(mma = .39, mmb = 0, pp = 1, dep = .30, type = "non-uniform", DIFcntry = "Med"),
  cc94 = list(mma = .91, mmb = 0, pp = 1, dep = .30, type = "non-uniform", DIFcntry = "Med"),
  cc95 = list(mma = .10, mmb = 0, pp = 2, dep = .30, type = "non-uniform", DIFcntry = "Med"),
  cc96 = list(mma = .39, mmb = 0, pp = 2, dep = .30, type = "non-uniform", DIFcntry = "Med"),
  cc97 = list(mma = .91, mmb = 0, pp = 2, dep = .30, type = "non-uniform", DIFcntry = "Med"),
  cc98 = list(mma = .10, mmb = 0, pp = 3, dep = .30, type = "non-uniform", DIFcntry = "Med"),
  cc99 = list(mma = .39, mmb = 0, pp = 3, dep = .30, type = "non-uniform", DIFcntry = "Med"),
  cc100 = list(mma = .91, mmb = 0, pp = 3, dep = .30, type = "non-uniform", DIFcntry = "Med"),

  cc101 = list(mma = .10, mmb = 0, pp = 1, dep = .50, type = "non-uniform", DIFcntry = "Med"),
  cc102 = list(mma = .39, mmb = 0, pp = 1, dep = .50, type = "non-uniform", DIFcntry = "Med"),
  cc103 = list(mma = .91, mmb = 0, pp = 1, dep = .50, type = "non-uniform", DIFcntry = "Med"),
  cc104 = list(mma = .10, mmb = 0, pp = 2, dep = .50, type = "non-uniform", DIFcntry = "Med"),
  cc105 = list(mma = .39, mmb = 0, pp = 2, dep = .50, type = "non-uniform", DIFcntry = "Med"),
  cc106 = list(mma = .91, mmb = 0, pp = 2, dep = .50, type = "non-uniform", DIFcntry = "Med"),
  cc107 = list(mma = .10, mmb = 0, pp = 3, dep = .50, type = "non-uniform", DIFcntry = "Med"),
  cc108 = list(mma = .39, mmb = 0, pp = 3, dep = .50, type = "non-uniform", DIFcntry = "Med"),
  cc109 = list(mma = .91, mmb = 0, pp = 3, dep = .50, type = "non-uniform", DIFcntry = "Med"),
  
  cc110 = list(mma = 0, mmb = 0,   pp = 0, dep = .10,   type = "No DIF",  DIFcntry = "No DIF"),
  cc111 = list(mma = 0, mmb = 0,   pp = 0, dep = .30,   type = "No DIF",  DIFcntry = "No DIF"),
  cc112 = list(mma = 0, mmb = 0,   pp = 0, dep = .50,   type = "No DIF",  DIFcntry = "No DIF")
  # 
  # cc137 = list(mma = .10, mmb = 0, pp = 1, dep = .10, type = "non-uniform", DIFcntry = "High"),
  # cc138 = list(mma = .39, mmb = 0, pp = 1, dep = .10, type = "non-uniform", DIFcntry = "High"),
  # cc139 = list(mma = .91, mmb = 0, pp = 1, dep = .10, type = "non-uniform", DIFcntry = "High"),
  # cc140 = list(mma = .10, mmb = 0, pp = 2, dep = .10, type = "non-uniform", DIFcntry = "High"),
  # cc141 = list(mma = .39, mmb = 0, pp = 2, dep = .10, type = "non-uniform", DIFcntry = "High"),
  # cc142 = list(mma = .91, mmb = 0, pp = 2, dep = .10, type = "non-uniform", DIFcntry = "High"),
  # cc143 = list(mma = .10, mmb = 0, pp = 3, dep = .10, type = "non-uniform", DIFcntry = "High"),
  # cc144 = list(mma = .39, mmb = 0, pp = 3, dep = .10, type = "non-uniform", DIFcntry = "High"),
  # cc145 = list(mma = .91, mmb = 0, pp = 3, dep = .10, type = "non-uniform", DIFcntry = "High"),
  # 
  # cc146 = list(mma = .10, mmb = 0, pp = 1, dep = .30, type = "non-uniform", DIFcntry = "High"),
  # cc147 = list(mma = .39, mmb = 0, pp = 1, dep = .30, type = "non-uniform", DIFcntry = "High"),
  # cc148 = list(mma = .91, mmb = 0, pp = 1, dep = .30, type = "non-uniform", DIFcntry = "High"),
  # cc149 = list(mma = .10, mmb = 0, pp = 2, dep = .30, type = "non-uniform", DIFcntry = "High"),
  # cc150 = list(mma = .39, mmb = 0, pp = 2, dep = .30, type = "non-uniform", DIFcntry = "High"),
  # cc151 = list(mma = .91, mmb = 0, pp = 2, dep = .30, type = "non-uniform", DIFcntry = "High"),
  # cc152 = list(mma = .10, mmb = 0, pp = 3, dep = .30, type = "non-uniform", DIFcntry = "High"),
  # cc153 = list(mma = .39, mmb = 0, pp = 3, dep = .30, type = "non-uniform", DIFcntry = "High"),
  # cc154 = list(mma = .91, mmb = 0, pp = 3, dep = .30, type = "non-uniform", DIFcntry = "High"),
  # 
  # cc155 = list(mma = .10, mmb = 0, pp = 1, dep = .50, type = "non-uniform", DIFcntry = "High"),
  # cc156 = list(mma = .39, mmb = 0, pp = 1, dep = .50, type = "non-uniform", DIFcntry = "High"),
  # cc157 = list(mma = .91, mmb = 0, pp = 1, dep = .50, type = "non-uniform", DIFcntry = "High"),
  # cc158 = list(mma = .10, mmb = 0, pp = 2, dep = .50, type = "non-uniform", DIFcntry = "High"),
  # cc159 = list(mma = .39, mmb = 0, pp = 2, dep = .50, type = "non-uniform", DIFcntry = "High"),
  # cc160 = list(mma = .91, mmb = 0, pp = 2, dep = .50, type = "non-uniform", DIFcntry = "High"),
  # cc161 = list(mma = .10, mmb = 0, pp = 3, dep = .50, type = "non-uniform", DIFcntry = "High"),
  # cc162 = list(mma = .39, mmb = 0, pp = 3, dep = .50, type = "non-uniform", DIFcntry = "High"),
  # cc163 = list(mma = .91, mmb = 0, pp = 3, dep = .50, type = "non-uniform", DIFcntry = "High")
)
  
#=============================================================
## Fixed Factors
#=============================================================
 

#fixed factors
N <- 500 # number of persons
G <- 61  # number of groups
I <- 18 # number of items per testlet
TT <- 2 # number of testlets
ITT <- I*TT #total number of items

#=============================================================
#directory to save data
#=============================================================
path = "/N/project/Ar00157/Paper_3_MVM/Results/"

#=============================================================
## Functions
#=============================================================
 
twopl_trt <- function(a,b,theta, dep){
  sigma = rnorm(1, mean = 0, sd = dep)
  prob=exp(1.7*(a*(theta-b-sigma)))/(1+exp(1.7*(a*(theta-b-sigma))))
  return(prob)
}

twopl <- function(a,b,theta){
  prob=exp(1.7*(a*(theta-b)))/(1+exp(1.7*(a*(theta-b))))
  return(prob)
}

#=============================================================
#=============================================================
#                    Data Generation Function
#=============================================================
#=============================================================

data_generation = function(ii, cc){
#Calling conditions
  pp <- DIF_conditions[[cc]][["pp"]]
  type <- DIF_conditions[[cc]][["type"]]
  DIFcntry <- DIF_conditions[[cc]][["DIFcntry"]]
  dep <- DIF_conditions[[cc]][["dep"]]
  mmb <- DIF_conditions[[cc]][["mmb"]]
  mma <- DIF_conditions[[cc]][["mma"]]
  if(DIFcntry == "high") mmb = mmb*(-1) #if high DIF is negative, item is easier   
#set seed
set.seed(1992*ii)
#=====================================================================
#item parameter generation 
#=====================================================================
b <- round( stats::runif(ITT, -2.5, 1.14), 2) # difficulty 
a <- round( stats::runif(ITT,  0.4 , 2.28), 2) # discrimination



#Conditionals for DIF type 

if(type == "uniform"){
  # Difficulty parameter
  b1 = b
  #for first testlet
  b1[1:pp] = b[1:pp]+mmb
  #for second testlet
  b1[19:(19+(pp-1))] = b[19:(19+(pp-1))]+mmb
  # Discrimination parameter
  a1 = a
}else if(type == "non-uniform"){
  # Difficulty parameter
  b1 = b
  # Discrimination parameter
  a1 = a
  #for first testlet
  a1[1:pp] = a[1:pp]+mma
  #for second testlet
  a1[19:(19+(pp-1))] = a[19:(19+(pp-1))]+mma
}else if(type == "No DIF"){
  a1 = a
  b1 = b
}

# No DIF item parameters
bM <- matrix( b, nrow=N, ncol=ITT, byrow=TRUE )
aM <- matrix( a, nrow=N, ncol=ITT, byrow=TRUE )
# DIF item parameters
bMDIF <- matrix( b1, nrow=N, ncol=ITT, byrow=TRUE )
aMDIF <- matrix( a1, nrow=N, ncol=ITT, byrow=TRUE )


#=====================================================================
# simulate testlet's 
#=====================================================================
ut <- matrix(0,nrow=N, ncol=TT )
for (tt in 1:TT){
  ut[,tt] <- stats::rnorm( N, mean = 0 , sd=dep )
}
ut1 <- ut[, rep(1:TT,each=I) ]
#=====================================================================
# person parameter generation 
#=====================================================================
## generating overall means and sds per group
# generate group means (randomly generate group means)
theta <- matrix(0,nr=N,nc=G)
genmu <- matrix(0, nr = G, nc=2)
# for means  
for(g in 2:G){
  if(g==2 & DIFcntry == "Low"){genmu[g,1] <- -1.80
  }else if(g==2 & DIFcntry == "Med"){genmu[g,1] <- .11   
  }else if(g==2 & DIFcntry == "High"){genmu[g,1] <- 0.81}
  else genmu[g,1] <- runif(1,min=-.82,max=.81)
}

#for standard deviations
for (g in 1:G) {
  if(g==2) genmu[g,2] <- 1
  else genmu[g,2] <- runif(1, min= .60, max = 1.24)
}

# random generation of individual proficiency per group
for(g in 1:G){
  theta[,g] <- rnorm(N,genmu[g,1],genmu[g,2])   
}

#=====================================================================
# calculate response probability
#=====================================================================
resp = NULL
for (g in 1:G) {
  if(g == 2){
    prob <- matrix( stats::pnorm( aMDIF*theta[,g] + ut1 + bMDIF ), N, ITT)
    Y <- (matrix( stats::runif(N*ITT), N, ITT) < prob )*1  
  }else{
  prob <- matrix( stats::pnorm( aM*theta[,g] + ut1 + bM ), N, ITT)
  Y <- (matrix( stats::runif(N*ITT), N, ITT) < prob )*1
  colMeans(Y)}
  resp = rbind(resp, Y)
}

# Generated data responses
data_resp = data.frame(rep = ii,
                       cond = cc,
                       DIFgrp = DIFcntry,
                       DIFmmb = mmb,
                       DIFmma = mma,
                       DIFper = pp,
                       DIFtype = type,
                       effectT = dep,
                       group = rep(1:G, each = N),
                       resp)


# Generated item responses
All_items = data.frame(rep = ii,
                       cond = cc,
                       DIFgrp = DIFcntry,
                       DIFmmb = mmb,
                       DIFmma = mma,
                       DIFper = pp,
                       DIFtype = type,
                       effectT = dep,
                       item = paste0("X", 1:36),
                       a = a, 
                       b = b,
                       a_DIF = aMDIF[1,],
                       b_DIF = bMDIF[1,])

R.utils::saveObject(data_resp, file = paste0(path, "Data_generation/Responses_cond", cc,"_r", ii, ".Rbin"))
R.utils::saveObject(genmu, file = paste0(path, "Data_generation/GenGroup_cond", cc,"_r", ii, ".Rbin"))

data_gen = list(data_resp = data_resp, genmu = genmu, theta = theta, All_items = All_items)
return(data_gen)

}
  



#=============================================================
#=============================================================
#                           IRT analysis 
#=============================================================
#=============================================================

irt_analysis = function(ii, cc, data_resp, All_items, genmu){
  #Calling conditions
  pp <- DIF_conditions[[cc]][["pp"]]
  type <- DIF_conditions[[cc]][["type"]]
  DIFcntry <- DIF_conditions[[cc]][["DIFcntry"]]
  dep <- DIF_conditions[[cc]][["dep"]]
  mmb <- DIF_conditions[[cc]][["mmb"]]
  mma <- DIF_conditions[[cc]][["mma"]]
  if(DIFcntry == "high") mmb = mmb*(-1) #if high DIF is negative, item is easier   
  # group 
  grpmn <- genmu[,1]
  grpsd <- genmu[,2]
  
  # Multigroup estimation
  irt_int <- TAM::tam.mml.2pl(data_resp[, grepl("X", colnames(data_resp))],
                            irtmodel = "2PL",
                            est.variance = FALSE,
                            group = data_resp$group,
                            control = list(maxiter = 10000), verbose = FALSE)
  
  # rmsd
  fmod <- CDM::IRT.RMSD(irt_int)
  rmsd_tam = fmod$RMSD
  
  #save
  tam_results = data.frame(rep = ii,
                       cond = cc,
                       DIFgrp = DIFcntry,
                       DIFmmb = mmb,
                       DIFmma = mma,
                       DIFper = pp,
                       DIFtype = type,
                       effectT = dep,
                       rmsd_mod = "tam",
                       item = paste0("X", 1:36),
                       rmsd_tam)
  
  R.utils::saveObject(tam_results, file = paste0(path, "RMSD_IRT/rmsd_tam_cond", cc,"_r", ii, ".Rbin"))
  
  
  #save(irt_int, file = paste0(path, "Tam mod/tam_mod_int_cond", cc,"_r", ii, ".RData"))
  
  All_items = full_join(All_items, irt_int$item_irt, by = "item")

  # To obtain item parameter estimates per group
  
  for(g in 1:G){
    beta.fixed = irt_int$beta.fixed
    beta.fixed[3] <- 0
  
    #for the Items with maximum score of 0 error, change one response to be 1 
    item_Sum = colSums(data_resp[data_resp$group == g, grepl("X", colnames(data_resp))])
    #zero_sum =  (which(item_Sum == 0))
    #non_zero_sum = (which(item_Sum != 0))
    less_than_five_sum = (which(item_Sum < 5))
    five_or_more_sum = (which(item_Sum >= 5))
    
    irt_int_g <- TAM::tam.mml.2pl(data_resp[data_resp$group == g, names(five_or_more_sum)],
                             irtmodel = "2PL",
                             est.variance = FALSE, 
                             beta.fixed = beta.fixed,
                             beta.inits = beta.fixed,
                             control = list(maxiter = 10000), verbose = F, 
                             item.elim = TRUE)  
    
    #save irt estimates by group
    irt_items_g = as.data.frame(irt_int_g$item_irt)
    
    if(g == 2){
    ## mean - sigma
    est_a = irt_items_g$alpha
    est_b = irt_items_g$beta

    u <- mean(est_a)/mean(All_items$a_DIF)
    v <- mean(All_items$b_DIF) - u*mean(est_b)

    irt_items_g$rs_est_a <- est_a/u          #rescaled item discrimination
    irt_items_g$rs_est_b <- u*est_b + v      #rescaled item difficulty  
    }else{
      ## mean - sigma
    est_a = irt_items_g$alpha
    est_b = irt_items_g$beta

    u <- mean(est_a)/mean(All_items$a)
    v <- mean(All_items$b) - u*mean(est_b)

    irt_items_g$rs_est_a <- est_a/u          #rescaled item discrimination
    irt_items_g$rs_est_b <- u*est_b + v      #rescaled item difficulty
    }
    
    
    names(irt_items_g)[2:5] <- c(paste0("est_a_tam_G", g), paste0("est_b_tam_G", g),
                                 paste0("rsc_a_tam_G", g),paste0("rsc_b_tam_G", g))
    All_items[, paste0("N.responses_G", g)] = item_Sum
    All_items = dplyr::left_join(All_items, irt_items_g, by = "item")
    
   # save(irt_int_g, file = paste0(path, "Tam mod/tam_mod_g", g, "_cond", cc,"_r", ii, ".RData"))
    
  }
  
  R.utils::saveObject(All_items, file = paste0(path, "IRT_estimates/Items_cond", cc,"_r", ii, ".Rbin"))
  
  #=================================================================================================
  # RMSD computation
  #=================================================================================================
  
  # Compute group weights
  # grpmn <- apply(theta, 2, mean)
  # grpsd <- apply(theta, 2, sd)
  

  # defining the weights
  qt <- seq(-4,4,by=.2)
  n.qt <- length(qt) 
  grpwt <- matrix(NA,nr=n.qt,nc=G)
      
  for(g in 1:G){
    for(q in 1:n.qt){
        grpwt[q,g] <- dnorm(qt[q],grpmn[g],grpsd[g])/sum(dnorm(qt,grpmn[g],grpsd[g]))
          }
        }
  # Compute posterior and observed icc 
  obs <- array(NA,dim=c(n.qt,ITT,G))
  for(g in 1:G){
  ipar <- All_items %>%
    select(item,paste0("rsc_a_tam_G",g), paste0("rsc_b_tam_G",g)) %>%
    mutate(Cats = 2)
  names(ipar)[2:3]<- c("a", "b")
  gdata <- data_resp[data_resp$group == g, grepl("X", colnames(data_resp))]
  
  itemtrace <- matrix(NA,nr=n.qt,nc=ITT)
  itemtrace0 <- matrix(NA,nr=n.qt,nc=ITT)
  itemtrace1 <- matrix(NA,nr=n.qt,nc=ITT)
  itemtrace2 <- matrix(NA,nr=n.qt,nc=ITT)
  pseudo <- matrix(0,nr=n.qt,nc=ITT)
  deno <- matrix(0,nr=n.qt,nc=ITT)
  for(i in 1:ITT){
    for(q in 1:n.qt){
      if(ipar$Cats[i]==2){
        itemtrace[q,i] <- twopl(ipar$a[i],ipar$b[i],qt[q])
      }
     
    }
  }
  n.ex <- nrow(gdata)
  for(p in 1:n.ex){
    likelihood <- rep(1,n.qt)
    for(i in 1:ITT){
      if(ipar$Cats[i]==2){
        if(gdata[p,i]==1) likelihood <- likelihood*itemtrace[,i]
        else likelihood <- likelihood*(1-itemtrace[,i])
      }
      if(ipar$Cats[i]>2){
        if(gdata[p,i]==2) likelihood <- likelihood*itemtrace2[,i]
        if(gdata[p,i]==1) likelihood <- likelihood*itemtrace1[,i]
        if(gdata[p,i]==0) likelihood <- likelihood*itemtrace0[,i]
      }
    }
    posterior <- likelihood*grpwt[,g]
    
    # normalize posterior 
    expd <- sum(posterior)
    posterior <- (posterior/expd)
    
    # put this response pattern (pseudo counts)
    for(i in 1:ITT){
      deno[,i] <- deno[,i] + posterior
      if(ipar$Cats[i]==2){
        if(gdata[p,i]==1) pseudo[,i] <- pseudo[,i] + posterior
      }
      if(ipar$Cats[i]>2){
        if(gdata[p,i]==2) pseudo[,i] <- pseudo[,i] + posterior
        if(gdata[p,i]==1) pseudo[,i] <- pseudo[,i] + 0.5*posterior
      }
    }
   }
  obs[,,g] <- pseudo/deno
 }

# Compute expected icc
 exp <- array(NA,dim=c(n.qt,ITT,G))
 for(g in 1:G){
  ipar <- All_items %>%
    select(item,paste0("rsc_a_tam_G",g), paste0("rsc_b_tam_G",g)) %>%
    mutate(Cats = 2)
  names(ipar)[2:3]<- c("Slope", "Difficulty")
  for(i in 1:ITT){
    if(ipar$Cats[i]==2){
      for(q in 1:n.qt){
        exp[q,i,g] <- twopl(ipar$Slope[i],ipar$Difficulty[i],qt[q])
      }
    }
    # if(ipar$Cats[i]>2){
    #   for(q in 1:n.qt){
    #     tmp.p <- gpcm(ipar$Slope[i],ipar$Difficulty[i],c(ipar$IRT_Step1[i],ipar$IRT_Step2[i]),qt[q],ipar$Cats[i])
    #     exp[q,i,g] <- 0.5*tmp.p[2] + tmp.p[3]
    #   }
    # }
  }
 }

 # Compute MD and RMSDs
 # md <- matrix(NA,nr=ITT,nc=G)
 rmsd_irt <- matrix(NA,nr=ITT,nc=G)
 for(g in 1:G){
  #ipar <- itempar[itempar[,"Group"]==100+g,]
  grp.range.ind <- which(round(grpwt[,g],2) %in% 0.01)
  grp.range.length <- length(grp.range.ind)
  minind <- grp.range.ind[1]
  maxind <- grp.range.ind[grp.range.length]
  for(i in 1:ITT){
    #md[i,g] <- sum((obs[minind:maxind,i,g]-exp[minind:maxind,i,g])*grpwt[minind:maxind,g])
   rmsd_irt[i,g] <- sqrt(sum((obs[minind:maxind,i,g]-exp[minind:maxind,i,g])^2*grpwt[minind:maxind,g]))
    
  }
}   
  
  irt_results = data.frame(rep = ii,
                       cond = cc,
                       DIFgrp = DIFcntry,
                       DIFmmb = mmb,
                       DIFmma = mma,
                       DIFper = pp,
                       DIFtype = type,
                       effectT = dep,
                       rmsd_mod = "irt",
                       item = paste0("X", 1:36),
                       rmsd_irt)
  
  R.utils::saveObject(irt_results, file = paste0(path, "RMSD_IRT/rmsd_cond", cc,"_r", ii, ".Rbin"))
}

#=============================================================
#=============================================================
#                      TRT analysis 
#=============================================================
#=============================================================

trt_analysis = function(ii, cc, data_resp, All_items, genmu, dep){
  #Calling conditions
  pp <- DIF_conditions[[cc]][["pp"]]
  type <- DIF_conditions[[cc]][["type"]]
  DIFcntry <- DIF_conditions[[cc]][["DIFcntry"]]
  dep <- DIF_conditions[[cc]][["dep"]]
  mmb <- DIF_conditions[[cc]][["mmb"]]
  mma <- DIF_conditions[[cc]][["mma"]]
  if(DIFcntry == "high") mmb = mmb*(-1) #if high DIF is negative, item is easier   
  #Group item parameter estimation
  for(g in 1:G){
  trt_int_g <- estimates <- trt_items_g <- NULL
  
  resp <- data_resp[data_resp$group == g, grepl("X", colnames(data_resp))]
  I<-dim(resp)[1]
  J<-dim(resp)[2]
  n_testlet<-2
  ind_testlet<-rep(seq(1,2),each=18)
  data_trt<-list(n_student =I, n_item=J,
               n_testlet = n_testlet, ind_testlet = ind_testlet, Y=t(resp))
  
  
  if(g == 2 & DIFcntry == "Low"){
  trt_int_g <- rstan::sampling(object = twoTRT_code_DIFlow, data = data_trt, iter =300,chains = 1)  
  }else{
  trt_int_g <- rstan::sampling(object = twoTRT_code , data = data_trt, iter =300,chains = 1)
  }
  #load(paste0(path, "Rstan mod/trt_mod_g1_cond1_r1.RData"))
  #load(paste0(path, "Rstan mod/trt_mod_g1_cond1_r2.RData"))
  #save results
  #save(trt_int_g, file = paste0(path, "Rstan mod/trt_mod_g", g, "_cond", cc,"_r", ii, ".RData"))
  
  estimates <- rstan::summary(trt_int_g,par=c("alpha","beta"))[[1]]
  estimates <- as.data.frame(estimates)
  estimates$par <-  rownames(estimates)
  
  
  trt_items_g = data.frame(a = estimates[grepl("alpha", estimates$par), "mean"],
                           b = estimates[grepl("beta", estimates$par), "mean"],
                           a_Rhat = estimates[grepl("alpha", estimates$par), "Rhat"],
                           b_Rhat = estimates[grepl("beta", estimates$par), "Rhat"],
                           a_neff = estimates[grepl("alpha", estimates$par), "n_eff"],
                           b_neff = estimates[grepl("beta", estimates$par), "n_eff"])
  
  
  if(g == 2){
    ## mean - sigma
    est_a = trt_items_g$a
    est_b = trt_items_g$b

    u <- mean(est_a)/mean(All_items$a_DIF)
    v <- mean(All_items$b_DIF) - u*mean(est_b)

    trt_items_g$rs_est_a <- est_a/u          #rescaled item discrimination
    trt_items_g$rs_est_b <- u*est_b + v      #rescaled item difficulty  
    }else{
      ## mean - sigma
    est_a = trt_items_g$a
    est_b = trt_items_g$b

    u <- mean(est_a)/mean(All_items$a)
    v <- mean(All_items$b) - u*mean(est_b)

    trt_items_g$rs_est_a <- est_a/u          #rescaled item discrimination
    trt_items_g$rs_est_b <- u*est_b + v      #rescaled item difficulty
    }
  
 
  old_names <- c( "a",  "b",  "a_Rhat", "b_Rhat", "a_neff", "b_neff", "rs_est_a", "rs_est_b")
  new_names <- c( paste0("est_a_trt_G", g), paste0("est_b_trt_G", g),
                  paste0("a_Rhat_g", g), paste0("b_Rhat_g", g), 
                  paste0("a_neff_g", g), paste0("b_neff_g", g),
                  paste0("rsc_a_trt_G", g), paste0("rsc_b_trt_G", g)) #rsc means rescaled
         
  data.table::setnames(trt_items_g, old = old_names, new = new_names)
  
  
  All_items = cbind(All_items, trt_items_g)
  
  
  }
  
  R.utils::saveObject(All_items, file = paste0(path, "TRT_estimates/Items_cond", cc,"_r", ii, ".Rbin"))
  #=================================================================================================
  # RMSD computation
  #=================================================================================================
  grpmn <- genmu[,1]
  grpsd <- genmu[,2]
  

  # defining the weights
  qt <- seq(-4,4,by=.2)
  n.qt <- length(qt) 
  grpwt <- matrix(NA,nr=n.qt,nc=G)
      
  for(g in 1:G){
    for(q in 1:n.qt){
        grpwt[q,g] <- dnorm(qt[q],grpmn[g],grpsd[g])/sum(dnorm(qt,grpmn[g],grpsd[g]))
          }
        }
  # Compute posterior and observed icc 
  obs <- array(NA,dim=c(n.qt,ITT,G))
  for(g in 1:G){
  ipar <- All_items %>%
    select(item,paste0("rsc_a_trt_G",g), paste0("rsc_b_trt_G",g)) %>%
    mutate(Cats = 2)
  names(ipar)[2:3]<- c("a", "b")
  gdata <- data_resp[data_resp$group == g, grepl("X", colnames(data_resp))]
  
  itemtrace <- matrix(NA,nr=n.qt,nc=ITT)
  itemtrace0 <- matrix(NA,nr=n.qt,nc=ITT)
  itemtrace1 <- matrix(NA,nr=n.qt,nc=ITT)
  itemtrace2 <- matrix(NA,nr=n.qt,nc=ITT)
  pseudo <- matrix(0,nr=n.qt,nc=ITT)
  deno <- matrix(0,nr=n.qt,nc=ITT)
  for(i in 1:ITT){
    for(q in 1:n.qt){
      if(ipar$Cats[i]==2){
        itemtrace[q,i] <- twopl_trt(ipar$a[i],ipar$b[i],qt[q], dep=dep)
      }
     
    }
  }
  n.ex <- nrow(gdata)
  for(p in 1:n.ex){
    likelihood <- rep(1,n.qt)
    for(i in 1:ITT){
      if(ipar$Cats[i]==2){
        if(gdata[p,i]==1) likelihood <- likelihood*itemtrace[,i]
        else likelihood <- likelihood*(1-itemtrace[,i])
      }
      if(ipar$Cats[i]>2){
        if(gdata[p,i]==2) likelihood <- likelihood*itemtrace2[,i]
        if(gdata[p,i]==1) likelihood <- likelihood*itemtrace1[,i]
        if(gdata[p,i]==0) likelihood <- likelihood*itemtrace0[,i]
      }
    }
    posterior <- likelihood*grpwt[,g]
    
    # normalize posterior 
    expd <- sum(posterior)
    posterior <- (posterior/expd)
    
    # put this response pattern (pseudo counts)
    for(i in 1:ITT){
      deno[,i] <- deno[,i] + posterior
      if(ipar$Cats[i]==2){
        if(gdata[p,i]==1) pseudo[,i] <- pseudo[,i] + posterior
      }
      if(ipar$Cats[i]>2){
        if(gdata[p,i]==2) pseudo[,i] <- pseudo[,i] + posterior
        if(gdata[p,i]==1) pseudo[,i] <- pseudo[,i] + 0.5*posterior
      }
    }
   }
  obs[,,g] <- pseudo/deno
 }

# Compute expected icc
 exp <- array(NA,dim=c(n.qt,ITT,G))
 for(g in 1:G){
  ipar <- All_items %>%
    select(item,paste0("rsc_a_trt_G",g), paste0("rsc_b_trt_G",g)) %>%
    mutate(Cats = 2)
  names(ipar)[2:3]<- c("Slope", "Difficulty")
  for(i in 1:ITT){
    if(ipar$Cats[i]==2){
      for(q in 1:n.qt){
        exp[q,i,g] <- twopl_trt(ipar$Slope[i],ipar$Difficulty[i],qt[q],dep=dep)
      }
    }
    # if(ipar$Cats[i]>2){
    #   for(q in 1:n.qt){
    #     tmp.p <- gpcm(ipar$Slope[i],ipar$Difficulty[i],c(ipar$IRT_Step1[i],ipar$IRT_Step2[i]),qt[q],ipar$Cats[i])
    #     exp[q,i,g] <- 0.5*tmp.p[2] + tmp.p[3]
    #   }
    # }
  }
 }

 # Compute MD and RMSDs
 # md <- matrix(NA,nr=ITT,nc=G)
 rmsd_trt <- matrix(NA,nr=ITT,nc=G)
 for(g in 1:G){
  #ipar <- itempar[itempar[,"Group"]==100+g,]
  grp.range.ind <- which(round(grpwt[,g],2) %in% 0.01)
  grp.range.length <- length(grp.range.ind)
  minind <- grp.range.ind[1]
  maxind <- grp.range.ind[grp.range.length]
  for(i in 1:ITT){
    #md[i,g] <- sum((obs[minind:maxind,i,g]-exp[minind:maxind,i,g])*grpwt[minind:maxind,g])
   rmsd_trt[i,g] <- sqrt(sum((obs[minind:maxind,i,g]-exp[minind:maxind,i,g])^2*grpwt[minind:maxind,g]))
    
  }
}   
  
  trt_results = data.frame(rep = ii,
                       cond = cc,
                       DIFgrp = DIFcntry,
                       DIFmmb = mmb,
                       DIFmma = mma,
                       DIFper = pp,
                       DIFtype = type,
                       effectT = dep,
                       rmsd_mod = "trt",
                       item = paste0("X", 1:36),
                       rmsd_trt)
  
  
  R.utils::saveObject(trt_results, file = paste0(path, "RMSD_TRT/rmsd_trt_cond", cc,"_r", ii, ".Rbin"))   
          
     
}
  

#=============================================================
#=============================================================
#               Function to Paralelize Simulation
#=============================================================
#=============================================================         
data_gen_analysis = function(cc, ii){
  #for (ii in 1:200) {
  
  # run data generation function
  data_gen = data_generation(ii = ii,cc = cc)
  
  # extract objects useful for analyzing the data
  data_resp = data_gen[["data_resp"]]
  theta = data_gen[["theta"]]
  genmu = data_gen[["genmu"]]
  All_items = data_gen[["All_items"]]
  
  #IRT analysis
  irt_analysis(ii = ii, cc = cc, data_resp = data_resp, All_items = All_items, genmu = genmu)
  
  #TRT analysis
  trt_analysis(ii = ii, cc = cc,data_resp = data_resp, All_items = All_items, genmu = genmu, dep = dep)
  
  
  #} #end for loop through replications
} #end of function

#=============================================================
#=============================================================
#                         Parallelize
#=============================================================
#=============================================================
# no_cores <- detectCores() 
# registerDoParallel(no_cores) 
mclapply(1:112, function(cc) data_gen_analysis(ii = 121, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 122, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 123, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 124, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 125, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 126, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 127, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 128, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 129, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 130, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 131, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 132, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 133, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 134, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 135, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 136, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 137, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 138, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 139, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 140, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 141, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 142, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 143, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 144, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 145, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 146, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 147, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 148, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 149, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 150, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 151, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 152, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 153, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 154, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 155, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 156, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 157, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 158, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 159, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 160, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 161, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 162, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 163, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 164, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 165, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 166, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 167, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 168, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 169, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 170, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 171, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 172, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 173, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 174, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 175, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 176, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 177, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 178, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 179, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 180, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 181, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 182, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 183, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 184, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 185, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 186, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 187, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 188, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 189, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 190, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 191, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 192, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 193, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 194, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 195, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 196, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 197, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 198, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 199, cc), mc.cores = 128)
mclapply(1:112, function(cc) data_gen_analysis(ii = 200, cc), mc.cores = 128)

# 
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 2), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 3), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 4), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 5), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 6), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 7), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 8), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 9), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 10), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 11), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 12), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 13), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 14), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 15), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 16), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 17), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 18), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 19), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 20), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 21), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 22), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 23), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 24), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 25), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 26), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 27), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 28), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 29), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 30), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 31), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 32), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 33), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 34), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 35), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 36), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 37), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 38), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 39), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 40), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 41), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 42), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 43), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 44), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 45), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 46), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 47), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 48), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 49), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 50), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 51), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 52), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 53), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 54), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 55), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 56), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 57), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 58), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 59), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 60), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 61), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 62), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 63), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 64), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 65), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 66), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 67), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 68), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 69), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 70), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 71), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 72), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 73), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 74), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 75), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 76), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 77), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 78), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 79), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 80), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 81), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 82), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 83), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 84), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 85), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 86), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 87), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 88), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 89), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 90), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 91), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 92), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 93), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 94), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 95), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 96), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 97), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 98), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 99), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 100), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 101), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 102), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 103), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 104), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 105), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 106), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 107), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 108), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 109), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 110), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 111), mc.cores = 128)
# mclapply(1:20, function(ii) data_gen_analysis(ii, cc = 112), mc.cores = 128)
