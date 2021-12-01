##   Copyright 2021 Neil Ferguson, Imperial College London
##
##   Licensed under the Apache License, Version 2.0 (the "License");
##   you may not use this file except in compliance with the License.
##   You may obtain a copy of the License at
##
##     http://www.apache.org/licenses/LICENSE-2.0
##
##   Unless required by applicable law or agreed to in writing, software
##   distributed under the License is distributed on an "AS IS" BASIS,
##   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
##   See the License for the specific language governing permissions and
##   limitations under the License.

source("libs.r")
source("plots.r")
source("read_data.r")

options(mc.cores=12)
options(max.print=10000)

#for reproducibility
set.seed(610981)

# data sets included in each set to be fitted. Only 2 and 4 used for ATTACCC
study_id_set <- list(1,2,c(1,2),3,4)
# input data sets
data_files <- c("ct_dat_refined.csv","20210920_ATACCC_L.csv","20210920_ATACCC_E_L.csv","SG_Ct_GT5.csv")
study_id <- c(1,2,3,4)
data_names <- c("kis","ataN","kis_ata","ataE","sg")
# different grouping types in hierachical model - only voc used here
group_names <- c("_none","_voc","_age","_study")
# parameters (by dataset) for translation of Ct values to ln(VL)
data_scale <- c(1.5677, 1.418, 1.394, 1.418)
data_offset <- c(6.119, 3.435, 3.145, 3.435)
# filename tags for different inclusion criteria for data included in fit. Only F used here
red_names <- c("A","R","G","M","F")
# values of hyperpriors for growth (first entry is "uninformative", second "informative")
hp_v_g <- c(1.0,0.9)
hp_v_sd_g <- c(1.4,0.2)
# file name tags for uninformative vs informative priors on growth
inform_prior_fn=c("t2","p1")

# prams common for all runs
group_by_which_cat <- 1  # 0 for no grouping, 1 for voc
red_to_full_prof <- 4  # see read_data.r - 4 removes subjects without full followup unless they are incident cases 
fit_model <- TRUE # set to FALSE to load previous results
inc_vacc <- 2  # 0 for exclude vacc, 1 for include, 2 for include as separate VOC category
t_sd=4.0 # prior on SD of timing of peak VL

# set up vectors to define params which vary between runs
runs=12
ng.runs=c(4,1,4,4,1,1,4,1,4,4,1,1) # number of groups to fit
corr.runs=c(1,1,0,1,0,1,1,1,0,1,0,1) # flag specifying if correlation structure fitted
inform_prior.runs=c(1,1,1,1,0,0,1,1,1,1,0,0) # flag specifying if informative priors on growth to be used
age_slope.runs=c(1,0,0,0,0,0,1,0,0,0,0,0) # flag specifying on whether age dependence of peak VL fitted
data.runs=c(2,2,2,2,2,2,4,4,4,4,4,4) # data set to use - 2=ORF1ab, 4=E

#
# run fits
# much better distributed over multiple servers
# each fit requires ~50GB+ of RAM, takes 4-10h (including loo_cv) on 10 cores
#
for(i in runs) {
  num_groups <- ng.runs[i]
  data_set_num <- data.runs[i] # 1=kis, 2=ata, 3=both, 4=ata_E, 5=Singapore
  
#file name for outputs
    model_name <- paste0("f_",data_names[data_set_num],"_4",
                       inform_prior_fn[inform_prior.runs[run]+1],
                       red_names[red_to_full_prof+1],
                       "_",num_groups,"group", 
                       group_names[group_by_which_cat+1],
                       "_incvacc",inc_vacc,"_tsd",t_sd,
                       "_corr",fit_corr)
  if(age_slope.runs[i]==1) model_name <- paste0(model_name,"_age1")
  
  data_list=read_data(data_files[study_id_set[[data_set_num]]],red_to_full_prof,inc_vacc,study_id[study_id_set[[data_set_num]]])
  
  data.stan=list(N=data_list$N,
                 subj_voc=data_list$subjects$voc,
                 subj_age_cat=data_list$subjects$age_cat,
                 subj_age=data_list$subjects$age,
                 NS=4,
                 subj_study=data_list$subjects$study_id,
                 M=data_list$M,
                 obs_id=data_list$data$id,
                 obs_day=data_list$data$TestDateIndex,
                 obs_ct=data_list$data$invCT,
                 pr_v_s=3.0,        # mean of prior for SD of CT measurements
                 pr_v_s_sd=3.0,     # SD of prior for SD of CT measurements
                 pr_fp=5.0,         # mean of prior for error proportion
                 pr_fp_sd=2.0,      # SD of prior of error proportion
                 pr_fp_v=3.0,       # prior for mean of error CT distribution
                 pr_fp_v_s=3.0,     # prior for SD of error CT distribution
                 pr_age_slope=0.5, # width of prior of slope of peak VL with age
                 K=3,
                 v_min=c(0,0,0),
                 hp_v=c(15,hp_v_g[inform_prior.runs[run]+1],0.5), # hyperprior on kinetic parameter group means
                 hp_v_sd=c(15,hp_v_sd_g[inform_prior.runs[run]+1],1.4), # hyperprior on kinetic parameter group SDs
                 hp_v_sd_sd=c(10,1.0,1.0), # hyper prior on within-group variation
                 hp_t_max_sd=t_sd,  # hyperprior for SD of t_max
                 eta=1,
                 ct_to_vl_sc=data_scale,
                 ct_to_vl_offset=data_offset,
                 fit_t_max_sd=0,
                 fit_corr=corr.runs[i],
                 fit_age_slope=age_slope.runs[i],
                 do_simp_err=0, 
                 do_simp_vl=0, # if this is set to 1, fit uses piece wise linear log VL model
                 NG=num_groups, # 2 if v[] fitted separately for VOC and non-VOC, 1 otherwise
                 NV=4,  # num of VOC categories
                 group_by_which_cat=group_by_which_cat, # 1 = voc, 2 = age, 3 = study
                 Nsamp=100)
  
# initial values for stan chains - necessary to constrain some to avoid numerical issues
    
  initf<- function(chain_id) {
    list(l_fp = 3+(runif(1)-0.5), v_s=3+(runif(1)-0.5), age_slope=0.1+0.3*runif(1), 
         fp_ct_sd=3+(runif(1)-0.5), fp_ct_mean=0+(runif(1)-0.5), 
         v=array(c(15,1,0.5),dim=c(3,num_groups)),v_sd=c(2,0.5,0.5))
  }
  
  
  if(fit_model==TRUE) {
 
# run stan model - note high adapt_delta, low stepsize, fairly long chains
    fit = stan(file="model.stan", chains=10, cores=10, data=data.stan,
               iter=7000, init=initf, init_r=0.1, warmup=3000, seed=497101, thin=2, 
               control = list(adapt_delta = 0.994, max_treedepth = 13, stepsize=0.01))
    saveRDS(fit,paste0(model_name,".rds"))
    check_hmc_diagnostics(fit)
# need to have moment_match=TRUE to get valuely decent loo_cv diagnostics    
    loo_res=loo(fit,cores=10,moment_match=TRUE)
    saveRDS(loo_res,paste0(model_name,"_loo.rds"))
    loo_res
    
  } else {
    
    fit=readRDS(paste0(model_name,".rds"))
    loo_res=readRDS(paste0(model_name,"_loo.rds"))
    
  }
# generate plot of fit to data
  
  plot_graphs(data_list, fit, model_name, num_groups, inc_vacc, group_names=c("pre-Alpha","Alpha","Delta","Delta(vacc)"))
}

