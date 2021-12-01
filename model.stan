//   Copyright 2021 Neil Ferguson, Imperial College London
//
//   Licensed under the Apache License, Version 2.0 (the "License");
//   you may not use this file except in compliance with the License.
//   You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
//   Unless required by applicable law or agreed to in writing, software
//   distributed under the License is distributed on an "AS IS" BASIS,
//   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//   See the License for the specific language governing permissions and
//   limitations under the License.

functions {

  real logVL(real t, row_vector vp, real peak_mult, int K, int do_simp_vl) {
    real vl;
    real x;
    real y;
// only K=3 used for paper - exponential growth, then exponential decay
    if(K==5) {
      y=1.0/(1.0+peak_mult*vp[1]/vp[5]);  
      x=(1.0-y)*vp[3]+y*vp[4];
      vl=log(peak_mult*vp[1])+log(vp[2]+x)-log_sum_exp(log(x)-vp[2]*t,log(vp[2])-log_sum_exp(log(1.0-y)-vp[3]*t,log(y)-vp[4]*t));
    } else if(K==4) {
        vl=log(peak_mult*vp[1])+log(vp[2]+vp[3]+vp[4])-log_sum_exp(log_sum_exp(log(vp[3])-vp[2]*t,log(vp[2])+vp[3]*t),log(vp[4]));       
    } else {
        if(do_simp_vl==1) {
// piecewise linear
          if(t<=0)
            vl=log(peak_mult*vp[1])+vp[2]*t;
        else
            vl=log(peak_mult*vp[1])-vp[3]*t;
        } else
// smooth function which still gives exponential growth and decay
            vl=log(peak_mult*vp[1])+log(vp[2]+vp[3])-log_sum_exp(log(vp[3])-vp[2]*t,log(vp[2])+vp[3]*t); 
    }
 
    return vl;
  }

}

data {
  int<lower=1> N;             // number of subjects
  int<lower=0> subj_voc[N];   // 0/1 marker for VOC
  int<lower=0> subj_age_cat[N];// 0/1 marker for age group
  vector<lower=0>[N] subj_age;  // age in years
  int<lower=1> NS;            // number of studies included
  int<lower=1> subj_study[N]; // study index for each subject
  int<lower=1> M;             // number of observations
  int<lower=0> obs_id[M];     // subject id for observation
  vector[M] obs_day;          // observation day
  vector<lower=0>[M] obs_ct;  // 40-Ct (where 40 is LOD)
  
// values for priors params which don't vary across subjects
  real<lower=0> pr_v_s;        // prior mean of SD of PCR measurement of CT
  real<lower=0> pr_v_s_sd;     // SD of prior for SD of PCR measurement of CT
  real<lower=0> pr_fp;         // prior mean of prob of false positive
  real<lower=0> pr_fp_sd;      // SD of prior for mean of prob of false positive
  real<lower=0> pr_fp_v_s;     // prior for SD of "error" CT distribution 
  real<lower=0> pr_fp_v;       // width of prior for "error" CT distribution mean and SD
  real<lower=0> pr_age_slope;  // width of prior of slope of peak VL with age
  
// values for hyperpriors
  int<lower=1> K; // number of correlated viral kinetic parameters
  vector<lower=0>[K] v_min;     // min value of each param
  vector[K] hp_v; // hyperprior on kinetic parameter group means
  vector<lower=0>[K] hp_v_sd; // hyperprior on kinetic parameter group SDs
  vector<lower=0>[K] hp_v_sd_sd; // hyperprior on within-group SD
  real<lower=0> hp_t_max_sd; // sigma for t_max
  real<lower=0> eta; // prior on correl matrix

  real<lower=0> ct_to_vl_sc[NS]; // scaling from CT to natural log VL
  real<lower=0> ct_to_vl_offset[NS]; // offset from CT to natural log VL
  int<lower=0,upper=1> fit_t_max_sd;
  int<lower=0,upper=1> fit_corr;
  int<lower=0,upper=1> fit_age_slope;
  int<lower=0,upper=1> do_simp_err;
  int<lower=0,upper=1> do_simp_vl;
  int<lower=1> NG; // number of VOC groups to be fitted (1 for only 1 group)
  int<lower=1> NV; // number of VOC types 
  int<lower=1,upper=3> group_by_which_cat; // 1 = voc, 2 = age, 3 = study
  int<lower=0> Nsamp; // number of samples from posterior per iteration
}

transformed data {
 
  vector[N] subj_age_offset=log(subj_age)-log(50);  // age in years
  vector[K] zeros = rep_vector(0, K);
  int<lower=1> group_index[N];
  if(NG==1 || group_by_which_cat==0) {
  for(n in 1:N)
    group_index[n]=1;     
  } else if(group_by_which_cat==1) {
  for(n in 1:N)
    group_index[n]=subj_voc[n]+1; 
  } else if (group_by_which_cat==2) {
  for(n in 1:N)
    group_index[n]=subj_age_cat[n]+1;    
  } else {
  for(n in 1:N)
    group_index[n]=subj_study[n]; 
  }
}


parameters {
// hyperparameters
  matrix[K,NG] v;        // viral growth rate
  vector<lower=0>[K] v_sd;     // coeff of var

  cholesky_factor_corr[K] Lc;
  
  real<lower=0> t_max_sd;    // standard deviation of tmax
  
// parameters  
  real<lower=1> v_s;          // sigma for VL
  vector[N] p_t_max;          // infection time for each subject
  matrix[N, K] n_v;            

  real l_fp;           // false negative probability
  real fp_ct_mean;
  real<lower=1> fp_ct_sd;
  
  real age_slope;  // slope of peak VL dependence on age
}

transformed parameters {
  real ll_total;              // log likelihood for all subjects
  vector[M] log_lik;          // log likelihood for each observation
  vector[N] age_slope_mult;   // subject specific age-dependent multiplier for PeakVL
  matrix<lower=0>[N, K] p_v;
  real t_sd;
  real fp;
  
  if(fit_age_slope==0)
    age_slope_mult=rep_vector(1,N);
  else
    age_slope_mult=exp(age_slope * subj_age_offset);

// enclose rest in a block to avoid variables being output   
  {
    real vl;            // modelled log10 viral titre
    real prob;          // probability test is postive or negative
    real l_fpm;
    real e_v_s;
    int id;
    real n_voc;
    real n_nvoc;

    if(fit_t_max_sd==1)
      t_sd=t_max_sd;
    else
      t_sd=hp_t_max_sd;

    fp=exp(-l_fp);
    l_fpm=log(1.0-fp);

		for(n in 1:N) 
		  p_v[n,] = v_min' + exp(v[,group_index[n]]' + v_sd' .* n_v[n,]);
// don't fit variation in peakVL by group when fitting age slope
		if(fit_age_slope==1 && NG>1) 
		  p_v[,1] = v_min[1] + exp(v[1,1] + v_sd[1] * n_v[,1]);

 // calculate the likelihood for each possible days of infection for each subject 
    if(do_simp_err==1)
      for (m in 1:M) {
          id=obs_id[m];
          vl=logVL(obs_day[m]-p_t_max[id],p_v[id,],age_slope_mult[id],K,do_simp_vl);
          if(obs_ct[m] <= 0) { // undetected
            log_lik[m] = log_sum_exp(l_fpm+normal_lcdf(0 | ct_to_vl_sc[subj_study[id]]*(vl-ct_to_vl_offset[subj_study[id]]), v_s), -l_fp);  // undetectable
          } else {  // known Ct obs above detection thre
            log_lik[m] = l_fpm+normal_lpdf(obs_ct[m] | ct_to_vl_sc[subj_study[id]]*(vl-ct_to_vl_offset[subj_study[id]]), v_s );  // detectable
          }
        }    
    else
      for (m in 1:M) {
        id=obs_id[m];
        vl=logVL(obs_day[m]-p_t_max[id],p_v[id,],age_slope_mult[id],K,do_simp_vl);
        if(obs_ct[m] <= 0) { // undetected
          log_lik[m] = log_sum_exp(l_fpm+normal_lcdf(0 | ct_to_vl_sc[subj_study[id]]*(vl-ct_to_vl_offset[subj_study[id]]), v_s), -l_fp+normal_lcdf(0 | fp_ct_mean, fp_ct_sd));  // undetectable
        } else {  // known Ct obs above detection thre
          log_lik[m] = log_sum_exp(l_fpm+normal_lpdf(obs_ct[m] | ct_to_vl_sc[subj_study[id]]*(vl-ct_to_vl_offset[subj_study[id]]), v_s ), -l_fp+normal_lpdf(obs_ct[m] | fp_ct_mean, fp_ct_sd) );  // detectable
        }
      }
// sums the likelihood (not log likelihood) over all possible days of infection for each subject
    ll_total=sum(log_lik);
  }

}


model {

//hyperpriors 
  for(i in 1:NG)
    v[,i]~normal(hp_v,hp_v_sd);
  v_sd~normal(0, hp_v_sd_sd);
  Lc ~ lkj_corr_cholesky(eta);
  
// per subject params - non-centred representation
  if(fit_corr==1) {
    for(n in 1:N)
      n_v[n,]~multi_normal_cholesky(zeros, Lc);
  } else {
      for(n in 1:N)
      n_v[n,]~normal(0, 1);
  }

  t_max_sd~normal(0,hp_t_max_sd); # hyper prior on t_sd - not used for paper
  
  p_t_max~normal(0,t_sd); // time of peak VL (max observed CT in data at day=0)


// params which don't vary across subjects

  v_s~normal(pr_v_s,pr_v_s_sd); // sd of CT measurements 
  
  l_fp~normal(pr_fp,pr_fp_sd); //  prior for "error" probability
  fp_ct_mean~normal(0,pr_fp_v);
  fp_ct_sd~normal(pr_fp_v_s,pr_fp_v);
  
  age_slope~normal(0,pr_age_slope); // slope of peak VL dependence on age

  target+=ll_total;

}

generated quantities {
  
  matrix[K,K] Omega; //correlation matrix
  matrix[K,NG] v_mean; // mean peakVL, growth/decay rates
  vector[K] v_cv; // Coeff of var of params
  vector[N] p_ln_v_max; // peak ln(VL) for each subject
  vector[N] p_ln_v_delta; // not used for paper
// variables calculated as averages over voc groups
  vector[NV+1] v_a_voc; // growth rate
  vector[NV+1] v_b_voc; // decay rate
  vector[NV+1] v_c_voc; // not used for paper
  vector[NV+1] v_max_voc; // peak VL
  vector[NV+1] v_delta_voc; //not used for paper
  vector[NV+1] t_a_voc; // time to peak (truncated at 14 days) - not used for paper
  vector[NV+1] t_b_voc; // time from peak (truncated at 28 days) - not used for paper
  vector[NV+1] t2_a_voc; // doubling time - not used for paper
  vector[NV+1] t2_b_voc; // halving time - not used for paper
  vector[NV+1] t2_c_voc; // not used for paper
  vector[NV+1] auc_voc; // area under ln(VL) curve - not used for paper
  vector<lower=0>[42*N] pred_ct; // VL trajectory for each subject
  vector<lower=0>[N] prop_inf; // not used for paper
  vector<lower=0>[N] AUC; // not used for paper
  vector<lower=0>[168] inf; // not used for paper
  vector<lower=0>[168] tot_inf; // not used for paper
  real infT_L; // not used for paper
  real infT_U; // not used for paper
  real infT_R; // not used for paper
  vector[N] p_infT_L; // not used for paper
  vector[N] p_infT_U; // not used for paper
  vector[N] p_infT_R; // not used for paper
  vector[N] p_vT_L; // not used for paper
  vector[N] p_vT_U; // not used for paper
  
  matrix[Nsamp,K] p_v_s[NG]; // sampled param values - not used for paper
  
  {
    vector[NV+1] n_voc;
    int id;
    int i;
    real tpi;
    real vl;
    real last_inf;
    real last_VL;
    real tot_vl;
    real lod;

     // correlation matrix
    
    Omega = multiply_lower_tri_self_transpose(Lc); 
    for(n in 1:K)
      Omega[n,n]=Omega[n,n]+0.01*normal_rng(0,1);  // to avoid NaN Rhat for diagonal elements

    
     // calculate log v_max and v_delta
     
     for(n in 1:N) {
       p_ln_v_max[n]=log(p_v[n,1]);
       p_ln_v_delta[n]=0;
     }
     
     if(K>3)
        for(n in 1:N) {
          p_ln_v_delta[n]=log(p_v[n,4]);
        }    
     // calculate means and CVs of params
     
     for(j in 1:NG)
        v_mean[,j]=exp(v[,j] + v_sd .* v_sd/2);
     v_cv=sqrt(exp(v_sd .* v_sd)-1); 

     
     // calculate predicted values and infectious profile
     
    i=0;
    for(x in 1:168) {
      tot_inf[x]=0;
    }
    for (n in 1:N) {
      lod=ct_to_vl_offset[subj_study[n]]; // limit of detection
      last_inf=0.0;
      last_VL=-1.0;
      tot_vl=0.0;
      p_vT_L[n]=-14;
      p_vT_U[n]=28;
      for(x in 1:168) {
        tpi=(x-1)*0.25-14;
        vl=logVL(tpi,p_v[n,],age_slope_mult[n],K,do_simp_vl)-lod;
        if(vl>=0) tot_vl += vl*0.25;  // calc AUC above lod
        inf[x]=last_inf+exp(vl+lod); // calc inf for abs VL
        last_inf=inf[x];
        tot_inf[x]+=last_inf;
        if(last_VL<=0 && vl>0) p_vT_L[n]=-tpi; // calc vT_L above lod
        if(vl>0) p_vT_U[n]=tpi; // calc vT_U above lod
        last_VL=vl;
        if(floor(tpi)==tpi) {
          i += 1;
          vl=logVL(tpi-p_t_max[n],p_v[n,],age_slope_mult[n],K,do_simp_vl)-lod;
          if(vl < 0 || is_inf(vl)) vl=0;  // report VL below detection limit as 0
          pred_ct[i] = ct_to_vl_sc[subj_study[n]]*vl;
        }
      }
      prop_inf[n]=last_inf;
      AUC[n]=tot_vl;
      for(x in 1:168) inf[x]/=last_inf;
      for(x in 1:167) {
        if(inf[x+1]>0.025 && inf[x]<=0.025) p_infT_L[n]=-(x*0.25-14);
        if(inf[x+1]>0.975 && inf[x]<=0.975) p_infT_U[n]=(x-1)*0.25-14;      
      }
      p_infT_R[n]=p_infT_U[n]+p_infT_L[n];
    }
    last_inf=tot_inf[168];
    for(n in 1:N) prop_inf[n]/=last_inf;
    for(x in 1:168) tot_inf[x]/=last_inf;
    for(x in 1:167) {
      if(tot_inf[x+1]>0.025 && tot_inf[x]<=0.025) infT_L=-(x*0.25-14);
      if(tot_inf[x+1]>0.975 && tot_inf[x]<=0.975) infT_U=(x-1)*0.25-14;      
    }
    infT_R=infT_U+infT_L;
    
    for(n in 1:NV+1) {
      n_voc[n]=0;
      auc_voc[n]=0;
      v_a_voc[n]=0;
      v_b_voc[n]=0;
      v_c_voc[n]=0;
      v_max_voc[n]=0;
      v_delta_voc[n]=0;
      t_a_voc[n]=0;
      t_b_voc[n]=0;
      t2_a_voc[n]=0;
      t2_b_voc[n]=0;
      t2_c_voc[n]=0;
    }
    
    for(n in 1:N) {
      i=subj_voc[n]+1;
      n_voc[i]+=1;
      auc_voc[i]+=AUC[n];
      v_max_voc[i]+=p_ln_v_max[n];
      v_a_voc[i]+=p_v[n,2];
      v_b_voc[i]+=p_v[n,3];
      if(K>3) v_delta_voc[i]+=p_ln_v_delta[n];
      if(K==5) {
        v_c_voc[i]+=p_v[n,4];
        t2_c_voc[i]+=(log(10)/p_v[n,4]);
      }
      t_a_voc[i]+=p_vT_L[n];
      t_b_voc[i]+=p_vT_U[n];
      t2_a_voc[i]+=(log(10)/p_v[n,2]);
      t2_b_voc[i]+=(log(10)/p_v[n,3]);

    }
    for(n in 1:NV) {
      n_voc[NV+1]=n_voc[NV+1]+n_voc[n];
      auc_voc[NV+1]=auc_voc[NV+1]+auc_voc[n];
      v_max_voc[NV+1] =v_max_voc[NV+1]+v_max_voc[n];
      v_delta_voc[NV+1] = v_delta_voc[NV+1]+v_delta_voc[n];
      v_a_voc[NV+1] =v_a_voc[NV+1]+v_a_voc[n];
      v_b_voc[NV+1] =v_b_voc[NV+1]+v_b_voc[n];
      v_c_voc[NV+1] =v_c_voc[NV+1]+v_c_voc[n];
      t_a_voc[NV+1] =t_a_voc[NV+1]+t_a_voc[n];
      t_b_voc[NV+1] =t_b_voc[NV+1]+t_b_voc[n];
      t2_a_voc[NV+1] =t2_a_voc[NV+1]+t2_a_voc[n];
      t2_b_voc[NV+1] =t2_b_voc[NV+1]+t2_b_voc[n];      
    }
 
    for(n in 1:NV+1) {
      auc_voc[n] /= n_voc[n];
      v_max_voc[n] /= n_voc[n];
      v_delta_voc[n] /= n_voc[n];
      v_a_voc[n] /= n_voc[n];
      v_b_voc[n] /= n_voc[n];
      v_c_voc[n] /= n_voc[n];
      t_a_voc[n] /= n_voc[n];
      t_b_voc[n] /= n_voc[n];
      t2_a_voc[n] /= n_voc[n];
      t2_b_voc[n] /= n_voc[n];
    }
  }    
// this generates simulated VL trajectories for each group - not used in paper  
  {
    vector[K] n_v_s;
    
    if(fit_corr==1) {
      for(i in 1:NG)
        for(n in 1:Nsamp)
        {
          n_v_s=multi_normal_cholesky_rng(zeros, Lc);
    		  p_v_s[i,n,]=v_min' + exp(v[,i]' + v_sd' .* n_v_s');
        }
    } else {
      vector[K] ones = rep_vector(1, K);      
      for(i in 1:NG)
        for(n in 1:Nsamp)
        {
          n_v_s=to_vector(normal_rng(zeros, ones));
    		  p_v_s[i,n,]=v_min' + exp(v[,i]' + v_sd' .* n_v_s');
        }
    }
  }
}

