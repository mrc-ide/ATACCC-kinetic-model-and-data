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
library(gridExtra)

options(mc.cores=12)
options(max.print=10000)



### compare ORF1ab data data fits ("N" is a old tag for "New")
loo_res.N4pca=readRDS("f_ataN_4p1F_4group_voc_incvacc2_tsd4_corr1_age1_loo.rds")
loo_res.N1pc=readRDS("f_ataN_4p1F_1group_voc_incvacc2_tsd4_corr1_loo.rds")
loo_res.N4p=readRDS("f_ataN_4p1F_4group_voc_incvacc2_tsd4_corr0_loo.rds")
loo_res.N4pc=readRDS("f_ataN_4p1F_4group_voc_incvacc2_tsd4_corr1_loo.rds")
loo_res.N1t=readRDS("f_ataN_4t2F_1group_voc_incvacc2_tsd4_corr0_loo.rds")
loo_res.N1tc=readRDS("f_ataN_4t2F_1group_voc_incvacc2_tsd4_corr1_loo.rds")
loo_comp.N=loo_compare(loo_res.N4pca,loo_res.N4pc,loo_res.N1tc,loo_res.N1pc,loo_res.N4p,loo_res.N1t)
print(loo_comp.N, simplify = FALSE, digits = 3)

### compare E gene data fits
loo_res.E4pca=readRDS("f_ataE_4p1F_4group_voc_incvacc2_tsd4_corr1_age1_loo.rds")
loo_res.E1pc=readRDS("f_ataE_4p1F_1group_voc_incvacc2_tsd4_corr1_loo.rds")
loo_res.E4p=readRDS("f_ataE_4p1F_4group_voc_incvacc2_tsd4_corr0_loo.rds")
loo_res.E4pc=readRDS("f_ataE_4p1F_4group_voc_incvacc2_tsd4_corr1_loo.rds")
loo_res.E1t=readRDS("f_ataE_4t2F_1group_voc_incvacc2_tsd4_corr0_loo.rds")
loo_res.E1tc=readRDS("f_ataE_4t2F_1group_voc_incvacc2_tsd4_corr1_loo.rds")
loo_comp.E=loo_compare(loo_res.E4pca,loo_res.E4pc,loo_res.E1pc,loo_res.E4p,loo_res.E1tc,loo_res.E1t)
print(loo_comp.E, simplify = FALSE, digits = 3)

### outputs for Table S4
loo_res.N4pca
loo_res.N4pc
loo_res.N1tc
loo_res.N1pc
loo_res.N4p
loo_res.N1t

loo_res.E4pca
loo_res.E4pc
loo_res.E1pc
loo_res.E4p
loo_res.E1tc
loo_res.E1t


### load best fit models

fit.N4pca=readRDS("f_ataN_4p1F_4group_voc_incvacc2_tsd4_corr1_age1.rds")
fit.E4pca=readRDS("f_ataE_4p1F_4group_voc_incvacc2_tsd4_corr1_age1.rds")

fit.N4pc=readRDS("f_ataN_4p1F_4group_voc_incvacc2_tsd4_corr1_age1.rds")
fit.E4pc=readRDS("f_ataE_4p1F_4group_voc_incvacc2_tsd4_corr1_age1.rds")

### ORF1ab gene results
fit=fit.N4pca
# fit=fit.E4pca  # for E-gene results

#extract posterior samples
posterior <- rstan::extract(fit, inc_warmup = FALSE, permuted = TRUE)

### looking at individual-level parameter estimates (not included in paper)

# assemble dataframe of individual level mean posterior estimates and 95% CrI
group_names=c("pre-Alpha","Alpha","Delta-unvacc","Delta-vacc")
PersonID=data_list$subjects$PersonID
voc=1+data_list$subjects$voc
age=data_list$subjects$age

# means
PeakVL=colMeans(posterior$p_ln_v_max)/log(10)
VLg=colMeans(posterior$p_v[,,2])/log(10)
VLd=colMeans(posterior$p_v[,,3])/log(10)

# CrI for Peak VL
tmp=as.data.frame(posterior$p_ln_v_max)/log(10)
PeakVL_r <-  tmp %>%
  summarise(across(.cols = everything(),~quantile(.x, probs = c(0.025,0.975))))
PeakVL_r=as.data.frame(t(PeakVL_r))
names(PeakVL_r)=c("PeakVL_0.025","PeakVL_0.975")
# CrI for growth rate
tmp=as.data.frame(posterior$p_v[,,2])/log(10)
VLg_r <-  tmp %>%
  summarise(across(.cols = everything(),~quantile(.x, probs = c(0.025,0.975))))
# CrI for decline rate
tmp=as.data.frame(posterior$p_v[,,3])/log(10)
VLg_r=as.data.frame(t(VLg_r))
names(VLg_r)=c("VLg_0.025","VLg_0.975")
VLd_r <-  tmp %>%
  summarise(across(.cols = everything(),~quantile(.x, probs = c(0.025,0.975))))
VLd_r=as.data.frame(t(VLd_r))
names(VLd_r)=c("VLd_0.025","VLd_0.975")
# put in single dataframe
results <- as.data.frame(cbind(PersonID,age,voc,PeakVL,PeakVL_r,VLg,VLg_r,VLd,VLd_r))
results$voc_group <- as.factor(group_names[results$voc])
results$voc_group <- factor(results$voc_group,levels=group_names)

# now calculate averages across VOC groups
results.sum=NULL
for(i in 1:4) {
  PeakVL=mean(posterior$v_max_voc[,i])/log(10)
  PeakVL_0.025=quantile(posterior$v_max_voc[,i],probs=0.025)/log(10)
  PeakVL_0.975=quantile(posterior$v_max_voc[,i],probs=0.975)/log(10)
  VLg=mean(posterior$v_a_voc[,i])/log(10)
  VLg_0.025=quantile(posterior$v_a_voc[,i],probs=0.025)/log(10)
  VLg_0.975=quantile(posterior$v_a_voc[,i],probs=0.975)/log(10)
  VLd=mean(posterior$v_b_voc[,i])/log(10)
  VLd_0.025=quantile(posterior$v_b_voc[,i],probs=0.025)/log(10)
  VLd_0.975=quantile(posterior$v_b_voc[,i],probs=0.975)/log(10)
  voc=i
  tmp=as.data.frame(t(c(voc,PeakVL,PeakVL_0.025,PeakVL_0.975,VLg,VLg_0.025,VLg_0.975,VLd,VLd_0.025,VLd_0.975)))
  names(tmp)=c("voc","PeakVL","PeakVL_0.025","PeakVL_0.975","VLg","VLg_0.025","VLg_0.975","VLd","VLd_0.025","VLd_0.975")
  results.sum=rbind(results.sum,tmp)
}
results.sum$voc_group <- as.factor(group_names[results.sum$voc])
results.sum$voc_group <- factor(results.sum$voc_group,levels=group_names)


#plot dependence of peak VL on age

res.f=results[results$PeakVL>1.4,]
res.f$PeakVL=res.f$PeakVL-results.sum$PeakVL[res.f$voc]
res.f$PeakVL_0.025=res.f$PeakVL_0.025-results.sum$PeakVL[res.f$voc]
res.f$PeakVL_0.975=res.f$PeakVL_0.975-results.sum$PeakVL[res.f$voc]
res.f$log_age=log(res.f$age)
res.f$g="regression"
p <- ggplot(res.f, aes(x=log_age, y=PeakVL, ymin=PeakVL_0.025, ymax=PeakVL_0.975, colour=voc_group)) +
  geom_point() +
#  geom_smooth(data=res.f, method='lm', formula= y~x) +
  geom_smooth(data=res.f, method='lm', formula= y~x, aes(x=log_age, y=PeakVL, colour=g)) +
#  geom_pointrange(size = 0.3,position="jitter") + 
  xlab("Age") +ylab("log10(Peak VL)") + # coord_cartesian(ylim=c(6,9.5)) +
  theme_classic() + theme(legend.position="right") + 
  theme(strip.background = element_blank()) 

p
ggsave(paste0("figAge.pdf"), p, width=8, height = 4.5)

cor.test(res.f$PeakVL,res.f$log_age,method = "pearson")

# nice plot of group level and individual level estimates

p1<-ggplot(results, aes(x=1, y=PeakVL, ymin=PeakVL_0.025, ymax=PeakVL_0.975, colour=voc_group)) + 
  facet_wrap(~voc_group) + 
  geom_boxplot(data=results.sum,stat = "identity",aes(lower=PeakVL_0.025,middle=PeakVL,upper=PeakVL_0.975)) +
  geom_pointrange(size = 0.3, position = "jitter") + 
  xlab("") +ylab("log10(Peak VL)") + coord_cartesian(ylim=c(6,9.5)) +
  theme_classic() + theme(legend.position="none") + 
  theme(strip.background = element_blank(), axis.line.x = element_blank()) +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())

p2<-ggplot(results, aes(x=1, y=VLg, ymin=VLg_0.025, ymax=VLg_0.975, colour=voc_group)) + 
  facet_wrap(~voc_group) + 
  geom_boxplot(data=results.sum,stat = "identity",aes(lower=VLg_0.025,middle=VLg,upper=VLg_0.975)) +
  geom_pointrange(size = 0.3, position = "jitter") +
  xlab("") +ylab("VL growth rate (log10 units/day)") + coord_cartesian(ylim=c(0,15)) +
  theme_classic() + theme(legend.position="none") + 
  theme(strip.background = element_blank(), axis.line.x = element_blank()) +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())

p3<-ggplot(results, aes(x=1, y=VLd, ymin=VLd_0.025, ymax=VLd_0.975, colour=voc_group)) + 
  facet_wrap(~voc_group) + 
  geom_boxplot(data=results.sum,stat = "identity",aes(lower=VLd_0.025,middle=VLd,upper=VLd_0.975)) +
  geom_pointrange(size = 0.3, position = "jitter") +
  xlab("") +ylab("VL decline rate (log10 units/day)") + coord_cartesian(ylim=c(0,2.5)) +
  theme_classic() + theme(legend.position="none") + 
  theme(strip.background = element_blank(), axis.line.x = element_blank()) +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())

pl1 = grid.arrange(p1,p2,p3,ncol=3)
pl1

ggsave(paste0("figNew.svg"), pl1, width=8, height = 4.5)


## calculate posterior probabilities of kinetic parameters differing between groups - Tables S9/S10

# .p variables are group level, .mp are within-sample averages

post.v_max.p=array(rep(0,16),dim=c(4,4))
post.v_a.p=array(rep(0,16),dim=c(4,4))
post.v_b.p=array(rep(0,16),dim=c(4,4))

post.v_max.mp=array(rep(0,16),dim=c(4,4))
post.v_a.mp=array(rep(0,16),dim=c(4,4))
post.v_b.mp=array(rep(0,16),dim=c(4,4))


for(i in 1:4)
  for(j in 1:4) {
    tmp=ifelse(posterior$v[,1,j] >= posterior$v[,1,i],1,0)
    post.v_max.p[i,j]=sum(tmp)/length(tmp)
    tmp=ifelse(posterior$v_mean[,2,j] >= posterior$v_mean[,2,i],1,0)
    post.v_a.p[i,j]=sum(tmp)/length(tmp)
    tmp=ifelse(posterior$v_mean[,3,j] >= posterior$v_mean[,3,i],1,0)
    post.v_b.p[i,j]=sum(tmp)/length(tmp)
	
    tmp=ifelse(posterior$v_max_voc[,j] >= posterior$v_max_voc[,i],1,0)
    post.v_max.mp[i,j]=sum(tmp)/length(tmp)
    tmp=ifelse(posterior$v_a_voc[,j] >= posterior$v_a_voc[,i],1,0)
    post.v_a.mp[i,j]=sum(tmp)/length(tmp)
    tmp=ifelse(posterior$v_b_voc[,j] >= posterior$v_b_voc[,i],1,0)
    post.v_b.mp[i,j]=sum(tmp)/length(tmp)
    tmp=ifelse(posterior$t_a_voc[,j] >= posterior$t_a_voc[,i],1,0)
  }
post.v_max.p
post.v_a.p
post.v_b.p
post.v_max.mp
post.v_a.mp
post.v_b.mp

# posterior estimates for Tables S5/S6
# peak VL, growth rate and decline rate have been divided by ln(10) in Tables 3, S5 & S6 to put results on log10 VL scale

# within sample
print(fit, probs=c(0.025, 0.5, 0.975), pars = c("v_max_voc","v_a_voc","v_b_voc"), digits = 3)

# group level (note v[1,j] is peak VL on ln(VL) scale, v_mean[2,j] is growth rate on linear scale, v[3,j] is decline rate])


print(fit, probs=c(0.025, 0.5, 0.975), pars = c("v","v_mean"), digits = 3)

# estimates for tables S7/S8 

print(fit, probs=c(0.025, 0.5, 0.975), pars = c("age_slope","v", "v_sd","Omega","v_s","fp","fp_ct_mean",
                                                "fp_ct_sd"), digits = 3)

