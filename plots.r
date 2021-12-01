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

plot_graphs <- function(data_list, fit, fn, NG, inc_vacc, group_names=c("non-B.1.1.7","B.1.1.7")) {
# only including code to plot fit of model to data here 
  data <- data_list$data
  subjects <- data_list$subjects
  N <- data_list$N
  posterior <- extract(fit, inc_warmup = FALSE, permuted = TRUE)
  pred_ct_samples <- as.data.frame(posterior$pred_ct)
  
  pred_ct_mean <- pred_ct_samples %>%
    summarise(across(.cols = everything(),.fns=mean)) %>% t(.)
  pred_ct_0.5 <- pred_ct_samples %>%
    summarise(across(.cols = everything(),~quantile(.x, probs = 0.5))) %>% t(.) 
  pred_ct_0.025 <- pred_ct_samples %>%
    summarise(across(.cols = everything(),~quantile(.x, probs = 0.025))) %>% t(.) 
  pred_ct_0.975 <- pred_ct_samples %>%
    summarise(across(.cols = everything(),~quantile(.x, probs = 0.975))) %>% t(.) 
  pred_id <- rep(1:N, 1, each=42)
  pred_t <- rep(-14:27, N)
  
  data3 <- data
  
  x=log(10) # all VL plotted on log10 rather than ln scale
  
  data3$study_id = subjects$study_id[data3$id]
  data3$data_scale = data_scale[data3$study_id]
  data3$data_offset = data_offset[data3$study_id]
  data3$vl = (data3$data_offset+data3$invCT/data3$data_scale)/x

  data3$B117Status <- as.factor(group_names[data3$B117Status+1])
  data3$B117Status <- factor(data3$B117Status,levels=group_names)
  
  data4 <- data %>% group_by(id) %>% summarise(B117Status=mean(B117Status),PersonID=mean(PersonID),Symptomatic=mean(Symptomatic))

  B117Status <- rep(data4$B117Status, 1, each=42)
  PersonID <- rep(data4$PersonID, 1, each=42)
  Symptomatic <- rep(data4$Symptomatic, 1, each=42)
  data2 <- as.data.frame(cbind(pred_id,PersonID,Symptomatic,pred_t,B117Status,pred_ct_mean,pred_ct_0.5,pred_ct_0.025,pred_ct_0.975))
  names(data2) <- c("id","PersonID","Symptomatic","t","B117Status","invCT","ct_0.5","ct_0.025","ct_0.975")

  
  data2$study_id=subjects$study_id[data2$id]
  data2$data_scale=data_scale[data2$study_id]
  data2$data_offset=data_offset[data2$study_id]
  data2$vl = (data2$data_offset+data2$ct_0.5/data2$data_scale)/x
  data2$vl_0.025 = (data2$data_offset+data2$ct_0.025/data2$data_scale)/x
  data2$vl_0.975 = (data2$data_offset+data2$ct_0.975/data2$data_scale)/x
  
  data2$B117Status <- as.factor(group_names[data2$B117Status+1])
  data2$B117Status <- factor(data2$B117Status,levels=group_names)
 
  data3$symp[is.na(data3$Symptomatic)]="U"
  data3$symp[data3$Symptomatic==1]="S"
  data3$symp[data3$Symptomatic==0]="A"
  data2$symp[is.na(data2$Symptomatic)]="U"
  data2$symp[data2$Symptomatic==1]="S"
  data2$symp[data2$Symptomatic==0]="A"
  data3$pID2=paste0(data3$PersonID,":",data3$symp)
  data2$pID2=paste0(data2$PersonID,":",data2$symp)

  pl0 <- ggplot(data=data3,aes(x=TestDateIndex,y=vl,colour=B117Status)) + 
    labs(colour="Variant") +
    geom_point(colour="black") +
    geom_ribbon(data=data2,aes(x=t,y=vl,ymin=vl_0.025,ymax=vl_0.975, linetype=NA),alpha=0.5,show.legend = FALSE) + 
    geom_line(data=data2,aes(x=t,y=vl),size=1) + 
    facet_wrap(~pID2) +theme_classic() + theme(legend.position="bottom") + 
    scale_x_continuous("Day relative to peak", limits = c(-14,28), breaks = c(-14,0,14)) +
    scale_y_continuous("log10(viral copies/ml)", limits = c(1,10), breaks = seq(2,10,by=2)) +
    theme(panel.grid.major = element_line(colour="lightgrey",linetype="dashed", size = 0.1),
          strip.background = element_blank(),strip.text.x = element_blank())

  ggsave(paste0(fn, "_plot0.pdf"), pl0, width=7, height = 10)
  ggsave(paste0(fn, "_plot0.svg"), pl0, width=7, height = 10)
  
 
}
