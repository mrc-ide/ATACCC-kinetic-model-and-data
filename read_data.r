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

read_data = function(fn,full_profile,inc_vacc=1,study_id){
  
  n <- length(study_id)
  id_start <- 1
  for(i in 1:n) {
      data0 <- read.csv(fn[i])
      data0 <-  data0 %>% select(PersonID,Vaccinated,VaccInclude,B117Status,Age,AgeGT35,VacVOC,TestDateIndex, 
                                 FirstSwabNotMax,FullFollowup,fullProfile,Symptomatic,FirstSwabHalfMaxVL,CtT1)
      data <- data0[data0$TestDateIndex>=-14 & data0$TestDateIndex<28,]
      if(full_profile==1) 
        data <- data[data$fullProfile==1,]    # removes subjects missing a full VL profile
      else if(full_profile==2) 
        data <- data[data$FirstSwabHalfMaxVL==1,] # removes subjects where first VL measure is >50% of max VL
      else if(full_profile==3) 
        data <- data[data$FirstSwabNotMax==1,] # removes subjects where first VL measure is within 2 of max VL
      else if(full_profile==4) 
        data <- data[(data$FullFollowup==1 | data$fullProfile==1) & data$VaccInclude==1,] # removes subjects without full followup unless they are incident cases

      if(inc_vacc==0)
      {
        data <- data[!is.na(data$Vaccinated) & data$Vaccinated==0 ,]
      }
      M <- nrow(data)
      data$invCT <- 40.0-data$CtT1
      data$PersonID <- data$PersonID+10000*(i-1)
      if(full_profile!= 5 && inc_vacc==2) data$B117Status=data$VacVOC
      subjects <- data %>% group_by(PersonID) %>% summarise(voc=mean(B117Status),age=mean(Age),age_cat=mean(AgeGT35),symp=mean(Symptomatic))
      subjects$symp[is.na(subjects$symp)]="unknown symptoms"
      subjects$symp[subjects$symp==1]="symptomatic"
      subjects$symp[subjects$symp==0]="asymptomatic"
      subjects$study_id <- study_id[i]
      N <- nrow(subjects)
      data_id <- data %>% group_by(PersonID) %>% summarise(PersonID=mean(PersonID),id=mean(PersonID))
      data_id$id <- id_start:(N+id_start-1)
      id_start <- id_start+N
      data <- left_join(data,data_id)
      if(i==1) {
        N_full <- N
        M_full <- M
        subjects_full <- subjects
        data_full <- data
      }
      else {
        N_full <- N_full+N
        M_full <- M_full+M
        subjects_full <- rbind(subjects_full,subjects)
        data_full <- rbind(data_full,data)
      }
  }
  list(N=N_full, M=M_full, NS=n, subjects=subjects_full, data=data_full)
}
