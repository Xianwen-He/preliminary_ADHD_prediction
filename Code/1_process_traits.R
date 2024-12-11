### This script aims to select and process the traits in interest.
### It also process the data.frame for modeling with traits.
### The selection and processing are based on the EDA.
### @Xianwen He, Nov 25, 2025

library(dplyr)

### load full data set
alldata <- readRDS("./Data/alldata.rds")

### select traits in interest
traits.in.interest <- c('id', "PMAT24_A_CR", "ReadEng_AgeAdj", "PicVocab_AgeAdj",
                        "IWRD_TOT", "ProcSpeed_AgeAdj", "VSPLOT_TC", "SCPT_TPRT",
                        "ListSort_AgeAdj", "PicSeq_AgeAdj",
                        "SSAGA_ChildhoodConduct", "SSAGA_Alc_12_Frq_Drk", "SSAGA_TB_Reg_CPD",
                        'SSAGA_Times_Used_Illicits', 'SSAGA_Times_Used_Stimulants',
                        'Height', 'BMI', 'SSAGA_Income', 'SSAGA_Educ', 'age', 'gender',
                        'ASR_Attn_Raw', 'ASR_Totp_Raw', 'DSM_Adh_Raw')
demographic.data <- alldata %>% select(all_of(traits.in.interest))


### process the traits

# whether reported problems
demographic.data$SSAGA_ChildhoodConduct_binary <- factor(demographic.data$SSAGA_ChildhoodConduct > 0, levels=c(F, T))
table(demographic.data$SSAGA_ChildhoodConduct_binary)

# drunk frequency
# the frequency of being drunk
drk_ordered <-  rep(NA, nrow(demographic.data))
drk_ordered[demographic.data$SSAGA_Alc_12_Frq_Drk==1] <- 1  # male-low-freq
drk_ordered[demographic.data$SSAGA_Alc_12_Frq_Drk==2] <- 1  # female-low-freq
drk_ordered[demographic.data$SSAGA_Alc_12_Frq_Drk==3] <- 2  # high-freq
drk_ordered[demographic.data$SSAGA_Alc_12_Frq_Drk==4] <- 0  # never
demographic.data$SSAGA_Alc_12_Frq_Drk_ordered <- factor(drk_ordered, levels=c(0, 1, 2), ordered=TRUE)
# whether being drunk
demographic.data$SSAGA_Alc_12_Frq_Drk_binary <- factor(demographic.data$SSAGA_Alc_12_Frq_Drk<4, levels=c(F, T))

# whether used illicit drugs
demographic.data$SSAGA_Times_Used_Illicits_binary <- factor(demographic.data$SSAGA_Times_Used_Illicits > 0, levels=c(F, T))

# whether used stimulants
demographic.data$SSAGA_Times_Used_Stimulants_binary <- factor(demographic.data$SSAGA_Times_Used_Stimulants > 0, levels=c(F, T))

# process ordinal demographic variables
demographic.data$SSAGA_Income_ordered <- factor(demographic.data$SSAGA_Income, 
                                                levels = 1:8, ordered=T)
demographic.data$SSAGA_Educ_ordered <- factor(demographic.data$SSAGA_Educ,
                                              levels=11:17, ordered=T)
demographic.data$gender_factor <- factor(demographic.data$gender-1,
                                         levels=c(0, 1), ordered=FALSE)

# logarithm for DSM scores
demographic.data$DSM_Adh_Raw_log <- log(demographic.data$DSM_Adh_Raw+1)

# save the data
save(demographic.data, file='./Data/demographic.data.RData')


### process the data for modeling
# this is based on the results of EDA

# the response variable is log(DSM)
traits.log.dat <- demographic.data %>% select('id', "PMAT24_A_CR", "ReadEng_AgeAdj", "PicVocab_AgeAdj",
                                              "IWRD_TOT", "ProcSpeed_AgeAdj", "VSPLOT_TC", "SCPT_TPRT",
                                              "ListSort_AgeAdj", "PicSeq_AgeAdj",
                                              "SSAGA_ChildhoodConduct_binary", "SSAGA_Alc_12_Frq_Drk_binary",
                                              'SSAGA_Times_Used_Illicits_binary', 'SSAGA_Times_Used_Stimulants_binary',
                                              'BMI', 'gender_factor',
                                              'DSM_Adh_Raw_log')
# the response variable is DSM
traits.dat <- demographic.data %>% select('id', "PMAT24_A_CR", "ReadEng_AgeAdj", "PicVocab_AgeAdj",
                                          "IWRD_TOT", "ProcSpeed_AgeAdj", "VSPLOT_TC", "SCPT_TPRT",
                                          "ListSort_AgeAdj", "PicSeq_AgeAdj",
                                          "SSAGA_ChildhoodConduct_binary", "SSAGA_Alc_12_Frq_Drk_binary",
                                          'SSAGA_Times_Used_Illicits_binary', 'SSAGA_Times_Used_Stimulants_binary',
                                          'BMI', 'gender_factor',
                                          'DSM_Adh_Raw')

# remove NAs
traits.dat <- na.omit(traits.dat)
traits.log.dat <- na.omit(traits.log.dat)

save(traits.log.dat, file='./Data/traits.log.dat.RData')
save(traits.dat, file='./Data/traits.dat.RData')

