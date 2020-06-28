# Appendix Statistical Codes 

#############################################################################
## Implementation of “cloning, censoring, and weighting” with R (selected) ##
#############################################################################

# R version 3.5.2 with the following packages
library(dplyr) # version 0.8.3
library(sandwich) # version 2.4-0
library(lmtest) # version 0.9-36
# Data setup, “data” is the long format dataset with a maximum 26 weeks follow-up for each injection. 
# Variables: EMPI_cum_dmab, ID; week, an indicator of follow-up week; delay, days delayed for this injection; injection, received a subsequent injection (1) or not (0); injection_l1, lag one week for the variable injection; p_injection, probability of receiving a subsequent injection for each injection-week; 
# For the data structure, please see https://github.com/houchenlyu/DMAb/variables.
data <- readRDS("data.RDS") 

########################################################################
## Step 1: Estimate Probability for Artificial Censoring Weights      ##
########################################################################
# In this study, artificial censoring is determined by when the person received the subsequent injection; we estimated the probability of being artificially censored from the following two models in the original dataset: The first model was restricted to the first 4 weeks of follow-up, and the second model to subsequent weeks. 
# first model was restricted to the first 4 weeks of follow-up;
data_restrict_first <- data %>% filter(injection_l1==0 & week<=4)
dFit <- glm(injection ~ age + sex + cci_score + f4_score + b_fracture4 + t_antidepressant+ t_anycancer + t_asthmacopd + t_corticosteroids + t_cvd + t_dementia + t_epilepsy2 + t_falls + t_parkinsons + t_ra_sle + t_renal+ t_type2 
+ cum_dmab + week + weeksq, data= data_restrict_first, family=binomial(link = "logit"))
data_restrict_first$p_injection <- predict(dFit, newdata=data_restrict_first, type='response')
# second model was restricted to the left weeks of follow-up;
data_restrict_second <- data %>% filter(injection_l1==0 & week>=5)
dFit <- glm(injection ~ age + sex + cci_score + f4_score + b_fracture4 + t_antidepressant+ t_anycancer + t_asthmacopd + t_corticosteroids + t_cvd + t_dementia + t_epilepsy2 + t_falls + t_parkinsons + t_ra_sle + t_renal+ t_type2 
+ cum_dmab + week + weeksq, data= data_restrict_second, family=binomial(link = "logit"))
data_restrict_second$p_injection <- predict(dFit, newdata=data_restrict_second, type='response')
# combind the probability
data_restrict <- bind_rows(data_restrict_first, data_restrict_second)
# merge the probability to main dataset
data <- left_join(data, data_restrict[,c("EMPI_cum_dmab","week","p_injection")], by= c("EMPI_cum_dmab"="EMPI_cum_dmab", "week"="week"))
# set the probability to 1 for weeks that were already on injection
data <- data %>% mutate(p_injection = ifelse(is.na(p_injection), 1, p_injection))
########################################################################
## Step 2: Performing cloning and censoring                           ##
########################################################################
## Sep 2.1, cloning: make 3 copies of this dataset
# copy 1 indexed by g1 (on time); censor, a indicator of natural censoring (1 or 0);
copy_1 <- data %>% mutate(group=1) %>% filter(censor ==0 )
# copy 2 indexed by g2 (short delay).
copy_2 <- data %>% mutate(group=2) %>% filter(censor ==0 )
# copy 3 indexed by g3 (long delay).
copy_3 <- data %>% mutate(group=3) %>% filter(censor ==0 )
## Step 2.2, censoring. In each g-specific copy, creating an indicator of artificial censoring: an individual will be artificially censored in the first interval when his data becomes inconsistent with g. In this copy, the time when he is artificially censored should be the last week. Tips: the variable "delay": how many days have been delayed for each subsequent injection; censor_g1, censor_g2, and censor_g3 are indicators of artificial censoring for g1, g2, and g3.
## for g1: rules, patients received their subsequent injections within 4 weeks; 
temp_1 <- copy_1 %>% group_by(EMPI_cum_dmab) %>% filter(delay <= 28) %>% mutate(censor_g1 = 0)
temp_2 <- copy_1 %>% group_by(EMPI_cum_dmab) %>% filter(delay > 28) %>% filter (week<= 4) %>% mutate(censor_g1= ifelse(week==4, 1, 0))
copy_1 <- bind_rows(temp_1, temp_2)
## for g2: rules, patients received their subsequent injections between 5 week or 16 weeks;
temp_1 <- copy_2 %>% group_by(EMPI_cum_dmab) %>% filter(delay > 28 & delay <= 112) %>% mutate(censor_g2 = 0)
# There are two patterns: <= 28 days or > 112 days.
temp_2 <- copy_2 %>% group_by(EMPI_cum_dmab) %>% filter(delay <= 28) %>% filter (day_t_start <= delay) %>% arrange(EMPI_cum_dmab, week) %>% mutate(censor_g2= ifelse(row_number()==n(), 1, 0))
temp_3 <- copy_2 %>% group_by(EMPI_cum_dmab) %>% filter(delay >112) %>% filter (week<=16 ) %>% mutate (censor_g2= ifelse(week==16, 1, 0))
copy_2 <- bind_rows(temp_1, temp_2, temp_3)
## for g3: rules, patients received their subsequent injections after 16 weeks;
temp_1 <- copy_3 %>% group_by(EMPI_cum_dmab) %>% filter(delay > 112) %>% mutate(censor_g3 = 0)
temp_2 <- copy_3 %>% group_by(EMPI_cum_dmab) %>% filter(delay <= 112) %>% filter (day_t_start <= delay) %>% arrange(EMPI_cum_dmab, week) %>% mutate(censor_g3= ifelse(row_number()==n(), 1, 0))
copy_3 <- bind_rows(temp_1, temp_2)
#########################################################################
## Step 3: Construct Inverse Probability Artificial Censoring Weights  ##
#########################################################################
# in each copy, estimate the probability of NOT being artificially censored at each time k, conditional on (i): past treatment and confounder history, week (k), and (ii): previous survival and previously free of all censoring (artificial and natural). Once a patients received treatment, the future weights would be 1 based on deterministic knowledge.
# for g1
copy_1 <- copy_1 %>% 
mutate(p_censor_g1= ifelse(injection==1, 0, p_injection), p_uncensor_g1= 1- p_censor_g1)
copy_1$g1_cw_denCont <- with(copy_1, censor_g1*p_censor_g1 + (1-censor_g1)*(1-p_censor_g1))
copy_1$g1_cw_deno <- with(copy_1, ave(g1_cw_denCont, EMPI_cum_dmab, FUN=cumprod))
copy_1$artificial_cw_t <- 1/ copy_1$g1_cw_deno
# for g2
copy_2 <- copy_2 %>% 
  mutate(p_censor_g2 = ifelse(injection==1, 0,ifelse(week>=5 & week <=16, 0, p_injection)), p_uncensor_g2 = ifelse (week == 16 & !is.na(p_injection), p_injection, 1- p_censor_g2))
copy_2 <- copy_2 %>% arrange(EMPI_cum_dmab, week) 
copy_2$g2_cw_deno <- with(copy_2, ave(p_uncensor_g2, EMPI_cum_dmab, FUN=cumprod))
copy_2$artificial_cw_t <- 1/ copy_2$g2_cw_deno
# for g3
copy_3 <- copy_3 %>% 
  mutate(p_censor_g3 = ifelse(injection==1 | week>=17, 0, p_injection), p_uncensor_g3= 1 - p_censor_g3)
copy_3 <- copy_3 %>% arrange(EMPI_cum_dmab, week) 
copy_3$g3_cw_deno <- with(copy_3, ave(p_uncensor_g3, EMPI_cum_dmab, FUN=cumprod))
copy_3$artificial_cw_t <- 1/ copy_3$g3_cw_deno
# Combine copies to create the final dataset
combined_final_data <- bind_rows(copy_1, copy_2, copy_3)
# Truncate the weights at 99.5%;
threshold <- quantile(combined_final_data$artificial_cw_t, 0.995) 
combined_final_data$artificial_cw_t[combined_final_data$artificial_cw_t> threshold] <- threshold
# Create the final primary weights: stab_cw is natural censoring weights (see the full analytical codes)
combined_final_data$artificial_cw_t <- combined_final_data$artificial_cw_t * combined_final_data$stab_cw 
########################################################################
## Step 4: Primary analysis with final weights                       ##
########################################################################
## Full adjusted model, estimate the HRs and 95% CI;
res.cox <- glm(composite_fx_final ~ factor(group) + week+ weeksq + age + sex + cci_score + major_fx_history +Oral_BP_length + f4_score + cum_dmab , data = combined_final_data, weights = artificial_cw_t, family = "quasibinomial")
# to compute conservative 95% CIs with robust standard error for HR estimates. 
exp(coef(res.cox)[['factor(group)2']]) # to get Hazard Ratio
exp(coefci(res.cox, parm = 'factor(group)2', vcov=vcovHC(res.cox, type="HC1")))
exp(coef(res.cox)[['factor(group)3']]) # to get Hazard Ratio
exp(coefci(res.cox, parm = 'factor(group)3', vcov=vcovHC(res.cox, type="HC1")))
# p for trend
res.cox <- glm(composite_fx_final ~ group + week+ weeksq + age + sex + cci_score + major_fx_history +Oral_BP_length + f4_score + cum_dmab , data = combined_final_data, weights = artificial_cw_t, family = "quasibinomial")
coeftest(res.cox, vcov=vcovHC(res.cox, type="HC1")) 

