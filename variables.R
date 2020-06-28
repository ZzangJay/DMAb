#### Variable list: 
### 1. Common variables:
# EMPI, subjecti id;
# EMPI_cum_dmab, injection id; 
# cum_dmab, number of denosumab doses received at baseline;
# day0, time zero;
# week, an indicator of follow-up week; 
# weeksq, square of week;
# censor, a binary indicator of naturally censored (1) or no censored(0); 
# delay, days delayed for this injection; 
# injection, a binary variable indicating whether received a subsequent injection (1) or not (0); 
# injection_l1, lag one week for the variable injection; 
# p_injection, probability of receiving a subsequent injection for each injection-week; 
# artificial_cw_t, artifical censoring weights; 
# stab_cw, natural censoring weights;
# group, assigned group, on time (1), short delay (2), and long delay (3);
# composite_fx_final, a binary endpoint, fracture happened (1) or not (0)
# censor_g1, a binary variable for artifical censoring of on time (1, censored)
# censor_g2, a binary variable for artifical censoring of short dely (1, censored)
# censor_g3, a binary variable for artifical censoring of long delay (1, censored)
### 2. Baseline covariates: 
# age, a continous variable for baseline age;
# sex, a binary variable for gender (male: , female: )
# cci_score, Charlson comorbidity index;
# f4_score, risk of having any osteoporotic (i.e. hip, wrist, shoulder or spine) fracture within the next 10 years;
# b_fracture4, history of any osteoporotic (i.e. hip, wrist, shoulder or spine) fracture at baseline;
# Oral_BP_length, oral bisphosphoante treatment length (year). 
### 3. Time-varying covarites: 
# t_antidepressant, current antidepressants;
# t_anycancer, diagnosis of cancer;
# t_asthmacopd, diagnosis of asthma or chronic obstructive pulmonary disease;
# t_corticosteroids, current corticosteroids;
# t_cvd, diagnosis of cardiovascular diseases;
# t_dementia, diagnosis of dementia;
# t_epilepsy2, diagnosis of epilepsy or current use of anticonvulsants;
# t_falls, history of falls;
# t_parkinsons, diagnosis of Parkinson's diseases;
# t_ra_sle, diagnosis of rheumatoid arthritis or systemic lupus erythematosus;
# t_renal, diagnosis of renal disease;
# t_type2, diagnosis of type 2 diabetes;
