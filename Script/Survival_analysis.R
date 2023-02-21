# load packages -----
library(dplyr)
library(purrr)
library(coxme)
library(survminer)
library(survival)
# load data ---- 

surv_weight <- readRDS("./Data/surv_weigt.rds")
surv <-readRDS("./Data/surv.rds")

# Survival effectif table --- 
eff_table <- surv_weight%>%
  group_by(Traitement)%>%
  dplyr::summarise(n= n())

# Survival analysis for each lot ------

# function that performs survival analyses for each lot
ana_surv <- function(df,...){
  df.mod <- df %>% mutate_if(is.factor, fct_drop)
  if (nrow(df.mod) > 0){
    cox <-  try(coxme(Surv(Temps,Statut) ~ Traitement + (1|ID_bovin/Gobelet), data=df.mod)) # mixed-effect cox model
    emm <- try(summary(emmeans(cox, revpairwise  ~ Traitement), type = "response", infer = TRUE)) # marginal means and multiple comparisons
    surv <- try(survfit(Surv(Temps,Statut)~Traitement, data=df.mod)) # survival data
    plot <- try(ggsurvplot(surv, data = df.mod, legend="right")) # curves
    med <- try(surv_median(surv)) # median survival time
    return(list(cox=cox,emm=emm,surv=surv,plot=plot,med=med)) # return a list with all the created objects
  } else {
    return(NULL)
  }
  
}

## Application to data of each lot ----

tib_survi <- surv_weight%>%group_by(Lot)%>%nest()
tib_survi <- tib_survi%>%mutate(res_surv=map(data, ana_surv))
saveRDS(tib_survi, "./Output/tib_survie.rds")

##############################
# Dose response analysis ----- 

# Here are presented results for model selected in script DRCmodel
# Plot follows model fit 

## Calcul of mortality proportion in control group per lot 

surv_control<- surv %>% dplyr::filter(Traitement=="Temoin")%>%
  mutate(m7 = ifelse(Temps <= 7, 1, 0), m10 = ifelse(Temps <= 10, 1, 0), m13=ifelse(Temps <= 13, 1, 0))%>% 
  mutate_if(is.factor, fct_drop)

## Estimating basal mortality proportion 
n_mort_7d<- surv_control%>%filter(m7=="1")%>%summarise(n1_7d=n())#438/5384
n_mort_10d<- surv_control%>%filter(m10=="1")%>%summarise(n1_10d=n())# 743/5384
n_mort_13d<- surv_control%>%filter(m13=="1")%>%summarise(n1_13d=n())#1474/5384

# subset of L0 (before injection) and temoin group 
swp_sub <- subset(surv_w_pk, Lot!="L0"& Traitement!="Temoin")%>%droplevels()

# creation of a binomial variable, 1=dead at 7 day, 0=alive at 7 day post feeding
swp_sub_binom <- swp_sub %>% 
  mutate(m7 = ifelse(Temps <= 7, 1, 0), m10 = ifelse(Temps <= 10, 1, 0), m13=ifelse(Temps <= 13, 1, 0))%>% 
  mutate_if(is.factor, fct_drop)

# check if pb with data 
L5_IVOMEC <- swp_sub_binom%>%filter(Lot=='L5'& Traitement=="IVOMEC_D")
n17IVOMEC <- swp_sub_binom%>%group_by(Plasmatic_concentration, Lot)%>%filter(m7==1&Traitement=="IVOMEC_D")%>%summarise(n1_7=n())%>%ungroup()
n17mdc <- swp_sub_binom%>%group_by(Plasmatic_concentration, Lot)%>%filter(m7==1&Traitement=="Medincell_F3_2")%>%summarise(n1_7=n())%>%ungroup()

# pb of mortality for L12 and L8 for IVomec, no pb of higher mortality in mdc or control groups
# I retrieve this lot for IVOMEC 

IVOMECtot <- surv_w_pk %>% filter(Traitement=="IVOMEC_D"& Lot!="L8"& Lot!="L12")
MDCtot <- surv_w_pk %>% filter(Traitement=="Medincell_F3_2")
data_drc_tot <- rbind(IVOMECtot, MDCtot)
data_drc_tot <- data_drc_tot %>% 
  mutate(m7 = ifelse(Temps <= 7, 1, 0), m10 = ifelse(Temps <= 10, 1, 0), m13=ifelse(Temps <= 13, 1, 0))%>% 
  mutate_if(is.factor, fct_drop)


## MODEL APPLICATION for 7, 10 and 14 day effect (mortality) post feeding ------

# 7 days 
drc.mod7.LL4.2 <-drm(m7~Plasmatic_concentration, IVM_Formulation, data=data_drc_tot, 
                     fct=LL.4(fixed=c(NA, 0.08, 1, NA),
                              names = c("Slope", "Lower Limit", "Upper Limit", "ED50")), type="binomial") 
summary(drc.mod7.LL4.2)

# 7 days with L0 (for comparison purpose)
# drc.mod7.LL4.2.tot <-drm(m7~Plasmatic_concentration, IVM_Formulation, data=data_drc_tot, 
#                      fct=LL.4(fixed=c(NA, 0.08, 1, NA),
#                               names = c("Slope", "Lower Limit", "Upper Limit", "ED50")), type="binomial") 
# summary(drc.mod7.LL4.2.tot)
# not useful as we have already calculated the baseline mortality and fixed it in the model 

# PLOT 
plot(drc.mod7.LL4.2, log="x", legendPos = c(9, 1), 
     xlab= "Ivermectin plasma concentration (log) ng/mL", 
     ylab= "7 days cumulative mortality", 
     col=c("#FE6100", "#037153"), 
     legendText = c("IVOMEC-D®", "IVM-BEPO®"))


# Comparison of ED50 
ED(drc.mod7.LL4.2, c(50, 90), interval = "delta")# Estimated ED50 
compParm(drc.mod7.LL4.2, "ED50", "/")# no difference detected for ED50 
compParm(drc.mod7.LL4.2, "Slope", "/")# Signif difference between both formulation for slope, impact on time to reach efficacy 
EDcomp(drc.mod7.LL4.2, c(90, 90))# comparaison of LC90 between both formulation 

# 10 days 
drc.mod10.LL4.2 <-drm(m10~Plasmatic_concentration, IVM_Formulation, data=data_drc_tot, 
                      fct=LL.4(fixed=c(NA, 0.13, 1, NA),
                               names = c("Slope", "Lower Limit", "Upper Limit", "ED50")), type="binomial") 
summary(drc.mod10.LL4.2)

# PLOT 
plot(drc.mod10.LL4.2, log="x", legendPos = c(9, 1), 
     xlab= "Ivermectin plasma concentration (log) ng/mL", 
     ylab= "10 days cumulative mortality",
     col=c("#FE6100", "#037153"), 
     legendText = c("IVOMEC-D®", "IVM-BEPO®"))

# Comparison of ED50 
ED(drc.mod10.LL4.2, c(50, 95), interval = "delta")# Estimated ED50 
compParm(drc.mod10.LL4.2, "ED50", "/")# no difference detected for ED50 
compParm(drc.mod10.LL4.2, "Slope", "/")# Signif difference between both formulation for slope, impact on time to reach efficacy 
EDcomp(drc.mod10.LL4.2, c(90, 90))# comparaison of LC90 between both formulation 

# 13 days 
drc.mod13.LL4.2 <-drm(m13~Plasmatic_concentration, IVM_Formulation, data=data_drc_tot, 
                      fct=LL.4(fixed=c(NA, 0.27, 1, NA),
                               names = c("Slope", "Lower Limit", "Upper Limit", "ED50")), type="binomial") 
summary(drc.mod13.LL4.2)

# PLOT 
plot(drc.mod13.LL4.2, log="x", legendPos = c(9, 1), 
     xlab= "Ivermectin plasma concentration (log) ng/mL", 
     ylab= "13 days cumulative mortality",
     col=c("#FE6100", "#037153"), 
     legendText = c("IVOMEC-D®", "IVM-BEPO®"))

# Comparison of ED50 
ED(drc.mod13.LL4.2, c(50, 95), interval = "delta")# Estimated ED50 
compParm(drc.mod13.LL4.2, "ED50", "/")# no difference detected for ED50 
compParm(drc.mod13.LL4.2, "Slope", "/")# Signif difference between both formulation for slope, impact on time to reach efficacy 
EDcomp(drc.mod13.LL4.2, c(90, 90))# comparaison of LC90 between both formulation 



