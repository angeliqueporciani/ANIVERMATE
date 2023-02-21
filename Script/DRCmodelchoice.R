# Script for DRC model choice (ANIVERMATE project)
# Load packages ----
library(drc)
library(dplyr)

# Load data ---
surv_weight <- readRDS("./Data/surv_weigt.rds")
surv <-readRDS("./Data/surv.rds")

# IVOMECtot <- surv_w_pk %>% filter(Traitement=="IVOMEC_D"& Lot!="L8"& Lot!="L12")
MDCtot <- surv_w_pk %>% filter(Traitement=="Medincell_F3_2")
data_drc_tot <- rbind(IVOMECtot, MDCtot)
data_drc_tot <- data_drc_tot %>% 
  mutate(m7 = ifelse(Temps <= 7, 1, 0), m10 = ifelse(Temps <= 10, 1, 0), m13=ifelse(Temps <= 13, 1, 0))%>% 
  mutate_if(is.factor, fct_drop)

# MODEL CHOICE -----

# Fitting of different dose response curve model : log-logistic with adjustement of basal mortality (observed) and weibull 

# simpler model : all parameter allow to vary for both formulation 
drc.mod1.LL4.1 <-drm(m7~Plasmatic_concentration, IVM_Formulation, data=data_drc_tot, fct=LL.4(), 
                     type="binomial") 
summary(drc.mod1.LL4.1)
plot(drc.mod1.LL4.1)

# Fixed model :  lower and upper limit fixed, slope and ED50 allow to vary between formulation
drc.mod1.LL4.2 <-drm(m7~Plasmatic_concentration, IVM_Formulation, data=data_drc_tot, 
                     fct=LL.4(fixed=c(NA, 0.08, 1, NA),
                              names = c("Slope", "Lower Limit", "Upper Limit", "ED50")), type="binomial") 
summary(drc.mod1.LL4.2)
plot(drc.mod1.LL4.2)
drc.mod1.LL4.2$fct
## Different param only for intercept and final 
drc.mod1.LL4.3 <-drm(m7~Plasmatic_concentration, IVM_Formulation, data=data_drc_tot, 
                     fct=LL.4(fixed=c(NA, NA, 1, NA),
                              names = c("Slope", "Lower Limit", "Upper Limit", "ED50")), type="binomial",
                     pmodels=list(~1,~IVM_Formulation-1,~1, ~IVM_Formulation)) 
summary(drc.mod1.LL4.3)
plot(drc.mod1.LL4.3)

# comparison of best model between one with parameter slope and ED 50 to vary between formulation 
AIC(drc.mod1.LL4.1,drc.mod1.LL4.2,drc.mod1.LL4.3)# drc.mod1.LL4.2 with lower AIC 

compParm(drc.mod1.LL4.2, "ED50", "/")

# LL3 model with upper limit fixed at 1 and lower limit fixed at 0.2 (need to be adapted to observed mortality for same day in control)
#same as LL4.2
drc.mod1.LL3u.1 <-drm(m7~Plasmatic_concentration, IVM_Formulation, data=data_drc, 
                      fct=LL.3u(fixed=c(NA, 0.08, NA),
                                names = c("Slope", "Lower Limit", "ED50")), type="binomial")

summary(drc.mod1.LL3u.1)

drc.mod1.LL3u.2 <-drm(m7~Plasmatic_concentration, IVM_Formulation, data=data_drc, 
                      fct=LL.3u(fixed=c(NA, NA, NA),
                                names = c("Slope", "Lower Limit", "ED50")), type="binomial")
summary(drc.mod1.LL3u.2)

# selection of different model type (function mselect compare different type of model)
mselect(drc.mod1.LL4.2, list(LL.2(),LL.3u(),LL.3u(fixed=c(NA, 0.08, NA)), W1.4(fixed=c(NA, 0.08, 1, NA))), 
        linreg=F, icfct=AIC, sorted="IC")

# Results of mselect
# logLik       IC Lack of fit
# LL.4  -1691.623 3391.246   0.9995684
# LL.3u -1691.623 3391.246   0.9995684
# LL.3u -1690.789 3393.577   0.9978017
# W1.4  -1706.185 3420.371   0.9999983
# LL.2  -2072.347 4152.694   0.9966621

# validation of the model 
plot(fitted(drc.mod1.LL4.2), residuals(drc.mod1.LL4.2))# not so bad 
modelFit(drc.mod1.LL4.2)# seems ok 

# plot of the best model selected 
plot(drc.mod1.LL4.2, log="x", col = c(2,5), legendPos = c(4, 1), 
     xlab= "Ivermectin Plasmatic concentration (log) ng/mL", ylab= "7 days mortality probability")

# using some tools for model validation (residuals checking and goodness of fit test plus graphics observation), 
# we can conclude that the model selected LL4.2 is good enough to be used for estimate ED50 and Slope. 
