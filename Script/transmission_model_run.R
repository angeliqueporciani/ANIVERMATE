# Transmission model run 

# load packages 
library(tidyverse)
library(coxme)
library(gamm4)
library(emmeans)

source("./Script/fun_mod_ANIVM.R")
source("./src/data_loading.R")
PK_dataW <- readRDS("./Data/PKdata.rds")
surv_w_pk <- readRDS("./Data/surv_w_pk.rds")

# Mdc simulation 
## Initialisation of cox and gam model_MDC -----
PK_MDC <- PK_dataW %>%
  filter(IVM_Formulation=="Medincell")%>%# selection que du traitement MDC
  droplevels()%>%
  mutate(day=DAI)# ajout d'une colone day pour coller au nomenclature des fonctions du modèle.

modPK <- mgcv::gam(Plasmatic_concentration ~ s(day)+s(ID_bovin, bs="re"), data = PK_MDC)
plot(modPK)

# emm_pred_pk.mdc <- emmeans(modPK, ~day, at=list(day=seq(1,165,1)),type="response") %>% as.data.frame
# ggplot(data=emm_pred_pk.mdc, aes(x=day, y=emmean, ymin=lower.CL, ymax=upper.CL)) +
# 	geom_line() +
# 	geom_ribbon(alpha=0.5) +
# 	ylim(0,30)


swp_sub_mdc <- subset(surv_w_pk, Lot!="L0"& IVM_Formulation=="Medincell")%>%droplevels()
cox.mdc <- coxme(Surv(Temps,Statut) ~ Plasmatic_concentration + (1|DAI)+(1|Bovin/Gobelet), data=swp_sub_mdc)
modSurv <- cox.mdc 

## define scenarios----- 
l.LLIN <- c(0.2,0.5,0.8)		# LLIN coverage (reported use)
l.pref <- c(0.2,0.5,0.8) # vector's preference for humans
l.ratio <- c(0.5,1,2) # cattle:human ratio
p.c.treated <- c(0,1)

## Recup EMM ----
#modPK <-c(emm_pred_pk.mdc$emmean)# pour IVOMEC ? 
#modPK <- gam.PK.mdc

## Simulations -----
# this could take 10min to run
df_res_mdc <- expand_grid(l.LLIN,l.pref,l.ratio,p.c.treated)
df_res_mdc <- df_res_mdc %>% mutate(SIR = pmap(list(l.LLIN,l.pref,l.ratio,p.c.treated), SIR_IVM)) # la fonction SIR_IVM retourne une liste de plusieurs vecteurs

## Save results 
save(df_res_mdc, file = "./Output/df_res_mdc_3g_211.RData")

# IVOMEC D simulation ----

PK_IVOM <- PK_dataW %>%
  filter(IVM_Formulation=="IVOMEC D")%>%# selection que du traitement MDC
  droplevels()%>%
  mutate(day=DAI)# ajout d'une colone day pour coller au nomenclature des fonctions du modèle.
# #!!!! WARNING ici si day = DPI on considère qu'il n'y a eu qu'une injection. !!!!!!!
# 

# gam.PK.ivomec <- mgcv::gam(Plasmatic_concentration ~ s(day)+s(ID_bovin, bs="re"), data = PK_IVOM, family=gaussian(link = "log"))
# ggemmeans(gam.PK.ivomec,terms=c("day")) %>% plot()
# summary(gam.PK.ivomec )
# # calcul of estimated marginal means
# emm_pred_pk.ivomec <- emmeans(gam.PK.ivomec, ~day, at=list(day=seq(0,200,1)),type="response") %>% as.data.frame
# #plot of EMM
# ggplot(data=emm_pred_pk.ivomec, aes(x=day, y=response, ymin=lower.CL, ymax=upper.CL)) +
# 	geom_line() +
# 	geom_ribbon(alpha=0.5) 

## simulation PK data
library(linpk)
t.obs <- seq(0, 211, 1)
y1 <- pkprofile(t.obs, cl=1.5, vc=11, ka=2, dose=list(t.dose=seq(0,211, 28),amt=1250))
plot(y1, col="red")

swp_sub_ivm <- subset(surv_w_pk, Lot!="L0"& IVM_Formulation=="IVOMEC D")%>%droplevels()
cox.ivm <- coxme(Surv(Temps,Statut) ~ Plasmatic_concentration+ (1|DAI)+(1|Bovin/Gobelet), data=swp_sub_ivm)

### define scenarios----
l.LLIN <- c(0.2,0.5,0.8)		# LLIN coverage (reported use)
l.pref <- c(0.2,0.5,0.8) # vector's preference for humans
l.ratio <- c(0.5,1,2) # cattle:human ratio
p.c.treated <- c(0,1)

# model to use for PK data anc Cox ---- 
modPK <- y1
#modPK <- gam.PK.mdc
modSurv <- cox.ivm

# Simulations Ivomec -----
# this could take 10min to run

df_res_IVO <- expand_grid(l.LLIN,l.pref,l.ratio,p.c.treated)
df_res_IVO <- df_res_IVO %>% mutate(SIR = pmap(list(l.LLIN,l.pref,l.ratio,p.c.treated), SIR_IVM)) # la fonction SIR_IVM retourne une liste de plusieurs vecteurs

save(df_res_IVO, file = "./Output/df_res_ivmD_3g211.RData")



