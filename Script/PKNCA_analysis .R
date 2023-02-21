# load package----
library(PKNCA)
library(dplyr)
library(forcats)
library(kableExtra)

# load data ----
PK_dataW <- readRDS("./Data/PKdata.rds")
poids <- readRDS("./Data/poidsbovins.rds")
IVMdungW <- readRDS("./Data/PKdung.rds")

substi<-function(x) {gsub("[,]",".",x) } 

# PKNCA Analysis plasma IVM ----

# Subset for IVOMEC and MDC
df.PK.ID <- PK_dataW%>% filter(IVM_Formulation=="IVOMEC D")%>% mutate_if(is.factor, fct_drop)%>%mutate(HAI=DPI*24)
df.PK.ID_24 <- df.PK.ID[c(1:24),]
df.PK.MDC <- PK_dataW%>% filter(IVM_Formulation=="Medincell", DAI<=211)%>% mutate_if(is.factor, fct_drop)

# Dose per bovin calculation 

#IVOMEC-D
d_dose_ID <- df.PK.ID_24 %>% filter(DPI==0)
poids_D0 <- poids%>%filter(Date_pesee=='2019-04-30')%>%mutate(Poids=as.numeric(substi(Poids)))
d_dose_ID <- d_dose_ID%>% inner_join(poids_D0, by =c("ID_bovin"))
d_dose_ID <- d_dose_ID%>% mutate(Dose=0.4*Poids)

# MDC 
d_dose_MDC <- df.PK.MDC %>% filter(DAI==0)
d_dose_MDC <- d_dose_MDC%>% inner_join(poids_D0, by =c("ID_bovin"))%>%mutate(Poids=as.numeric(substi(Poids)))
d_dose_MDC <- d_dose_MDC%>% mutate(Dose=2.4*Poids)


# Creation of concentration object (for IVOMEC D only for the first month of injection (0-28day))
conc_obj_ID <- PKNCAconc(df.PK.ID_24, Plasmatic_concentration~DAI|Bovin)
conc_obj_MDC <- PKNCAconc(df.PK.MDC, Plasmatic_concentration~DAI|Bovin)

# dose object
dose_obj_ID <-	PKNCAdose(d_dose_ID,	Dose~DAI|Bovin)
dose_obj_MDC <-	PKNCAdose(d_dose_MDC,	Dose~DAI|Bovin)


# PKNCA.set.summary(
#   name=c("auclast", "cmax", "tmax", "half.life", "aucinf.obs"),
#   point=business.geomean,
#   spread = business.geocv,
#   description="geometric mean and CV"
# )

# Parameters estimation for IVOMEC D----

# Interval manual (0-28 days)

intervals_manual_ID <- data.frame(start=0,
                                  end=28,
                                  half.life=TRUE,
                                  cmax=TRUE,
                                  tmax=TRUE,
                                  aucinf.obs=TRUE,
                                  auclast=TRUE)

data_obj_manual <- PKNCAdata(conc_obj_ID, dose_obj_ID,
                             intervals=intervals_manual_ID)

results_obj_manual <- pk.nca(data_obj_manual)
kable(summary(results_obj_manual), "markdown")#Estimated parameters 

# Run this code to have sd as spread estimation for cmax and not CV(%)
# PKNCA.set.summary(
#   name= "cmax",
#   description= "Cmax and geo sd",
#   point= geomean,
#   spread=geosd,
#   rounding = list(signif = 3),
#   reset = FALSE
# )
# results_obj_manual <- pk.nca(data_obj_manual)
# kable(summary(results_obj_manual), "markdown")# On garde ces paramètres là pour la publi.
# 


# Parameters estimation for MDC -----
intervals_manual_MDC <- data.frame(start=0,
                               end=211,
                               cmax=TRUE,
                               tmax=TRUE,
                               half.life=TRUE,
                               aucinf.obs=TRUE,
                               auclast=TRUE)


data_obj_manual_MDC <- PKNCAdata(conc_obj_MDC, dose_obj_MDC,
                             intervals=intervals_manual_MDC)

results_obj_manual_MDC <- pk.nca(data_obj_manual_MDC)
kable(summary(results_obj_manual_MDC), "markdown")# Estimated parameters for MDC 

##############################
# PK parameters for dung ----
##############################

df.PK.MDC_dung <- IVMdungW%>% filter(IVM_Formulation=="Medincell", DAI<=211)%>% mutate_if(is.factor, fct_drop)
df.PK.ivomec_dung <- IVMdungW%>% filter(IVM_Formulation=="IVOMEC D", DAI<30)%>% mutate_if(is.factor, fct_drop)

# faire un df dose avec T0 et poids des bovins 
d.dose.dung.mdc <-  df.PK.MDC_dung%>%slice(4:7)%>%mutate(DAI=rep(0), Dung_concentration=rep(0))%>%
  inner_join(poids_D0, by =c("ID_bovin"))%>%
  mutate(Dose=2.4*Poids)
  
d.dose.dung.id <-  df.PK.ivomec_dung%>%slice(1:4)%>%mutate(DAI=rep(0), Dung_concentration=rep(0))%>%
  inner_join(poids_D0, by =c("ID_bovin"))%>%
  mutate(Dose=0.4*Poids)

# conc object dung -----
conc_obj_ID_dung <- PKNCAconc(df.PK.ivomec_dung, Dung_concentration~DAI|ID_bovin)
conc_obj_MDC_dung <- PKNCAconc(df.PK.MDC_dung, Dung_concentration~DAI|ID_bovin)
#plot(conc_obj_MDC_dung$data$DAI, conc_obj_MDC_dung$data$Dung_concentration)

# dose object dung ----
dose_obj_ID_dung <-	PKNCAdose(d.dose.dung.id,	Dose~DAI|ID_bovin)
dose_obj_MDC_dung <-	PKNCAdose(d.dose.dung.mdc ,	Dose~DAI|ID_bovin)

# Parameters estimation for MDC dung -----
intervals_manual_MDC <- data.frame(start=0,
                                   end=204,
                                   cmax=TRUE,
                                   tmax=TRUE,
                                   half.life=TRUE,
                                   aucinf.obs=TRUE,
                                   auclast=TRUE)

data_obj_manual_MDC_dung <- PKNCAdata(conc_obj_MDC_dung, dose_obj_MDC_dung,
                                 intervals=intervals_manual_MDC)
res_MDC_dung <- pk.nca(data_obj_manual_MDC_dung)
# Estimated parameters for MDC dung 
kable(summary(res_MDC_dung), "markdown")

# Parameters estimation for IVOMEC D dung -----

data_obj_ID <- PKNCAdata(conc_obj_ID_dung, dose_obj_ID_dung ,
                             intervals=intervals_manual_ID)

res_ID_dung <- pk.nca(data_obj_ID)
# Estimated parameters for IVOMEC D dung 
kable(summary(res_ID_dung), "markdown")# On garde ces paramètres là pour la publi. 

