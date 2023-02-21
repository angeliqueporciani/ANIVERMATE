# Blood feeding rate analysis 
# load packages -----
library(tidyverse)
library(emmeans)
library(kableExtra)
library(ggplot2)
library(glmmTMB)
library(stringr)
library(car)
library(DHARMA)

# source data -----
source("./src/data_loading.R")
source("./src/fun_ana_surv.R")
All_results <- readRDS("./Output/all_results.rds")

Gorgement <- read.csv("./Data/Gorgement.csv", sep=";")


## Data management ----

Gorgement <- Gorgement%>%mutate(Total= Gorgés+Non_gorgés, Taux_gorgement2=Gorgés/Total, Num_Boucle=as.factor(Num_Boucle), Traitement=as.factor(Traitement))
Gorgement$Traitement[Gorgement$Traitement == "IVOMEC"] <- "IVOMEC D"
Gorgement <- Gorgement %>%  mutate_if(is.character,as.factor) %>%
  mutate_if(is.factor,droplevels)

# Effectif --- 

eff <- Gorgement%>%group_by(Lot, Traitement)%>%drop_na()%>%
  dplyr::summarise(n= sum(Total))
kable(eff,"markdown")

eff <- Gorgement%>%group_by(Traitement)%>%drop_na()%>%
  dplyr::summarise(n= sum(Total))

kable(eff,"markdown")

## Binomial model 
glm.binom1 <- glmmTMB(Taux_gorgement2~ Traitement+(1|Lot)+(1|Num_Boucle), data=Gorgement, weights=Total, family=binomial(link = "logit"))
summary(glm.binom1)
Anova(glm.binom1)
plotResiduals(glm.binom1)#model ok 
emmeans(glm.binom1, ~Traitement, type="response")
