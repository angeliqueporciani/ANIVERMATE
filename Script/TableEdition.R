# load package----
# load data ---- 
tib_surv <- readRDS("./Output/tib_survie.rds")

# Median survival time -----

surv_med <- list()
for (i in 1:nrow(tib_surv)){
  if(is.null(tib_surv$res_surv[[i]]$med)){
    surv_med[[i]] <- NULL }
  else{
    surv_med[[i]] <- tib_surv$res_surv[[i]]$med
    surv_med[[i]]$Lot <- rep(tib_surv$Lot[[i]])
    surv_med[[i]]$DAI <- rep(unique(tib_surv$data[[i]]$DAI))
  }}

df_med <- data.frame(do.call(rbind, surv_med))

#df_med <- subset(df_med, Lot!="L0")
df_med$strata <- factor(df_med$strata ,levels = c("Traitement=Temoin", "Traitement=IVOMEC_D", "Traitement=Medincell_F3_2"))
#write_csv(HR_df, "./Output/df_med.csv")
#saveRDS(df_med, "./Output/df_med.rds")

# Hazard ratio between treatment group for each lot ---- 
HR_list <- list()
for (i in 1:nrow(tib_surv))
{ if(is.null(tib_surv$res_surv[[i]]$emm)){HR_list[[i]] <- NULL}
  else{
    HR_list[[i]] <- as.data.frame(tib_surv$res_surv[[i]]$emm$contrasts)%>%
      mutate(Lot=tib_surv$Lot[[i]], DAI=unique(tib_surv$data[[i]]$DAI))}
}
HR_df <- do.call(rbind, HR_list)

library(stringr)
round_perso <- function(x){
  round(x, digits = 3)
}
HR_df <- HR_df%>%mutate(Comp = as.character(contrast),DAI2 = as.numeric(as.character(DAI)))%>%
  dplyr::filter(str_detect(Comp, ".*Temoin"))%>%mutate_if(is.numeric, round_perso)


HR_df$DAI[HR_df$DAI==-18]<-0
#write_csv(HR_df, "./Output/HR_tab.csv")
#saveRDS(HR_df, "./Output/HR_df.rds")
