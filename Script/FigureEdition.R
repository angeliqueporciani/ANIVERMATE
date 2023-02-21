# Code for figure edition for the paper entitled : An environmentally conscious One Health approach to control residual malaria transmission.
# Heinrich et al., submited

# Plot Median of survival and ±95%CI (Fig3a)-----

df_med <- readRDS("./Output/df_med.rds")

data.rect <- data.frame(xmin=-Inf, xmax=+Inf, ymin=7, ymax=13)# pour avoir un rectangle gris dans le graph suivant
df_med <- df_med %>% drop_na(median)
df_med$DAI[df_med$DAI==0]<- c(-18)

ggplot() +
  geom_rect(
    data = data.rect,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax
    ),
    colour = "grey",
    alpha = 0.2
  ) +
  geom_point(data = df_med, aes(x = DAI, y = median, colour = strata)) +
  #    geom_segment(data=injection_date_plot, aes(x=DAI,y=position,yend=0,xend=DAI), color='blue', size=0.3, arrow=arrow(length=unit(0.3, "cm")))+
  # scale_y_continuous(expand=c(0,0))+
  # coord_cartesian(ylim = c(0, 35), clip="off")+
  geom_errorbar(
    data = df_med,
    aes(
      x = DAI,
      y = median,
      ymin = lower,
      ymax = upper,
      color = strata
    ),
    alpha = 0.2
  ) +
  ylab("Survival Median Time (Days)") +
  xlab("Days from the start of the experiment") +
  #geom_vline(xintercept = c(60, 90), linetype="dashed", size=0.5)+
  scale_x_continuous(breaks = c(-18, seq(0, 220, 15)), limits = c(-18, 220)) +
  scale_y_continuous(breaks = seq(0, 30, 5), limits = c(0, 30)) +
  theme_light() + labs(colour = "Group") + scale_color_manual(
    values = c("#648FFF", "#FE6100", "#037153"),
    labels = c("Control", "IVOMEC-D®", "IVM-BEPO®")
  ) + theme(
    text = element_text(
      size = 12,
      color = "black",
      family = "Arial"
    ),
    axis.text.x = element_text(angle = 90),
    panel.grid.minor = element_blank()
  )

ggsave("./Figure/MedtimevsDAI.jpeg", width = 20, height = 10, units = "cm")

# Kaplan Meier Curves ----
surv_weight <- readRDS("./Data/surv_weigt.rds")

## Before injection 
surv_T0 <- surv_weight%>% filter(Lot=="L0")%>%droplevels()
surv_cattle <- survfit(Surv(Temps,Statut)~Bovin, data=surv_T0)
ggsurvplot(surv_cattle, data = surv_T0, legend="right", title="Initial survival curve")

## After injection 
surv_weight <- surv_weight%>% mutate(DAI = str_pad(DAI, 3, pad = "0"))  # data managment for having DAI in good order in graph 
surv_weight <- surv_weight %>% dplyr::filter(DAI==c("000", "002"," 007", "014", "021","028"))

surv_DAI <- survfit(Surv(Temps,Statut)~DAI, data=surv_weight) # survival data
ggsurvplot_facet(surv_DAI, data = surv_weight, facet.by = "Traitement", legend.title = "DAI")

surv_trait <- survfit(Surv(Temps,Statut)~Traitement, data=surv_weight) # survival data
ggsurvplot_facet(surv_trait, data = surv_weight, facet.by = "DAI", legend.title = "Treatment", legend.labs = c("IVOMEC-D®", "IVM-BEPO®", "Control"), palette = c("#FE6100","#037153", "#648FFF"))# plot not in good order even releveled factor ... 



