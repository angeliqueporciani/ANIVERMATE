hIVM_mod <- function(modPK, modSurv, day){	
  if(is.numeric(modPK))
    {
    predPK <- modPK[day]
    if (day<211){	
      param <- c(predPK)
      coef <- coef(modSurv)
      h <- exp(sum(param*coef))				# predict h for a given IVM concentration
    } else {
      h <- 1
  }}else { 
  # extract h (hazard) of mortality predicted by the cox model for a day post IVM injection
  predPK <- emmeans(modPK, ~day, at=list(day=day),type="response") %>% as.data.frame %>% .$emmean				# predict IVM concentration according to PK gam model
  if (predPK > 0){						# exclude PK predictions of negative IVM concentration
    param <- c(predPK)
    coef <- coef(modSurv)
    h <- exp(sum(param*coef))				# predict h for a given IVM concentration
  } else {
    h <- 1
  }} 
  return(h)
}


SIR_IVM <- function(Ch = 0.5, pref.h = 0.5, ratio = 1, p.c.treated = 1){
  
  ### cattle / human pop data and LLIN protection
  #Ch <- 0.5			 				# Net coverage (reported use) 
  pii <- 0.9							# Proportion of exposure to bite during which ITN is in use
  #Nc <- Nc 							# Number of cattle
  #Nh <- Pop <- Nh					# Number of humans
  #ratio <- Nc/Nh
  #ratio <- 1					# ratio Ncattle:Nhumans
  #Nh.u <- Nh * (1 - Ch * pii)		# Mean number of unprotected humans
  #Nh.p <- Nh * Ch * pii				# Mean number of protected humans
  
  # Permanet 2, Kou Valley
  UA <- c(185,90) 	# Unfed Alive (24h) mosquitoes in control hut and intervention hut
  UD <- c(36,123) 	# Unfed Dead (24h) mosquitoes in control hut and intervention hut
  FA <- c(679,93) 	# Fed Alive (24h) mosquitoes in control hut and intervention hut
  FD <- c(8,23) 		# Fed Dead (24h) mosquitoes in control hut and intervention hut
  T <- UA+UD+FA+FD
  
  ### Taking into account the deterence of the LLIN (NEW!!!)
  #pref.h.u <- round( T[1]/(T[1]+T[2]), 2)		# preference for unprotected human (against protected human) = Nb of vector in the control hut / (Nb in the contral hut + Nb in the treated hut)
  #pref.h.p <- 1 - pref.h.u						# preference for protected human
  
  #p.ent.h.u <- dhyper(2, (1-Ch*pii)*Pop, Ch*pii*Pop, 2, log = FALSE) + dhyper(1, (1-Ch*pii)*Pop, Ch*pii*Pop, 2, log = FALSE) * (1 - pref.h.p)		# proba to enter a house with an unprotected people
  #p.ent.h.p <- dhyper(0, (1-Ch*pii)*Pop, Ch*pii*Pop, 2, log = FALSE) + dhyper(1, (1-Ch*pii)*Pop, Ch*pii*Pop, 2, log = FALSE) * pref.h.p			# proba to enter a house with a protected people
  # giving a mosquito randomly released in the village, it have the choice between two humans (either protected vs unprotected, unprot. vs unprot., or prot. vs prot.),
  # what is the chance that it choose an unprotected or protected human ?
  #p.ent.h.p <- Ch*pii			# simplification, proba to encounter a protected human (if anthropophiliq)
  #p.ent.h.u <- 1-p.ent.h.p	# simplification, proba to encounter an unprotected human (if anthropophiliq)
  
  
  
  μ_h.u <- (FD[1]+UD[1]) / T[1] # proba of death when entering a hut with an unprotected host
  μ_h.p <- (FD[2]+UD[2]) / T[2] # proba of death when entering a hut with a protected host
  
  #μ_h <- μ_h.u * p.ent.h.u + μ_h.p * p.ent.h.p # proba of death when feeding on humans
  #μ_h <- (μ_h.p - μ_h.u)/(1-μ_h.u) * p.ent.h.p # proba of death due to the insecticide when feeding on humans (i.e corrected mortality)
  
  
  # HBI calculation
  #feed.h.u.hbi <- (FD[1]+FA[1])/T[1]		# Successful feeding probability when entering a hut with an unprotected human (proportion of Fed in the control hut, for HBI calculation)
  #feed.h.p.hbi <- (FD[2]+FA[2])/T[2]		# Successful feeding probability when entering a hut with a protected human (proportion of Fed in the treated hut, for HBI calculation)
  feed.h.u.hbi <- FA[1]/(FA[1]+UA[1])		# Successful feeding probability when surviving in a hut with an unprotected human (proportion of Fed alive among alive mosquitoes in the control hut, for HBI calculation)
  feed.h.p.hbi <- FA[2]/(FA[1]+UA[1])		# Successful feeding probability when surviving in hut with a protected human (proportion of Fed alive among alive mosquitoes in the treated hut, for HBI calculation)
  
  #f.h.u.hbi <- p.ent.h.u * feed.h.u.hbi		# probability of feeding on an unprotected human while searching, for HBI calculation
  #f.h.p.hbi <- p.ent.h.p * feed.h.p.hbi		# probability of feeding on a protected human while searching, for HBI calculation
  #f.h.hbi <- 	f.h.u.hbi + f.h.p.hbi				# probability of feeding while searching (when feeding on human), for HBI calculation
  #f.c.hbi <- feed.h.u.hbi									# probability of feeding while searching (when feeding on cattle), for HBI calculation
  
  ### To find the HBI (human blood index) from cattle/human pop. data and olfactometer data
  #pref.h <- 0.5						# preference for human (against cattle) = proportion of vectors that choose the human odor in a dual-choice olfactometer (Lefevre et al.)
  pref.c <- (1 - pref.h)			# preference for cattle
  #odds.pref <- pref.c / pref.h			# odds of choosing a cattle (based on olfactometer data)
  
  
  
  #ph.u <- Nh.u / (Nh+Nc) 				# proportion of hosts that are unprotected human
  ph <- 1/(1+ratio)						# proportion of host that are humans
  ph.u <- ph * (1 - Ch * pii)				# proportion of hosts that are unprotected human
  #ph.p <- Nh.p / (Nh+Nc)					# proportion of hosts that are protected human
  ph.p <- ph * Ch * pii					# proportion of hosts that are protected human
  #pc	<- Nc / (Nh+Nc)						# proportion of hosts that are cattle
  pc <- 1 - ph							# proportion of hosts that are cattle
  Ph <- (ph.u * feed.h.u.hbi + ph.p * feed.h.p.hbi) * pref.h		# proba that a vector feed on a human
  Pc <- pc * feed.h.u.hbi * pref.c 								# proba that a vector feed on a cattle
  HBI <- Ph / (Ph + Pc) 					# deduced HBI (post-feeding HBI; after taking into account the post-feeding mortality)
  
  # mortality probability due to insecticide
  μ_h <- ph.p * pref.h * (μ_h.p - μ_h.u)/(1-μ_h.u) # proba that a vector die due to LLIN while searching for a blood meal
  
  
  # Slater et al 2014 modified
  
  day1 <- 200													# time of treatment
  dur.eff <- 211											# duration of ivermectin efficacy
  t <- day1 + dur.eff											# duration of obs (days)
  
  S_v <- vector(mode = "numeric", length = t)
  E_v <- vector(mode = "numeric", length = t)
  I_v <- vector(mode = "numeric", length = t)
  S_v_ivm <- vector(mode = "numeric", length = t)
  E_v_ivm <- vector(mode = "numeric", length = t)
  I_v_ivm <- vector(mode = "numeric", length = t)
  
  ### Pf transmission parameters
  #HBI <- 0.5
  prev <- 0.5													# Pf prevalence in the human population
  k <- 0.1														# infectiousness: probability that a vector become infectious while taking a blood meal on an infectious host (included gemotocyte to sporozoite proba)
  p.infectiousbite <- k * prev * HBI  # probability that a vector become infectious (exposed and then infectious) while taking a bloodmeal
  #g <- 2															# Duration of the gonotrophiq cycle (= mean time between two consecutive bloodmeal)
  g <- 3
  br <- 1/g														# biting rate
  n <- 11															# duration of extrinsic incubation = length of sporogony (in days)
  
  
  ### ivermectin treatment parameters
  #p.h.treated <- 0											# proportion of human treated with ivermectin
  #p.c.treated <- 1								# proportion of cow treated with ivermectin
  #p.ive.h <- HBI * p.h.treated 					# probability that a vector take a bloodmeal on an ivermectin treated human
  p.ive.c <- (1-HBI) * p.c.treated				# probability that a vector take IVM on cattle while taking a bloodmeal 
  
  
  ###
  Em <- 300										# number of daily emergences
  S_v[1] <- Em
  E_v[1] <- 0
  I_v[1] <- 0
  
  s <- day1											# time of first host IVM ingestion
  l <- dur.eff										# duration of IVM efficacy
  
  
  μ_v0 <- 0.1										# baseline daily mortality rate (bloodmeal without IVM)
  μ_v_d <- rep(μ_v0,t)							# vector of mortality rates for vectors biting on a given day on IVM treated cattle (before and after injection)
  for (i in (s+1):(s+l)){							# fill μ_v_d (Need a constant survival rate over mosquito ages)
    # h <- hIVM_mod(modPK, modSurv,i-s,form)
    h <- hIVM_mod(modPK, modSurv,i-s) # predicted mortality
    if (h*μ_v0 > μ_v0){	# if predicted mortality is higher than baseline mortality
      μ_v_d[i] <- h*μ_v0	# fill μ_v_d with predicted mortality (hazard x μ_v0)
    } else {μ_v_d[i] <- μ_v0}					# 
  }
  
  move <- function(t, to, d){
    if ((t >= to) & (t < (to + d))){
      return (1)
    } else {
      return (0)
    }
  }
  
  for (i in 1:t){
    S_v[i+1] <- eval(S_v[i] 											# vector of Pf susceptibles no IVM
                     + Em 															                # new emergences
                     - μ_v0 * S_v[i] 												            # baseline (daily) mortality for non IVM anopheles
                     - br *  S_v[i] * (                                 # blood searching mosquitoes
                       μ_h                                                   # death due to LLIN
                        + (1 - μ_h) * p.infectiousbite                        # become exposed
                        + (1 - μ_h) * p.ive.c * move(i,s,l)                   # take IVM
                     ))
    
    S_v_ivm[i+1] <- eval(S_v_ivm[i] 									# vector of Pf susceptibles with IVM
                         + move(i,s,l) * (1 - μ_h) * br * p.ive.c * S_v[i] 				# new Pf susceptible with IVM (from S_v)
                         - (1 - μ_h) * br * p.infectiousbite * S_v_ivm[i] 				# exposed to Pf (move to E_v_ivm)
                         - μ_v_d[i] * S_v_ivm[i] 										# mortality for IVM anopheles
                         - br * μ_h * S_v_ivm[i])									# mortality due to LLIN
    
    E_v[i+1] <- eval(E_v[i] 											# vector of Pf exposed no IVM
                     + (1 - μ_h) * br * p.infectiousbite * S_v[i] 					# new exposed to Pf (from S_v)
                     - (1/n) * E_v[i] 												# become infectious (move to I_v)
                     - μ_v0 * E_v[i] 												# mortality for non IVM anopheles
                     - move(i,s,l) * (1 - μ_h) * br * p.ive.c * E_v[i] 							# take IVM (move to E_v_ivm)
                     - br * μ_h * E_v[i])										# mortality due to LLIN
    
    
    E_v_ivm[i+1] <- eval(E_v_ivm[i] 									# vector of Pf exposed with IVM
                         + move(i,s,l) * (1 - μ_h) * br * p.ive.c * E_v[i] 							# new Pf exposed with IVM (from E_v)
                         + (1 - μ_h) * br * p.infectiousbite * S_v_ivm[i]  				# new Pf exposed with IVM (from S_v_ivm)
                         - (1/n) * E_v_ivm[i] 											# became infectious (move to I_v_ivm)
                         - μ_v_d[i] * E_v_ivm[i] 										# mortality for IVM anopheles
                         - br * μ_h * E_v_ivm[i])									# mortality due to LLIN
    
    I_v[i+1] <- eval(I_v[i] 											# vector of Pf infectious no IVM
                     + (1/n) * E_v[i] 												# new Pf infectious (from E_v)
                     - μ_v0 * I_v[i] 												# mortality for non IVM anopheles
                     - move(i,s,l) * (1 - μ_h) * br * p.ive.c * I_v[i] 							# take IVM (move to I_v_ivm)
                     - br * μ_h * I_v[i])										# mortality due to LLIN
    
    I_v_ivm[i+1] <- eval(I_v_ivm[i] 									# vector of Pf infectious with IVM
                         + move(i,s,l) * (1 - μ_h) * br * p.ive.c * I_v[i] 							# new Pf infectious with IVM (from I_v)
                         + (1/n) * E_v_ivm[i] 											# new Pf infectious with IVM (from E_v_ivm)
                         - μ_v_d[i] * I_v_ivm[i] 										# mortality for IVM anopheles
                         - br * μ_h * I_v_ivm[i])									# mortality due to LLIN
  }
  
  pop <- S_v + S_v_ivm + E_v + E_v_ivm + I_v + I_v_ivm
  inf <- I_v + I_v_ivm
  SR <- inf / pop
  return(list(pop=pop, inf=inf, SR=SR, HBI=HBI, μ_v_d=μ_v_d))
}

