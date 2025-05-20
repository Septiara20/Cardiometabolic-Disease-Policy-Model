library('eha')
library('data.table')
library('flexsurv')
library(ggplot2)
source('SurvivalAnalysis.R')
##########################################################################################
##========================================================================================
## START OF SEMIMARKOV FUNCTION
##========================================================================================
##########################################################################################

##========================================================================================
## THIS FUNCTION BUILDS A SEMI-MARKOV MULTI-STATE MODEL AND RETURNS THE ASSOCIATED STATE
## OCCUPANCY PROBABILITIES (PROBABILITIES OVER TIME OF BEING IN EACH STATE)

## MODIFIBLE ARGUMENTS:
## ntrans       NUMBER OF TRANSITIONS IN THE MODELLING
## ncovs        NUMBER OF COVARIATE PARAMETERS IN THE MODEL (EXCLUDING INTERCEPT) FOR EACH
##              TRANSITION I.E. ONE FOR EACH BINARY AND CONTINUOUS VARIABLE AND K-1 FOR  
##              EACH CATEGORICAL VARIABLE, WHERE K= NO OF CATEGORIES
## covs         VARIABLE NAMES FOR THE COVARIATES
## coveval      VALUE AT WHICH TO EVALUATE EACH COVARIATE
## dist         DISTRIBUTION TO USE FOR EACH TRANSITION. OPTIONS ARE:
##              wei FOR WEIBULL, exp FOR EXPONENTIAL, gom FOR GOMPERTZ,
##              logl FOR LOGLOGISTIC,logn FOR LOGNORMAL AND gam FOR GENERALISED GAMMA.
##              IF FITTING THE SAME MODEL OVER THE OBSERVED AND EXTRAPOLATED PERIOD 
##              DISTRIBUTION WILL SPAN THE WHOLE TIME HORIZON. IF FITTING A DIFFERENT 
##              MODEL OVER THE OBSERVED PERIOD THAN FOR THE EXTRAPOLATION, DISTRIBUTION
##              WILL SPAN THE OBSERVED PERIOD ONLY.
## dist2        ONLY APPLICABLE IF FITTING A DIFFERENT MODEL OVER THE EXTRAPOLATION PERIOD
##              THAN THE OBSERVED PERIOD.DISTRIBUTION TO USE FOR EACH TRANSITION OVER THE
##              EXTRAPOLATION PERIOD. OPTIONS ARE:
##              wei FOR WEIBULL, exp FOR EXPONENTIAL, gom FOR GOMPERTZ,
##              logl FOR LOGLOGISTIC,logn FOR LOGNORMAL AND gam FOR GENERALISED GAMMA.
## timeseq      THE TIME POINTS TO USE FOR PREDICTIONS OVER THE OBSERVED PERIOD OF THE STUDY.
##              THE FIRST ARGUMENT OF seq SHOULD BE THE START TIME, THE SECOND ARGUMENT THE 
##              END TIME AND THE THIRD ARGUMENT THE TIME INCREMENT. 
## timeseq_ext  THE TIME POINTS TO USE FOR PREDICTIONS OVER THE EXTRAPOLATION PERIOD.
##              THE FIRST ARGUMENT OF seq SHOULD BE THE START TIME, THE SECOND ARGUMENT THE 
##              END TIME AND THE THIRD ARGUMENT THE TIME INCREMENT. 
## data         DATASET TO USE FOR MODELLING
## seedno       NUMBER TO USE TO SET THE RANDOM NUMBER GENERATOR SO THAT SIMULATIONS CAN
##              BE REPLICATED EXACTLY. IF NOT REQUIRED SET TO NULL
## M            NUMBER OF SIMULATIONS USED TO CALCULATE PROBABILITIES 
## trans        TRANSITION MATRIX
## predinitial  PREDICT FROM INITIAL STATE? EITHER TRUE OR FALSE
## predfrom     IF NOT PREDICTING FROM INITIAL STATE, NUMBER OF STATE IN WHICH TO START
##              PREDICTION 
##========================================================================================
semiMarkov<-function(ntrans=3, ncovs=c(1,1,2), 
                     covs=rbind("covariate1", "covariate1",c("covariate1", "covariate2")),
                     coveval=rbind(0,0,c(0,1)),  
                     dist=cbind("wei", "wei", "wei"),
                     dist2=cbind(NA, NA, NA),
                           timeseq=seq(0,4,1/12), 
                          timeseq_ext=seq(49/12,15,1/12),
              data=msmcancer, seedno=12345, M=100,
              trans=tmat,predinitial=TRUE, predfrom=2){
  #### set up required lists 
  models<-vector("list", ntrans)
  fmla<-vector("list", ntrans)
  covars<-vector("list", ntrans)
  datasub<-vector("list", ntrans)
  lp<-vector("list", ntrans)
  coeffs<-vector("list", ntrans)
  x<-vector("list", ntrans)
  cumHaz<-vector("list", ntrans) 
  cumHaz_ext<-vector("list", ntrans) 
  kappa<-vector("list", ntrans) 
  gamma<-vector("list", ntrans) 
  mu<-vector("list", ntrans) 
  sigma<-vector("list", ntrans) 
  z<-vector("list", ntrans) 
  u<-vector("list", ntrans) 
  
  #### create the timepoints
  tt2<-timeseq_ext
  for (i in 1:ntrans) {
    
    if (is.na(dist2[i])==TRUE) tt<-c(timeseq,timeseq_ext)
    if (is.na(dist2[i])==FALSE) tt<-timeseq
  
  #### coefficients from modelling of each transition

  covars[[i]]<-covs[i,1:ncovs[i]]
  fmla[[i]]<-as.formula(paste("Surv(time,status)~ ",paste(covars[[i]],
                          collapse= "+"))) 
  datasub[[i]]<-subset(data,trans==i) 
  x[[i]]<-coveval[i,1:ncovs[i]]
    if (dist[i]=="wei") {
      print(paste("Running Weibull model for transition", i))
      models[[i]]<-phreg(fmla[[i]],dist="weibull", data=datasub[[i]])$coeff
      coeffs[[i]]<-models[[i]][1:ncovs[i]] 
      lp[[i]]<-sum(coeffs[[i]]* x[[i]] )
      #### cumulative hazards for each transition
      cumHaz[[i]]<-exp(-exp(models[[i]][ncovs[i]+2])*models[[i]][ncovs[i]+1]+lp[[i]])*
      tt^exp(models[[i]][ncovs[i]+2])    
    }
  if (is.na(dist2[i])==FALSE & dist2[i] =="wei") {
    print(paste("Running Weibull model for transition", i))
    models[[i]]<-phreg(fmla[[i]],dist="weibull", data=datasub[[i]])$coeff
    coeffs[[i]]<-models[[i]][1:ncovs[i]] 
    lp[[i]]<-sum(coeffs[[i]]* x[[i]] )
    #### cumulative hazards for each transition
    cumHaz_ext[[i]]<-exp(-exp(models[[i]][ncovs[i]+2])*models[[i]][ncovs[i]+1]+lp[[i]])*
      tt2^exp(models[[i]][ncovs[i]+2])  
  }
    if (dist[i]=="exp") {
      print(paste("Running exponential model for transition", i))
      models[[i]]<-phreg(fmla[[i]],dist="weibull", shape=1, data=datasub[[i]])$coeff
      coeffs[[i]]<-models[[i]][1:ncovs[i]] 
      lp[[i]]<-sum(coeffs[[i]]* x[[i]] )
      #### cumulative hazards for each transition
      cumHaz[[i]]<-exp(-models[[i]][ncovs[i]+1]+lp[[i]])*tt    
    } 
  if (is.na(dist2[i])==FALSE & dist2[i] =="exp") {
    print(paste("Running exponential model for transition", i))
    models[[i]]<-phreg(fmla[[i]],dist="weibull",shape=1, data=datasub[[i]])$coeff
    coeffs[[i]]<-models[[i]][1:ncovs[i]] 
    lp[[i]]<-sum(coeffs[[i]]* x[[i]] )
    #### cumulative hazards for each transition
    cumHaz_ext[[i]]<-exp(-models[[i]][ncovs[i]+1]+lp[[i]])*tt2    
  }
  
    if (dist[i]=="gom") {
      print(paste("Running Gompertz model for transition", i))
      models[[i]]<-phreg(fmla[[i]],dist="gompertz", param="rate", data=datasub[[i]])$coeff
      coeffs[[i]]<-models[[i]][1:ncovs[i]] 
      lp[[i]]<-sum(coeffs[[i]]* x[[i]] )
      #### cumulative hazards for each transition
      cumHaz[[i]]<-exp(models[[i]][ncovs[i]+2]+lp[[i]])*(1/models[[i]][ncovs[i]+1])*
        (exp(models[[i]][ncovs[i]+1]*tt) -1)   
    }
  if (is.na(dist2[i])==FALSE & dist2[i] =="gom") {
    print(paste("Running Gompertz model for transition", i))
    models[[i]]<-phreg(fmla[[i]],dist="gompertz", param="rate", data=datasub[[i]])$coeff
    coeffs[[i]]<-models[[i]][1:ncovs[i]] 
    lp[[i]]<-sum(coeffs[[i]]* x[[i]] )
    #### cumulative hazards for each transition
    cumHaz_ext[[i]]<-exp(models[[i]][ncovs[i]+2]+lp[[i]])*(1/models[[i]][ncovs[i]+1])*
      (exp(models[[i]][ncovs[i]+1]*tt2) -1)  
  }
    if (dist[i]=="logl") {
      print(paste("Running loglogistic model for transition", i))
      models[[i]]<-aftreg(fmla[[i]],dist="loglogistic",  data=datasub[[i]])$coeff
      coeffs[[i]]<-models[[i]][1:ncovs[i]] 
      lp[[i]]<-sum(coeffs[[i]]* x[[i]] )
      #### cumulative hazards for each transition
      cumHaz[[i]]<- -log(1/(1+(exp(-(models[[i]][ncovs[i]+1]-lp[[i]]))*tt)^
                              (1/(exp(-models[[i]][ncovs[i]+2])))))
    }
  if (is.na(dist2[i])==FALSE & dist2[i] =="logl") {
    print(paste("Running loglogistic model for transition", i))
    models[[i]]<-aftreg(fmla[[i]],dist="loglogistic",  data=datasub[[i]])$coeff
    coeffs[[i]]<-models[[i]][1:ncovs[i]] 
    lp[[i]]<-sum(coeffs[[i]]* x[[i]] )
    #### cumulative hazards for each transition
    cumHaz_ext[[i]]<--log(1/(1+(exp(-(models[[i]][ncovs[i]+1]-lp[[i]]))*tt2)^
                               (1/(exp(-models[[i]][ncovs[i]+2])))))
  }
    if (dist[i]=="logn") {
      print(paste("Running lognormal model for transition", i))
      models[[i]]<-aftreg(fmla[[i]],dist="lognormal",  data=datasub[[i]])$coeff
      coeffs[[i]]<-models[[i]][1:ncovs[i]] 
      lp[[i]]<-sum(coeffs[[i]]* x[[i]] )
      #### cumulative hazards for each transition
      cumHaz[[i]]<- -log(1-pnorm((log(tt)-(models[[i]][ncovs[i]+1]-lp[[i]]))/
                                   (exp(-models[[i]][ncovs[i]+2]))))  
    }
  if (is.na(dist2[i])==FALSE & dist2[i] =="logn") {
    print(paste("Running lognormal model for transition", i))
    models[[i]]<-aftreg(fmla[[i]],dist="lognormal",  data=datasub[[i]])$coeff
    coeffs[[i]]<-models[[i]][1:ncovs[i]] 
    lp[[i]]<-sum(coeffs[[i]]* x[[i]] )
    #### cumulative hazards for each transition
    cumHaz_ext[[i]]<--log(1-pnorm((log(tt2)-(models[[i]][ncovs[i]+1]-lp[[i]]))/
                                    (exp(-models[[i]][ncovs[i]+2])))) 
  }
  
    if (dist[i] == "gamma") {
      print(paste("Running gamma model for transition", i))
      models[[i]] <- flexsurvreg(fmla[[i]], dist = "gamma", data = datasub[[i]])
      shape <- models[[i]]$res[1, 1]
      rate <- models[[i]]$res[2, 1]
      coeffs[[i]] <- models[[i]]$coefficients[1:ncovs[i]]
      lp[[i]] <- sum(coeffs[[i]] * x[[i]])
      scale <- 1 / rate
      cumHaz[[i]] <- -log(1 - pgamma(tt, shape = shape, rate = rate * exp(-lp[[i]])))
    }
  
    if (!is.na(dist2[i]) && dist2[i] == "gamma") {
      print(paste("Running gamma model for transition", i))
      models[[i]] <- flexsurvreg(fmla[[i]], dist = "gamma", data = datasub[[i]])
      shape <- models[[i]]$res[1, 1]
      rate <- models[[i]]$res[2, 1]
      coeffs[[i]] <- models[[i]]$coefficients[1:ncovs[i]]
      lp[[i]] <- sum(coeffs[[i]] * x[[i]])
      cumHaz_ext[[i]] <- -log(1 - pgamma(tt2, shape = shape, rate = rate * exp(-lp[[i]])))
    }

    
    if (dist[i]=="gam") {
      print(paste("Running generalised gamma model for transition", i))
      models[[i]]<-flexsurvreg(fmla[[i]],dist="gengamma",data=datasub[[i]])$res
      kappa[[i]]<- models[[i]][3,1]
      gamma[[i]]<-(abs(kappa[[i]]))^(-2)
      coeffs[[i]]<-models[[i]][4:(ncovs[i]+3),1]
      lp[[i]]<-sum(coeffs[[i]]* x[[i]] )
      mu[[i]]<- models[[i]][1,1]  +lp[[i]] 
      sigma[[i]]<-  models[[i]][2,1]
      z[[i]]<-rep(0,length(tt))
      z[[i]]<- sign(kappa[[i]])*((log(tt)-mu[[i]])/sigma[[i]])
      u[[i]]<-gamma[[i]]*exp((abs(kappa[[i]]))*z[[i]])
      if(kappa[[i]]>0){
        cumHaz[[i]]<--log(1-pgamma(u[[i]],gamma[[i]]))
      }
      if(kappa[[i]]==0){
        cumHaz[[i]]<--log(1-pnorm(z[[i]]))
      }
      if(kappa[[i]]<0){
        cumHaz[[i]]<--log(pgamma(u[[i]],gamma[[i]]))
      }
    }
  
  if (is.na(dist2[i])==FALSE & dist2[i] =="gam") {
    print(paste("Running generalised gamma model for transition", i))
    models[[i]]<-flexsurvreg(fmla[[i]],dist="gengamma",data=datasub[[i]])$res
    kappa[[i]]<- models[[i]][3,1]
    gamma[[i]]<-(abs(kappa[[i]]))^(-2)
    coeffs[[i]]<-models[[i]][4:(ncovs[i]+3),1]
    lp[[i]]<-sum(coeffs[[i]]* x[[i]] )
    mu[[i]]<- models[[i]][1,1]  +lp[[i]] 
    sigma[[i]]<-  models[[i]][2,1]
    z[[i]]<-rep(0,length(tt2))
    z[[i]]<- sign(kappa[[i]])*((log(tt2)-mu[[i]])/sigma[[i]])
    u[[i]]<-gamma[[i]]*exp((abs(kappa[[i]]))*z[[i]])
    if(kappa[[i]]>0){
      cumHaz_ext[[i]]<--log(1-pgamma(u[[i]],gamma[[i]]))
    }
    if(kappa[[i]]==0){
      cumHaz_ext[[i]]<--log(1-pnorm(z[[i]]))
    }
    if(kappa[[i]]<0){
      cumHaz_ext[[i]]<--log(pgamma(u[[i]],gamma[[i]]))
    }
  }
  
      ## Royston-Parmar Spline
    if (dist[i] == "rps1") {
      print(paste("Running RPS model for transition", i))
      models[[i]] <- flexsurvspline(fmla[[i]], data = datasub[[i]], scale = "hazard", k = 1)
      coeffs[[i]] <- models[[i]]$coefficients[1:ncovs[i]]
      lp[[i]] <- sum(coeffs[[i]] * x[[i]])
      surv_pred <- summary(models[[i]], t = tt, newdata = as.data.frame(t(x[[i]])))
      S_t <- sapply(surv_pred, function(x) x$est)
      cumHaz[[i]] <- -log(S_t)
    }

    if (!is.na(dist2[i]) && dist2[i] == "rps1") {
      print(paste("Running RPS model for transition", i))
      models[[i]] <- flexsurvspline(fmla[[i]], data = datasub[[i]], scale = "hazard", k = 1)
      coeffs[[i]] <- models[[i]]$coefficients[1:ncovs[i]]
      lp[[i]] <- sum(coeffs[[i]] * x[[i]])
      surv_pred_ext <- summary(models[[i]], t = tt2, newdata = as.data.frame(t(x[[i]])))
      S_t_ext <- sapply(surv_pred_ext, function(x) x$est)
      cumHaz_ext[[i]] <- -log(S_t_ext)
    }
  
        ## Royston-Parmar Spline
    if (dist[i] == "rps2") {
      print(paste("Running RPS model for transition", i))
      models[[i]] <- flexsurvspline(fmla[[i]], data = datasub[[i]], scale = "hazard", k = 2)
      coeffs[[i]] <- models[[i]]$coefficients[1:ncovs[i]]
      lp[[i]] <- sum(coeffs[[i]] * x[[i]])
      surv_pred <- summary(models[[i]], t = tt, newdata = as.data.frame(t(x[[i]])))
      S_t <- sapply(surv_pred, function(x) x$est)
      cumHaz[[i]] <- -log(S_t)
    }

    if (!is.na(dist2[i]) && dist2[i] == "rps2") {
      print(paste("Running RPS model for transition", i))
      models[[i]] <- flexsurvspline(fmla[[i]], data = datasub[[i]], scale = "hazard", k = 2)
      coeffs[[i]] <- models[[i]]$coefficients[1:ncovs[i]]
      lp[[i]] <- sum(coeffs[[i]] * x[[i]])
      surv_pred_ext <- summary(models[[i]], t = tt2, newdata = as.data.frame(t(x[[i]])))
      S_t_ext <- sapply(surv_pred_ext, function(x) x$est)
      cumHaz_ext[[i]] <- -log(S_t_ext)
    }
  
        ## Royston-Parmar Spline
    if (dist[i] == "rps3") {
      print(paste("Running RPS model for transition", i))
      models[[i]] <- flexsurvspline(fmla[[i]], data = datasub[[i]], scale = "hazard", k = 3)
      coeffs[[i]] <- models[[i]]$coefficients[1:ncovs[i]]
      lp[[i]] <- sum(coeffs[[i]] * x[[i]])
      surv_pred <- summary(models[[i]], t = tt, newdata = as.data.frame(t(x[[i]])))
      S_t <- sapply(surv_pred, function(x) x$est)
      cumHaz[[i]] <- -log(S_t)
    }

    if (!is.na(dist2[i]) && dist2[i] == "rps3") {
      print(paste("Running RPS model for transition", i))
      models[[i]] <- flexsurvspline(fmla[[i]], data = datasub[[i]], scale = "hazard", k = 3)
      coeffs[[i]] <- models[[i]]$coefficients[1:ncovs[i]]
      lp[[i]] <- sum(coeffs[[i]] * x[[i]])
      surv_pred_ext <- summary(models[[i]], t = tt2, newdata = as.data.frame(t(x[[i]])))
      S_t_ext <- sapply(surv_pred_ext, function(x) x$est)
      cumHaz_ext[[i]] <- -log(S_t_ext)
    }
  if (is.na(dist2[i])==FALSE) cumHaz[[i]]<-c(cumHaz[[i]], cumHaz_ext[[i]])
  } 
    Haz<-unlist(cumHaz)
    tt<-c(timeseq,timeseq_ext)
    newtrans<-rep(1:ntrans,each=length(tt))
    time<-rep(tt,ntrans)
    Haz<-cbind(time=as.vector(time),Haz=as.vector(Haz),trans=as.vector(newtrans))
    Haz<-as.data.frame(Haz)
  #### state occupancy probabilities
    set.seed(seedno)
    if (predinitial==TRUE) {
    stateprobs <- mssample(Haz=Haz,trans=trans,tvec=tt,clock="reset", M=M)
    }
    if (predinitial==FALSE) {
      stateprobs<-mssample(Haz=Haz,trans=trans, tvec=tt,clock="reset", M=M,
                       history=list(state=predfrom,time=0,tstate=NULL))  
    }  
return(stateprobs) 
}  
##########################################################################################
##========================================================================================
## END OF SEMIMARKOV FUNCTION
##========================================================================================
##########################################################################################



doMarkov <- function(msdataContinuous) {
  # This section sets up the covariates in the experiment!
  covariatesInExperiment <- c("startingAge", "gender", "ethnicity",
                              "deprivationIndex", "cvdFH", "diabetesFH",
                              "atrialFib", "hypertension", "hyperlipidaemia", "latestSmokingStatus", "alcoholStatus",
                              "bmi", "hdl", "ldl", "triglycerides", 
                              "cholesterol", "glucose", "sbp", "dbp")
  
  # This section just does prep
####################################################################################
  # Define parameters for the semiMarkov function
  markovData <- copy(msdataContinuous)
  
  # Add the starting age - make startingAge non time-dependent.
  markovData <- markovData %>%
  mutate(patid = as.numeric(as.character(patid)))
  
  markovData <- markovData %>%
  left_join(wideTransitionTableContSub %>% select(patid, age) %>% rename(startingAge = age), by = "patid")
  
  ntrans <- 13  # Number of transitions

  # Repeat covariatesInExperiment for each state - required by Markoc function.
  covs <- do.call(rbind, replicate(ntrans, covariatesInExperiment, simplify = FALSE))

  # Automatically calculate ncovs from covs
  ncovs <- rep(length(covariatesInExperiment), ntrans)
  
  deprivationMap <- c("Most privileged" = 1, "Privileged" = 2, "Standard" = 3, "Deprived" = 4, "Most deprived" = 5)
  markovData[, deprivationIndex := names(deprivationMap)[match(deprivationIndex, deprivationMap)]]
  
  alcoholNamesMap <- c("AlcoholConsumptionLevel1" = "Safe alcohol", "AlcoholConsumptionLevel2" = "High alcohol", "AlcoholConsumptionLevel3" = "Alcoholic")
  markovData[, alcoholStatus := alcoholNamesMap[alcoholStatus]]
  
  markovData$gender <- ifelse(markovData$gender == "Female", 0, 1)
################################################################################
  smokingMap <- c("Non smoker" = 0,  "Ex smoker" = 1, "Active smoker" = 2)
  deprivationMap <- c("Most privileged" = 1, "Privileged" = 2, "Standard" = 3, "Deprived" = 4, "Most deprived" = 5)
  alcoholMap <- c("Safe alcohol" = 1, "High alcohol" = 2, "Alcoholic" = 3)
  genderMap <- c("Female" = 0, "Male" = 1)
  ethnicityMap <- c("White" = 0, "Black" = 1, "Asian" = 2, "Other" = 3)
################################################################################
  markovData$cvdFH <- as.integer(markovData$cvdFH)
  markovData$diabetesFH <- as.integer(markovData$diabetesFH)
  markovData$atrialFib <- as.integer(markovData$atrialFib)
  markovData$hypertension <- as.integer(markovData$hypertension)
  markovData$hyperlipidaemia <- as.integer(markovData$hyperlipidaemia)
  markovData$latestSmokingStatus <- smokingMap[markovData$latestSmokingStatus]
  markovData$deprivationIndex <- deprivationMap[markovData$deprivationIndex]
  markovData$alcoholStatus <- alcoholMap[markovData$alcoholStatus]
  markovData$ethnicity <- ethnicityMap[markovData$ethnicity]

###########################################################################################
  
#### IMPORTANT PART: Your values for the experiment, and which states they apply to
  
    #covariatesInExperiment <- c("startingAge", "gender", "ethnicity",
    #                          "deprivationIndex", "cvdFH", "diabetesFH",
    #                          "atrialFib", "hypertension", "hyperlipidaemia", "latestSmokingStatus", "alcoholStatus",
    #                          "bmi", "hdl", "ldl", "triglycerides", 
    #                          "cholesterol", "glucose", "sbp", "dbp")
  
  initialCovariateValues <- c(
    startingAge = 34,
    gender = unname(genderMap["Female"]),       
    ethnicity = unname(ethnicityMap["Asian"]),    
    deprivationIndex = 1,                         
    cvdFH = 0,
    diabetesFH = 1,
    atrialFib = 0,
    hypertension = 0,
    hyperlipidaemia = 0,
    latestSmokingStatus = unname(smokingMap["Non smoker"]),
    alcoholStatus = unname(alcoholMap["Safe alcohol"]),
    bmi = 22,   # Healthy weight range
    hdl = 1.8,  # Ideal HDL >1.6
    ldl = 2.5,  # Ideal LDL <3.0
    triglycerides = 1.2,  # Ideal <1.7
    cholesterol = 4.5,  # Ideal <5.0
    glucose = 4.8,  # Ideal fasting glucose 4.5-5.0
    sbp = 115,  # Ideal SBP <120
    dbp = 75   # Ideal DBP <80
  )

  interventionCovariateValues <- c(
    startingAge = 34,
    gender = unname(genderMap["Female"]),
    ethnicity = unname(ethnicityMap["Asian"]),
    deprivationIndex = 1,
    cvdFH = 0,
    diabetesFH = 1,
    atrialFib = 0,
    hypertension = 0,
    hyperlipidaemia = 0,
    latestSmokingStatus = unname(smokingMap["Non smoker"]),
    alcoholStatus = unname(alcoholMap["Safe alcohol"]),
    bmi = 22,   # Healthy weight range maintained
    hdl = 1.8,  # Ideal HDL >1.6
    ldl = 2.5,  # Ideal LDL <3.0
    triglycerides = 1.2,  # Ideal <1.7
    cholesterol = 4.5,  # Ideal <5.0
    glucose = 4.8,  # Ideal fasting glucose 4.5-5.0
    sbp = 115,  # Ideal SBP <120
    dbp = 75   # Ideal DBP <80
  )
  
  covEval <- rbind(
    initialCovariateValues, # Disease-free -> Diabetes
    initialCovariateValues, # Disease-free -> MI
    initialCovariateValues, # Disease-free -> Stroke
    initialCovariateValues, # Disease-free -> Death
    interventionCovariateValues, # Diabetes -> MI
    interventionCovariateValues, # Diabetes -> Stroke
    interventionCovariateValues, # Diabetes -> Death
    interventionCovariateValues, # MI -> Post-MI
    interventionCovariateValues, # MI -> Death
    interventionCovariateValues, # Post-MI -> Death
    interventionCovariateValues, # Stroke -> Post-Stroke
    interventionCovariateValues, # Stroke -> Death
    interventionCovariateValues # Post-Stroke -> Death
  )

# This section sets up the distributions and times
#####################################################################
  
dist <- cbind(
    "rps3", "rps3", "rps3", "gam", # Disease-free → T2DM
    "rps1", "rps2", "rps3",         # T2DM → Transitions
    "rps3", "rps3", "rps3",         # MI → Post-MI, MI → Death, Post-MI → Death
    "rps3", "rps2", "rps3"          # Stroke → Post-Stroke, Stroke → Death, Post-Stroke → Death
)  #dist2 <- cbind("wei", "wei", "wei", "wei", "wei", "wei", "wei", "wei", "wei", "wei", "wei", "wei", "wei")
  timeSeq <- seq(0, 10950, 1)  

  timeSeqExt = seq(10950, 36500, 1)
  #data <- subset(markovData, time > 0)
  markovData$time[markovData$time == 0] <- 1
  data <- markovData
  tmat <- transitionMatrix
#####################################################################
message("Running model... please wait.")
  # Run the semiMarkov function
  results <- semiMarkov(
    ntrans = ntrans,
    ncovs = ncovs,
    covs = covs,
    coveval = covEval,
    dist = dist,
    #dist2 = dist2,
    timeseq = timeSeq,
    timeseq_ext = timeSeqExt,
    data = data,
    seedno = 133,
    M = 100,
    trans = tmat,
    predinitial = TRUE
  )
  results <- applyLifeTableCorrection(results, initialCovariateValues["startingAge"], initialCovariateValues["gender"])
  return(results)
}

applyLifeTableCorrection <- function(markovOutput, startingAge, gender) {
    setDT(markovOutput)

    # Load life tables
    male_life_table <- fread("MaleLifeTable.csv")
    female_life_table <- fread("FemaleLifeTable.csv")

    # Compute patient age in years
    markovOutput[, age := startingAge + floor(time / 365)]

    # Ensure age is within a valid range for merging
    # Rename probability columns to actual state names
    stateNamesMarkov <- c("Disease-free", "T2DM", "MI", "Stroke", "Post-MI", "Post-Stroke", "Death")
    setnames(markovOutput, old = paste0("pstate", 1:7), new = stateNamesMarkov)

    # Select correct life table based on gender
    life_table <- if (gender == 1) male_life_table else female_life_table

    # Explicitly merge qx using base R merge()
    markovOutput <- merge(markovOutput, life_table[, .(age, qx)], by = "age", all.x = TRUE)
    markovOutput[, qx := ifelse(age > 100, qx[age == 100], qx)]
    # Debug: Check qx values after merging
    print("Summary of qx values after merging:")
    print(summary(markovOutput$qx))

    # ** Apply life table mortality consistently throughout the entire timeline**
    markovOutput[, Death_before := Death]
    
    
# Convert annual qx to daily probability
markovOutput[, scaled_qx := 1 - (1 - qx)^(1 / 365)]

# Apply cumulative product approach for compounding risk
markovOutput[, Death := 1 - cumprod(1 - scaled_qx)]

# Ensure probabilities stay in bounds
    markovOutput[, Death := pmin(Death, 1)]  

    # Adjust non-death states proportionally to keep total = 1**
    markovOutput[, (setdiff(stateNamesMarkov, "Death")) := 
                     lapply(.SD, function(x) x * (1 - Death) / (1 - Death_before)), 
                 .SDcols = setdiff(stateNamesMarkov, "Death")]

    # Debug: Check changes in Death probability
    print("Summary of Death probabilities after life table application:")
    print(summary(markovOutput$Death))

    # Return adjusted Markov output
    return(markovOutput)
}


getMeanAge <- function(msdataContinuous) {
  markovData <- copy(msdataContinuous)
  # Filter the dataset to get the first record for each patid
starting_ages <- markovData[markovData$Tstart == 0 & markovData$from == 1, ]

# Extract the ages for each patid by only taking the first record for each unique patid
# 'starting_ages' now contains the relevant subset; 
# the next step depends on whether 'patid' is uniquely identifying rows.
if (!"patid" %in% colnames(starting_ages)) {
  stop("Column 'patid' is not found in msdataContinuous.")
}

# Calculate the mean starting age
mean_starting_age <- mean(starting_ages$age, na.rm = TRUE)
return(mean_starting_age)
}

makeGraph <- function(markovOutput) {
    simAllcov <- copy(markovOutput)
  setDT(simAllcov)
  
  stateNamesMarkov <- c("Disease-free", "T2DM", "MI", "Stroke", "Post-MI", "Post-Stroke", "Death")

  # Convert time 
  simAllcov[, time := time / 365]

  # Ensure total probability sums correctly
  simAllcov[, total_probability := rowSums(.SD), .SDcols = StateNamesMarkov]
  summary(simAllcov$total_probability)
  head(simAllcov[, .(time, total_probability)])

  # Convert data to long format for ggplot
  simAllcov_long <- melt(simAllcov, id.vars = "time", 
                         measure.vars = StateNamesMarkov, 
                         variable.name = "State", value.name = "Probability")

  # Create the stacked area plot
  ggplot(simAllcov_long, aes(x = time, y = Probability, fill = State)) +
    geom_area(alpha = 0.7) +
    labs(title = "State Occupancy Probabilities Over 100 Years",
         x = "Time (Years)", 
         y = "Probability") +
    scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 10)) + 
    scale_fill_brewer(palette = "Set1") +  
    theme_grey(base_size = 14) 
}

test <- function(msdata, covEvalRow = 4) {
  library(flexsurv)
  library(survival)

  covariatesInExperiment <- c("startingAge", "gender", "ethnicity",
                              "deprivationIndex", "cvdFH", "diabetesFH",
                              "atrialFib", "hypertension", "hyperlipidaemia", "latestSmokingStatus", "alcoholStatus",
                              "bmi", "hdl", "ldl", "triglycerides", 
                              "cholesterol", "glucose", "sbp", "dbp")
  
  cat("Extracting transition 4 data...\n")
  markovData <- copy(msdata)
  data4 <- subset(markovData, trans == 4)
  cat("Rows in data4:", nrow(data4), "\n")

  # Print uniqueness of each covariate
  cat("Checking unique values...\n")
  cov_counts <- sapply(data4[, covariatesInExperiment, drop = FALSE], function(x) length(unique(x)))
  print(cov_counts)

  cat("Fitting flexsurvspline model...\n")
  mod <- tryCatch({
    flexsurvspline(Surv(time, status) ~ ., data = data4, scale = "hazard", k = 3)
  }, error = function(e) {
    cat("⚠️ Model fitting failed: ", conditionMessage(e), "\n")
    return(NULL)
  })

  if (is.null(mod)) return(NULL)

  # Prepare prediction covariates
  x4 <- covEval[covEvalRow, ]
  newdata <- as.data.frame(t(x4))
  colnames(newdata) <- covariatesInExperiment

  # Generate prediction timepoints
  tt <- seq(3032, 11244, by = 100)

  cat("Running prediction...\n")
  pred <- tryCatch({
    summary(mod, t = tt, newdata = newdata)
  }, error = function(e) {
    cat("⚠️ Prediction failed: ", conditionMessage(e), "\n")
    return(NULL)
  })

  if (is.null(pred)) return(NULL)

  # Extract and check survival probabilities
  S_vals <- sapply(pred, function(x) x$est)
  cat("Survival range: ", range(S_vals), "\n")

  if (any(S_vals < 0 | S_vals > 1)) {
    cat("❌ Invalid survival values: ", S_vals[S_vals < 0 | S_vals > 1], "\n")
  } else {
    cat("✅ Survival values are valid.\n")
  }

  return(S_vals)
}
