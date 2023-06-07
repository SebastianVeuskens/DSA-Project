# Sebastian Veuskens, Jose Gabriel Escarraman, Desmond Reynolds, Max Schlake


# 1) SETUP ####
rm(list=ls())

# Load libraries
library(survival)
library(ggplot2)
library(rstudioapi)
library(DescTools)
library(survminer)
library(mfp)
library(flexsurv)

# Load data
setwd(file_path <- dirname(rstudioapi::getSourceEditorContext()$path))

heart_data <- read.csv("../data/S1Data.csv")

heart_data$Gender <- as.factor(heart_data$Gender)
heart_data$Smoking <- as.factor(heart_data$Smoking)
heart_data$Diabetes <- as.factor(heart_data$Diabetes)
heart_data$BP <- as.factor(heart_data$BP)
heart_data$Anaemia <- as.factor(heart_data$Anaemia)

attach(heart_data) 

heart_data$Age_group <- cut(Age, quantile(Age), include.lowest = TRUE)
# SOLVED: Should we continue with all the variables below or just leave them out 
# because it is too much? -> Continue with all of them 
heart_data$Ejection.Fraction_group <- cut(Ejection.Fraction, quantile(Ejection.Fraction), include.lowest = TRUE)
heart_data$Sodium_group <- cut(Sodium, quantile(Sodium), include.lowest = TRUE)
heart_data$Creatinine_group <- cut(Creatinine, quantile(Creatinine), include.lowest = TRUE)

heart_data$Pletelets_group <- cut(Pletelets, quantile(Pletelets), include.lowest = TRUE)
heart_data$CPK_group <- cut(CPK, quantile(CPK), include.lowest = TRUE)
# SOLVED: 'include.lowest = TRUE' was added as otherwise we get some NAs if the 
# value of the variable is equal to the lower boundary of the interval. To 
# compare, see 'summary(heart_data)' after changing the data without this option

detach(heart_data)
attach(heart_data)


# 2) EXPLORATORY ANALYSIS #### 
## 2.1) Chi-Squared tests and histograms ####
PercTable(Event, Gender, rfrq="001", margins = c(1,2))
chisq.test(table(Event, Gender), correct = FALSE)

PercTable(Event, Smoking, rfrq="001", margins = c(1,2))
chisq.test(table(Event, Smoking), correct = FALSE)

PercTable(Event, Diabetes, rfrq="001", margins = c(1,2))
chisq.test(table(Event, Diabetes), correct = FALSE)

PercTable(Event, Diabetes, rfrq="001", margins = c(1,2))
chisq.test(table(Event, BP), correct = FALSE)

PercTable(Event, Anaemia, rfrq="001", margins = c(1,2)) 
chisq.test(table(Event, Anaemia), correct = FALSE)
# Conclusion: No association between any of the variables and the patient's 
# death; lowest p_value for BP (blood pressure)

hist(Age, breaks = 15, freq = F)
# Conclusion: Age is not normally distributed

rate <- sum(Event==1) / sum(TIME) * 100 
rate 

summary(heart_data)
####### COMMENT 
# - Event: 32% of the patients died; 
# - TIME: follow-up time is between 4 and 285 days;
# - Gender: 194 men (= 1) / 105 women (= 0)
# - Smoking: 96 yes (= 1) / 203 no (= 0)
# - Diabetes: 125 yes (= 1) / 174 no (= 0)
# - BP: 105 yes (= 1) / 194 no (= 0)
# - Anaemia: patients with haematocrit < 36 (= 1) / >= 36 (= 0)
# - Age: average = 60.83
# - Ejection.Fraction: average = 38.08
# - Sodium: average = 136.6
# - Creatinine: average = 1.394
# - Pletelets: average = 263358
# - CPK: average = 581.8

## 2.2) Kaplan-Meier estimator ####
surv_obj <- Surv(TIME, Event) 
result_simple <- survfit(surv_obj~ 1, conf.type = "log-log") 
names(result_simple) 
summary(result_simple)

plot(result_simple)

jpeg("../images/surv_curve_overall.jpeg", quality = 75)
ggsurvplot(result_simple, xlab = "Days", ylab = "Overall survival probability",
            data = heart_data, legend = "none")
dev.off()

## 2.3) Additional estimates ####
# Cumulative hazard 
simple_cumhaz <- round(result_simple$cumhaz, 4)
plot(simple_cumhaz, type="l")

# Fleming-Harrington (Nelson-Aalen) estimator 
result_simple_na <- survfit(Surv(TIME, Event) ~ 1, conf.type = "log-log") 
plot(result_simple_na) 


# 3.) UNIVARIATE ANALYSIS #### 
## 3.1) Survival curves ####
jpeg("../images/surv_curve_Gender.jpeg", quality = 75)
ggsurvplot(survfit(surv_obj~Gender,), conf.int = T, xlab = "Days", ylab = "Overall survival probability", data = heart_data)
dev.off()

jpeg("../images/surv_curve_Smoking.jpeg", quality = 75)
ggsurvplot(survfit(surv_obj~Smoking,), conf.int = T, xlab = "Days", ylab = "Overall survival probability", data = heart_data)
dev.off()

ggsurvplot(survfit(surv_obj~Diabetes,), conf.int = T, xlab = "Days", ylab = "Overall survival probability", data = heart_data)

jpeg("../images/surv_curve_BP.jpeg", quality = 75)
ggsurvplot(survfit(surv_obj~BP,), conf.int = T, xlab = "Days", ylab = "Overall survival probability", data = heart_data)
dev.off()

jpeg("../images/surv_curve_Anaemia.jpeg", quality = 75)
ggsurvplot(survfit(surv_obj~Anaemia,), conf.int = T, xlab = "Days", ylab = "Overall survival probability", data = heart_data)
dev.off()
####### COMMENT
# Already here we can see that BP and Anaemia have well separated survival 
# curves, which means that there might be a univariate effect of these variables 
# on the survival probability

ggsurvplot(survfit(surv_obj~Age_group,), conf.int = T, xlab = "Days", ylab = "Overall survival probability", data = heart_data)
ggsurvplot(survfit(surv_obj~Ejection.Fraction_group,), conf.int = T, xlab = "Days", ylab = "Overall survival probability", data = heart_data) # seems time dependent 
ggsurvplot(survfit(surv_obj~Sodium_group,), conf.int = T, xlab = "Days", ylab = "Overall survival probability", data = heart_data) # seems partially time dependent 

jpeg("../images/surv_curve_Creatinine_group.jpeg", quality = 75)
ggsurvplot(survfit(surv_obj~Creatinine_group,), conf.int = T, xlab = "Days", ylab = "Overall survival probability", data = heart_data, legend = 'top')
dev.off()

jpeg("../images/surv_curve_Pletelets.jpeg", quality = 75)
ggsurvplot(survfit(surv_obj~Pletelets_group,), conf.int = T, xlab = "Days", ylab = "Overall survival probability", data = heart_data) # seems not important
dev.off()

ggsurvplot(survfit(surv_obj~CPK_group,), conf.int = T, xlab = "Days", ylab = "Overall survival probability", data = heart_data) # seems not important 
####### COMMENT
# No clear univariate effect for Pletelets and CPK. For all the other variables, 
# we have at least one level with a clearly separated survival curve. For 
# Creatinine, we could consider merging the two groups in the middle, 
# ([0.9, 1.1] and [1.1, 1.4]).
# SOLVED: We leave Creatinine with 4 levels.

## 3.2) Cloglog visual for PH-assumption ####  
# SOLVED: Legend location (removed legend box with 'bty' argument)
for(i in c(3:7, (ncol(heart_data) - 5):ncol(heart_data))){
    temp_name = colnames(heart_data)[i]
    file_path = paste("../images/cloglog_transform_", temp_name, ".jpeg", sep = '')
    jpeg(file_path, quality = 75)
    plot(survfit(surv_obj~heart_data[,i]), fun="cloglog", yscale=-1, col=1:nlevels(heart_data[,i]), xlab="Days", ylab="Estimated -log(-log S(t))")
    legend("topright", bty = "n", title= colnames(heart_data)[i], legend = levels(heart_data[,i]), col=1:nlevels(heart_data[,i]), lty=1)
    dev.off()
}

####### COMMENT
# PH may
# - stand for: BP, Anaemia, Creatinine_group
# - not stand for: Gender, Smoking, Diabetes, Age_group, Ejection.Fraction_group
# Sodium_group, Pletelets_group, CPK_group

## 3.3) Logrank and Wilcoxon tests ####
# SOLVED: Add grouped variables in both tests -> Or we leave it as it is

####### COMMENT
# We should consider the logrank test only for those variables for which PH 
# holds (BP, Anaemia, Creatinine_group). For the rest, the Wilcoxon is preferred

# Logrank test 
survdiff(surv_obj~Gender)
survdiff(surv_obj~Smoking)
survdiff(surv_obj~Diabetes)
survdiff(surv_obj~BP) # significant
survdiff(surv_obj~Anaemia) # almost significant
survdiff(surv_obj~Age_group) # significant
survdiff(surv_obj~Ejection.Fraction_group) # significant
survdiff(surv_obj~Sodium_group) # significant
survdiff(surv_obj~Creatinine_group) # significant
survdiff(surv_obj~Pletelets_group)
survdiff(surv_obj~CPK_group)

# Wilcoxon test (Prentice correction)
survdiff(surv_obj~Gender, rho = 1)
survdiff(surv_obj~Smoking, rho = 1)
survdiff(surv_obj~Diabetes, rho = 1)
survdiff(surv_obj~BP, rho = 1) # significant
survdiff(surv_obj~Anaemia, rho = 1) # almost significant
survdiff(surv_obj~Age_group, rho = 1) # significant
survdiff(surv_obj~Ejection.Fraction_group, rho = 1) # significant
survdiff(surv_obj~Sodium_group, rho = 1) # significant
survdiff(surv_obj~Creatinine_group, rho = 1) # significant
survdiff(surv_obj~Pletelets_group, rho = 1)
survdiff(surv_obj~CPK_group, rho = 1)

## 3.4) KM versus COX fit #### 
# SOLVED: We could add some more variables with the same analysis here -> But I 
# think it is enough this way 
jpeg("../images/KM_vs_Cox_BP.jpeg", quality = 75)
plot(survfit(surv_obj ~ BP, data=heart_data), col=1:nlevels(BP)) # KM plot
lines(survfit(coxph(surv_obj ~ BP, data=heart_data), newdata=data.frame(BP=c(0, 1)),se.fit=F), col=1:nlevels(BP), lty=2) # Cox predicted
legend("bottomright", title= "BP status", legend = levels(BP), col=1:nlevels(BP), lty=1)
dev.off()
####### COMMENT
# The Cox model seems to fit the KM estimate quite well, no systematic over- or 
# underestimation is visible (except maybe for BP status 1 at the end)

jpeg("../images/KM_vs_Cox_Anaemia.jpeg", quality = 75)
plot(survfit(surv_obj ~ Anaemia, data=heart_data), col=1:nlevels(Anaemia)) # KM plot
lines(survfit(coxph(surv_obj ~ Anaemia, data=heart_data), newdata=data.frame(Anaemia=c(0, 1)),se.fit=F), col=1:nlevels(Anaemia), lty=2) # Cox predicted
legend("bottomright", title= "Anaemia status", legend = levels(Anaemia), col=1:nlevels(Anaemia), lty=1)
dev.off()
####### COMMENT
# The Cox model seems to fit the KM estimate quite well, no systematic over- or 
# underestimation is visible

jpeg("../images/KM_vs_Cox_Creatinine_group.jpeg", quality = 75)
plot(survfit(surv_obj ~ Creatinine_group, data=heart_data), col=1:nlevels(Creatinine_group)) # KM plot
lines(survfit(coxph(surv_obj ~ Creatinine_group, data=heart_data), newdata=data.frame(Creatinine_group=levels(heart_data$Creatinine_group),se.fit=F)), col=1:nlevels(Creatinine_group), lty=2) # Cox predicted
legend("bottomleft", title= "Creatinine_group status", legend = levels(Creatinine_group), col=1:nlevels(Creatinine_group), lty=1)
dev.off()
####### COMMENT
# The Cox model seems to fit the KM estimate well for the first 3 groups.
# However, there is under- and overstimation for the last (blue) group

jpeg("../images/KM_vs_Cox_CPK_group.jpeg", quality = 75)
plot(survfit(surv_obj ~ CPK_group, data=heart_data), col=1:nlevels(CPK_group)) # KM plot
lines(survfit(coxph(surv_obj ~ CPK_group, data=heart_data), newdata=data.frame(CPK_group=levels(heart_data$CPK_group),se.fit=F)), col=1:nlevels(CPK_group), lty=2) # Cox predicted
legend("bottomleft", title= "CPK_group status", legend = levels(CPK_group), col=1:nlevels(CPK_group), lty=1)
dev.off()

jpeg("../images/KM_vs_Cox_Pletelets_group.jpeg", quality = 75)
plot(survfit(surv_obj ~ Pletelets_group, data=heart_data), col=1:nlevels(Pletelets_group)) # KM plot
lines(survfit(coxph(surv_obj ~ Pletelets_group, data=heart_data), newdata=data.frame(Pletelets_group=levels(heart_data$Pletelets_group),se.fit=F)), col=1:nlevels(Pletelets_group), lty=2) # Cox predicted
legend("bottomleft", title= "Pletelets_group status", legend = levels(Pletelets_group), col=1:nlevels(Pletelets_group), lty=1)
dev.off()

jpeg("../images/KM_vs_Cox_Ejection.Fraction_group.jpeg", quality = 75)
plot(survfit(surv_obj ~ Ejection.Fraction_group, data=heart_data), col=1:nlevels(Ejection.Fraction_group)) # KM plot
lines(survfit(coxph(surv_obj ~ Ejection.Fraction_group, data=heart_data), newdata=data.frame(Ejection.Fraction_group=levels(heart_data$Ejection.Fraction_group),se.fit=F)), col=1:nlevels(Ejection.Fraction_group), lty=2) # Cox predicted
legend("bottomleft", title= "Ejection.Fraction_group status", legend = levels(Ejection.Fraction_group), col=1:nlevels(Ejection.Fraction_group), lty=1)
dev.off()

# 4.) MULTIVARIATE COX MODELS #### 
## 4.1) Base model: cox_all ####
cox_all <- coxph((surv_obj ~ Age + Ejection.Fraction + Sodium + 
                    Creatinine + Pletelets + CPK + Gender + 
                    Smoking + Diabetes + BP + Anaemia), method = "breslow") # breslow is Wilcoxon test
cox_all 

jpeg("../images/ggforest_hazard_ratio.jpeg", quality = 100)
ggforest(cox_all, data = heart_data)
dev.off()

####### COMMENT
# Significant risk factors: Age, Creatinine, BP, Anaemia
# Significant protective factors: Ejection.Fraction
# Borderline cases: Sodium, CPK 

jpeg("../images/cph_model.jpeg", quality = 75)
plot(survfit(cox_all), main = "cph model", xlab="Days", )
lines(result_simple, col="green")
legend("bottomleft", legend = c("Cox", "KM"), col = c("black", "green"), lty = c(1, 1))
dev.off()

## 4.2) Reduced model: cox_reduced ####
update(cox_all, .~. - Pletelets)
update(cox_all, .~. - Pletelets - Smoking)
update(cox_all, .~. - Pletelets - Smoking - Diabetes)
update(cox_all, .~. - Pletelets - Smoking - Diabetes - Gender)

## 4.3) Testing the PH assumption for the reduced model ####
cox_reduced <- update(cox_all, .~. - Pletelets - Smoking - Diabetes - Gender)

# Checking: H0 = PH assumption holds for all covariates
# SOLVED: Should we do the PH test on the cox_all or the cox_reduced data set 
# (makes sense on all, but analysis outline says reduced) -> We can choose, he 
# recommended both
test_ph_reduced1 <- cox.zph(cox_reduced, transform="identity", terms = F) 
test_ph_reduced1
test_ph_reduced2 <- cox.zph(cox_reduced, transform="rank", terms = F) 
test_ph_reduced2
test_ph_reduced3 <- cox.zph(cox_reduced, transform="log", terms = F) 
test_ph_reduced3
test_ph_reduced4 <- cox.zph(cox_reduced, transform="km", terms = F) 
test_ph_reduced4
####### COMMENT
# Globally, the PH assumption stands in all tests. Only Ejection.Fraction shows
# inconsistent behaviour for the PH assumption under 'identity' and 'km' 
# transformations

for (i in 1:length(test_ph_reduced1)) {
  temp_name = rownames(test_ph_reduced1$table)[i]
  file_path = paste("../images/schoenfeld_res_identity_", temp_name, ".jpeg", sep = '')
  jpeg(file_path, quality = 75)  
  plot(test_ph_reduced1[i])
  dev.off()
}
####### COMMENT
# Age declines at the end
# Ejection.Fraction steadily declines (little)
# Sodium rises and then declines 
# Creatinine rising
# CPK rising then falling 
# BP shaped as wave 
# Anaemia slump at beginning, otherwise straight 


## 4.4) Cox model with time varying coefficients  ####
# SOLVED: Which variables should we include here? -> All, but only 
# Ejection.Fraction shows changing behaviour over time 

# Time dependent model
cox_time = coxph(surv_obj ~ Age + Sodium + Creatinine + Pletelets + CPK + Gender
                 + Ejection.Fraction + Smoking + Diabetes + BP + Anaemia 
                 + tt(Ejection.Fraction), method = "breslow", data = heart_data,
                 tt = function(x, t, ...) x * log(t))

cox_time

# Reduce time-dependent model
update(cox_time, .~. - Smoking)
update(cox_time, .~. - Smoking - Pletelets)
update(cox_time, .~. - Smoking - Pletelets - Diabetes)
update(cox_time, .~. - Smoking - Pletelets - Diabetes - Ejection.Fraction)
update(cox_time, .~. - Smoking - Pletelets - Diabetes - Ejection.Fraction - Gender)
update(cox_time, .~. - Smoking - Pletelets - Diabetes - Ejection.Fraction - Gender - Sodium)
update(cox_time, .~. - Smoking - Pletelets - Diabetes - Ejection.Fraction - Gender - Sodium - Anaemia)
update(cox_time, .~. - Smoking - Pletelets - Diabetes - Ejection.Fraction - Gender - Sodium - Anaemia - CPK)

cox_time_reduced <- update(cox_time, .~. - Smoking - Pletelets - Diabetes - Ejection.Fraction - Gender - Sodium - Anaemia - CPK)


# 5.) TEST LINEARITY ####

## 5.1) Martingale residuals ####
# TODO: They seem all quite non-linear, use fixes as described in lecture 5.5 
# SOLVE: Using fractional polynomials
jpeg("../images/martingale_res_grid_plot.jpeg", quality = 75, width = 1200, height = 600)
par(mfrow = c(2, 3), mar = c(4, 4, 2, 1))
mart_res_reduced <- residuals(cox_reduced, type = "martingale")
cex.axis = 1.2
scatter.smooth(mart_res_reduced ~ Age, xlab = "Age", cex.lab = cex.axis)
scatter.smooth(mart_res_reduced ~ Ejection.Fraction, xlab = "Ejection Fraction", cex.lab = cex.axis)
scatter.smooth(mart_res_reduced ~ Sodium, xlab = "Sodium", cex.lab = cex.axis)
scatter.smooth(mart_res_reduced ~ Creatinine, xlab = "Creatinine", cex.lab = cex.axis)
scatter.smooth(mart_res_reduced ~ CPK, xlab = "CPK", cex.lab = cex.axis)
dev.off()

## 5.2) Deviance residuals #### 
# (Can help to identify outliers (subjects with poor fit))
### Ejection.Fraction, Sodium and Creatine have a strange behaviour

jpeg("../images/deviance_res_grid_plot.jpeg", quality = 75, width = 1200, height = 600)
par(mfrow = c(2, 3), mar = c(4, 4, 2, 1))
dev_res_reduced <- residuals(cox_reduced, type = "deviance")
cex.axis = 1.2
scatter.smooth(dev_res_reduced ~ Age, xlab = "Age", cex.lab = cex.axis)
scatter.smooth(dev_res_reduced ~ Ejection.Fraction, xlab = "Ejection Fraction", cex.lab = cex.axis)
scatter.smooth(dev_res_reduced ~ Sodium, xlab = "Sodium", cex.lab = cex.axis)
scatter.smooth(dev_res_reduced ~ Creatinine, xlab = "Creatinine", cex.lab = cex.axis)
scatter.smooth(dev_res_reduced ~ CPK, xlab = "CPK", cex.lab = cex.axis)
dev.off()

jpeg("../images/dev_res_vs_lin_pred_reduced.jpeg", quality = 75, width = 700, height = 450)
dev_res_reduced <- residuals(cox_reduced, type = "deviance") 
lin_pred_reduced <- cox_reduced$linear.predictors
scatter.smooth(dev_res_reduced~lin_pred_reduced)
abline(h=2, lty=3)
abline(h=-2, lty=3)
dev.off()

summary(dev_res_reduced) 

####### COMMENT
# We do not seem to have any outliers, as all deviance residuals are < 3 in 
# absolute terms (note: deviance residuals are standard normal). The mean is 
# close to zero.
# However the residuals are not evenly scattered and we are seeing a slightly 
# decreasing trend. Also, we can recognize two clusters in the data.

## 5.3) Fractional Polynomials ####
# Solving non-linearity in Age, Ejection.Fraction, Sodium, Creatinine and CPK 
# using Fractional polynomials:

cox_all_mfp = mfp(surv_obj ~ fp(Age) + fp(Ejection.Fraction) + fp (Sodium) 
                 + fp(Creatinine) + fp(CPK) + Gender + Smoking + Diabetes + BP 
                 + Anaemia, family=cox, method="breslow",verbose=T, select=1, 
                 alpha=0.05, data=heart_data)
cox_all_mfp

####### COMMENT
# The results tells that the variables Ejection.Fraction should be added to the 
# model with a negative power of -2 and the variable Creatinine with a negative 
# power of -1, the other variables are okay with just the linear form.


# 6) COMPARE THE MODELS AND CHOOSE THE BEST #### 

# # Loglikelihood
# log_reduced <- cox_reduced$loglik 
# log_time_reduced <- cox_time_reduced$loglik
# log_all_mfp <- cox_all_mfp$loglik
# 
# # Likelihood ratio test statistics # -> should not be used to compare because 
# # models are not nested
# -2*diff(log_reduced)
# -2*diff(log_time_reduced)
# -2*diff(log_all_mfp)

# AIC
aic_reduced = AIC(cox_reduced)
aic_time_reduced = AIC(cox_time_reduced)
aic_all_mfp = AIC(cox_all_mfp)


# 7.) CLASSICAL PARAMETRIC MODELS ####

# TODO: describe the hazard function of the selected model
# Exp has the lowest AIC

## 7.1) PH models
# Fit Weibull model
weibullph_model <- flexsurvreg(surv_obj ~ Age + Ejection.Fraction + Sodium + Creatinine + Pletelets + CPK + Gender + Smoking + Diabetes + BP + Anaemia, 
                               method = "Nelder-Mead", hessian = FALSE,
                               data = heart_data, dist = "weibullPH")
weibullph_model
plot(weibullph_model,type="hazard")

# Fit exponential model
expph_model <- flexsurvreg(surv_obj ~ Age + Ejection.Fraction + Sodium + Creatinine + Pletelets + CPK + Gender + Smoking + Diabetes + BP + Anaemia, data = heart_data,
                           method = "Nelder-Mead", hessian = FALSE,
                           # control = flexsurv:::control.flexsurvreg(btol=1e-6, stol=1e-6),
                           ,dist = "exponential")
expph_model
plot(expph_model,type="hazard")


# Fit GOMPERTZ model
gompertzph_model <- flexsurvreg(surv_obj ~ Age + Ejection.Fraction + Sodium + Creatinine + Pletelets + CPK + Gender + Smoking + Diabetes + BP + Anaemia, 
                                method = "Nelder-Mead", hessian = FALSE,
                                data = heart_data, dist = "gompertz")
gompertzph_model
plot(gompertzph_model,type="hazard")


# 8.) Parametric models with only AFT parametrization ####

# Fit log-normal model AFT
lognormalaft_model<- flexsurvreg(surv_obj ~ Age + Ejection.Fraction + Sodium + Creatinine + Pletelets + CPK + Gender + Smoking + Diabetes + BP + Anaemia, 
                                 method = "Nelder-Mead", hessian = FALSE,
                                 data = heart_data, dist = "lnorm")
lognormalaft_model
plot(lognormalaft_model,type="hazard")

# Fit Generalized gamma (stable) AFT
gamaaft_model <- flexsurvreg(surv_obj ~ Age + Ejection.Fraction + Sodium + Creatinine + Pletelets + CPK + Gender + Smoking + Diabetes + BP + Anaemia, 
                             method = "Nelder-Mead", hessian = FALSE,
                             data = heart_data, dist = "gengamma")
gamaaft_model
plot(gamaaft_model,type="hazard")

#	Fit Log-logistic AFT
loglogisticaft_model <- flexsurvreg(surv_obj ~ Age + Ejection.Fraction + Sodium + Creatinine + Pletelets + CPK + Gender + Smoking + Diabetes + BP + Anaemia, 
                                    method = "Nelder-Mead", hessian = FALSE,
                                    data = heart_data, dist = "llogis")
loglogisticaft_model
plot(loglogisticaft_model,type="hazard")


# 9.) Compare models using AIC and BIC ####
AIC(expph_model,weibullph_model, gompertzph_model,lognormalaft_model,
    gamaaft_model,loglogisticaft_model )
BIC(expph_model,weibullph_model, gompertzph_model,lognormalaft_model,
    gamaaft_model,loglogisticaft_model)

# Select the model with the lowest AIC
#TODO: best model is expph_model, describe hazard fuction

expph_model
plot(expph_model)
plot(expph_model,type="hazard")
plot(expph_model,type="cumhaz")


