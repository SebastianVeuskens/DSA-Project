# Max Schlake, Jose Gabriel Escarraman, Desmond Reynolds, Sebastian Veuskens
rm(list=ls())

# Load libraries
library(survival)
library(ggplot2)
library(rstudioapi)
library(DescTools)
library(survminer)

# Load data
setwd(file_path <- dirname(rstudioapi::getSourceEditorContext()$path))

heart_data <- read.csv("../data/S1Data.csv")


heart_data$Gender <- as.factor(heart_data$Gender)
heart_data$Smoking <- as.factor(heart_data$Smoking)
heart_data$Diabetes <- as.factor(heart_data$Diabetes)
heart_data$BP <- as.factor(heart_data$BP)
heart_data$Anaemia <- as.factor(heart_data$Anaemia)

attach(heart_data) 

heart_data$Age_group <- cut(Age, quantile(Age))
# SOLVED: Should we continue with all the variables below or just leave them out because it is too much? -> Continue with all of them 
heart_data$Ejection.Fraction_group <- cut(Ejection.Fraction, quantile(Ejection.Fraction))
heart_data$Sodium_group <- cut(Sodium, quantile(Sodium))
heart_data$Creatinine_group <- cut(Creatinine, quantile(Creatinine))
heart_data$Pletelets_group <- cut(Pletelets, quantile(Pletelets))
heart_data$CPK_group <- cut(CPK, quantile(CPK))

detach(heart_data)
attach(heart_data)


# EXPLORATORY ANALYSIS #### 
PercTable(Event, Gender, rfrq="001", margins = c(1,2))
chisq.test(table(Event, Gender), correct = FALSE)

PercTable(Event, Smoking, rfrq="001", margins = c(1,2))
chisq.test(table(Event, Smoking), correct = FALSE)

PercTable(Event, Diabetes, rfrq="001", margins = c(1,2))
chisq.test(table(Event, Diabetes), correct = FALSE)

PercTable(Event, Anaemia, rfrq="001", margins = c(1,2)) 
chisq.test(table(Event, Anaemia), correct = FALSE)

hist(Age, breaks = 15, freq = F)

rate <- sum(Event==1) / sum(TIME) * 100 
rate 

summary(heart_data)

# Kaplan-Meier estimator
surv_obj <- Surv(TIME, Event) 
result_simple <- survfit(surv_obj~ 1, conf.type = "log-log") 
names(result_simple) 
summary(result_simple)

plot(result_simple)
ggsurvplot(result_simple, xlab = "Days", ylab = "Overall survival probability",
            data = heart_data, legend = "none")

# Additional estimates ####
# Cumulative hazard 
simple_cumhaz <- round(result_simple$cumhaz, 4)
plot(simple_cumhaz, type="l")

# Fleming-Harrington (Nelson-Aalen) estimator 
result_simple_na <- survfit(Surv(TIME, Event) ~ 1, conf.type = "log-log") 
plot(result_simple_na) 

# Univariate Analysis #### 
ggsurvplot(survfit(surv_obj~Gender,), conf.int = T, xlab = "Days", ylab = "Overall survival probability", data = heart_data)
ggsurvplot(survfit(surv_obj~Smoking,), conf.int = T, xlab = "Days", ylab = "Overall survival probability", data = heart_data)
ggsurvplot(survfit(surv_obj~Diabetes,), conf.int = T, xlab = "Days", ylab = "Overall survival probability", data = heart_data)
ggsurvplot(survfit(surv_obj~BP,), conf.int = T, xlab = "Days", ylab = "Overall survival probability", data = heart_data)
ggsurvplot(survfit(surv_obj~Anaemia,), conf.int = T, xlab = "Days", ylab = "Overall survival probability", data = heart_data)

ggsurvplot(survfit(surv_obj~Age_group,), conf.int = T, xlab = "Days", ylab = "Overall survival probability", data = heart_data)
ggsurvplot(survfit(surv_obj~Ejection.Fraction_group,), conf.int = T, xlab = "Days", ylab = "Overall survival probability", data = heart_data) # seems time dependend 
ggsurvplot(survfit(surv_obj~Sodium_group,), conf.int = T, xlab = "Days", ylab = "Overall survival probability", data = heart_data) # seems partially time dependend 
ggsurvplot(survfit(surv_obj~Creatinine_group,), conf.int = T, xlab = "Days", ylab = "Overall survival probability", data = heart_data)
ggsurvplot(survfit(surv_obj~Pletelets_group,), conf.int = T, xlab = "Days", ylab = "Overall survival probability", data = heart_data) # seems not important
ggsurvplot(survfit(surv_obj~CPK_group,), conf.int = T, xlab = "Days", ylab = "Overall survival probability", data = heart_data) # seems not important 

# Cloglog visual for PH-assumption  
# FIXME: Legend location 
for(i in c(3:7, (ncol(heart_data) - 5):ncol(heart_data))){
    # dev.new()
    plot(survfit(surv_obj~heart_data[,i]), fun="cloglog", yscale=-1, col=1:nlevels(heart_data[,i]), xlab="Days", ylab="Estimated -log(-log S(t))")
    legend("topright", title= colnames(heart_data)[i], legend = levels(heart_data[,i]), col=1:nlevels(heart_data[,i]), lty=1)
}

# SOLVED: Add grouped variables in both tests -> Or we leave it as it is
# Logrank test ####
survdiff(surv_obj~Gender)
survdiff(surv_obj~Smoking)
survdiff(surv_obj~Diabetes)
survdiff(surv_obj~BP) # significant
survdiff(surv_obj~Anaemia) # almost significant

# Wilcoxon test (Prentice correction) #### 
survdiff(surv_obj~Gender, rho = 1)
survdiff(surv_obj~Smoking, rho = 1)
survdiff(surv_obj~Diabetes, rho = 1)
survdiff(surv_obj~BP, rho = 1) # significant
survdiff(surv_obj~Anaemia, rho = 1) # almost significant

####################################
####### COMMENT ####################
# All Covariates are not significant except BP, although Anaemia shows aberrant behaviour that is almost significant 

# KM versus COX fit #### 
# SOLVED: We could add some more variables with the same analysis here -> But I think it is enough this way 
plot(survfit(surv_obj ~ BP, data=heart_data), col=1:nlevels(BP)) # KM plot
lines(survfit(coxph(surv_obj ~ BP, data=heart_data), newdata=data.frame(BP=c(0, 1)),se.fit=F), col=1:nlevels(BP), lty=2) # Cox predicted
legend("bottomright", title= "BP status", legend = levels(BP), col=1:nlevels(BP), lty=1)
####################################
####### COMMENT ####################
# The Cox model seems to fit the KM estimate quite well, no systematic over- or underestimation is visible (except maybe for BP status 1 at the end)

plot(survfit(surv_obj ~ Anaemia, data=heart_data), col=1:nlevels(Anaemia)) # KM plot
lines(survfit(coxph(surv_obj ~ Anaemia, data=heart_data), newdata=data.frame(Anaemia=c(0, 1)),se.fit=F), col=1:nlevels(Anaemia), lty=2) # Cox predicted
legend("bottomright", title= "Anaemia status", legend = levels(Anaemia), col=1:nlevels(Anaemia), lty=1)
####################################
####### COMMENT ####################
# The Cox model seems to fit the KM estimate quite well, no systematic over- or underestimation is visible


# Multivariate Cox models #### 
cox_all <- coxph((surv_obj ~ Age + Ejection.Fraction + Sodium + Creatinine + Pletelets + CPK + Gender + Smoking + Diabetes + BP + Anaemia), method = "breslow") # breslow is Wilcoxon test
cox_all 

ggforest(cox_all, data = heart_data)

plot(survfit(cox_all), main = "cph model", xlab="Days")
lines(result_simple, col="grey")

update(cox_all, .~. - Pletelets)
update(cox_all, .~. - Pletelets - Smoking)
update(cox_all, .~. - Pletelets - Smoking - Diabetes)
update(cox_all, .~. - Pletelets - Smoking - Diabetes - Gender)

cox_reduced <- update(cox_all, .~. - Pletelets - Smoking - Diabetes - Gender)

# SOLVED: Should we do the PH test on the cox_all or the cox_reduced data set (makes sense on all, but analysis outline says reduced) -> We can choose, he recommended both 
test_ph_reduced <- cox.zph(cox_reduced, transform="identity", terms = F) 
test_ph_reduced # Only Ejection.Fraction shows inconsistent behaviour for the ph assumption

for (i in 1:length(test_ph_reduced)) {
    plot(test_ph_reduced[i])
}
####################################
####### COMMENT ####################
# Age declines at the end
# Ejection.Fraction steadily declines (little)
# Sodium rises and then declines 
# Creatinine rising
# CPK rising then falling 
# BP shaped as wave 
# Anaemia slump at beginning, otherwise straight 

test_ph_reduced <- cox.zph(cox_reduced, transform="log", terms = F)
test_ph_reduced # Only Ejection.Fraction shows inconsistent behaviour for the ph assumption, but it is not significant anymore 

for (i in 1:length(test_ph_reduced)) {
    plot(test_ph_reduced[i])
}


# COX model with time varying coefficients  
# SOLVED: Which variables should we include here? -> All, but only Ejection.Fraction shows changing behaviour over time 
cox_time <- coxph(surv_obj ~ Age + Gender + Smoking + Diabetes + BP + Anaemia +tt(Ejection.Fraction), method = "breslow", data = heart_data, tt = function(x, t, ...) x * t)
cox_time 

cox_time <- coxph(surv_obj ~ Age + Gender + Smoking + Diabetes + BP + Anaemia +tt(Ejection.Fraction), method = "breslow", data = heart_data, tt = function(x, t, ...) x * log(t))
cox_time 

# Reduce time-dependent model #### 
update(cox_time, .~. - Smoking)
update(cox_time, .~. - Smoking - Diabetes)
update(cox_time, .~. - Smoking - Diabetes - Gender)
update(cox_time, .~. - Smoking - Diabetes - Gender - Anaemia)

cox_time_reduced <- update(cox_time, .~. - Smoking - Diabetes - Gender - Anaemia)

# Test linearity ####

# Martingale residuals 
# TODO: They seem all quite non-linear, use fixes as described in lecture 5.5 
mart_res_reduced <- residuals(cox_reduced) 
scatter.smooth(mart_res_reduced ~ Age)
scatter.smooth(mart_res_reduced ~ Ejection.Fraction)
scatter.smooth(mart_res_reduced ~ Sodium)
scatter.smooth(mart_res_reduced ~ Creatinine)
scatter.smooth(mart_res_reduced ~ CPK)

# Compare models and choose best one #### 

# Loglikelihood
log_1 <- cox_reduced$loglik 
log_2 <- cox_time_reduced$loglik 

# Likelihood ratio test statistics
-2*diff(log_1)
-2*diff(log_2) 
