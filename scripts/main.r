# Max Schlake, Jose Gabriel Escarraman, Desmond Reynolds, Sebastian Veuskens
rm(list=ls())

# Load libraries
library(survival)
library(ggplot2)
library(rstudioapi)
library(DescTools)

# Load data
setwd(file_path <- dirname(rstudioapi::getSourceEditorContext()$path))

heart_data <- read.csv("../data/S1Data.csv")

attach(heart_data) 

heart_data$Event <- as.factor(heart_data$Event)

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
surv_obj <- Surv(TIME, Event==1) 
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

# FIXME: Strange plot, does not make sense yet
for(i in c(3:7, ncol(heart_data))){
    dev.new()
    plot(survfit(surv_obj~heart_data[,i]), fun="cloglog", yscale=-1, col=1:nlevels(heart_data[,i]), xlab="Days", ylab="Estimated -log(-log S(t))")
    # legend("bottomleft", title= colnames(heart_data)[i], legend = levels(heart_data[,i]), col=1:nlevels(heart_data[,i]), lty=1)
}

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

# Multivariate Cox models #### 
cox_all <- coxph((surv_obj ~ Gender + Smoking + Diabetes + BP + Anaemia), method = "breslow") # breslow is Wilcoxon test
cox_all 

ggforest(cox_all, data = heart_data)
# All Covariates are not significant except BP, although Anaemia shows aberrant behaviour that is almost significant 


