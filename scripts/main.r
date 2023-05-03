# Max Schlake, Jose Gabriel Escarraman, Desmond Reynolds, Sebastian Veuskens
# Load libraries
library(survival)
library(ggplot2)
library(tidyverse)

# Load data
setwd(file_path <- dirname(rstudioapi::getSourceEditorContext()$path))

heart_data <- read.csv("../data/S1Data.csv")

attach(heart_data) 

# EXPLORATORY ANALYSIS #### 
ggplot2(data = heart_data) +
    geom_bar(mapping = aes(x = Diabetes) + 
    geom_boxplot(mapping = (aes(x = Smoking))))

# Kaplan-Meier estimator
summary(heart_data)
result_simple <- survfit(Surv(TIME, Event) ~ 1, conf.type = "log-log") 
result_simple
summary(result_simple)

plot(result_simple)
# Cumulative hazard 
simple_cumhaz <- round(result_simple$cumhaz, 4)
plot(simple_cumhaz, type="l")

# Fleming-Harrington (Nelson-Aalen) estimator 
result_simple_na <- survfit(Surv(TIME, Event) ~ 1, conf.type = "log-log") 
plot(result_simple_na) 

