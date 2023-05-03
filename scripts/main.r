# Max Schlake, Jose Gabriel Escarraman, Desmond Reynolds, Sebastian Veuskens
setwd(file_path <- dirname(rstudioapi::getSourceEditorContext()$path))

heart_data <- read.csv("../data/S1Data.csv")
