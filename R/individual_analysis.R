rm(list = ls())
source("R/functions.R")
source("R/classes.R")       

jag_i <- individual_analysis$new(id = 114)
jag_i$compare_dispersal()
