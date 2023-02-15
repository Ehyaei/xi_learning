options(scipen = 999)
options(knitr.kable.NA = '')

##%######################################################%##
#                                                          #
####                    load library                    ####
#                                                          #
##%######################################################%##

library(tidyverse)
library(MPIThemes)
library(XICOR)
library(ggplot2)
library(latex2exp)
library(ggdag)
library(dagitty)
library(ggdag)
library(ggcorrplot)
library(ggpubr)
############################################################
#                                                          #
#                   Set Color and Theme                    #
#                                                          #
############################################################

set_color_theme()

############################################################
#                                                          #
#                      Load functions                      #
#                                                          #
############################################################

source("manuscript/chatterjee.R")
source("manuscript/linear_scm_generator.R")
