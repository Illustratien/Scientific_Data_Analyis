# install necessary packages

# combined raw data and remove outlier
source(file.path("scripts/pre-processing","BLUEs.R"))
# combined management files
source(file.path("scripts/pre-processing","BP.R"))
# visualization and overview
source(file.path("scripts/pre-processing","R2sma.R"))


#4
source(file.path("scripts/figure","R2sma_groupping.R"))

#7
source(file.path("scripts/figure","BP_plot.R"))