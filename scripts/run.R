# install necessary packages
dir.create("output",showWarnings = FALSE)
dir.create("figure",showWarnings = FALSE)
# calculate BLUEs
source(file.path("scripts/pre-processing","BLUEs.R"))
# calculate Breeding progress
source(file.path("scripts/pre-processing","BP.R"))
# calculate consistency
source(file.path("scripts/pre-processing","R2sma.R"))

# generate figures
figure.path <- 'scripts/figures'
eachfig <- list.files(figure.path)
for(i in eachfig){
  print(i)
  figure.path <- 'scripts/figures'
  source(file.path(figure.path,i))
}

figure.path <- 
source(file.path('scripts/figures','Fig9.R'))
