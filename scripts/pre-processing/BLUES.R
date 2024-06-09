rm(list=ls())
pacman::p_load(dplyr,purrr,foreach,lme4)
source("scripts/fun/BLUEs_fun.R")
# -------------------------------------------------------------------------
data <- read.csv2("data/BRIWECS_data_publication.csv") %>% 
  mutate(across(BBCH59:Protein_yield,as.numeric),
         tinfect=rowSums(across(Stripe_rust:Fusarium), na.rm=TRUE)
         ) %>% 
  tidyr::pivot_longer(BBCH59:tinfect,
                      names_to = "trait",values_to = "Trait")

data_list <- data %>% 
  group_by(Location,Treatment,Year,trait) %>% 
  group_split() %>% 
  .[!map_lgl(.,~{length(which(!is.na(.x$Trait)))==1})] 

# run parallel -------------------------------------------------------------------------
n.cores <- parallel::detectCores() - 1
#create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)

system.time(
  res <- foreach(
    i  = 1:length(data_list),
    dfj = data_list,
    .packages = c('dplyr',"lme4")
  ) %dopar% {
    Blues(i,dfj)
  }
)

# doParallel::stopImplicitCluster()
# saveRDS(res,"data/BLUEs/BLUE_original_list.RDS",compress = T)

resdf <- res %>% 
  .[!map_lgl(.,~{is.character(.x)})] %>% 
  map_dfr(.,~{.x})

saveRDS(resdf,"output/BLUE.RDS",compress = T)
