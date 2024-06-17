rm(list = ls())
pacman::p_load(dplyr,purrr,toolPhD,ggraph,igraph,toolPhD)
source("scripts/fun/Net_link.R")
# -------------------------------------------------------------------------
trait.vector <-c('BBCH59','Straw','Grain','Seedyield',
                 'BBCH87',"RIE",
                 'RUE',
                 "LAI",
                 "k",
                 'Crude_protein','TGW',"Harvest_Index_bio")

canopy_trait<- readRDS("data/additional_canopy_trait.RDS")%>% 
  dplyr::select(-Date) %>% 
  mutate(Sowing_date=as.numeric(Sowing_date))

raw_data <- read.csv2("data/BRIWECS_data_publication.csv") %>% 
  mutate(across(BBCH59:Protein_yield,as.numeric)
  ) %>% 
  tidyr::pivot_longer(BBCH59:Protein_yield,
                      names_to = "trait",values_to = "Trait") %>% 
  filter(
    Treatment=="HN_WF_RF"&Location%in%c("HAN","KIE"),#|Location=="KIE"
    trait%in%trait.vector,
    Year%in%c(2015:2017)) %>% 
  mutate(Trait=case_when(grepl("BBCH",trait)~360-Sowing_date+Trait,
                         T~Trait))%>%
  dplyr::select(-Sowing_date) %>%
  bind_rows(.,canopy_trait %>% dplyr::select(-Sowing_date)) %>% 
  # some traits has only one value in one Block of certain year, like BBCH87 in 2017 has only 25 values
  group_by(BRISONr,Year,Location,Treatment,trait) %>% 
  summarise(Trait=mean(Trait,na.rm = T)) %>% ungroup() %>% 
  mutate(env =interaction(Year,Location,Treatment)) %>% 
  tidyr::pivot_wider(names_from=trait,values_from = Trait)

# name of stability index
SIname <- 'Genotypic superiority measure'

SIshort <- 'Pi'

stable_table <- 
  raw_data %>%
  group_by(Location,Treatment) %>% 
  group_split() %>% 
  map_dfr(.,function(data){
    map(trait.vector,function(trait){
      # print(trait)
      rn <- paste0(SIshort,trait)
      
      sdf <- data%>%
        dplyr::select(c("BRISONr","Year",trait)) %>%
        na.omit(.)
      if(nrow(sdf)==0){
        
      }else{
        res <- 
          tryCatch(
            toolStability::genotypic_superiority_measure(sdf,
                                                         trait,
                                                         'BRISONr',c("Year"),unit.correct = T),
            warning=function(w) {w; print(paste(trait,data$Location[1],data$Treatment[1]))}
          )%>%
          rename_with(~paste0(SIshort,trait),starts_with("genotypic.superiority.measure")) %>% 
          mutate(Location=data$Location[1],
                 Treatment=data$Treatment[1])
        return(res)
      }
    })%>% .[map_lgl(.,~{!is.null(.x)})]%>% 
      Reduce("full_join",.) %>% suppressMessages() 
  })

# # reorder the column based on the order of first trait then SI
stable_table <- stable_table %>% 
  dplyr::select(matches("^[^MP]"), starts_with("M"), starts_with("P")) %>% 
  rename_with(~gsub("Mean.","",.x),starts_with("Mean"))
# -------------------------------------------------------------------------
thresh <- 0# not applying threshold
network_type <- "Trait_network"

re <- stable_table %>% 
  group_by(Treatment,Location) %>% 
  group_split() %>% 
  map_dfr(.,function(datat){
    
    if(network_type=='SI_network'){
      subtble <- datat %>% .[,c(1,grep('Pi',names(.)))]%>%
        select(where(~!all(is.na(.))))
    }else{
      subtble <- datat %>% .[,grep('^(?!Pi)',names(.),perl=T)]%>%
        select(where(~!all(is.na(.))))
    }
    nodename <- names(subtble)[4:ncol(subtble)]
    
    node0 <- data.frame(
      # node index
      Id=1:length(nodename),
      Mix=gsub('Pi','',nodename),
      trait=nodename,
      SI=ifelse(grepl('Pi',nodename),yes =1,no =2),
      group=case_when(grepl('BBCH',nodename)~"P",
                      nodename%in%c('Kernel','TKW','Seedyield','Crude_protein')~"Y",
                      T~"C"),
      Shape=case_when(grepl('BBCH',nodename)~"circle",
                      nodename%in%c('Kernel','TKW','Seedyield','Crude_protein')~"square",
                      T~"rectangle")
    )
    
    
    link <- Net_link(node0,
                     subtble %>% select(-c(Treatment,Location)),
                     thresh.pp = thresh)%>%
      mutate(sign=as.factor(sign))# positive or negative
    # for labelling purpose, select unique from and to 
    node <- node0[unique(c(link$from,link$to)),] %>% 
      .[order(.$Id),] %>% 
      mutate(group=case_when(group=='P'~'phenology',
                             group=='Y'~'yield component',
                             T~'canopy'),
             Name=case_when(grepl('Pi',trait)~ paste0('bolditalic(P)[bold(\"i,',' ',gsub('Pi','',trait),'\")]'),
                            T~trait)
      )
    # for exporting table
    nam_con <- node0$trait%>%sub('Pi','',.)
    names(nam_con) <- node0$Id
    link_tb <- link%>%
      mutate(from=nam_con[from%>%as.character()],
             to=nam_con[to%>%as.character()],
             type=network_type,
             Treatment=subtble$Treatment[1],
             Location=subtble$Location[1]
      )%>%
      rename("r"="line")%>%
      dplyr::select(from,to,r,Treatment,Location) 
    return(link_tb)
  })

write.csv(re,"output/network_link.csv",row.names = F)
