rm(list=ls())
pacman::p_load(dplyr,purrr,foreach,ggplot2,doParallel,toolPhD,ggpmisc)
Rtabl_raw<- function(env_ele,dataf){
  vec<-c("trait",env_ele) 
  env_ref <- c("Treatment","Location","Year")
  d_l <- dataf %>% 
    group_by_at(vec) %>% 
    group_split() %>% 
    # check each group has minimum 10 values
    map(.,~{
      a <- .x
      .x%>%ungroup() %>% 
        group_by_at(setdiff(env_ref,env_ele)) %>% 
        summarise(n=length(unique(Trait)),.groups = "drop") %>% 
        filter(n>10) %>% dplyr::select(-n) %>% left_join(a,by=setdiff(env_ref,env_ele))
    }) %>% 
    # check at least 2 unique group to create correlation 
    .[map_lgl(.,~{.x %>% .[,setdiff(env_ref,env_ele)] %>% distinct() %>% nrow()>1})] 
  
  cor_vec <- foreach(
    data = d_l,
    .packages = c('dplyr',"tidyr","purrr","toolPhD")
  ) %dopar% {
    wide_df<- data %>%dplyr::select(-Env) %>%
      tidyr::pivot_wider(names_from = setdiff(env_ref,env_ele),
                         values_from = "Trait",names_sep = "-") 
    # %>% 
    # removeColumnsWithNA()
    
    tar <- paste0("(BRISONr|trait|data|",paste(env_ele,collapse="|"),")")
    cor_vec <-names(wide_df)  %>% .[!grepl(tar,.)] %>% 
      combn(.,2) %>% matrix(ncol = 2,byrow = T) %>% as.data.frame() 
    
    smalist <- sma_cor(wide_df %>% dplyr::select(-c(BRISONr:trait)))
    
    cor_vec<-  cor_vec %>% 
      mutate(npoint= map_dbl(1:nrow(cor_vec),~{
        wide_df[,cor_vec[.x,] %>% unlist()] %>% na.omit() %>% nrow()
      }))%>%
      tidyr::unite("combi",V1:V2,sep='\n') %>% 
      mutate(
        # R2_pearson=cor(wide_df %>% select(-c(BRISONr:trait,data)), use="complete.obs") %>%
        #   .[upper.tri(.)] ,
        # R2_spearman=cor(wide_df %>% select(-c(BRISONr:trait,data)), use="complete.obs",method = "spearman") %>%
        #   .[upper.tri(.)] ,
        R2_sma= smalist[[1]] %>%.[upper.tri(.)],
        sma_sig=smalist[[2]] %>%.[upper.tri(.)],
        sma_slope=smalist[[3]] %>%.[upper.tri(.)],
        trait=data[["trait"]][1],
        type=paste(env_ele,collapse="-"))
    
    
    if(length(env_ele)==1){
      cor_vec <- cor_vec %>% mutate(
        factor=data[[env_ele]][1]
      )
    }else{
      
      f <- wide_df%>% .[1,] %>% 
        mutate(
          factor=across(all_of(env_ele)) %>% 
            reduce(.,paste, sep = "-")) %>% 
        .$factor %>% unlist()
      
      cor_vec <- cor_vec %>% 
        mutate(factor= f)
    }
    
    cor_vec%>%
      na.omit() %>% 
      tidyr::pivot_longer(starts_with("R2"),values_to = "R2",names_to = "R_type") %>% 
      mutate(factor=as.character(factor))
  } %>% purrr::map_dfr(.,~{.x})
  
  
  return(cor_vec)
}
# data --------------------------------------------------------------------
blue <- readRDS("output/BLUE.RDS") %>%
  dplyr::filter(Treatment%in%c("HN_WF_RF","LN_NF_RF","HN_NF_RF"),
                !Location%in%c('RHH'))

datadf<- 
  blue %>%
  filter(!Location%in%c("DKI"),Year<2018)%>%
  tidyr::unite("Env",c( 'Location',"Treatment","Year"),sep="-",remove = F) %>%
  dplyr::select(BRISONr:emmean,Env:trait,Treatment) %>%
  rename(Trait=emmean)

# R table -------------------------------------------------------------------------
# create grouping combination of environments
e_vec <-list("Treatment","Location","Year",
             c("Treatment","Location"),c("Treatment","Year"))
env_ref <- c("Treatment","Location","Year")

n.cores <- parallel::detectCores() - 1
#create the cluster
my.cluster <- parallel::makeCluster(
  n.cores,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)
# 
R_raw <-
  map_dfr(e_vec,~{
    Rtabl_raw(.x,datadf)
  }) %>%
  mutate(f1=case_when(grepl("\\-",factor)~strsplit(factor,'-') %>% map_chr(.,~{.x[1]}),
                      T~factor),
         f2=case_when(grepl("\\-",factor)~strsplit(factor,'-') %>% map_chr(.,~{.x[2]}),
                      T~NA),
         Rsign=ifelse(R2>0,"+","-"),
         R2=R2 %>% .^2) %>%
  mutate(across(c(f2,factor,f1),function(x)gsub("\\.0+",'',x)))
saveRDS(R_raw,"output/sma_Raw.RDS",compress = T)
