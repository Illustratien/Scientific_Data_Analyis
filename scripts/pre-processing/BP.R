rm(list = ls())
pacman::p_load(dplyr,ggplot2,scales,purrr,relaimpo,foreach)
source("scripts/fun/SW_fun.R")
# -------------------------------------------------------------------------
geno_info <-  read.csv("data/cultivar_info.csv")  %>% 
  relocate(RYear, .after = last_col()) %>% 
  dplyr::select(-c(quality:breeder)) %>% 
  # select only Kai's cultivars
  filter(kai==T)

env_gen_averaged <- readRDS("output/BLUE.RDS") %>% 
  merge(geno_info)%>%
  rename(Genotype=BRISONr) %>% 
  mutate_all(~ifelse(is.nan(.), NA, .)) %>% 
  mutate(Environment = paste(Year, Location, Treatment, sep="/")) %>% 
  filter(!trait=="RYear",
         !Location%in%c("RHH"),
         Year<2018,
         Treatment%in%c("HN_WF_RF","LN_NF_RF","HN_NF_RF"))

# slide window df and graph -------------------------------------------------------------------------
PointNr <- 10 ## number of point of which sliding w
## Do average and transform to long format

s_window_ls <- env_gen_averaged %>% 
  group_by(trait,Environment) %>% 
  arrange(Year,Location,Treatment,RYear) %>% 
  group_split() %>% 
  # rmove no value 
  .[!map_lgl(.,~{all(is.na(.x$emmean))})] %>% 
  .[!map_lgl(.,~{length(unique(.x$emmean))==1})] %>% 
  # remove only single value
  .[!map_lgl(.,~{length(which(!is.na(.x$emmean)))==1})] %>% 
  .[map_lgl(.,~{nrow(.x)>11})] # nrow(data) - window + 1 should > 1, otherwise error in line 26 SW_FUN

n.cores <- parallel::detectCores() - 1
#create the cluster
my.cluster <- parallel::makeCluster(
  n.cores,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)

system.time(
  res_ls <- foreach(
    y  = 1:length(s_window_ls),      # index
    x = s_window_ls,                 # content
    .packages = c("dplyr","purrr","ggplot2") # your required packages
  ) %dopar% {
    source("scripts/fun/SW_fun.R")
    ## To avoid NA in calculating mean
    trait <- x$trait[1]
    Environment <- x$Environment[1]
    # from tibble to dataframe
    x <- x %>% as.data.frame()
    RYear <- SlidingWindow("Namean", x$RYear, PointNr,1)
    slide_mean <- SlidingWindow("Namean", x$"emmean", PointNr,1)
    slide_sd <- SlidingWindow("Nasd", x$"emmean", PointNr,1)
    df_plot <-  data.frame(Parameter = trait, 
                           Environment=Environment,
                           RYear=RYear,
                           mean=slide_mean, 
                           sd=slide_sd)
    
    lr <- try(lm(slide_mean~RYear))   
    if(any(class(lr)=="try-error")){
      print(paste(trait, " in ",Environment, "has no results"))
      # statistics
      ## Ab_BP: absolute breeding progress; v1970: value at 1970; v2010: value at 2010
      df <-  data.frame(Parameter = trait, 
                        Environment=Environment, 
                        Ab_BP=NA,Ab_BPsd=NA ,Ab_BPCIl=NA,
                        Ab_BPCIu=NA,
                        v1970=NA, v2010=NA,
                        p=NA, R2=NA)
      res <- list(df_plot,df)
      
    }else{
      # # statistics
      # BP: slope, v: prediction at specific years.
      Ab_BP=summary(lr)$coefficients[2,1]
      bpsd=summary(lr)$coef[2,2]
      bp.ci.l=confint(lr,"RYear",level=.95)[1]
      bp.ci.u=confint(lr,"RYear",level=.95)[2]
      v1970=as.numeric(lr$coefficients[2])*1970+as.numeric(lr$coefficients[1])
      v2010=as.numeric(lr$coefficients[2])*2010+as.numeric(lr$coefficients[1])
      p=round(summary(lr)$coefficients[1,4], digits=4)
      R2=round(summary(lr)$r.square, digits = 2)
      
      df <-  data.frame(Parameter = trait,
                        Environment=Environment,
                        Ab_BP=Ab_BP, 
                        Ab_BPsd=bpsd,
                        Ab_BPCIl=bp.ci.l,
                        Ab_BPCIu=bp.ci.u,
                        v1970=v1970,
                        v2010=v2010,
                        p=p, R2=R2)
      # plot
      if(p<0.05){
        pl <-tryCatch({
          BP_plot(df_plot,trait,Environment,lr)
        },error=function(condi){
          message(y)
        })
        res <- list(df_plot,df,pp=BP_plot(df_plot,trait,Environment,lr))
        names(res)[3] <- with(x[1,],paste(trait,Year,Location,Treatment,sep="-"))
        
      }else{
        res <- list(df_plot,df)
      }
      
    }  
    res
  }
)

doParallel::stopImplicitCluster()
## Create the sliding window mean, do the graphes

s_window <- res_ls  %>% map_dfr(.,~{.x[[1]]})
s_window_statistics <- res_ls  %>% map_dfr(.,~{.x[[2]]})
s_window_plot <- res_ls  %>%
  .[map_lgl(.,~{length(.x)>2})] %>%
  map(.,~{
    return(.x[[3]])
  })

##write sliding window results
write.csv (s_window, file = "output/Slidingwindow_BLUES.csv",
           row.names = FALSE)

write.csv  (s_window_statistics,  
            file = "output/Slidingwindow_BLUES_statistics.csv",
            row.names = FALSE)