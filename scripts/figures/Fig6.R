rm(list=ls())
pacman::p_load(dplyr,toolPhD,ggplot2,purrr)

apsim_trait <- c("flowering","straw","grain_number","yield","maturity","lai","HI",
                 "grain_protein","grain_size",'y_extinct_coef','y_rue')

bri_trait <-c('BBCH59','Straw','Grain','Seedyield','BBCH87',
              'LAI',"Harvest_Index_bio",'Crude_protein','TGW','k','RUE')

abv <- c("FT","Straw","GN","GY","MT","LAI","HI","GP","TGW",'k','rue')
ltb<- data.frame(apsim_trait,bri_trait,id=1:length(apsim_trait),abv)

apsim <-  xlsx::read.xlsx("data/APSIM_network.xlsx",sheetIndex = 1) %>% 
  dplyr::select(from:r)%>% filter(group%in%c("TP","TT"))  

lk<- function(vec,coln,df){
  df$id[match(vec,df[[coln]])]
  # match(vec,df[[coln]])
}

ord <- function(df){
  df %>% rowwise() %>%
    mutate(f=min(from,to),t=max(from,to)) %>% 
    dplyr::select(-c(from,to))
}


lkr<- function(vec,coln,df){
  df$apsim_trait[match(vec,df[[coln]])]
}

# -------------------------------------------------------------------------
apsim_df <-  apsim %>% filter(group%in%c("TT",'TP'),
                              from%in%apsim_trait&to%in%apsim_trait) %>%
  mutate(across(from:to,function(x){lk(x,"apsim_trait",ltb)})) %>%
  dplyr::select(-group) %>%
  rename(r_apsim=r) %>% ord() %>% as_tibble()


# each treatment and location -------------------------------------------------------------------------
trt.list<- read.csv("data/network_link.csv")%>%
  filter(from%in%bri_trait&to%in%bri_trait) %>% 
  group_by(Treatment) %>% group_split() %>% 
  map_dfr(.,function(tdf){
    
    tt <- tdf %>% 
      group_by(Location) %>% 
      group_split() %>% 
      map(.,~{
        vnam <- paste0('r_',.x$Location[1]) %>% quo_name()
        .x %>% 
          mutate(across(from:to,function(x){lk(x,"bri_trait",ltb)}))%>%
          rename(!!vnam:=r) %>% ord() %>% dplyr::select(-c(Treatment,Location))
      })
    tt[[length(tt)+1]] <- apsim_df
    
    tt%>%.[c(length(tt),1:(length(tt)-1))] %>% 
      Reduce("left_join",.) %>%suppressMessages() %>% 
      mutate(across(f:t,function(x){lkr(x,"id",ltb)})) %>%
      relocate(f,t) %>% 
      mutate(Treatment=tdf$Treatment[1])
    # tt <- c(apsim_df,tt)
  }) %>%  ungroup() %>% 
  mutate(across(f:t,function(x)ltb$"abv"[match(x,ltb$apsim_trait)]),
         s=paste(f,t,sep="-"))

supp.labs <- grep("^r(?!_apsim)",names(trt.list),perl=T,value = T) %>% gsub("r_","",.)
names(supp.labs) <- grep("^r(?!_apsim)",names(trt.list),perl=T,value = T)

# -------------------------------------------------------------------------
ipnutdf <- trt.list %>%filter(Treatment=="HN_WF_RF") %>% 
  tidyr::pivot_longer(cols = matches("^r(?!_apsim)",perl=T),
                      names_to = "r_name",values_to = "r_value") %>%
  na.omit() %>% 
  filter(r_name%in%c("r_KIE","r_HAN")) %>% 
  rowwise() %>%
  mutate(rmse=sqrt(mean((r_apsim-r_value)^2)),
         k=case_when(rmse<0.09&(abs(r_value)>.5&abs(r_apsim)>.5)~"a",
                     rmse>=0.09&(abs(r_value)>.5&abs(r_apsim)>.5)~"b",
                     # abs(r_value)>.5|abs(r_apsim)>.5~"c",
                     T~"c"),
         t=case_when(rmse<0.09&(abs(r_value)>.5&abs(r_apsim)>.5)~"a",
                     rmse>=0.09&(abs(r_value)>.5&abs(r_apsim)>.5)~"b",
                     abs(r_value)>.5|abs(r_apsim)>.5~"c",
                     T~"d"),
  ) 

p <-
  ipnutdf %>% 
  ggplot()+
  geom_hline(yintercept = 0,color="darkgray",linetype=35)+
  geom_vline(xintercept = 0,color="darkgray",linetype=35)+
  geom_abline(slope=1,intercept = 0,linetype=8,
              # linewidth=2,
              color="darkgray",alpha=.5)+
  geom_point(aes(r_value,r_apsim,
                 # color=k,
                 # shape=t
                 ),shape=1,
             show.legend = F)+
  # scale_size_area()+
  # scale_shape_manual(values=c(a=2,b=2,c=1,d=0))+
  scale_y_continuous(limits=c(-1,1))+
  scale_x_continuous(limits=c(-1,1))+
  toolPhD::theme_phd_main(b=5,legend.position="bottom", legend.box="vertical")+
  facet_grid(~r_name,labeller = labeller(r_name=supp.labs,type=label_parsed))+
  xlab(parse(text="italic(r)~from~field~data ~2015-2017"))+
  ylab(parse(text="italic(r)~from~simulation~data"))+
  ggrepel::geom_text_repel(seed = 57921,
                           data=ipnutdf %>% filter(abs(r_value)<=.5&abs(r_apsim)<=.5),
                           mapping=aes(r_value,r_apsim,label=s),
                           max.overlaps = 999,
                           box.padding = .1,size=2)+
  ggrepel::geom_text_repel(seed = 57921,
                           data=ipnutdf %>% filter(abs(r_value)>.5|abs(r_apsim)>.5),
                           mapping=aes(r_value,r_apsim,label=s),max.overlaps = 999,
                           box.padding = .1,size=2.5,fontface="bold")

# p
tiff(filename="figure/Fig6.tiff",
     type="cairo",
     units="cm",
     compression = "lzw",
     width=19,
     height=12,
     pointsize=12,
     res=600,# dpi,
     family="Arial")
egg::tag_facet(p,open = "",close = "",tag_pool = LETTERS[1:2])
dev.off()
# table -------------------------------------------------------------------------
# 
# pair_filter<- function(a,b){
#   vec <- c(paste(a,b,sep="-"),
#            paste(b,a,sep="-")
#   )
#   ipnutdf %>% 
#     filter(s%in%vec) %>% 
#     arrange(r_name,abs(r_value),rmse)
# }
# pair_filter("GY","GN")
# pair_filter("GY","GP")
# pair_filter("TGW","GN")
# pair_filter("Straw","GY")
# pair_filter("HI","GY")
# pair_filter("Straw","GY")
# pair_filter("HI","GN")
# pair_filter("GP","Straw")
# 
# 
# rs<- ipnutdf %>% 
#   filter(abs(r_value)>.4) %>% 
#   arrange(r_name,abs(r_value),rmse)
# rs<- ipnutdf %>% 
#   filter(abs(r_apsim)>.5) %>% 
#   arrange(r_name,abs(r_value),rmse)