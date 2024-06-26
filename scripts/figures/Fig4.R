rm(list = ls())
pacman::p_load(dplyr,purrr,ggplot2,toolPhD,ggpmisc,ggh4x)
R_raw <- readRDS("output/sma_Raw.RDS") 
source("scripts/fun/Cor_fun.R")
R_raws <- R_raw%>%
  ungroup() %>% 
  tidyr::separate("combi",c("F1","F2"),sep="\n",)%>% 
  filter( R_type=="R2_sma")

# -------------------------------------------------------------------------
fillpalette <- viridis::viridis(3)
names(fillpalette) <- c("HN_WF_RF","HN_NF_RF","LN_NF_RF")

nice_break <- function(range.vec,nbk,digit,include=T){
  axis_breaks <- pretty(range.vec, nbk)
  
  axis_breaks <- round(axis_breaks,digit)
  if(include==T){
    axis_breaks <-  axis_breaks[axis_breaks >= min(range.vec) & axis_breaks <= max(range.vec)]
  }
  return(axis_breaks)
}

lsddf<- function(formu,dat){
  datinfo <- dat %>% dplyr::select(f2,factor,type,f1) %>% distinct()
  
  resl <- aov(as.formula(formu),data=dat) %>%
    agricolae::LSD.test(.,"f1", group=TRUE)
  if (is.null(resl)|length(unique(resl$groups$groups))==1){
    datinfo %>% mutate(
      R2=0,
      groups=''
    )
  }else{
    resl$means %>% .["Min"] %>%
      merge(.,resl$groups%>% .["groups"],by=0, all=TRUE) %>% 
      rename(R2=Min,f1=Row.names) %>% 
      left_join(datinfo,"f1")
  }
}

dfun <- function(x){
  x %>%
    .[,!sapply(x, function(x)all(is.na(x)))] %>% 
    as.data.frame() %>% 
    arrange({{typ}}) %>% 
    `rownames<-`(.[[typ]]) %>% 
    dplyr::select(-all_of(c("type","R_type","f2",typ))) %>%  
    rbind(rep(.8,ncol(.)) , rep(0,ncol(.)) , .)
}

# general view , regardless environmental group -------------------------------------------------------------------------
tr_vec<- c('Seedyield','Crude_protein','Straw','Harvest_Index_bio','Biomass',
           'TGW','Grain','Grain_per_spike_bio',"Spike_number_bio")
col_pal<- viridis::viridis(4)
names(col_pal) <- c("HN_NF_RF","HN_WF_RF","LN_NF_RF","LN_WF_RF")
shp2 <- c(1,2)
names(shp2) <- c("+","-")

sdf <- R_raw %>% dplyr::select(trait,sma_slope,sma_sig,R2,Rsign) %>% distinct() %>% 
  filter(trait%in%tr_vec) %>% 
  mutate(trait=factor(trait,levels=tr_vec)) %>%
  group_by(trait) %>%
  mutate(m=mean(R2),
         trait=case_when(trait=='Crude_protein'~'GP',
                         trait=='Seedyield'~'GY',
                         trait=='Grain'~'GN',
                         trait=='Harvest_Index_bio'~"HI",
                         trait=='Spike_number_bio'~'SN',
                         trait=='Grain_per_spike_bio'~'GpS',
                         trait=='Biomass'~'SDM',
                         T~trait
         ))

resl <- aov(as.formula("R2~trait"),data=sdf) %>%
  agricolae::LSD.test(.,"trait", group=TRUE) %>% .$groups %>% 
  tibble::rownames_to_column(var = "trait") %>% mutate(R2=-0.05)

plotdf <- sdf %>%
  summarise(
    r1=paste0(
      "M:",max(R2) %>% round(.,2),"\n",
      "A:",mean(R2)%>% round(.,2),"\n",
      "m:",min(R2) %>% round(.,2)),
    R2=1,
  ) %>% ungroup() 
p <- sdf%>% 
  ggplot(aes(reorder(trait, -m),R2))+
  ggplot2::geom_violin(alpha = 0.5, 
                       position = position_dodge(width = 1), size = 1, color = NA) + 
  ggbeeswarm::geom_quasirandom(size = 1.5, shape=1,alpha=.3,
                               dodge.width = 1) + 
  ggplot2::geom_boxplot(outlier.size = -1,  
                        position = position_dodge(width = 1), lwd = .6, color="darkblue",
                        width = 0.3, alpha = 0.05, show.legend = F) +
  toolPhD::theme_phd_talk(t=10,b=10) + 
  ggplot2::theme(panel.grid.major.x = element_blank(), 
                 legend.position = "bottom",
                 axis.title.x = element_blank())+
  geom_text(  data=plotdf,
              mapping=aes(trait,R2,label=r1),
              color="darkblue",size=4)+
  geom_text(  data=resl,
              mapping=aes(trait,R2,label=groups),
              color="darkred",size=7)+
  scale_y_continuous(limits = c(-0.1,1.1),breaks = c(0,.4,.8))+
  ylab(parse(text='italic(R)[sma]^2'))
p <- sdf%>%
  ggplot(aes(reorder(trait, -m),R2))+
  # ggplot2::scale_fill_viridis_d(option = "c") +
  ggplot2::geom_violin(alpha = 0.5,
                       position = position_dodge(width = 1), size = 1, color = NA) +
  ggbeeswarm::geom_quasirandom(size = 1.5, shape=1,alpha=.3,
                               dodge.width = 1,
                               aes(
                                 # shape=R_sign,
                                 # fill=sma_sig
                               )) +
  ggplot2::geom_boxplot(outlier.size = -1,
                        position = position_dodge(width = 1), lwd = .6, color="darkblue",
                        width = 0.3, alpha = 0.05, show.legend = F) +
  toolPhD::theme_phd_main(t=10,b=10) +
  ggplot2::theme(panel.grid.major.x = element_blank(),
                 legend.position = "bottom",
                 axis.title.x = element_blank())+
  geom_text(  data=plotdf,
              mapping=aes(trait,R2,label=r1),
              color="darkblue",size=3)+
  geom_text(  data=resl,
              mapping=aes(trait,R2,label=groups),
              color="darkred",size=5)+
  scale_y_continuous(limits = c(-0.1,1.1),breaks = c(0,.4,.8))+
  ylab(parse(text='italic(R)[sma]^2'))

tiff(filename="figure/Fig4.tiff",
     type="cairo",
     units="cm",
     compression = "lzw",
     width=16,
     height=10,
     pointsize=12,
     res=600,# dpi,
     family="Arial")
p
dev.off()

