##write sliding window results
rm(list=ls())
pacman::p_load(dplyr,purrr,ggplot2,hrbrthemes,ggh4x)
windowsFonts("Arial"=windowsFont("Arial"))
# -------------------------------------------------------------------------
traitlist <- c('Straw','Harvest_Index_bio',
               'TGW',
               "Grain_per_spike_bio",
               "Spike_number_bio",
               'Seedyield')
s_window <- read.csv (file = "output/Slidingwindow_BLUES_statistics.csv") %>% 
  dplyr::select(1:3) %>%
  tidyr::separate(Environment,c("Year","Location","Treatment"),"/",remove = F) %>% 
  filter(Parameter%in%traitlist,
         !Location%in%c("RHH"),
         Treatment%in%c("HN_WF_RF","LN_NF_RF","HN_NF_RF"),
         Year<2018) %>% 
  tidyr::pivot_wider(names_from = Parameter,values_from = Ab_BP)
# -------------------------------------------------------------------------
fit <- lm(Seedyield~TGW+Harvest_Index_bio+Straw+
            Grain_per_spike_bio+Spike_number_bio
          , data=s_window) 
out  <- relaimpo::calc.relimp(fit, type="lmg")

Si_reg <-  data.frame(explain=out$lmg) %>% 
  tibble::rownames_to_column("Trait") %>% 
  merge(data.frame(coef=coef(fit)) %>% 
          tibble::rownames_to_column("Trait")) %>%  
  mutate(lab="all",
         nobs=nrow(na.omit(s_window))
  )
# check colinearity
# PerformanceAnalytics::chart.Correlation(s_window %>% dplyr::select(Grain_per_spike_bio:TGW))

# # subgroup -------------------------------------------------------------------------
env_ref <- c("Treatment","Location","Year")
BP_wide <- s_window 
BP_regg<-  
  map_dfr(env_ref,function(envec){
    
    subdat <- BP_wide %>%
      split(.,BP_wide[[envec]])
    
    imap_dfr(subdat,function(edata,elevel){
      # print(elevel)
      fit3 <- lm(Seedyield~TGW+Harvest_Index_bio+Straw+
                   Grain_per_spike_bio+Spike_number_bio, data=na.omit(edata)) 
      
      # PerformanceAnalytics::chart.Correlation(na.omit(edata) %>% 
      #                                           .[,-c(1:4)]) %>% print()
      
      out  <- tryCatch({     relaimpo::calc.relimp(fit3, type="lmg")},
                       error=function(cond){
                         list(lmg=NA)
                       })
      data.frame(explain=out$lmg) %>% 
        tibble::rownames_to_column("Trait") %>% 
        merge(data.frame(coef=coef(fit3)) %>% 
                tibble::rownames_to_column("Trait")) %>% 
        mutate( 
          egp=envec,
          el=elevel,
          nobs=nrow(na.omit(edata))
        )
    })
    
  }) %>% 
  filter(!is.na(explain))

# -------------------------------------------------------------------------
mergedf <- BP_regg %>%
  rbind(.,Si_reg %>% rename(egp=lab) %>% mutate(el="all")) %>% 
  mutate(Trait=case_when(Trait=='Crude_protein'~'GP',
                         Trait=='Seedyield'~'GY',
                         Trait=='Harvest_Index_bio'~"HI",
                         Trait=='Grain_per_spike_bio'~'GpS',
                         Trait=='Spike_number_bio'~'SN',
                         T~Trait)) %>%
  filter(nobs>6)

table <- mergedf %>% dplyr::select(-coef) %>% 
  tidyr::pivot_wider(values_from = explain,names_from = Trait) %>% 
  rowwise() %>% 
  mutate(sumnumeric = sum(c_across(GpS:TGW), na.rm = T),
         across(where(is.numeric),~toolPhD::round_scale(.x)))


fil.pal <- c("#264653",'#287271',"#2A9D8F",
             "#8AB17D","#E9C46A","#F4A261",
             "#E76F51")
names(fil.pal) <-  c("GP",
                     "GpS",
                     "HI",
                     "SN",
                     "Straw",
                     "TGW")

p <-
  mergedf%>% 
  mutate(egp=gsub("all"," ",egp),
         el=paste0(el,"\n(",nobs,")"),
         explain=round(explain*100,0)
  ) %>% 
  ggplot(.,aes(fill=Trait, y=explain, x=interaction(el,egp))) + 
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values=fil.pal)+ 
  scale_x_discrete(guide = guide_axis_nested(check.overlap = T)) +
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        ggh4x.axis.nestline = element_line(),
        axis.title =  element_text(color="black" ,size=7),
        legend.title= element_text(color="black",size=6  ,face = 'bold'),
        strip.text =  element_text(color="black"  ,size=5.5,face = 'bold'),
        legend.text = element_text(color="black" ,size=5  ,face = 'bold'),
        axis.text =   element_text(color="black"   ,size=5  ,face = 'bold'),
        legend.position = "bottom"
  )+
  xlab("levels of single grouping for growing conditions (Y/L/M)")+
  ylab(
    parse(text='relative~importance~of~BP[trait]~to~BP[yield]~"(%)"'))+
  geom_text(aes(label=round(explain,2)),size =2.5,fontface="bold", 
            color="white",family="Arial",
            position = position_stack(vjust = 0.5))+
  guides(fill=guide_legend(title=parse(text="BP[trait]")))


tiff(filename="figure/Fig8.tiff",
    type="cairo",
    units="cm",
    compression = "lzw",
    width=12,
    height=10,
    pointsize=12,
    res=600,# dpi,
    family="Arial")
p %>% print()
dev.off()

# coef -------------------------------------------------------------------------
BP_coef<-  
  map_dfr(env_ref,function(envec){
    
    subdat <- BP_wide %>%
      split(.,BP_wide[[envec]])
    
    imap_dfr(subdat,function(edata,elevel){
      fit3 <- lm(
        Seedyield~TGW+Harvest_Index_bio+Straw+
          Grain_per_spike_bio+Spike_number_bio, data=na.omit(edata)) 
      
      data.frame(coef=coef(fit3)) %>% 
        tibble::rownames_to_column("Trait") %>% 
        mutate( 
          coef=paste(toolPhD::round_scale(coef),
                     toolPhD::sig_pvalue(summary(fit3)$coefficients[,4])),
          egp=envec,
          el=elevel,
          nobs=nrow(na.omit(edata))
        )
    })
    
  }) %>% filter(nobs>6)


tbla <- BP_coef %>%
  rbind(data.frame(coef=coef(fit)) %>% 
          tibble::rownames_to_column("Trait") %>% 
          mutate( 
            coef=paste(toolPhD::round_scale(coef),
                       toolPhD::sig_pvalue(summary(fit)$coefficients[,4])),
            egp="all",
            el="all",
            nobs=nrow(na.omit(s_window))
          )) %>% 
  dplyr::select(-nobs,-egp) %>%
  tidyr::pivot_wider(values_from = coef,names_from = el) %>% 
  mutate(Trait=case_when(Trait=='Crude_protein'~'GP',
                         Trait=='Seedyield'~'GY',
                         Trait=='Harvest_Index_bio'~"HI",
                         Trait=='Grain_per_spike_bio'~'GpS',
                         Trait=='Spike_number_bio'~'SN',
                         T~Trait)) %>%
  
  arrange(Trait) %>%
  mutate(Trait=paste("BP[",Trait,"]")) 
# %>% 
# gridExtra::tableGrob(.,theme=gridExtra::ttheme_minimal(
#   core = list(fg_params=list(cex = .85,parse=T)),
#   colhead = list(fg_params=list(cex =.85,parse=T)),
#   rowhead = list(fg_params=list(cex = .85)))) 
# grid::grid.draw(tbla)


xlsx::write.xlsx(tbla, "output/Table.xlsx",append=T,
                 sheetName = "BP")
