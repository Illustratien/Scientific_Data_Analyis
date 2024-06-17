# SI multilinear regression 
rm(list=ls())
pacman::p_load(dplyr,purrr,ggplot2,hrbrthemes,toolStability)
# -------------------------------------------------------------------------
v<- c('widetilde(c[i])',#adjusted coefficient variation
      "r[i]^2", #coefficient of regression
      "b[i]", # Coefficient of regression
      "s[di]^2",# deviation mean squares
      "W[i]",# ecovalence
      'S["xi"]^2',#environmental variance
      "D[i]^2",# genotypic stability
      "P[i]",# Genotypic superiority measure
      'sigma[i]^2',#stability variance
      "S[i]*4"#variance of rank
)
names(v)<- c("Adjusted.coefficient.of.variation",
             "Coefficient.of.determination",
             "Coefficient.of.regression",
             "Deviation.mean.squares",
             "Ecovalence",
             "Environmental.variance",
             "Genotypic.stability",
             "Genotypic.superiority.measure",
             "Stability.variance",
             "Variance.of.rank")

type <- c("x",
  "parametric",
  "parametric",
  "parametric",
  "parametric",
  "parametric",
  "parametric",
  "parametric",
  "non-parametric")

# c(
  # ,
  # "dynamic")

# -------------------------------------------------------------------------
blue <- readRDS("output/BLUE.RDS") %>%
  dplyr::filter(Treatment%in%c("HN_WF_RF","LN_NF_RF","HN_NF_RF"))
# -------------------------------------------------------------------------
m.gen <- 
  map_dfr(c('Seedyield','Straw','Harvest_Index_bio',
            'Biomass',"Crude_protein",
            'TGW','Grain','Grain_per_spike_bio',
            "Spike_number_bio","tinfect") ,function(tr){
              
              blue %>% 
                dplyr::filter(
                  trait==tr,
                  !Location%in%c("RHH"),
                  Treatment%in%c("HN_WF_RF","LN_NF_RF","HN_NF_RF"),
                  Year<2018
                ) %>% 
                group_by_at(c("Treatment","Location","Year","trait")) %>% 
                summarise(BRISONr=length(unique(BRISONr)),.groups="drop") 
            })
# make sure balance number of years for each location, treatment and trait 
env_list <- m.gen %>%
  group_by(Location,Treatment,trait) %>%
  summarise(n=n(),.groups="drop") %>% ungroup() %>% 
  dplyr::filter(n==3) %>% dplyr::select(-n)
# -------------------------------------------------------------------------
data_list<-   blue %>%
  filter(!Location%in%c("RHH"),
         Treatment%in%c("HN_WF_RF","LN_NF_RF","HN_NF_RF"),
         Year<2018) %>% 
  tidyr::unite("Env",c( 'Location',"Treatment","Year"),sep="-",remove = F) %>%
  dplyr::select(BRISONr:emmean,Env:trait,Treatment) %>%
  rename(Trait=emmean) %>%
  left_join(env_list,.) 

Si_df <- data_list%>%
  group_by(trait) %>% 
  group_split() %>% 
  map_dfr(.,~{
    table_stability(
      .x,
      genotype = "BRISONr",
      environment ="Env",
      trait="Trait",
      unit.correct = T,
      lambda = quantile(.x$Trait,.95)
    ) %>% 
      mutate(trait=.x$trait[1]) %>% 
      mutate_all(~ifelse(is.nan(.), NA, .))
  }) 

Si_noTFI<- Si_df %>% 
  dplyr::select(-c(Mean.Trait,Normality)) %>% 
  tidyr::pivot_longer(Safety.first.index:Ecovalence.modified,names_to="Si",
                      values_to="Si.value") %>%
  tidyr::pivot_wider(names_from="trait",values_from ="Si.value") %>%
  filter(!is.na(Seedyield),!Si%in%c("Safety.first.index","Ecovalence.modified")) %>% 
  dplyr::select(-tinfect)
Mean_TFI<- Si_df %>% 
  dplyr::select(-c(Normality)) %>% 
  filter(trait=="tinfect") %>%
  mutate(across(Safety.first.index:Ecovalence.modified,~Mean.Trait)) %>% 
  dplyr::select(-Mean.Trait) %>% 
  tidyr::pivot_longer(Safety.first.index:Ecovalence.modified,names_to="Si",
                      values_to="Si.value") %>%
  tidyr::pivot_wider(names_from="trait",values_from ="Si.value") %>%
  filter(!Si%in%c("Safety.first.index","Ecovalence.modified"))
Si <- merge(Si_noTFI,Mean_TFI) %>%  
  group_by(Si) %>% group_split()

# reliampo-------------------------------------------------------------------------
Si_reg<- map_dfr(Si[2:10],~{
  fit <- lm(Seedyield~TGW+Harvest_Index_bio+Straw+Crude_protein+
              Grain_per_spike_bio+Spike_number_bio+tinfect, data=.x) 
  out  <- relaimpo::calc.relimp(fit, type="lmg")
  
  data.frame(explain=out$lmg) %>% 
    tibble::rownames_to_column("Trait") %>% 
    mutate(Si=.x$Si[1]
           # nobs=nrow(.x)
           )
}) %>%  mutate(lab=match(Si,names(v)) %>% v[.]) 
# -------------------------------------------------------------------------
fil.pal <- c("#264653",'#287271',"#2A9D8F",
             "#8AB17D","#E9C46A","#F4A261",
             "#E76F51")
names(fil.pal) <-  c("SI[GP]",
                     "SI[GpS]",
                     "SI[HI]",
                     "SI[SN]",
                     "SI[Straw]",
                     "SI[TGW]",
                     "TFI")
p <- ggplot(Si_reg %>% 
              mutate(Trait=case_when(Trait=='Crude_protein'~"SI[GP]",
                                     Trait=='Seedyield'~"SI[GY]",
                                     Trait=='Harvest_Index_bio'~"SI[HI]",
                                     Trait=='Grain_per_spike_bio'~"SI[GpS]",
                                     Trait=='Spike_number_bio'~"SI[SN]",
                                     Trait=="tinfect"~"TFI",
                                     Trait=="Straw"~"SI[Straw]",
                                     Trait=="TGW"~"SI[TGW]",
                                     T~Trait))%>%
              arrange(lab) ,
            aes(fill=Trait, y=explain, x=lab)) + 
  geom_bar(position="stack", stat="identity")+
  scale_x_discrete(labels = function(l) parse(text=l))+
  scale_fill_manual(values=fil.pal,labels=scales::parse_format())+
  toolPhD::theme_phd_main()+
  ylab(parse(text="relative~importance~of~regressors~to~SI[yield]"))+
  theme(axis.text.x = element_text(vjust=.5),
        legend.text.align = 0)+
  guides(fill=guide_legend(title=parse(text="regressors")))+
  #legend.position = "bottom")+
  xlab(parse(text='stability~indices~(SI)'))+
  geom_text(aes(label=round(explain,2)),size = 3,fontface="bold", 
            color="white",
            # angle=90,
            position = position_stack(vjust = 0.5))
# p
table <- Si_reg %>% 
  dplyr::select(-lab) %>% 
  tidyr::pivot_wider(values_from = explain,names_from = Trait) %>% 
  rowwise() %>% 
  mutate(sumnumeric = sum(c_across(TGW:tinfect), na.rm = T)) %>% 
  mutate(across(where(is.numeric),~toolPhD::round_scale(.x)))

# coefficient-------------------------------------------------------------------------
Si_coef<- map_dfr(Si[2:10],~{
  fit <- lm(
    Seedyield~TGW+Harvest_Index_bio+Straw+Crude_protein+
      Grain_per_spike_bio+Spike_number_bio+tinfect, data=.x) 
  out  <- relaimpo::calc.relimp(fit, type="lmg")
  
  data.frame(coef=coef(fit)
             
  ) %>% 
    tibble::rownames_to_column("Trait") %>% 
    mutate(
      coef=paste(toolPhD::round_scale(coef),
                 '^"',
                 toolPhD::sig_pvalue(summary(fit)$coefficients[,4]),'"'),
      
      Si=.x$Si[1])
}) %>%
  mutate(lab=match(Si,names(v)) %>% v[.])

tbla <- Si_coef %>% dplyr::select(-Si) %>% 
  arrange(lab) %>% 
  tidyr::pivot_wider(values_from = coef,names_from = lab) %>% 
  mutate(Trait=case_when(Trait=='Crude_protein'~'GP',
                         Trait=='Seedyield'~'GY',
                         Trait=='Harvest_Index_bio'~"HI",
                         Trait=='Grain_per_spike_bio'~'GpS',
                         Trait=='Spike_number_bio'~'SN',
                         Trait=="tinfect"~"TFI",
                         T~Trait)) %>% 
  arrange(Trait) %>%
  mutate(Trait=paste("SI[",Trait,"]")) %>% 
  .[,c(1,7,2,3,8,4,9,5,10,6)] %>% 
  gridExtra::tableGrob(.,theme=gridExtra::ttheme_minimal(
    core = list(fg_params=list(cex = .85,parse=T)),
    colhead = list(fg_params=list(cex =.85,parse=T)),
    rowhead = list(fg_params=list(cex = .85)))) 
grid::grid.draw(tbla)


figdata <- cowplot::plot_grid(tbla, p, nrow = 2,rel_widths =  c(1, 1),labels = "AUTO")

# tiff(filename="output/Fig9.tiff",
#      type="cairo",
#      units="cm",
#      compression = "lzw",
#      width=18,
#      height=24,
#      pointsize=12,
#      res=600,# dpi,
#      family="Arial")
# figdata
# dev.off()


tiff(filename="figure/Fig9.tiff",
     type="cairo",
     units="cm",
     compression = "lzw",
     width=16,
     height=14,
     pointsize=12,
     res=600,# dpi,
     family="Arial")
p
dev.off()


# excel -------------------------------------------------------------------


Si_coef<- map_dfr(Si[2:10],~{
  fit <- lm(
    Seedyield~TGW+Harvest_Index_bio+Straw+Crude_protein+
      Grain_per_spike_bio+Spike_number_bio+tinfect, data=.x) 
  out  <- relaimpo::calc.relimp(fit, type="lmg")
  
  data.frame(coef=coef(fit)
             
  ) %>% 
    tibble::rownames_to_column("Trait") %>% 
    mutate(
      coef=paste(toolPhD::round_scale(coef),
                 toolPhD::sig_pvalue(summary(fit)$coefficients[,4])),
      
      Si=.x$Si[1])
}) %>%
  mutate(lab=match(Si,names(v)) %>% v[.])

tbla <- Si_coef %>% dplyr::select(-Si) %>% 
  arrange(lab) %>% 
  tidyr::pivot_wider(values_from = coef,names_from = lab) %>% 
  mutate(Trait=case_when(Trait=='Crude_protein'~'GP',
                         Trait=='Seedyield'~'GY',
                         Trait=='Harvest_Index_bio'~"HI",
                         Trait=='Grain_per_spike_bio'~'GpS',
                         Trait=='Spike_number_bio'~'SN',
                         Trait=="tinfect"~"TFI",
                         T~Trait)) %>% 
  arrange(Trait) %>%
  mutate(Trait=paste("SI[",Trait,"]")) %>% 
  .[,c(1,7,2,3,8,4,9,5,10,6)] 
xlsx::write.xlsx(tbla, "output/Table.xlsx",append = T,
                 sheetName = "SI")
