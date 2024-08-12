rm(list=ls())
pacman::p_load(dplyr,purrr,ggplot2)
level_order <- c('HN_WF_RF', 'HN_NF_RF', 'LN_NF_RF') 

traitlist <- c('Straw','Harvest_Index_bio','TGW',
               "Spike_number_bio",'Grain_per_spike_bio',
               'Biomass','Seedyield','Grain')

un <- xlsx::read.xlsx("data/Unit.xlsx",sheetIndex = 1) %>% 
  mutate(unit=gsub('\\#','Nbr ',unit),
         unit=gsub('\\*',' x ',unit),
         unit=ifelse(is.na(unit),' ',unit)) %>% 
  rename(Parameter=trait)

blue <- readRDS("output/BLUE.RDS") %>%
  dplyr::filter(Treatment%in%level_order,
                !Location%in%c('RHH'),
                Year<2018)
m.gen <- 
  map_dfr(c('Seedyield','Straw','Harvest_Index_bio','Biomass',
            'TGW','Grain','Grain_per_spike_bio',"Spike_number_bio") ,
          function(tr){
            blue %>% 
              dplyr::filter(
                trait==tr) %>% 
              group_by_at(c("Treatment","Location","Year","trait")) %>% 
              summarise(BRISONr=length(unique(BRISONr)),.groups="drop") 
          })

env_trait<-  m.gen %>%
  group_by(Location,Treatment,trait) %>%
  summarise(n=n(),.groups="drop") %>% ungroup() %>% 
  filter(n==3) %>% dplyr::select(-n)%>% rename(Parameter=trait)

s_window_statistic <- read.csv("output/Slidingwindow_BLUES_statistics.csv") %>% 
  mutate(across(Ab_BP:Ab_BP,as.numeric),
         Treatment=Environment %>% strsplit(.,"/") %>% 
           map_chr(.,~{.x[3]}) %>% as.factor(),
         Year=Environment %>% strsplit(.,"/") %>% 
           map_chr(.,~{.x[1]}),
         Location=Environment %>% strsplit(.,"/") %>% 
           map_chr(.,~{.x[2]}),
         Treatment_Location=paste(Treatment, Location, sep="-"),
         Treatment_Year=paste(Treatment, Year, sep="-")
  ) %>% 
  left_join(.,un,"Parameter") %>% 
  left_join(env_trait,.,by = join_by(Location, Treatment, Parameter)) %>% 
  filter(!is.na(Ab_BP)) %>% 
  group_by(Parameter) %>%
  mutate(m=mean(Ab_BP),
         M=max(Ab_BP) %>%toolPhD::round_scale() %>% as.numeric(),
         cv=(sd(Ab_BP)/mean(Ab_BP))%>% toolPhD::round_scale(),
         r1=paste0(
           "M:",M ,"\n",
           "A:",mean(Ab_BP)%>% toolPhD::round_scale(),"\n",
           "m:",min(Ab_BP) %>% toolPhD::round_scale()),
         unit=case_when(grepl("\\/",unit)~paste0("(",unit," year",")"),
                        T~paste0("(",unit," /year",")"))
  ) 

# cvdf <- s_window_statistic %>%
#   dplyr::select(abbrev,cv,unit) %>%
#   distinct() %>% arrange(desc(cv)) %>% 
#   mutate(cvnam=paste0('CV[',"BP[",abbrev,"]",']','*"="*',cv),
#          abbrev=factor(abbrev,levels=orderid))
# cvdf
mordid <- s_window_statistic %>%
  dplyr::select(abbrev,Ab_BP,unit,Parameter) %>% 
  group_by(Parameter,abbrev) %>% 
  summarise(
    m=length(which(Ab_BP<0)),
    m=case_when(m==0~-1*mean(Ab_BP),
                T~m),.groups="drop") %>%
  arrange(m)


orderid<- s_window_statistic %>%ungroup() %>% 
  dplyr::select(abbrev,cv,unit) %>%
  distinct() %>% 
  mutate(cv=as.numeric(cv)) %>% 
  arrange(desc(cv))%>%
  dplyr::select(abbrev) %>% 
  mutate(letter=letters[1:n()])

# levels(sdf$abbrev)
# sdf$unit %>% unique()
p <- s_window_statistic%>% 
  mutate(
    unit=gsub("\\s+","~",unit) %>% 
      gsub("\\/",'*"\\/"*',.) %>% 
      gsub("\\~\\*",'\\~',.),
    letter=letters[match(abbrev,mordid$abbrev)],
    cvnam=paste0(
      '~"("*',
      letter,
      '*")"~',
      "BP[",abbrev,"]",'~',unit),
    p=case_when(
      p < 0.001~ "***",
      p >= 0.001 & p < 0.01 ~ "**",
      p >= 0.01 & p < 0.05~ "*",
      T~"not-significant"),
    cvnam2=paste0('~CV[',"BP[",abbrev,"]",']','*"="*',cv),
    mnam=paste0('bar(BP[',abbrev,'])*"="*',toolPhD::round_scale(m)),
    abbrev=factor(abbrev,levels=mordid$abbrev)) %>% 
  ggplot(aes(abbrev,Ab_BP,group=abbrev))+
  guides(color=guide_legend(title=parse(text="italic(p)-value")))+
  # ggplot2::scale_color_viridis_d(option = "D",name=) + 
  ggplot2::geom_violin(alpha = 0.5, 
                       position = position_dodge(width = 1), 
                       linewidth = .5) + 
  ggbeeswarm::geom_quasirandom(size = 3, shape=1,stroke=1,aes(
    color=p),
    dodge.width = 1) + 
  ggplot2::geom_boxplot(outlier.size = -1,  
                        position = position_dodge(width = 1), 
                        lwd = .6, color="darkorange",
                        width = 0.3, alpha = 0.05, show.legend = F) +
  toolPhD::theme_phd_facet(t=10,b=30) + 
  ggplot2::theme(panel.grid.major.x = element_blank(), 
                 legend.position = c(.8,.2),
                 # strip
                 strip.text = element_text(size=.5),
                 axis.title.x = element_blank(),
                 axis.title.y = element_blank(),
                 axis.text.x=element_blank())+
  scale_y_continuous(scales::cut_short_scale())+
  geom_hline(yintercept=0,linetype=3,linewidth=1,color="darkgray")+
  facet_wrap(~
               cvnam+
               cvnam2+mnam,
             # abbrev+unit+cvnam,
             scales = "free",
             dir="h",
             labeller =  label_parsed
             # stickylabeller::label_glue('({.L}) {abbrev} {unit}\n{cvnam}')
  )
p
png(filename="figure/Fig7.png",
    type="cairo",
    units="cm",
    width=20,
    height=22,
    pointsize=5,
    res=650,# dpi,
    family="Arial")
p
dev.off()

