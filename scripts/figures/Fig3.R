rm(list=ls())
pacman::p_load(purrr,dplyr,toolPhD,ggplot2,scales)
# devtools::install_github("rensa/stickylabeller")
unit<- xlsx::read.xlsx("data/Unit.xlsx",sheetIndex = 1) %>% 
  mutate(unit=case_when(is.na(unit)~" ",
                        grepl("\\#",unit)~gsub("\\#","Nbr.",unit),
                        T~unit))
raw <- read.csv2("data/BRIWECS_data_publication.csv") %>% 
  mutate(across(BBCH59:Protein_yield,as.numeric))

col_pal<- c("#df436b","#eb1010","#6ac0c9", "#0860c9",
            "#136940",  "#f9912d","#00A337","#850242","#f9ca6c")


# preprocessing -------------------------------------------------------------------------
long <- raw  %>% 
  tidyr::pivot_longer(BBCH59:Protein_yield,
                      names_to="trait",values_to = "Trait") %>% 
  filter(!is.na(Trait)) %>% left_join(unit,"trait")
names(col_pal) <- long$Treatment %>% unique()

fig1_sub <- raw %>% 
  mutate(Environment = paste(Location, Year, sep = "_")) %>%
  dplyr::select(Environment,Treatment,Seedyield,Harvest_Index_bio,Grain,Straw) %>% 
  tidyr::pivot_longer(Seedyield:Straw,values_to = "trait",names_to="Trait")%>%
  left_join(unit %>% 
              rename(Trait=trait) %>% dplyr::select(-Full.name),by="Trait") %>% 
  mutate(unit=case_when(!is.na(unit)~paste0("(",unit,")"),
                        T~""),
         Nam=paste(Trait,"\n",unit)
  )
# data range density -------------------------------------------------------------------------
# range 
# fig1_sub %>% 
#   group_by(Trait) %>% summarise(m=min(trait,na.rm = T),
#                                 M=max(trait,na.rm = T))
# density plot
fig1 <- fig1_sub %>% rename(Management=Treatment) %>% 
  ggplot() +
  aes(x = trait, 
      y = Environment, fill = Management,color=Management) +
  # ggridges::theme_ridges()+
  theme_classic()+
  theme(legend.position  = "bottom",
        axis.title.x = element_blank(),
        strip.background = element_blank()) +
  ggridges::geom_density_ridges(
    alpha = 0.5,
    # size=.3,
    # linewidth=.2,
    scale=1,# height
    rel_min_height=0.005# width higher when value is small
  ) +
  ylab("Year x Location (Y/L)")+
  scale_fill_manual(values=col_pal)+
  scale_color_manual(values=col_pal)+
  # nord::scale_color_nord('aurora')+
  # nord::scale_fill_nord('aurora')+
  scale_y_discrete(drop=FALSE) +
  ggh4x::facet_nested(~Nam,nest_line=T, 
                      switch = "x",# place strip to bottom
                      scales = "free_x",
                      independent = "x")
# number of observation -------------------------------------------------------------------------
fig2 <- long %>% 
  group_by(trait,sample.source) %>% summarise(n=n(),.groups = "drop") %>% 
  group_by(sample.source) %>% 
  mutate(trait = forcats::fct_reorder(trait, n)) %>%
  ggplot( aes(x=trait, y=n)) +
  geom_segment( aes(xend=trait, yend=0)) +
  geom_point( size=4, color="orange") +
  coord_flip() +
  ggh4x::facet_nested(sample.source~.,
                      nest_line=T, 
                      switch = "y",# place strip to bottom
                      scales = "free_y",space ="free_y"
                      # independent = "y"
  )+
  ggtitle(sprintf("total number of observation: %s",nrow(long)))+
  xlab("")+ylab("number of observations")+
  ggrepel::geom_text_repel(aes(label=n),
                           size=2.7,
                           hjust=0,
                           box.padding = -.1,
                           point.padding = 0,
                           nudge_y=0,
                           # nudge_x=0,
                           direction="x"
  )+
  scale_y_continuous(
    labels =label_number(scale_cut = cut_short_scale()),
    limits=c(0,37000)
  )+
  toolPhD::theme_phd_facet(b=10,r=10,strp.txt.siz = 8,
                           plot.title = element_text(size=10))
# -------------------------------------------------------------------------
cp <- cowplot::plot_grid(fig1+
                           theme(legend.key.size = unit(.5,"line"),
                                 legend.position = "top",
                                 legend.text = element_text(size=4),
                                 legend.title=element_text(size=5),
                                 axis.title = element_text(size=6),
                                 strip.text = element_text(size=6),
                                 plot.margin = margin(r=0,l=0),
                                 axis.text=element_text(size=5))+
                           guides(colour = guide_legend(nrow = 2),
                                  fill = guide_legend(nrow = 2)),
                         fig2+
                           theme_classic() +
                           theme(
                             strip.background = element_blank(),
                             plot.title = element_text(size=8),
                             plot.margin = margin(t=20,r=3,l=0),
                             # strip.text = element_text(size=6),
                             axis.title = element_text(size=6),
                             axis.text.x=element_text(size=5),
                             axis.text.y=element_text(size=5)),
                         nrow=1,labels = c("a","b"),align = "hv")%>% 
  suppressWarnings() %>% suppressMessages()
tiff(filename="figure/Fig3.tiff",
    type="cairo",
    units="cm",
    compression = "lzw",
    width=20,
    height=12,
    pointsize=3,
    res=600,# dpi,
    family="Arial")
print(cp)
dev.off()
