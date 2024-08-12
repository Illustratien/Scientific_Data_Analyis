rm(list=ls())
pacman::p_load(purrr,dplyr,toolPhD,ggplot2,scales)
# devtools::install_github("rensa/stickylabeller")
unit<- xlsx::read.xlsx("data/Unit.xlsx",sheetIndex = 1) %>% 
  mutate(unit=case_when(is.na(unit)~" ",
                        grepl("\\#",unit)~gsub("\\#","Nbr.",unit),
                        T~unit))
raw <- read.csv2("data/BRIWECS_data_publication.csv") %>% 
  mutate(across(BBCH59:Protein_yield,as.numeric))
# -------------------------------------------------------------------------
s <- raw %>% 
  dplyr::select(Treatment,Year,Location) %>% 
  tidyr::separate(Treatment,into=c("nitrogen\nfertilizer",
                                   "fungicide\napplication",
                                   "water\navailability")) %>% 
  distinct() %>% 
  mutate(phase=case_when(Year<2018~"Phase I",
                         T~"Phase II") %>%
           factor(.,levels=c("Phase II","Phase I")),
         Year=as.character(Year))

tbla <-s %>% dplyr::select(1:3) %>%
  distinct() %>% 
  gridExtra::tableGrob(.,theme=gridExtra::ttheme_minimal(core = list(fg_params=list(cex = .70)),
                                                         colhead = list(fg_params=list(cex = .70)),
                                                         rowhead = list(fg_params=list(cex = .70)))) 
s <- raw %>% 
  dplyr::select(Treatment,Year,Location) %>% 
  tidyr::separate(Treatment,into=c("Nitrogen","Fungicide","Water_availability")) %>% 
  distinct() %>% 
  mutate(phase=case_when(Year<2018~"Phase I",
                         T~"Phase II") %>%
           factor(.,levels=c("Phase II","Phase I")),
         Year=as.character(Year))
mp <- s%>% 
  ggplot() +
  aes(x = interaction(Nitrogen,Fungicide), y = Year, color = interaction(Nitrogen,Fungicide)) +
  theme_classic() +
  ggh4x::facet_nested(phase~Water_availability+Location,nest_line = T,
                      scales="free",space = "free_x",switch = "both",
  )+
  scale_color_brewer(palette = "Set1") +
  geom_point(size=4,shape=15)+
  theme(legend.position = "bottom",
        axis.title=element_blank(),
        strip.text = element_text(size=12),
        axis.text.y=element_text(size=12),
        legend.text = element_text(size=10),
        # axis.line.x.bottom = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside"
  )
figdata <- cowplot::plot_grid(tbla, mp, nrow = 1,rel_widths =  c(1, 3))
tiff(filename="figure/Fig3.tiff",
    type="cairo",
    compress="lzw",
    units="cm",
    width=26,
    height=12,
    pointsize=3,
    res=650,# dpi,
    family="Arial")
figdata %>% print()
dev.off()