## ggsflabel has to be installed from github
# devtools::install_github("yutannihilation/ggsflabel")
rm(list=ls())
pacman::p_load(dplyr,readr,sf,ggthemes,rnaturalearth,ggrepel,
               ggsflabel,ggplot2,purrr,toolPhD,ggtern)
## Background map
ctrys <- ne_countries(continent = "europe", scale = 10, type = "countries", returnclass = "sf")
ctrys <- ctrys %>% 
  mutate("fill" = "grey90",
         fill = ifelse(sovereignt == "Germany", "darkgoldenrod1", fill))

## Locations
locations <- read.csv("data/BRIWECS_location_coordinates.csv") %>% 
  st_as_sf( coords = c("Long", "Lat"), crs = st_crs(ctrys)) 
names(locations)[1] <- "Location"
locations <- locations %>% 
  mutate("nudge" = case_when(Location == "KIE" ~ 2,
                             Location == "HAN" ~ -2,
                             Location == "QLB" ~ 2,
                             Location == "RHH" ~ 2,
                             Location == "KAL" ~ -2,
                             Location == "GGE" ~ 2))

locations_1 <- locations %>% 
  filter(Location %in% c("KIE", "QLB", "RHH", "GGE"))
locations_2 <- locations %>% 
  filter(Location %in% c("HAN", "KAL"))

## Map
fig1 <- 
  suppressMessages(
    ggplot() +
      theme_map() + 
      theme(
        panel.background = element_rect(fill = "dodgerblue3",
                                        colour = "transparent", linewidth = 0.25),
      ) +
      geom_sf(data = ctrys, aes(fill = fill), colour = "black", size = 0.25) +
      geom_sf_label_repel(data = locations_1, fill = "white", size = 2.5, aes(label = Location), nudge_x = 2, label.r = unit(0, "lines")) +
      geom_sf_label_repel(data = locations_2, fill = "white", size = 2.5, aes(label = Location), nudge_x = -2, label.r = unit(0, "lines")) +
      geom_sf(data = locations, shape = 21, fill = "red3", size = 2.5, stroke = 0.5) +
      coord_sf(xlim = c(-3, 23), ylim = c(43, 59)) +
      ggspatial::annotation_scale( location = "bl",
                                   # style="ticks",
                                   pad_y = unit(0.15, "in"),
                                   text_cex = 2.2)+
      scale_fill_identity()+
      ggspatial::annotation_north_arrow(location = "b", 
                                        pad_x = unit(0.1, "in"), 
                                        pad_y = unit(0.3, "in"),
                                        width = unit(.3, "in"),height = unit(.3, "in"),
                                        style = ggspatial::north_arrow_nautical(fill = c("grey40", "white"),
                                                                                line_col = "grey20")) )

# fig1
# -------------------------------------------------------------------------
soil<- xlsx::read.xlsx(file="data/soil.xlsx",stringsAsFactors=F,sheetIndex = 1) %>% 
  mutate(
    across(Clay:Sand,as.numeric),
    Sand=1-Clay-Silt,
    Location=gsub("DKI","KIE",Location)
  ) %>% distinct()
data(USDA)
USDA_text <- USDA  %>% group_by(Label) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE)

usd <- USDA_text %>% filter(Label%in%c('Sandy Loam','Sand','Silt Loam')) %>% 
  mutate(soiltype=c('Sandy loam','Sand','Silty loam'))

soil <- bind_rows(soil %>% filter(!is.na(Sand)) ,
                  soil %>% filter(is.na(Sand)) %>% dplyr::select(-c(Clay:Sand)) %>% 
                    left_join(.,usd,"soiltype")
)
# -------------------------------------------------------------------------
process_data <- function(df) {
  dfp <- df %>% 
    group_by(Location,Clay,Silt,Sand) %>%  
    group_split() %>% 
    map_dfr(.,function(dff){
      if(nrow(dff)>1){
        yf <- dff$Year %>% substr(start = 3,stop = 4) 
        if(all(yf %>% as.numeric() %>% diff() %>% unique()==1)){
          yf <- yf%>% range() %>% paste(.,collapse="-")
        }else{
          yf <- yf%>% paste(.,collapse="-")
        }
        dff[1,] %>% dplyr::select(Location,Clay,Silt,Sand) %>% mutate(g=paste(Location,yf),
                                                                      year=yf)
      }else{
        yf <- dff$Year %>% substr(start = 3,stop = 4) 
        dff[1,] %>% dplyr::select(Location,Clay,Silt,Sand) %>% mutate(g=paste(Location,yf),
                                                                      year=yf)
        
      }
    })
  return(dfp)
}
# ternary plot based on USDA standard -------------------------------------------------------------------------
#What was recorded in file about the soil type 
styp <- soil %>% group_by(Location) %>% dplyr::select(Location,soiltype) %>% distinct()
# check with plot using three ratio 
db <- soil %>%process_data() %>% 
  mutate(across(Clay:Sand,as.numeric)) %>% 
  na.omit()
dbm <- db %>% group_by(Location) %>% summarise(across(Clay:Sand,mean))
fig2 <- suppressMessages(
  ggplot(data = USDA, 
         aes(
           y = Clay,
           x = Sand,
           z = Silt)) +
    theme_bw()+
    coord_tern(L = "x", T = "y", R = "z") +
    geom_polygon(
      aes(fill = Label),
      alpha = 0.0,size = 0.5,
      color = "darkgray",show.legend = F) +
    ggalt::geom_encircle(data = db %>% filter(Location%in%c("KIE","GGE")),
                         mapping=aes(color=Location),size=1,alpha=.5, expand=.02,spread=0.001) +
    geom_text(data = USDA_text,
              mapping=aes( 
                label = Label),
              color = 'darkgray',alpha=.7,size = 1.75) +
    geom_point(
      data = db,
      size=1,shape=1,stroke=.5) +
    theme_showarrows() +
    labs(yarrow = "clay (%)",
         zarrow = "silt (%)",
         xarrow = "sand (%)") +
    theme_clockwise() +
    guides(fill=FALSE, color=FALSE)+
    geom_text(
      data = dbm,
      mapping=aes(label=Location,color=Location),
      size=1.5,
      hjust=+.2,vjust=-0.2,
      show.legend = F
    )+
    geom_text(
      data = db %>% filter(Location%in%c("KIE","GGE")),
      mapping=aes(label=year),
      size=1.5,show.legend = T,
      hjust=+.15,vjust=1.2
    )+
    theme(axis.title = element_blank()))
# grid.arrange(fig2)
# -------------------------------------------------------------------------
p <- cowplot::plot_grid(fig1,
                        grid.arrange(fig2),
                        nrow=1,labels = c("a","b"),rel_widths = c(.7,1.2),
                        align = "hv") %>% suppressWarnings()
tiff(filename="figure/Fig1.tiff",
     type="cairo",
     units="cm",
     width=16,
     height=8, compress="lzw",
     pointsize=3,
     res=600,# dpi,
     family="Arial")
p %>% print()
dev.off()

