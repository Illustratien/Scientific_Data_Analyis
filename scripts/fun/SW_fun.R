lr_fun<- function(x,tr,tt=FALSE){
  
  trait <- x$trait[1]
  Environment <- x$Environment[1]
  
  if(tt){
    if(tr=="emmean"){
      cgroup="Treatment"
    }else{
      cgroup="compare"
    }
    
    df_plot <- 
      x %>% group_by_at(cgroup) %>% group_split() %>% 
      map_dfr(.,~{
        x <- .x %>% data.frame()
        RYear <- SlidingWindow("Namean", x$RYear, PointNr,1)
        slide_mean <- SlidingWindow("Namean", x[[tr]], PointNr,1)
        slide_sd <- SlidingWindow("Nasd", x[[tr]], PointNr,1)
        lr <- try(lm(slide_mean~RYear)) 
        data.frame(Parameter = trait, 
                   Environment=Environment,
                   RYear=RYear,
                   mean=slide_mean, 
                   sd=slide_sd,
                   group2=.x[[cgroup]][1])
        
      })
    
    
    p <-BP_plot_gp(df_plot,trait,Environment)
  }else{
    x <- x %>% as.data.frame()
    
    RYear <- SlidingWindow("Namean", x$RYear, PointNr,1)
    slide_mean <- SlidingWindow("Namean", x[[tr]], PointNr,1)
    slide_sd <- SlidingWindow("Nasd", x[[tr]], PointNr,1)
    
    lr <- try(lm(slide_mean~RYear)) 
    df_plot <-  data.frame(Parameter = trait, 
                           Environment=Environment,
                           RYear=RYear,
                           mean=slide_mean, 
                           sd=slide_sd)
    p <- BP_plot(df_plot,trait,Environment,lr)
  }
  
  
  return(p)
}
fipal <- c("#1D6B51", "#51AA69","#C7EA56","#FF6A0E","#006896")
names(fipal) <- c("HN_WF","HN_NF","LN_NF","add_Fungicide","add_Nitrogen")
BP_plot_gp <- function(df_plot,trait,Environment){
  #plot breeding progress for each environment
  
  plt<- ggplot(df_plot %>% dplyr::filter(!is.na(RYear),!is.na(mean)),
               aes(x=RYear, y=mean,fill=group2,color=group2)) +
    # theme_bw() +
    labs(x="Year of release",
         y=trait, title=Environment %>% strsplit("-add") %>% unlist() %>% .[1]) + 
    geom_ribbon(aes(x = RYear, ymax = mean+sd,
                    ymin =mean-sd),show.legend = F, colour = NA,
                alpha = 0.1) +
    scale_fill_manual(values=fipal)+
    scale_color_manual(values=fipal)+
    geom_point(size=1.5,shape=1,alpha=.8) +
    stat_smooth(
      method="lm", se=FALSE, formula=y ~ x,linewidth=1,alpha=.3) +
    toolPhD::theme_phd_main(legend.position="bottom")+
    ggpmisc::stat_poly_eq(formula = y ~ x, 
                          aes(label = paste(
                            after_stat(eq.label),
                            after_stat(rr.label),
                            sep = "*\", \"*")),
                          label.x = 'right',
                          label.y='bottom',
                          size = 3.5)

  return(plt)
}

BP_plot <- function(df_plot,trait,Environment,lr){
  #plot breeding progress for each environment
  plt<- ggplot(df_plot %>% dplyr::filter(!is.na(RYear),!is.na(mean)),
               aes(x=RYear, y=mean)) +
    # theme_bw() +
    labs(x="Year of release",
         y=trait, title=Environment) + 
    geom_ribbon(aes(x = RYear, ymax = mean+sd,
                    ymin =mean-sd),fill="lightskyblue3",show.legend = F, 
                alpha = 0.4) +
    geom_point(size=1.5,shape=1) +
    stat_smooth(
      method="lm", se=FALSE, color="black",formula=y ~ x,linewidth=1) +
    toolPhD::theme_phd_main()+
    # theme(title = element_text(size=20),
    #       axis.title = element_text(size = 30), 
    #       axis.text = element_text(size=30))+
    ggpp::geom_text_npc(aes(npcx="right",npcy="bottom",
                            label =lm_eqn(lr)),parse = T,size=5)
  return(plt)
}

SlidingWindow <- function(FUN, data, window, step){
  # from EvobiR
  total <- length(data)
  spots <- seq(from = 1, to = (total - window + 1), by = step)
  vector(length = length(spots))
  result <-purrr::map_dbl(1:length(spots),function(i){
    match.fun(FUN)(data[spots[i]:(spots[i] + window - 1)])
  })
  # complete failure message
  if(!exists("result")) stop("Hmmm unknown error... Sorry")
  # 
  # return the result to the user
  return(result)
}

## select columns with combination of number and character
selec_col <- function(data, x, y){
  X <- as.data.frame(data[,x])
  Y <- as.data.frame(data[,y])
  if(length(y)==1)
  {names(Y) <- y}else{}
  out <- cbind(X, Y)
  return(out)
} 

### create equation in the graph
lm_eqn <- function(m){
  if (coef(m)[2] >= 0)
  {
    eq <-
      substitute(
        italic(y) == a + b %.%  italic(x) * "," ~  ~ italic(r) ^ 2 ~ "=" ~ r2,
        list(
          a = as.character(format(coef(m)[1], digits = 2)),
          b = as.character(format(coef(m)[2], digits = 2)),
          r2 = format(summary(m)$r.squared, digits = 2)
        )
      )
  } else{
    eq <-
      substitute(
        italic(y) == a - b %.%  italic(x) * "," ~  ~ italic(r) ^ 2 ~ "=" ~ r2,
        list(
          a = as.character(format(coef(m)[1], digits = 2)),
          b = as.character(format(abs(coef(m)[2]), digits = 2)),
          r2 = format(summary(m)$r.squared, digits = 2)
        )
      )
  }
  as.character(as.expression(eq));                 
}
Namean = function(x){mean(x, na.rm=TRUE)}
Nasd =  function(x){sd(x, na.rm=TRUE)}