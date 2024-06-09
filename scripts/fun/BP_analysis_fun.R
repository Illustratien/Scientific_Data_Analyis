col_pal <- ggsci::pal_locuszoom()(4)
names(col_pal) <- c("HN_NF","HN_WF","LN_NF","all")
x1 <- c(0,.5,0,0.5)
x2 <- c(0.5,1,0.5,1)
y1 <- c(0.5,.5,0,0)
y2 <- c(1,1,0.5,0.5)
v <- c(T,T,T,T)
thresh <- 0

add_label <- function(xfrac, yfrac, label, pos = 4, ...) {
  # for adding text in the graph *from web
  # https://seananderson.ca/2013/10/21/panel-letters/
  u <- par("usr")
  x <- u[1] + xfrac * (u[2] - u[1])
  y <- u[4] - yfrac * (u[4] - u[3])
  text(x, y, label, pos = pos,font=2,cex=1,...)
}

# color RampPalette always start from low value to high value
# for <0 blue to white
colfunc1 <- colorRampPalette(c('#1167b1',"#99efff"))(4)
# for >0 orange color start from white to orange
colfunc2<- colorRampPalette(c("#ffaf7a",'#c83200'))(4) 
# put the white in the middle
# from dark orange to dark blue
colfunc <- c(colfunc1,"#FFFFFF",colfunc2) %>%  rev(.)
# "#C83200" "#DA5B28" "#EC8551" "#FFAF7A" "#FFFFFF" "#99EFFF""#6BC1E4" "#3E94CB" "#1167B1"
# -------------------------------------------------------------------------

lm_prepare<- function(df,xname,yname){
  relative_col <- names(df)[grep(paste(xname,yname,sep="|"),names(df))]
  
  df <- df %>% dplyr::select(all_of(c("Treatment",relative_col))) %>% na.omit()
  if(nrow(df)>0){
    ## add se line. 
    combine_df <-rbind(df,df %>% mutate(Treatment='all')) %>% 
      rbind(.,df %>%filter(grepl("HN",Treatment)) %>%  mutate(Treatment='HN')) %>%
      rbind(.,df %>%filter(grepl("NF",Treatment)) %>%  mutate(Treatment='NF'))
    group_dat <- combine_df %>% group_by(Treatment) %>% group_split()
    
    lm.fit <- tryCatch({
      group_dat %>%  
        imap_dfr(.,~{
          # print(.y)
          # regression formular 
          formu <- paste0(yname,'~',xname) 
          fit <- lm(formu,data = .x,na.action=na.exclude)
          fs <- fit %>%summary()
          
          if(nrow(.x)<=2){
            data.frame( pvalue=FALSE,
                        r2= NA,
                        r=NA,
                        slope=NA,
                        slope_se=NA,
                        npoints=NA,
                        from=xname,
                        to=yname,
                        trait=paste(xname,yname,sep="_"),
                        Treatment=.x[["Treatment"]][1])
          }else{
            .x %>% summarise(
              pvalue=fs%>% coef() %>% .[xname,4]<0.05,
              r2= fs%>% .$r.squared,
              r=cor(.x[[xname]],.x[[yname]]),
              slope=fs$coefficients[2,1],
              slope_se=fs$coef[2,2],
              npoints=n(),
              from=xname,
              to=yname,
              trait=paste(xname,yname,sep="_"),
              Treatment=.x[["Treatment"]][1])
          }
          
          
        })})
    combine_df <- left_join(combine_df,lm.fit,"Treatment")
  }else{
    combine_df <- NULL
    lm.fit <- NULL
  }
  
  return(list(combine_df,lm.fit))
}
RBP_plot<- function(lmdf,xname,yname){
  # plot breeding progress 
  # xname x column name label
  # yname y column name labeil
  # add all group
  
  df <- lmdf %>% filter(Treatment%in%c("HN_NF","HN_WF","LN_NF"))
  xy.limits <- range( c(lmdf[[xname]],lmdf[[yname]]) )
  
  p <- ggplot(data=lmdf,
              aes_string(xname,yname,color='Treatment'))+
    geom_point(data=df,aes_string(xname,yname,color='Treatment'),size=2)+
    coord_fixed(ratio = 1)+
    scale_x_continuous(limits = xy.limits)+
    scale_y_continuous(limits = xy.limits)+
    geom_abline(slope=1,linetype=2,color="darkgrey",alpha=.8)+
    geom_hline(yintercept=0,linetype=2,color="darkgrey",alpha=.8)+
    geom_vline(xintercept=0,linetype=2,color="darkgrey",alpha=.8)+
    scale_color_manual(values=col_pal)+
    theme_test()+
    theme(legend.position = 'bottom')+
    labs(x=paste0('relative gentic gain of ',xname),
         y= paste0('relative gentic gain of ',yname))
  
  if(nrow(lmdf %>% filter(pvalue==T))>0){
    
    d <- lmdf%>% filter(pvalue==T,!Treatment%in%c("HN","NF"))
    p <- p+
      geom_line(data=d,stat="smooth",method = "lm", formula ='y~x',alpha = 0.5)+
      ggpmisc::stat_poly_eq(data=d,mapping=aes(
        label = paste(after_stat(eq.label),
                      after_stat(rr.label), sep = "*\", \"*")),
        label.x = 'right',
        label.y='bottom')
  }
  
  
  return(p)
}


ABP_plot<- function(lmdf,xname,yname){
  # plot breeding progress 
  # xname x column name label
  # yname y column name labeil
  # add all group
  
  df <- lmdf %>% filter(Treatment%in%c("HN_NF","HN_WF","LN_NF"))
  p <- ggplot(data=lmdf,
              aes_string(xname,yname,color='Treatment'))+
    geom_point(data=df,aes_string(xname,yname,color='Treatment'),size=2)+
    geom_hline(yintercept=0,linetype=2,color="darkgrey",alpha=.8)+
    geom_vline(xintercept=0,linetype=2,color="darkgrey",alpha=.8)+
    geom_errorbar(data=df,alpha=.3,mapping=aes_string(ymin = paste0('BPmin_',yname),ymax = paste0('BPmax_',yname))) + 
    geom_errorbarh(data=df,alpha=.3,mapping=aes_string(xmin = paste0('BPmin_',xname),xmax = paste0('BPmax_',xname)))+
    scale_color_manual(values=col_pal)+
    theme_test()+
    theme(legend.position = 'bottom')+
    labs(x=paste0('absolute gentic gain of ',xname),
         y= paste0('absolute gentic gain of ',yname))
  
  if(nrow(lmdf %>% filter(pvalue==T))>0){
    
    d <- lmdf%>% filter(pvalue==T,!Treatment%in%c("HN","NF"))
    p <- p+
      geom_line(data=d,stat="smooth",method = "lm", formula ='y~x',alpha = 0.5)+
      ggpmisc::stat_poly_eq(data=d,mapping=aes(
        label = paste(after_stat(eq.label),
                      after_stat(rr.label), sep = "*\", \"*")),
        label.x = 'right',
        label.y='bottom')
  }
  
  
  return(p)
}

pl<- list(c("HN_NF","HN_WF","LN_NF","all"),
     c("HN","NF","all"))

plot_net <- function(link,node,fname,labdf=NULL){
  for(pid in 1:length(pl)){
    
    # output  -------------------------------------------------------------------------
    tiff(filename=paste0(fname,"_",pid,".tiff"),
         units="cm",
         width=17.4,
         height=11.6,
         compression = "lzw",
         pointsize=12,
         res=500,# dpi,
         family="Arial"
    )
    plot.new()
    for ( i in 1:length(pl[[pid]])){
      if (i>1){  par(fig=c(x1[i], x2[i], y1[i],y2[i]),
                     mai = c(0.1, 0.1, 0.1, 0.1),
                     oma = c(.01, .01, 0.01, 0.01),
                     lwd=2,new=v[i])
      } else{  par(fig=c(x1[i], x2[i], y1[i],y2[i]),
                   mai = c(0.1, 0.1, 0.1, 0.1),
                   oma = c(.01, .01, 0.01, 0.01),
                   lwd=2)}
      sublink <- link %>% filter(Treatment==pl[[pid]][i],weight>thresh)
      subnodeid <- sublink%>% dplyr::select(from,to) %>% distinct() %>% unlist() %>% unique() 
      subnode <- node %>% filter(Id%in%subnodeid)
      if(!is.null(labdf)){
        subnode <- left_join(subnode,labdf %>% filter(Treatment==pl[[pid]][i]),by="Name") %>% 
          mutate(Name=paste0(Name,"\n(",m,")"))
      }
      
      net1 <- igraph::graph.data.frame(sublink,subnode, directed=F) 
      # igraph::V(net1)$color <-col_pal[3]
      igraph::V(net1)$frame.color <- NA
      igraph::E(net1)$width <- E(net1)$weight*3
      # igraph::E(net1)$width <- E(net1)$weight*8
      set.seed(777)   ## to make this reproducable
      l <- layout.fruchterman.reingold(net1)
      ii <- cut(sublink$r, breaks = c(seq(-1,-.33,len = 5),seq(.33, 1,len = 5)),
                include.lowest = TRUE)
      # remove unwanted text generated by cut()
      breaks <- c(-1,as.numeric(sub("[^,]*,([^]]*)\\]", "\\1", levels(ii))))
      ## Use bin indices, ii, to select color from vector of n-1 equally spaced colors
      edge_colors_trait<- colfunc[ii]
      V(net1)$label=V(net1)$Name %>% gsub("_","\n",.)
      E(net1)$label=E(net1)$npoints
      plot(net1,
           vertex.label.cex=0.31,
           vertex.label.font=2,
           vertex.shape="none",
           vertex.size=22,

           vertex.label.color='black',
           vertex.label.family='Arial',
           edge.color=edge_colors_trait,
           edge.label.color='black',
           edge.label.cex=.25,
           layout=l)
      
      add_label(0.02, 0.07, pl[[pid]][i])
      if(i==3){
        image.plot( legend.only=TRUE, zlim=c(-1,1),
                    col = colfunc,
                    breaks = breaks,
                    axis.args=list(cex.axis=.5,tck=1),
                    lab.breaks=breaks %>% round(.,2),
                    # addjust margin
                    smallplot=c(.05,.06, .08,.38)
        )
        # text(0.3,-.71,bquote(bolditalic(r)),cex = 1)
      }
    }
    
    
    dev.off()
  }
  
  
}
