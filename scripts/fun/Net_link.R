Net_link<- function(node0,
                    table,
                    thresh.pp,
                    node.condi=NULL){
  # for stable table 
  
  # node0: node dataframe including all the nodes
  # table : dataframe with each columns equals to one trait 
  # thresh : threshold to filter the link(correlation) below the threshold
  # node.condi : pass from rlang::quo(), a condition variable to filter the unwanted node
  
  # combine trait, stability indices and physiological parameters 
  uniform <- dplyr::select(table,-c("Genotype"))
  if(testit::has_warning(cor(uniform))){
    stop('uniform!')
  }
  # calculate corelation matrix
  cormatrx <-cor(uniform,use = 'complete.obs') 
  #copy and replace the lower triangle to NA
  up <- cormatrx
  up[lower.tri(up,diag = T)] <- NA
  # transform the upper triangle into linear by row
  line <- as.data.frame(na.omit(as.vector(t(up))))
  colnames(line) <- 'weight'
  # label correlation coefficient by sign
  line$sign <- ifelse(line$weight>0,1,-1)
  #generate from-to index 
  mat <- as.data.frame(t(combn(dim(node0)[1],2)))
  colnames(mat) <- c('from','to')
  # combine from-to index and correlation 
  link <- cbind(mat,line)%>%
    dplyr::mutate(
      # from link
      fl=node0[match(from,node0$Id),"group"],
      # to link
      tl= node0[match(to,node0$Id),"group"],
      # mix from to label
      group = paste0(fl,tl),
      # turn the weight into positive
      weight = abs(weight),
      line=weight*sign,
      fs=node0[match(from,node0$Id),"SI"],
      ts=node0[match(to,node0$Id),"SI"],
      SI=paste0(fs,ts)
    )
  # only select the corresponding relation ship between trait and stability indices
  if (!is.null(node.condi)){
    sub.link <- dplyr::filter(link,
                              !!node.condi)
  }else{sub.link <- link}
  
  condi.vec <- with(sub.link,which(weight>thresh.pp))
  sub.link <- sub.link[condi.vec,]
  sub.link%<>%mutate(group=case_when(grepl("YB",group)~"BY",
                                     grepl("YC",group)~"CY",
                                     grepl("YN",group)~"NY",
                                     T~group#stay the same if not match
  ),
  SI=case_when(grepl("21",SI)~"12",
               T~SI#stay the same if not match
  )
  )%>%dplyr::select(-c(fl,tl,fs,ts))
  
  return(sub.link)
}


Net_link_cor<- function(node0,
                        group_column,
                        table,
                        thresh.pp,
                        node.condi=NULL){
  # for trait correlation
  # node0: node dataframe including all the nodes
  # grouop_column: group column to excluded for the cor()
  # table : dataframe with each columns equals to one trait 
  # thresh : threshold to filter the link(correlation) below the threshold
  # node.condi : pass from rlang::quo(), a condition variable to filter the unwanted node
  
  nodename <- node0$ori
  # remove complete NA columns
  df <- table[,table %>% sapply(., function(x)all(!is.na(x)))] 
  
  # combine trait, stability indices and physiological parameters 
  # remove non-numeric column
  uniform <- df%>% dplyr::select(-!!group_column) 
  # check the uniformity between columns to avoid duplicates
  if(testit::has_warning(cor(uniform))){
    stop('uniform!')
  }
  # calculate corelation matrix, avoid NA in the columns for each traits
  cormatrx <-cor(uniform,use="pairwise.complete.obs") 
  #copy and replace the lower triangle to NA
  up <- cormatrx
  up[lower.tri(up,diag = T)] <- NA
  # transform the upper triangle into linear by row
  line <- up %>% t() %>% as.vector() %>% na.omit() %>% as.data.frame()
  # as.data.frame(na.omit(as.vector(t(up)))) # oneliner version
  
  colnames(line) <- 'weight'
  # label correlation coefficient by sign
  line$sign <- ifelse(line$weight>0,1,-1)
  #generate from-to index 
  mat <-uniform %>%  names() %>%
    paste(.,collapse="|") %>% 
    # match exactly the provided string only
    #https://stackoverflow.com/questions/69258011/r-exact-match-for-multiple-patterns
    paste0("\\b(",.,")\\b") %>% 
    grep(.,nodename) %>%  combn(.,2) %>% t() %>% as.data.frame() 
  
  colnames(mat) <- c('from','to')
  # combine from-to index and correlation 
  link <- cbind(mat,line)%>%
    dplyr::mutate(
      # from link
      fl=node0[match(from,node0$Id),"group"],
      # to link
      tl= node0[match(to,node0$Id),"group"],
      # mix from to label
      group = paste0(fl,tl),
      # turn the weight into positive
      weight = abs(weight),
      line=weight*sign
    )
  # only select the corresponding relation ship between trait and stability indices
  if (!is.null(node.condi)){
    sub.link <- dplyr::filter(link,!!node.condi)
  }else{sub.link <- link}
  
  condi.vec <- with(sub.link,which(weight>thresh.pp))
  sub.link <- sub.link[condi.vec,]%>%
    mutate(group=case_when(grepl("YN",group)~"NY",
                           grepl("YC",group)~"CY",
                           grepl("YS",group)~"SY",
                           grepl("NS",group)~"SN",
                           grepl("CS",group)~"SC",
                           grepl("CN",group)~"NC",
                           T~group),
           g=table[[group_column]][1]#stay the same if not match)
    )
  
  return(sub.link)
}

signif<-function(vec){
  # significance transformer for display
  purrr::map_chr(vec,function(x){
    if(is.na(x)) {y=NA
    }else if(x<0.001){y="***"	
    }else if(x>=0.001&x<0.01){y="**"
    }else if(x>=0.01&x<0.05){y="*"
    }else if(x>=1){y="error"
    }else{y=""}
    # if(x>=0.05&x<0.1){y="."}
    # if(x>=0.1&x<1){y="ns"}
    return(y)
  })
}

get_triangle<- function(df){
  # assign lower triangle with a number larger than 1 
  df[lower.tri(df,diag = T)] <- 20
  # transform the upper triangle into linear by row
  line <- df %>% t() %>% as.vector() %>% .[!.==20]
  # .[!is.na(.)] 
  # as.data.frame(na.omit(as.vector(t(up)))) # oneliner version
  return(line)
}

Edge_df<- function(node0,
                   group_column,
                   table,
                   thresh.pp,
                   node.condi=NULL,
                   method="pearson"){
  # add pvalue for each pair of correlation using Hmisc::rcorr()
  
  # Input
  # node0: node dataframe including all the nodes
  # grouop_column: group column to excluded for the cor()
  # table : dataframe with each columns equals to one trait 
  # thresh : threshold to filter the link(correlation) below the threshold
  # node.condi : pass from rlang::quo(), a condition variable to filter the unwanted node
  
  # Output
  # sublink: long format edge dataframe
  
  nodename <- node0$ori
  # remove complete NA columns
  df <- table[,table %>% sapply(., function(x) !all(is.na(x)))] 
  
  # combine trait, stability indices and physiological parameters 
  # remove non-numeric column
  uniform <- df%>% dplyr::select(-!!group_column) 
  # check the uniformity between columns to avoid duplicates
  if(testit::has_warning(cor(uniform))){
    stop('uniform!')
  }
  # only this function can provide the p-value of each 
  pmatrx <-Hmisc::rcorr(as.matrix(uniform),type=method)
  
  line <- get_triangle(pmatrx$r)%>% 
    data.frame(weight=.) %>% 
    mutate(
      # label correlation coefficient by sign
      sign = ifelse(weight>0,1,-1),
      p=get_triangle(pmatrx$P),
      plabel=signif(p),
      psig=ifelse(p<.05,1,-1),
    )
  #generate from-to index 
  mat <-uniform %>%  names() %>%
    paste(.,collapse="|") %>% 
    # match exactly the provided string only
    #https://stackoverflow.com/questions/69258011/r-exact-match-for-multiple-patterns
    paste0("\\b(",.,")\\b") %>% 
    grep(.,nodename) %>% 
    combn(.,2) %>% t() %>% as.data.frame() 
  
  colnames(mat) <- c('from','to')
  # combine from-to index and correlation 
  link <- cbind(mat,line)%>%
    dplyr::mutate(
      line=weight,
      weight = abs(weight)
    ) %>% dplyr::filter(!from==to)
  # only select the corresponding relation ship between trait and stability indices
  if (!is.null(node.condi)){
    sub.link <- dplyr::filter(link,!!node.condi)
  }else{sub.link <-link }
  
  condi.vec <- with(sub.link,which(weight>thresh.pp))
  sub.link <- sub.link[condi.vec,]%>%
    mutate(
      g=table[[group_column]][1],#stay the same if not match),
      method=method
    )
  
  return(sub.link)
}


tar_unify <- function(df,tar_char){
  # tar may appear either in from or to 
  # force it move to "to"
  # Input
  #df dataframe
  #tar_char  character of name of target
  
  df %>% 
    dplyr::filter(from==tar_char|to==tar_char) %>% 
    dplyr::mutate(from=case_when(from==tar_char~to,
                                 T~from),
                  to=tar_char)
}