
lsm<- function(model){
  # extract lsmeans from lmer model
  # rg    <-  lsmeans::ref.grid(model,pbkrtest.limit = 99999)
  # m.lsm <-  lsmeans::lsmeans(rg, pairwise ~ BRISONr, adjust = "tukey")
  # res   <-  data.table::as.data.table(m.lsm$lsmeans) 
  res <- data.table::as.data.table(summary( emmeans::emmeans(model,specs=~BRISONr,mode="asymptotic")))
  return(res)
}

add_info<- function(lstb,df){
  # add dataframe info for blue results
  sdf <- df %>% 
    dplyr::select(Location,Treatment,Year,trait) %>% 
    dplyr::distinct()
  
  cbind(lstb,sdf[rep(1,nrow(lstb)),])
}

rm_outlier <- function(x){
  x_m <- mean(x, na.rm=TRUE)
  x_4sd <- 4*sd(x, na.rm=TRUE)
  x <- ifelse(x >x_m+x_4sd| x <x_m-x_4sd,
              outlier_value, x) 
  return(x)
}

combiname <- function(x){
  # find name of trt combination
  x %>%
    ungroup() %>% 
    dplyr::select(Location,Treatment,Year,trait) %>% 
    .[1,] %>%unlist() %>%
    paste(.,collapse = "_")
} 

Blues <- function(i,df){
  # wrapper of three functions
  model <-  tryCatch(
    {
        lme4::lmer( Trait~ BRISONr+ (1|Row) +(1|Column),data = df,
                    control = lmerControl(
                      calc.derivs = FALSE,
                      optimizer = "bobyqa"
                    ) )
    },error=function(cond){
      print(i);print(cond)
      return("model failed")
    })
  if(!is.character(model)){
    lt <- tryCatch({
      lsm(model)
    },error=function(cond){
      return("lsm failed")
    })
    if(!is.character(lt)){
      res <- tryCatch({add_info(lt,df)},
                      error=function(cond){
                        return(NULL)
                      })
    }else{res <- lt}
  }else{res <- model}
  return(res)
}
