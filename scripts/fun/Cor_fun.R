# 
# cor.mtest <- function(mat, ...) {
#   mat <- as.matrix(mat)
#   n <- ncol(mat)
#   p.mat<- matrix(NA, n, n)
#   diag(p.mat) <- 0
#   for (i in 1:(n - 1)) {
#     for (j in (i + 1):n) {
#       tmp <- cor.test(mat[, i], mat[, j], ...)
#       p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
#     }
#   }
#   colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
#   p.mat
# }

cut_paste <- function(v1,v2,cmb){
  # for aesthetic pasting of vector 1 and vector 2
  v1nchr <- purrr::map_dbl(v1,~{nchar(.x)}) 
  v1.pad <- max(v1nchr)-v1nchr
  v2nchr <- purrr::map_dbl(v2,~{nchar(.x)})
  v2.pad <- max(v2nchr)-v2nchr
  purrr::map_chr(1:length(v1),~{
    paste0(v1[.x],
           paste(rep(" ",v1.pad[.x]),collapse = ""),
           cmb,v2[.x],
           paste(rep(" ",v2.pad[.x]),collapse = ""))
  })
}

out_fin <- function(vec){
  l <- mean(vec,na.rm=T)-2.5*sd(vec,na.rm=T)
  u <- mean(vec,na.rm=T)+2.5*sd(vec,na.rm=T)
  res <- rep(F,length(vec))
  res[vec<l|vec>u] <- T
  return(res)
}

# round_scale <-function(vec){
#   map_dbl(vec,~{
#     if(is.na(.x)){
#       .x
#     }else if (.x<.01) {
#       round(.x,3)
#     } else{
#       round(.x,2)
#     }
#   })
# }
# 
# sma_cor<- function(mat){
#   
#   nam_mat<- names(mat)
#   n <- ncol(mat)
#   names(mat) <- paste0("V",1:n)
#   
#   p.mat<- matrix(NA, n, n)
#   r.mat<- matrix(NA, n, n)
#   s.mat <- matrix(NA, n, n)
#   
#   diag(p.mat) <- 0
#   diag(r.mat) <- 1
#   diag(s.mat) <- 1
#   amat <- as.matrix(mat)
#   
#   for (i in 1:(n - 1)) {
#     for (j in (i + 1):n) {
#       
#       tmp <- smatr::sma(paste0(names(mat)[j],"~",names(mat)[i]), method = "SMA",
#                         data = mat%>% .[complete.cases(amat[,i],amat[,j]),]) %>% 
#         .$groupsummary
#       
#       p.mat[i, j] <- p.mat[j, i] <- tmp$pval
#       r.mat[i, j] <- r.mat[j, i] <- tmp$r2 %>% sqrt()*ifelse(tmp$Slope>0,1,-1)
#       s.mat[i, j] <- s.mat[j, i] <- tmp$Slope
#       
#     }
#   }
#   
#   
#   colnames(p.mat) <- rownames(p.mat) <- nam_mat
#   colnames(r.mat) <- rownames(r.mat) <- nam_mat
#   colnames(s.mat) <- rownames(s.mat) <- nam_mat
#   return(list(r.mat,p.mat,s.mat))
#   
# }
