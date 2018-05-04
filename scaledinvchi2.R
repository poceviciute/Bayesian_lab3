#n = number of draws
#df = degrees of freedom
#scale = scale
#If X~scale-inv-chi2(v,s) then (X/(v*s))~inv-chi2(v) and (v*s)/X ~ chi2(v)

scaled_inv_chi2 <- function(n, df, scale){
  df*scale/rchisq(n, df)
}