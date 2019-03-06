
# Clopper-Pearson Exact Interval (Lower Bound)
exact_two_sided_lb = function(x,n,alpha=0.05){
  
  # Beta distribution formulation
  return(qbeta((alpha/2),x,n-x+1))

  # Numeric approximation using the binomial cumulative distribution function
  # acc = 0.001
  # p = seq(0,1,acc)
  # return(p[which((1-pbinom(x-1,n,p))>(alpha/2))[1]-1])
  
  # library(binom)
  # return(subset(binom.confint(x,n,conf.level=1-alpha), method=='exact')$lower)
}

# Wilson Score Lower Bound
wilson_score_lb = function(x,n,alpha=0.05){
  z = qnorm(1-(alpha/2), 0, 1)
  return(((x+z*z/2)/(n + z*z)) - (z/(n + z*z))*sqrt((x*(n-x)/n) + z*z/4))

  # library(binom)
  # return(subset(binom.confint(x,n,conf.level=1-alpha), method=='wilson')$lower)
}

# Continuity Corrected Wilson Score Lower Bound
corrected_wilson_score_lb = function(x,n,alpha=0.05){
  z = qnorm(1-(alpha/2), 0, 1)
  lb = (2*n*p + z*z - z*sqrt(z*z - (1/n) + 4*n*p*(1-p) + 4*p - 2) - 1)/(2*(n + z*z))
  lb[lb < 0] = 0
  return(lb)
}

# Agresti-Coull
agresti_coull_lb = function(x,n,alpha=0.05){
  z = qnorm(1-(alpha/2), 0, 1)
  ptilde = (x+z*z/2)/(n + z*z)
  ntilde = n + z*z
  return(ptilde-z*sqrt(ptilde*(1-ptilde)/ntilde))
  
  # library(binom)
  # return(subset(binom.confint(164,200,conf.level=0.95), method=='agresti-coull')$lower)
}

# Normal Approximation (Wald)
wald_lb = function(x,n,alpha=0.05){
  z = qnorm(1-(alpha/2), 0, 1)
  p = x/n
  return(p-z*sqrt(p*(1-p)/n))  
}


# Adjusted Wald
adj_wald_lb = function(x,n,alpha=0.05){
  x = x + 2
  n = n + 4
  z = qnorm(1-(alpha/2), 0, 1)
  p = x/n
  return(p-z*sqrt(p*(1-p)/n))  
}


# Jeffreys
jeffreys_lb = function(x,n,alpha=0.05){
  return(qbeta((alpha/2),x + 0.5,n-x+0.5))
  
  # library(binom)
  # return(subset(binom.confint(164,200,conf.level=0.95), method=='agresti-coull')$lower)
}



df = data.frame(n = 1:1000)
for(i in 1:nrow(df)){
  n = df[i, 'n']
  x = 0:n
  p = x/n
  
  # Exact (Two-sided)
  df[i, 'Exact'] = sum(dbinom(x[exact_two_sided_lb(x,n) > 0.96],n,0.99))
  
  # Wilson Score
  df[i, 'Wilson'] = sum(dbinom(x[wilson_score_lb(x,n) > 0.96],n,0.99))

  # Continuity Corrected Wilson Score
  df[i, 'cWilson'] = sum(dbinom(x[corrected_wilson_score_lb(x,n) > 0.96],n,0.99))

  # Agresti-Coull
  df[i, 'Agresti-Coull'] = sum(dbinom(x[agresti_coull_lb(x,n) > 0.96],n,0.99))

  # Wald
  df[i, 'Wald'] = sum(dbinom(x[wald_lb(x,n) > 0.96],n,0.99))
  
  # Adjusted Wald
  df[i, 'adjWald'] = sum(dbinom(x[adj_wald_lb(x,n) > 0.96],n,0.99))
  
  # Jeffreys
  df[i, 'Jeffreys'] = sum(dbinom(x[jeffreys_lb(x,n) > 0.96],n,0.99))
  
}

# Melt df
df.m = melt(df, id.vars = 'n', variable.name = "Estimate", value.name = "Probability")


# # Exploring Results
# 
# library(ggplot2)
# library(reshape2)
# 
# default_plot = function(df){
#   p = ggplot(df, aes(x = n, y = Probability))
#   p = p + geom_point(alpha=0.4)
#   p = p + geom_hline(aes(yintercept=0.9))
#   p = p + xlab("Sample Size")
#   p = p + ylab("Probability that Lower Bound Estimate is > 0.96")
#   p
# }
# 
# 
# # Exact
est = 'Exact'
df.tmp = subset(df.m, Estimate == est)
df.tmp = df.tmp[order(df.tmp$n), ]
n = nrow(df.tmp)
p = df.tmp$Probability
df.exact.low = df.tmp[which(p[1:(n-1)] < p[2:n]), ]
df.exact.high = df.tmp[which(p[1:(n-1)] < p[2:n])+1, ]
# 
# p = default_plot(df.tmp)
# p
# p = p + geom_line(data = df.exact.low)
# p = p + geom_line(data = df.exact.high)
# p
# 
# # Lower Bound
# p = default_plot(df.exact.low) + geom_line()
# p
# p = p + xlim(c(250,450))
# p = p + xlim(c(320,370))
# p
# 
# 
# 
# # Wilson
est = 'Wilson'
df.tmp = subset(df.m, Estimate == est)
df.tmp = df.tmp[order(df.tmp$n), ]
n = nrow(df.tmp)
p = df.tmp$Probability
df.wilson.low = df.tmp[which(p[1:(n-1)] < p[2:n]), ]
df.wilson.high = df.tmp[which(p[1:(n-1)] < p[2:n])+1, ]
# 
# p = default_plot(df.tmp)
# p
# p = p + geom_line(data = df.wilson.low)
# p = p + geom_line(data = df.wilson.high)
# p
# 
# # Lower Bound
# p = default_plot(df.wilson.low) + geom_line()
# p
# p = p + xlim(c(250,450))
# p = p + xlim(c(320,370))
# p
# 
# 
# 
# # Corrected Wilson
est = 'cWilson'
df.tmp = subset(df.m, Estimate == est)
df.tmp = df.tmp[order(df.tmp$n), ]
n = nrow(df.tmp)
p = df.tmp$Probability
df.cwilson.low = df.tmp[which(p[1:(n-1)] < p[2:n]), ]
df.cwilson.high = df.tmp[which(p[1:(n-1)] < p[2:n])+1, ]
# 
# p = default_plot(df.tmp)
# p
# p = p + geom_line(data = df.cwilson.low)
# p = p + geom_line(data = df.cwilson.high)
# p
# 
# # Lower Bound
# p = default_plot(df.cwilson.low) + geom_line()
# p
# p = p + xlim(c(300,500))
# p = p + xlim(c(330,380))
# p
# 
# # Agresti-Coull
est = 'Agresti-Coull'
df.tmp = subset(df.m, Estimate == est)
df.tmp = df.tmp[order(df.tmp$n), ]
n = nrow(df.tmp)
p = df.tmp$Probability
df.agresti.low = df.tmp[which(p[1:(n-1)] < p[2:n]), ]
df.agresti.high = df.tmp[which(p[1:(n-1)] < p[2:n])+1, ]
# 
# p = default_plot(df.tmp)
# p
# p = p + geom_line(data = df.agresti.low)
# p = p + geom_line(data = df.agresti.high)
# p
# 
# # Lower Bound
# p = default_plot(df.agresti.low) + geom_line()
# p
# p = p + xlim(c(300,500))
# p = p + xlim(c(320,370))
# p
# 
# 
# # Wald
est = 'Wald'
df.tmp = subset(df.m, Estimate == est)
df.tmp = df.tmp[order(df.tmp$n), ]
n = nrow(df.tmp)
p = df.tmp$Probability
df.wald.low = df.tmp[which(p[1:(n-1)] < p[2:n]), ]
df.wald.high = df.tmp[which(p[1:(n-1)] < p[2:n])+1, ]
# 
# p = default_plot(df.tmp)
# p
# p = p + geom_line(data = df.wald.low)
# p = p + geom_line(data = df.wald.high)
# p
# 
# # Lower Bound
# p = default_plot(df.wald.low) + geom_line()
# p
# # p = p + xlim(c(300,500))
# # p = p + xlim(c(320,370))
# p
# 
# 
# # Adjusted Wald
est = 'adjWald'
df.tmp = subset(df.m, Estimate == est)
df.tmp = df.tmp[order(df.tmp$n), ]
n = nrow(df.tmp)
p = df.tmp$Probability
df.adjwald.low = df.tmp[which(p[1:(n-1)] < p[2:n]), ]
df.adjwald.high = df.tmp[which(p[1:(n-1)] < p[2:n])+1, ]
# 
# p = default_plot(df.tmp)
# p
# p = p + geom_line(data = df.adjwald.low)
# p = p + geom_line(data = df.adjwald.high)
# p
# 
# # Lower Bound
# p = default_plot(df.adjwald.low) + geom_line()
# p
# # p = p + xlim(c(300,500))
# # p = p + xlim(c(320,370))
# p
# 
# 
# 
# 
# # Jeffreys
est = 'Jeffreys'
df.tmp = subset(df.m, Estimate == est)
df.tmp = df.tmp[order(df.tmp$n), ]
n = nrow(df.tmp)
p = df.tmp$Probability
df.jeffreys.low = df.tmp[which(p[1:(n-1)] < p[2:n]), ]
df.jeffreys.high = df.tmp[which(p[1:(n-1)] < p[2:n])+1, ]
# 
# p = default_plot(df.tmp)
# p
# p = p + geom_line(data = df.jeffreys.low)
# p = p + geom_line(data = df.jeffreys.high)
# p
# 
# # Lower Bound
# p = default_plot(df.jeffreys.low) + geom_line()
# p
# # p = p + xlim(c(300,500))
# # p = p + xlim(c(320,370))
# p


