
df = data.frame(n = 1:1000)

# Maximum Likelihood Estimate
df = transform(df, MLE = 1-pbinom(0.98*n, n, 0.99))
# phat = x/n

# # Maximum Information (Shannon Entropy) Estimate
# MIE = function(x,n,acc=0.001){
#   p = seq(0,1,acc)
#   tmp = log(p/(1-p)) - (x/p) + (n-x)/(1-p)
#   return(p[which.min(abs(tmp))])
# }
# df$MIE = NA
# for (i in 1:nrow(df)){
#   n = df[i, 'n']
#   x = 0:n
#   phat = sapply(x, function(x){return(MIE(x,n))})
#   df[i, 'MIE'] = 1-pbinom(x[which(phat > 0.98)[1]-1],n,0.99)
# }

# Mean of the posterior given uniform prior Beta(0,0)
# yields Laplace's Law of Succession ("adding one success and one failure")
df = transform(df, Laplace = 1-pbinom(0.98*n + 0.96, n, 0.99))
# phat = (x+1)/(n+2)

# Mean of the posterior given Jeffreys prior Beta(-0.5,-0.5) or Beta(x,n-x)
df = transform(df, Jeffreys = 1-pbinom(0.98*n + 0.48, n, 0.99))
# phat = (x+0.5)/(n+1)

## Wilson Estimator
z = qnorm(0.975,0,1)
df = transform(df, Wilson = 1-pbinom((0.98*(n+z*z) - (z*z/2)), n, 0.99))

# Mean of the posterior with prior Beta(2,2)
# Also called Bayes estimate, and also the centre of the adjusted Wald CI)
z = qnorm(0.975,0,1)
df = transform(df, Bayes = 1-pbinom(0.98*n + 1.92, n, 0.99))

# This estimate appears to be dramatically more conservative than the other estimates, requiring a sample size on the order of n = 4000
# # Game Theoretic Steinhaus estimate
# df = transform(df, Steinhaus = 1-pbinom(0.98*(n + sqrt(n)) - sqrt(n)/2, n, 0.99))
# # phat = (x+sqrt(n)/2)/(n+sqrt(n))


# Melt df
df.m = melt(df, id.vars = 'n', variable.name = "Estimate", value.name = "Probability")




# # Exploring Results
# 
# library(ggplot2)
# library(reshape2)
# default_plot = function(df){
#   p = ggplot(df, aes(x = n, y = Probability))
#   p = p + geom_point(alpha=0.4)
#   p = p + geom_hline(aes(yintercept=0.9))
#   p = p + xlab("Sample Size")
#   p = p + ylab("Probability that Point Estimate is > 0.98")
#   p
# }
# 
# MLE
est = 'MLE'
df.tmp = subset(df.m, Estimate == est)
df.tmp = df.tmp[order(df.tmp$n), ]
n = nrow(df.tmp)
p = df.tmp$Probability
df.mle.low = df.tmp[which(p[1:(n-1)] < p[2:n]), ]
# p = default_plot(df.tmp)
# p
# p = p + geom_line(data = df.mle.low)
# p
# 
# # Lower Bound
# p = default_plot(df.mle.low) + geom_line()
# p
# p = p + xlim(c(200,400))
# p = p + xlim(c(250,300))
# p
# 
# 
# 
# # Jeffreys
est = 'Jeffreys'
df.tmp = subset(df.m, Estimate == est)
# p = default_plot(df.tmp)
# p
# 
# # Lower Bound
df.tmp = df.tmp[order(df.tmp$n), ]
n = nrow(df.tmp)
p = df.tmp$Probability
df.jeffreys.low = df.tmp[which(p[1:(n-1)] < p[2:n]), ]
df.jeffreys.high = df.tmp[which(p[1:(n-1)] < p[2:n]) + 1, ]
# p = default_plot(df.jeffreys.low) + geom_line()
# p
# p = p + xlim(c(300,500))
# p = p + xlim(c(300,380))
# p
# 
# 
# # Laplace
est = 'Laplace'
df.tmp = subset(df.m, Estimate == est)
# p = default_plot(df.tmp)
# p
# 
# # Lower Bound
df.tmp = df.tmp[order(df.tmp$n), ]
n = nrow(df.tmp)
p = df.tmp$Probability
df.laplace.low = df.tmp[which(p[1:(n-1)] < p[2:n]), ]
df.laplace.high = df.tmp[which(p[1:(n-1)] < p[2:n]) + 1, ]
# p = default_plot(df.laplace.low) + geom_line()
# p
# p = p + xlim(c(300,500))
# p = p + xlim(c(380,450))
# p
# 
# # Wilson
est = 'Wilson'
df.tmp = subset(df.m, Estimate == est)
# p = default_plot(df.tmp)
# p
# 
# # Lower Bound
df.tmp = df.tmp[order(df.tmp$n), ]
n = nrow(df.tmp)
p = df.tmp$Probability
df.wilson.low = df.tmp[which(p[1:(n-1)] < p[2:n]), ]
df.wilson.high = df.tmp[which(p[1:(n-1)] < p[2:n]) + 1, ]
# p = default_plot(df.wilson.low) + geom_line()
# p
# p = p + xlim(c(450,650))
# p = p + xlim(c(470,550))
# p
# 
# 
# 
# # Bayes
est = 'Bayes'
df.tmp = subset(df.m, Estimate == est)
# p = default_plot(df.tmp)
# p
# 
# # Lower Bound
df.tmp = df.tmp[order(df.tmp$n), ]
n = nrow(df.tmp)
p = df.tmp$Probability
df.bayes.low = df.tmp[which(p[1:(n-1)] < p[2:n]), ]
df.bayes.high = df.tmp[which(p[1:(n-1)] < p[2:n]) + 1, ]
# p = default_plot(df.bayes.low) + geom_line()
# p
# p = p + xlim(c(450,650))
# p = p + xlim(c(470,550))
# p













