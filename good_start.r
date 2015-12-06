# LEVELS
# Blocks (NOT ORDINAL)
#   1, 2, 3, 4, 5, 6
# Rows (NOT ORDINAL)
#   1, 2,..., 30
# Plots (NOT ORDINAL)
#   1, 2,..., 180
# Varieties (NOT ORDINAL)
#   1, 2, 3, 4, 5
# Nitrogen (ORDERED)
#   1 (early), 2 (mid), 3 (late)
# Inoc (NOT ORDINAL)
#   1 (not inoculated), 2 (inoculated)

# ANOVA
# block                 5
# -----------------------
# var                   4
# block:var            20
# -----------------------
# nitro                 2
# inoc                  1
# nitro:inoc            2
# nitro:var             8
# nitro:block          10
# inoc:var              4
# inoc:block            5
# nitro:inoc:var        8
# nitro:inoc:block     10
# nitro:var:block      40
# inoc:var:block       20
# nitro:inoc:var:block 40

logit <- function(x){return(log(x/(1-x)))}
expit <- function(x){return(exp(x)/(1+exp(x)))}

design <- data.frame('block' = as.numeric(gl(6, 30)),
                     'row' = as.numeric(gl(30, 6)),
                     'plot' = 1:180)
design <- cbind(design,
                'variety' = numeric(180),
                'inoc' = numeric(180),
                'nitrogen' = numeric(180))
set.seed(6823)
# Each column is a permutation of variety indices
permute.v <- matrix(replicate(6, sample(5, 5, replace = FALSE)), ncol = 6)
for(i in 1:180){
  design$variety[i] <- permute.v[(design$row[i]-1)%%5+1, design$block[i]]
}
ni <- matrix(c(1, 1, 1, 2, 2, 2, 1:3, 1:3), ncol = 2)
# Each column is a permutation of inoc:nitrogen indices
permute.ni <- matrix(replicate(30, sample(6, 6, replace = FALSE)), ncol = 30)
for(i in 1:30){
  design[((i-1)*6+1):(i*6), c('inoc', 'nitrogen')] <- ni[permute.ni[,i],]
}
#design






# Let's simulate some stuff
mu <- -3 # Intercept: Because logit(0.1) is about -2.2

set.seed(93)
# six blocks
b.eff <- rnorm(6,0,0.01)
# five varieties
v.eff <- rnorm(5,0,0.3)
# three nitrogen levels
n.eff <- sort(rnorm(3,0,0.4))
# two inoculation levels
i.eff <- c(-0.5,0.5)
#### kenny's dumb interaction effects
## two-way interactions
# block and variety
bv.eff <- rnorm(30,0,0.00001)
# block and nitrogen
bn.eff <- rnorm(18,0,0.00001)
# block and inoculation
bi.eff <- rnorm(12,0,0.00001)
# variety and nitrogen
vn.eff <- rnorm(15,0,0.01)
# variety and inoculation
vi.eff <- rnorm(10,0,0.03)
# nitrogen and inoculation
ni.eff <- rnorm(6,0,0.03)
## three-way interactions
# block, variety, and nitrogen
bvn.eff <- rnorm(90,0,0.00001)
# block, variety, and inoculation
bvi.eff <- rnorm(60,0,0.00001)
# block, nitrogen, and inoculation
bni.eff <- rnorm(36,0,0.00001)
# variety, nitrogen, and inoculation
vni.eff <- rnorm(30,0,0.001)
## four-way interaction
# bvni.eff <- rnorm(180,0,0.000001) # really is residual error since we do not have replication

#### update design (really data) matrix with the levels of each interaction
head(design)
## 2s
bv.int <- as.numeric(interaction(design[,1],design[,4]))
bn.int <- as.numeric(interaction(design[,1],design[,6]))
bi.int <- as.numeric(interaction(design[,1],design[,5]))
vn.int <- as.numeric(interaction(design[,4],design[,6]))
vi.int <- as.numeric(interaction(design[,4],design[,5]))
ni.int <- as.numeric(interaction(design[,6],design[,5]))
## 3s
bvn.int <- as.numeric(interaction(design[,1],design[,4],design[,6]))
bvi.int <- as.numeric(interaction(design[,4],design[,6]))
bni.int <- as.numeric(interaction(design[,4],design[,5]))
vni.int <- as.numeric(interaction(design[,6],design[,5]))

## update data frame
od <- data.frame(
  "blk" = factor(design[,1]),
  "vty" = factor(design[,4]),
  "nit" = factor(design[,6]),
  "ino" = factor(design[,5]),
  "bv.int" = factor(bv.int),
  "bn.int" = factor(bn.int),
  "bi.int" = factor(bi.int),
  "vn.int" = factor(vn.int),
  "vi.int" = factor(vi.int),
  "ni.int" = factor(ni.int),
  "bvn.int" = factor(bvn.int),
  "bvi.int" = factor(bvi.int),
  "bni.int" = factor(bni.int),
  "vni.int" = factor(vni.int)
)

# generate probablities
logit.pi.vec <- mu + b.eff[od[,1]] + v.eff[od[,2]] + n.eff[od[,3]] + i.eff[od[,4]] +
  bv.eff[od[,5]] + bn.eff[od[,6]] + bi.eff[od[,7]] +
  vn.eff[od[,8]] + vi.eff[od[,9]] + ni.eff[od[,10]] +
  bvn.eff[od[,11]] + bvi.eff[od[,12]] + bni.eff[od[,13]] + vni.eff[od[,14]]

# transform back to pi scale
pi.vec <- apply(cbind(logit.pi.vec), 1, expit) # pi = expit(X * beta)
# simulate some y's (leaves - out of 30 - that are contaminated)
set.seed(324)
y.vec <- apply(cbind(pi.vec), 1, function(p) rbinom(1,30,p)) # y_i ~ Bin(30,pi_i)

# attach to data set
od$"infected" <- y.vec
od$"uninfected" <- 30 - y.vec
# write.csv(od,"od.csv",row.names = FALSE)
summary(od)
names(od)

# glmer
require(lme4)
# start with glm to see if it runs
glm.mod <- glmer(cbind(infected,uninfected) ~ (nit+ino)^2 + (1|blk/vty),
                 data=od, family=binomial(link = "logit"),
                 control = glmerControl(optimizer = "bobyqa", optCtrl = list("maxfun" = 2000)))
# glm.mod.dmb <- glm(cbind(infected,uninfected) ~ factor(blk)+factor(vty)+factor(nit)+factor(ino)
#                +factor(bv.int)+factor(bn.int)+factor(bi.int)+factor(vn.int)+factor(vi.int)+factor(ni.int)
#                +factor(bvn.int)+factor(bvi.int)+factor(bni.int)+factor(vni.int),
#                data=od, family=binomial(link = "logit"))
summary(glm.mod)
anova(glm.mod, test="Chisq")
# summary(glm.mod.dmb)
# anova(glm.mod.dmb)



#### create some separation
y2.vec <- rep(0,length(y.vec))
## start with inoc
set.seed(324)
ino.rbin.1 <- rbinom(length(y.vec)*mean(od$ino==1),1,p=0.8) ## if not inoc
ino.rbin.2 <- rbinom(length(y.vec)*mean(od$ino==2),1,p=0.95) ## if inoc
# Brandon's trick
y2.vec[od$ino==1][ino.rbin.1==1] <- y.vec[od$ino==1][ino.rbin.1==1]
y2.vec[od$ino==2][ino.rbin.2==1] <- y.vec[od$ino==2][ino.rbin.2==1]
## now nitrogen
set.seed(833)
nit.rbin.1 <- rbinom(length(y.vec)*mean(od$nit==1),1,p=0.8) ## if early
nit.rbin.2 <- rbinom(length(y.vec)*mean(od$nit==2),1,p=0.95) ## if middle
nit.rbin.3 <- rbinom(length(y.vec)*mean(od$nit==3),1,p=0.95) ## if late
# Brandon's Trick
y2.vec[od$nit==1][nit.rbin.1==0] <- 0
y2.vec[od$nit==2][nit.rbin.2==0] <- 0
y2.vec[od$nit==3][nit.rbin.3==0] <- 0
## now variety, 2 and 5 are gonna be resistant (low probability)
set.seed(2841)
var.rbin.1 <- rbinom(length(y.vec)*mean(od$vty==1),1,p=0.95) ##
var.rbin.2 <- rbinom(length(y.vec)*mean(od$vty==2),1,p=0.8) ## very resistant
var.rbin.3 <- rbinom(length(y.vec)*mean(od$vty==3),1,p=0.95) ## if late
var.rbin.4 <- rbinom(length(y.vec)*mean(od$vty==4),1,p=0.95) ## if early
var.rbin.5 <- rbinom(length(y.vec)*mean(od$vty==5),1,p=0.8) ## very resistant
# The trick
y2.vec[od$vty==1][var.rbin.1==0] <- 0
y2.vec[od$vty==2][var.rbin.2==0] <- 0
y2.vec[od$vty==3][var.rbin.3==0] <- 0
y2.vec[od$vty==4][var.rbin.4==0] <- 0
y2.vec[od$vty==5][var.rbin.5==0] <- 0
## take a look
plot(table(y2.vec))
names(od)

# attach to data set
od$"infected2" <- y2.vec
od$"uninfected2" <- 30 - y2.vec
# write.csv(od,"od.csv",row.names = FALSE)
summary(od)
names(od)

# start with glm to see if it runs
glm.mod2 <- glm(cbind(infected2,uninfected2) ~ (blk+vty+nit+ino)^4, data=od, family=binomial(link = "logit"))
summary(glm.mod2)
anova(glm.mod2, test="Chisq")
# require(effects)
# plot(Effect(focal.predictors = c("vty","nit","ino"), mod = glm.mod2))



#### STANislaw
require(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
## Doug's Working Directory
# setwd("C:/Users/dsand_000/Desktop/Stats/S532/gitProj/Bayes-Group")
### set up data to send in to STAN
# sample sizes
n.o <- 180
m <- 30
cauchyscale <- 3
n.b <- 6
n.v <- 5
n.n <- 3
n.i <- 2
n.bv <- n.b*n.v
n.bn <- n.b*n.n
n.bi <- n.b*n.i
n.vn <- n.v*n.n
n.vi <- n.v*n.i
n.ni <- n.n*n.i
n.bvn <- n.b*n.v*n.n
n.bvi <- n.b*n.v*n.i
n.bni <- n.b*n.n*n.i
n.vni <- n.v*n.n*n.i

# vectors of levels
b.vec <- od$blk
v.vec <- od$vty
n.vec <- od$nit
i.vec <- od$ino
bv.vec <- od$bv.int
bn.vec <- od$bn.int
bi.vec <- od$bi.int
vn.vec <- od$vn.int
vi.vec <- od$vi.int
ni.vec <- od$ni.int
bvn.vec <- od$bvn.int
bvi.vec <- od$bvi.int
bni.vec <- od$bni.int
vni.vec <- od$vni.int

# infected
y <- od$infected




#### not at all useless
od.stan.data <- list(num_obs=n.o,
                     m=m,
                     cauchysd=cauchyscale,
                     num_b=n.b,
                     num_v=n.v,
                     num_n=n.n,
                     num_i=n.i,
                     num_bv=n.bv,
                     num_bn=n.bn,
                     num_bi=n.bi,
                     num_vn=n.vn,
                     num_vi=n.vi,
                     num_ni=n.ni,
                     num_bvn=n.bvn,
                     num_bvi=n.bvi,
                     num_bni=n.bni,
                     num_vni=n.vni,
                     y=y,
                     blk=as.numeric(b.vec),
                     vty=as.numeric(v.vec),
                     nit=as.numeric(n.vec),
                     ino=as.numeric(i.vec),
                     bv_int=as.numeric(bv.vec),
                     bn_int=as.numeric(bn.vec),
                     bi_int=as.numeric(bi.vec),
                     vn_int=as.numeric(vn.vec),
                     vi_int=as.numeric(vi.vec),
                     ni_int=as.numeric(ni.vec),
                     bvn_int=as.numeric(bvn.vec),
                     bvi_int=as.numeric(bvi.vec),
                     bni_int=as.numeric(bni.vec),
                     vni_int=as.numeric(vni.vec))
#stan_herpes <- stan_model(file = "./your_mom.stan", model_name = "ni")
#require(rstudioapi)
#options(mc.cores = 1)
#your_sister <- sampling(stan_herpes, chains = 4, iter = 500, data = weiner_logs, sample_file="trial")
# t1 <- read_stan_csv("trial_1.csv", col_major = TRUE)
# t2 <- read_stan_csv("trial_2.csv", col_major = TRUE)
# t3 <- read_stan_csv("trial_3.csv", col_major = TRUE)
# t4 <- read_stan_csv("trial_4.csv", col_major = TRUE)
trial5000 <- read_stan_csv(c("trial5000_1.csv","trial5000_2.csv","trial5000_3.csv","trial5000_4.csv"), col_major=TRUE)
# trial <- read_stan_csv(c("trial_1.csv","trial_2.csv","trial_3.csv","trial_4.csv"), col_major=TRUE)
# plot(your_sister, pars=c("sigma_b", "sigma_v", "sigma_n", "sigma_i", "sigma_bv", "sigma_bn","sigma_bi", "sigma_vn", "sigma_vi", "sigma_ni", "sigma_bvn", "sigma_bvi", "sigma_bni", "sigma_vni"), ci_level=0.5, outer_level=0.95, point_est="median")
plot(trial5000, pars=c("sigma_b", "sigma_v", "sigma_n", "sigma_i", "sigma_bv", "sigma_bn","sigma_bi", "sigma_vn", "sigma_vi", "sigma_ni", "sigma_bvn", "sigma_bvi", "sigma_bni", "sigma_vni"), ci_level=0.5, outer_level=0.95, point_est="mean")

od.cauchy  <- stan_model(file = 'cauchy.stan', model_name = 'cauchyprior')
od.cauchy.out <- sampling(od.cauchy, chains = 4, iter = 500, data = od.stan.data, sample_file = 'trial')

plot(od.cauchy.out, pars=c("sigma_b", "sigma_v", "sigma_n", "sigma_i", "sigma_bv", "sigma_bn","sigma_bi", "sigma_vn", "sigma_vi", "sigma_ni", "sigma_bvn", "sigma_bvi", "sigma_bni", "sigma_vni"), ci_level=0.5, outer_level=0.95, point_est="mean")

# A simple sensitivity analysis
scales <- c(1, 2, 3, 5, 7, 11)
sens <- lapply(scales, function(scale){
    stan.data <- list(num_obs = n.o, m = m, cauchysd = scale,
                      num_b = n.b, num_v = n.v, num_n = n.n, num_i = n.i,
                      num_bv = n.bv, num_bn = n.bn, num_bi = n.bi,
                      num_vn = n.vn, num_vi = n.vi, num_ni = n.ni,
                      num_bvn = n.bvn, num_bvi = n.bvi, num_bni = n.bni,
                      num_vni = n.vni, y = y, blk = as.numeric(b.vec),
                      vty = as.numeric(v.vec), nit = as.numeric(n.vec),
                      ino = as.numeric(i.vec), bv_int = as.numeric(bv.vec),
                      bn_int = as.numeric(bn.vec), bi_int = as.numeric(bi.vec),
                      vn_int = as.numeric(vn.vec), vi_int = as.numeric(vi.vec),
                      ni_int = as.numeric(ni.vec),
                      bvn_int = as.numeric(bvn.vec),
                      bvi_int = as.numeric(bvi.vec),
                      bni_int = as.numeric(bni.vec),
                      vni_int = as.numeric(vni.vec))
    return(sampling(od.cauchy, chains = 4, iter = 1000, data = stan.data))
  })
#save(sens, file = 'sens.rdata')
load('sens.rdata')

sensitivity.plot <- function(stan.list, pars, vals, xmin = 0, xmax = 1,
                             inner = 0.5, outer = 0.95,
                             ylab = 'Hyperparameter Value', ...){
# Creates a caterpillar plot comparing posterior intervals for each of
# several parameters for several values of one hyperparameter.
#  stan.list  A list of stan_fit objects
#  pars       Character vector of the parameters being examined
#  vals       Numeric vector of hyperparameter values that were used
#  xmin       Numeric (vector) of lower xlim values for each plot
#  xmax       Numeric (vector) of upper xlim values for each plot
#  inner      Confidence lever for outer interval
#  outer      Confidence level for outer interval
  if(length(xmin)<length(pars)) xmin <- rep(xmin, length(pars))
  if(length(xmax)<length(pars)) xmax <- rep(xmax, length(pars))
  for(j in 1:length(pars)){
    plot(NULL, xlim = c(xmin[j], xmax[j]), ylim = c(length(sens), 1),
         main = paste('Sensitivity Plot of', params[j]),
         xlab = params[j], ylab = ylab, yaxt = 'n', ...)
    axis(2, at = 1:length(sens), labels = scales, las = 1)
    for(i in 1:length(sens)){
      current <- unlist(extract(sens[[i]], pars = params[j]))
      segments(x0 = quantile(current, c((1-outer)/2, (1-inner)/2)),
               x1 = quantile(current, 1-c((1-outer)/2, (1-inner)/2)),
               y0 = i, lwd = c(1, 3))
      points(x = median(current), y = i, pch = 20)
    }
  }
}

params <- c("sigma_b", "sigma_v", "sigma_n", "sigma_i", "sigma_bv",
            "sigma_bn","sigma_bi", "sigma_vn", "sigma_vi", "sigma_ni",
            "sigma_bvn", "sigma_bvi", "sigma_bni", "sigma_vni")
xlims <- c(1, 1, 8, 20, 1, 1, 1, 1, 1, 3, 1, 1, 1, 5)
x11()
par(mfrow = c(3, 5))
sensitivity.plot(sens, params, scales, xmax = xlims, ylab = 'Prior Scale')

#for(j in 1:length(params)){
#  plot(NULL, xlim = c(0, xlims[j]), ylim = c(length(sens), 1),
#       main = paste('Comparison of', params[j], 'Prior Scales'),
#       xlab = params[j], ylab = 'Prior Scale', yaxt = 'n')
#  axis(2, at = 1:length(sens), labels = scales, las = 1)
#  for(i in 1:length(sens)){
#    sig.current <- unlist(extract(sens[[i]], pars = params[j]))
#    segments(x0 = quantile(sig.current, c(0.025, 0.25)),
#             x1 = quantile(sig.current, c(0.975, 0.75)),
#             y0 = i, lwd = c(1, 3))
#    points(x = median(sig.current), y = i, pch = 20)
#  }
#}

# Get the one Stan fit where we used a scale of 3
od.stan3 <- sens[[3]]

# Do some posterior predictin'
p.stan3 <- extract(od.stan3, pars = 'p')$p
y.pred3 <- apply(p.stan3, c(1, 2), function(p){return(rbinom(1, m, p))})

n.zero.pred3 <- apply(y.pred3, 1, function(x){return(sum(x==0))})
plot(table(n.zero.pred3)/sum(table(n.zero.pred3)),
     lwd = 3, xlim = c(15, 55), xaxt = 'n',
     main = 'Predicted Zero Counts',
     xlab = 'Number of Zeros', ylab = '')
abline(v = sum(y==0), lty = 'dashed', lwd = 2, col = 'red')
axis(1, seq(15, 55, 5))

max.pred3 <- apply(y.pred3, 1, max)
plot(table(max.pred3), xlim = c(6, 22), lwd = 3, xaxt = 'n',
     main = 'Predicted Maxima',
     xlab = 'Maximum', ylab = '')
abline(v = max(y), lty = 'dashed', lwd = 2, col = 'red')
axis(1, 6:22)

plot(x=y, y=y, type="n")


# Gelman ANOVA plot for superpop sds
batches3 <- c('block', 'variety', 'nitrogen', 'inoculation',
              'block * variety', 'block * nitrogen', 'block * inoculation',
              'variety * nitrogen', 'variety * inoculation',
              'nitrogen * inoculation',
              'block * variety * nitrogen', 'block * variety * inoculation',
              'block * nitrogen * inoculation',
              'variety * nitrogen * inoculation')
superpop.sds3 <- c("sigma_b", "sigma_v", "sigma_n", "sigma_i", "sigma_bv",
                   "sigma_bn","sigma_bi", "sigma_vn", "sigma_vi", "sigma_ni",
                   "sigma_bvn", "sigma_bvi", "sigma_bni", "sigma_vni")
dfs3 <- c(n.b, n.v, n.n, n.i, n.bv, n.bn, n.bi, n.vn, n.vi,
          n.ni, n.bvn, n.bvi, n.bni, n.vni) - 1
anova.superpop <- function(stan.fit, sources, sds, dfs, xlim = c(0, 1),
                           inner = 0.5, outer = 0.95,
                           width = c(6, 3, 25)){
# Creates a caterpillar plot comparing posterior intervals for each of
# several parameters for several values of one hyperparameter.
#  stan.fit  A stan_fit objects
#  sources   Character vector of batches being analyzed
#  sds       Character vector of the scale parameters being analyzed
#  dfs       Numeric vector of degrees of freedom for each batch
#  inner     Confidence lever for outer interval
#  outer     Confidence level for outer interval
#  width     Numeric vector of widths to send to layout()
  layout(matrix(1:3, nrow = 1), width = width)
  # Sources
  plot(NULL, xlim = c(0, 1), ylim = c(length(sources), 1), main = 'Sources',
       xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', frame.plot = FALSE,
       mar = c(4.1, 1.1, 4.1, 0.1))
  text(x = 1, y = 1:length(sources), labels = sources, pos = 2, offset = 0)
  # dfs
  plot(NULL, xlim = c(0, 1), ylim = c(length(sources), 1), main = 'df',
       xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', frame.plot = FALSE,
       mar = c(4.1, 0, 4.1, 0.1))
  text(x = 1, y = 1:length(sources), labels = dfs, pos = 2, offset = 0)
  # SDs
  plot(NULL, xlim = xlim, ylim = c(length(sds), 1),
       main = 'Estimated sd of Effects', xlab = '', ylab = '', yaxt = 'n',
       mar = c(4.1, 0, 4.1, 1.1))
  axis(3)
  for(j in 1:length(sources)){
    current <- unlist(extract(stan.fit, pars = sds[j]))
    segments(x0 = quantile(current, c((1-outer)/2, (1-inner)/2)),
             x1 = quantile(current, 1-c((1-outer)/2, (1-inner)/2)),
             y0 = j, lwd = c(1, 3))
    points(x = median(current), y = i, pch = 20)
  }
}

### git commit -am "message"
### git pull
### git push
