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
mu <- -2.2 # Intercept: Because logit(0.1) is about -2.2

set.seed(93)
# six blocks
b.eff <- rnorm(6,0,0.2)
# five varieties
v.eff <- rnorm(5,0,0.05)
# three nitrogen levels
n.eff <- sort(rnorm(3,0,0.7))
# two inoculation levels
i.eff <- c(-0.5,0.5)

# generate probablities
logit.pi.vec <- rep(0,180) # logit(pi) = X * beta
# block is col 1, variety is col 4, inoc is col 5, nit is col 6
for(i in 1:180){
  logit.pi.vec[i] <- mu + b.eff[design[i,1]] + v.eff[design[i,4]] + # 4 effects
                          n.eff[design[i,6]] + i.eff[design[i,5]] + 
                           b.eff[design[i,1]]*v.eff[design[i,4]] + # 6 two-way interactions
                           b.eff[design[i,1]]*n.eff[design[i,6]] +
                           b.eff[design[i,1]]*i.eff[design[i,5]] +
                           v.eff[design[i,4]]*n.eff[design[i,6]] +
                           v.eff[design[i,4]]*i.eff[design[i,5]] +
                           n.eff[design[i,6]]*i.eff[design[i,5]] +
                      b.eff[design[i,1]]*v.eff[design[i,4]]*n.eff[design[i,6]] + # 4 three-way interactions
                      b.eff[design[i,1]]*v.eff[design[i,4]]*i.eff[design[i,5]] +
                      b.eff[design[i,1]]*n.eff[design[i,6]]*i.eff[design[i,5]] +
                      v.eff[design[i,4]]*n.eff[design[i,6]]*i.eff[design[i,5]] +
     b.eff[design[i,1]]*v.eff[design[i,4]]*n.eff[design[i,6]]*i.eff[design[i,5]] # 1 four-way interaction
}
# transform back to pi scale
pi.vec <- apply(cbind(logit.pi.vec), 1, expit) # pi = expit(X * beta)
# simulate some y's (leaves - out of 30 - that are contaminated)
set.seed(324)
y.vec <- apply(cbind(pi.vec), 1, function(p) rbinom(1,30,p)) # y_i ~ Bin(30,pi_i)

# attach to data set
design$"infected" <- y.vec
design$"total" <- rep(30,180)
summary(design)
names(design)

# glmer
require(lme4)
# start with glm to see if it runs
glm.mod <- glm(cbind(infected,total-infected) ~ (block + variety + inoc + nitrogen)^4,
                data=design, family=binomial(link = "logit"))
summary(glm.mod)
