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
b.eff <- rnorm(6,0,0.01)
# five varieties
v.eff <- rnorm(5,0,0.05)
# three nitrogen levels
n.eff <- sort(rnorm(3,0,0.2))
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
glm.mod <- glm(cbind(infected,total-infected) ~ (block + variety + inoc + nitrogen)^4,
                data=design, family=binomial(link = "logit"))
summary(glm.mod)


### git commit -a -m "message"
### git pull
### git push
