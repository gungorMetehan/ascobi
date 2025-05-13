##################################
# LOAD PACKAGES ##################
##################################

library(fGarch)
library(catR)
library(psychometric)
library(tidyverse)
library(patchwork)
library(data.table)
library(psych)
library(mirt)
library(scirt)
library(lavaan)
library(cirt)

# ANTIALIASING FOR PLOTS (ONLY FOR WINDOWS)
trace(grDevices:::png, quote({
  if (missing(type) && missing(antialias)) {
    type <- "cairo-png"
    antialias <- "subpixel"
  }
}), print = FALSE)

##################################
# OPTIONS ########################
##################################

plotName <- "01"
fileName <- "x01.csv"
fitFileName <- "01-fit.csv"
irtFileName <- "01-irt.csv"
rep <- 100
sampleSize <- 100000
itemSize <- 20
bmin <- -3
bmax <- 3
amin <- 0.5
amax <- 1.5
gseed <- 123
thetamean <- 0
thetasd <- 1
D <- 1.702
dist <- 1 #1: Normal, 2: Positive Skewed, 3: Negative Skewed
model1 <- "F1 =~ i1 + i2 + i3 + i4 + i5 + i6 + i7 + i8 + i9 + i10 +
          i11 + i12 + i13 + i14 + i15 + i16 + i17 + i18 + i19 + i20"
model2 <- "F2 =~ i1 + i2 + i3 + i4 + i5 + i6 + i7 + i8 + i9 + i10 +
          i11 + i12 + i13 + i14 + i15 + i16 + i17 + i18 + i19 + i20 +
          i21 + i22 + i23 + i24 + i25 + i26 + i27 + i28 + i29 + i30 +
          i31 + i32 + i33 + i34 + i35 + i36 + i37 + i38 + i39 + i40"
model <- model1 #model1: itemSize = 20; model2: itemSize = 40

##################################
# PARAMETER GENERATION ###########
##################################
set.seed(gseed)
a <- runif(itemSize, amin, amax)
set.seed(gseed)
b <- runif(itemSize, bmin, bmax)
c <- rep(0.25, itemSize)
d <- rep(1, itemSize)
itemPar <- cbind(a, b, c, d)

##################################
# THETA GENERATION ###############
##################################
set.seed(gseed)
thseed <- round(runif(rep, min = 1, max = 10000))

theta <- list()
theta <- vector(mode = "list", length = rep)
system.time(
  for(i in 1:rep){
    set.seed(thseed[i])
    if(dist == 1){
      theta[[i]] <- rnorm(sampleSize, thetamean, thetasd)
    } else if(dist == 2){
      theta[[i]] <- rsnorm(sampleSize, thetamean, thetasd, xi = 1.5)
    } else if(dist == 3){
      theta[[i]] <- rsnorm(sampleSize, thetamean, thetasd, xi = 0.5)
    }
  })

##################################
# GROUP GENERATION ###############
##################################
# Generating 10 different groups
groups <- data.frame(group = as.numeric(rep(1:10, each = sampleSize / 10)))

##################################
# RESPONSE PATTERN GENERATION ####
##################################
response <- list()
response <- vector(mode = "list", length = rep)
system.time(
  for(i in 1:rep){
    response[[i]] <- genPattern(
      theta[[i]],
      itemPar,
      model = NULL,
      D = D,
      seed = thseed[i]
    )
  })

dimensionality <- list()
uniDim <- list()
system.time(
  for(i in 1:rep){
    colnames(response[[i]]) <- paste0("i",1:itemSize)
    dimensionality[[i]] <- cfa(model = model, data = response[[i]],
                               estimator = "WLSMV")
    tempfit <- fitMeasures(dimensionality[[i]])
    uniDim[[i]] <- data.frame(cfi = tempfit["cfi"],
                              rmsea = tempfit["rmsea"],
                              rmseal = tempfit["rmsea.ci.lower"],
                              rmseau = tempfit["rmsea.ci.upper"],
                              tli = tempfit["tli"]) 
  }
)

fitInd <- matrix(unlist(uniDim), nrow = rep, ncol = 5, byrow = T)
fitIndx <- as.data.frame(t(colMeans(fitInd)))
colnames(fitIndx) <- c("cfi", "rmsea", "rmsea-l", "rmsea-u", "tli")
write.csv2(fitIndx, fitFileName)

##################################
# RELIABILITY ####################
##################################

reliability <- list()
system.time(
  for(i in 1:rep){
    reliability[[i]] <- alpha(response[[i]])$total$std.alpha
  })

rel <- mean(unlist(reliability))

##################################
# IRT CALIBRATION ################
##################################
IRTa <- NULL # IRT a parameter
IRTa2 <- NULL # IRT a' parameter (2nd derivative)
IRTb <- NULL # IRT b parameter
mirtCluster(5) # Parallel computing
system.time(
  for(i in 1:rep){
    colnames(response[[i]]) <- paste0("i", 1:itemSize)
    irt.tmp <- mirt(response[[i]][,1:itemSize], 
                    1, 
                    itemtype = "2PL")
    theta.tmp <- fscores(irt.tmp, method = "MAP")
    theta.tmp <- theta.tmp[order(theta.tmp)]
    theta.tmp <- theta.tmp[sampleSize*0.9]
    irt.tmp.a <- coef(irt.tmp, IRTpars = TRUE, simplify = TRUE)[["items"]][,1]
    irt.tmp.b <- coef(irt.tmp, IRTpars = TRUE, simplify = TRUE)[["items"]][,2]
    irt.tmp.par <- data.frame(alpha = irt.tmp.a, beta = irt.tmp.b)
    irt.tmp.a2 <- dIRF(irt.tmp.par, theta.tmp) / 0.25
    IRTa <- rbind(IRTa, irt.tmp.a)
    IRTb <- rbind(IRTb, irt.tmp.b)
    IRTa2 <- rbind(IRTa2, irt.tmp.a2)
    print(i)
  })
mirtCluster(remove = TRUE)

##################################
# ITEM ANALYSIS ##################
##################################
iAnalysis <- list()
system.time(
  for(i in 1:rep){
    iAnalysis[[i]] <- item.exam(response[[i]], discrim = TRUE)
  })

##################################
# POST-PROCESSING ################
##################################
for(i in 1:rep){
  response[[i]] <- cbind(response[[i]],
                         rowSums(response[[i]])) 
}

for(i in 1:rep){
  response[[i]] <- response[[i]][order(response[[i]][,itemSize+1]),]
}

for(i in 1:rep){
  response[[i]] <- cbind(response[[i]], groups)
}

for(i in 1:rep){
  response[[i]] <- mutate(response[[i]], group = as.factor(group))
}

for(i in 1:rep){
  response[[i]] <- split(response[[i]], response[[i]]$group)
}

##################################
# CALCULATE THE PROBABILITIES ####
##################################
mx <- NULL
my <- NULL
mz <- NULL
for(h in 1:itemSize){
  for(i in 1:rep){
    for(j in 1:10){
      mx <- rbind(mx, mean(response[[i]][[j]][[h]]))
    }
    my <- rbind(my, mx)
    mx <- NULL
  }
  mz <- cbind(mz, my)
  my <- NULL
}

mz <- cbind(mz, rep(c(1:10), rep))
mz <- as.data.frame(mz)
colnames(mz) <- c(1:itemSize, "group")

##################################
# PLOTS ##########################
##################################
g <- list()
for(i in 1:itemSize){
  g[[i]] <-ggplot(data = data.frame(
    ort = get(paste0("mz"))[,i],
    group = mz[,itemSize+1])) +
    aes(group, ort) +
    geom_point() + 
    geom_smooth(method ="loess") +
    scale_x_continuous(breaks = c(1:10)) + 
    ylim(c(0,1)) +
    ggtitle(paste0("Item",i)) +
    xlab("Groups") +
    ylab("Probability")
}

(g[[1]] | g[[2]]) / (g[[3]] | g[[4]]) / (g[[5]] | g[[6]]) / (g[[7]] | g[[8]]) / (g[[9]] | g[[10]])
ggsave(paste0(plotName, "-01.png"), type = "cairo", scale = 1,
       width = 13, height = 20, units = "cm", dpi = 300)
(g[[11]] | g[[12]]) / (g[[13]] | g[[14]]) / (g[[15]] | g[[16]]) / (g[[17]] | g[[18]]) / (g[[19]] | g[[20]])
ggsave(paste0(plotName, "-02.png"), type = "cairo", scale = 1,
       width = 13, height = 20, units = "cm", dpi = 300)
(g[[21]] | g[[22]]) / (g[[23]] | g[[24]]) / (g[[25]] | g[[26]]) / (g[[27]] | g[[28]]) / (g[[29]] | g[[30]])
ggsave(paste0(plotName, "-03.png"), type = "cairo", scale = 1,
       width = 13, height = 20, units = "cm", dpi = 300)
(g[[31]] | g[[32]]) / (g[[33]] | g[[34]]) / (g[[35]] | g[[36]]) / (g[[37]] | g[[38]]) / (g[[39]] | g[[40]])
ggsave(paste0(plotName, "-04.png"), type = "cairo", scale = 1,
       width = 13, height = 20, units = "cm", dpi = 300)

# (g[[1]] + ggtitle("Type 1") | g[[2]] + ggtitle("Type 2")) / (g[[7]] + ggtitle("Type 3") + plot_spacer())

# ggsave("rjx/type_03.eps", device = "eps", scale = 1,
#        width = 19, height = 13,
#        units = "cm")
# 
# ggsave("rjx/type_03.png", type = "cairo", scale = 1,
#        width = 19, height = 13, units = "cm", dpi = 300)

# Extension
upper10 <- colMeans(mz[which(mz$group == 10),])
lower90 <- colMeans(mz[-which(mz$group == 10),])
newCoef <- upper10 - lower90

# Phi 10-90
phiList <- NULL
phiItem <- NULL
phiItem2List <- NULL
for(i in 1:rep){
  for(j in 1:itemSize){
    phi.a <- sum(mz[10*i,j]*(sampleSize/10))
    phi.b <- 10000 - phi.a
    phi.c <- sum(mz[((i*10)-9):((i*10)-1),j]*(sampleSize/10))
    phi.d <- 90000 - phi.c
    phiItem <- phi(matrix(c(phi.a, phi.b, phi.c, phi.d), 
                          nrow = 2, ncol = 2, byrow = T))
    phiItem2List <- cbind(phiItem2List, phiItem)
    phiItem <- NULL
  }
  colnames(phiItem2List) <- paste0("i",1:itemSize)
  phiList <- rbind(phiList, phiItem2List)
  phiItem2List <- NULL
}

phix <- colMeans(phiList)

# Combine Item Analyzes
tmpRIR <- NULL
tmpRIT <- NULL
tmpULI <- NULL
RIR <- NULL
RIT <- NULL
ULI <- NULL

for(j in 1:itemSize){
  for(i in 1:rep){
    tmpRIR <- rbind(tmpRIR,
                    iAnalysis[[i]]$Item.Tot.woi[j])
    tmpRIT <- rbind(tmpRIT,
                    iAnalysis[[i]]$Item.total[j])
    tmpULI <- rbind(tmpULI,
                    iAnalysis[[i]]$Discrimination[j])
  }
  RIR <- rbind(RIR, mean(tmpRIR))
  RIT <- rbind(RIT, mean(tmpRIT))
  ULI <- rbind(ULI, mean(tmpULI))
  tmpRIR <- NULL
  tmpRIT <- NULL
  tmpULI <- NULL
}

itemAN <- cbind(RIT, RIR, ULI, 
                colMeans(phiList), 
                colMeans(IRTa), 
                colMeans(IRTa2),
                newCoef[-(itemSize+1)])
colnames(itemAN) <- c("RIT", "RIR", "ULI", "PHI", "IRT", "IRT2", "EXT")
write.csv2(itemAN, fileName)
