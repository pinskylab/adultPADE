#### Let's first plot where these fish are coming from ####
setwd("/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation")

library(maps)
library(mapdata)
library(geosphere)

nvs_locs <- read.csv("allpops_combo.csv", header=TRUE) #all 241-9=232 fish divided into 5 populations; northvssouth.csv = lat/long for 32 fish from 5 populations
nvs_locs

#### Instead of lat & long, add a column for greater-circle distance from a southern point ####
longlat <- nvs_locs[,c("long", "lat")] # subset data to only longitude and latitude
longlat <- rbind(longlat, c(-80.546297, 28.384993)) # add a southern-most point
geodist <- distm(longlat, fun=distCosine) #matrix is in meters
hist(geodist[lower.tri(geodist)], nclass = 20)

# Distance matrix in km
geodistkm <- geodist * 0.001

# Add greater circle distance to data, minus the south reference point
nvs_locs$dist <- geodistkm[233,1:232]

#### Plotting covariance between five environmental variables ####
envs <- nvs_locs[-c(41, 96, 138, 199), c("lat", "long", "depth", "b_temp", "b_salin")]
dim(envs)
cors <- cor(envs, method = 'pearson') # calculate correlations
round(cors, 2) # coefficient of correlation
round(cors^2, 2) # coefficient of determination (pearson's coefficient of correlation squared)
cor.test(envs[,1], envs[,2])

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y)
  txt <- format(c(r^2, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = 1.5)
}

png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation/envs_covariance.png", width=7.5, height=7.5, res=300, units="in")
pairs(envs, upper.panel = panel.cor)
dev.off()

# Map of BayEnv populations
png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation/allpops_combo.png", width=5.5, height=7.5, res=300, units="in")

par(
  mar=c(5, 4, 4, 2), # panel magin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=18, # point size, which is the font size
  bg=NA
)

map("worldHires", "us", xlim=c(-85,-66), ylim=c(23,47.5), col="gray90", fill=TRUE) #plots the region of the USA that I want

points(nvs_locs$lon[1:40], nvs_locs$lat[1:40], pch=15, col="violet", cex=0.8) #plots sample sites
points(nvs_locs$lon[42:95], nvs_locs$lat[42:95], pch=15, col="blue", cex=0.8)
points(nvs_locs$lon[97:137], nvs_locs$lat[97:137], pch=15, col="green", cex=0.8)
points(nvs_locs$lon[139:198], nvs_locs$lat[139:198], pch=15, col="gold", cex=0.8) 
points(nvs_locs$lon[200:236], nvs_locs$lat[200:236], pch=15, col="tomato", cex=0.8) 

map.axes() #plots axes
title(xlab = "Longitude (°W)", ylab = "Latitude (°N)")

# title("Summer Flounder Sampling Sites")

dev.off()

#### PCA of these 232 individuals ####
library(ade4)
library(adegenet)
library(devtools)
library(hierfstat)
library(pegas)

# Reading in SNP data file containing only the first SNP at each locus
nvs_adults <- read.structure("structure_input_232forbayenv.str",
                         n.ind = 232, n.loc = 1137, col.lab = 1, col.pop = 2, row.marknames = 1, 
                         onerowperind = FALSE)

is.genind(nvs_adults)
head(indNames(nvs_adults),10)
locNames(nvs_adults)
sum <- summary(nvs_adults)
plot(sum$Hobs ~ sum$Hexp) # more heterozygosity observed than expected, what would cause this? outbreeding? introgression of one pop into another? balancing selection?
abline(a=0, b=1, col = "red")
barplot(sum$Hobs-sum$Hexp, main="Heterozygosity: observed-expected",
        ylab="Hobs - Hexp", ylim=c(-0.1,0.5))

bartlett.test(list(sum$Hobs, sum$Hexp)) #variances are not homogenous
t.test(sum$Hobs, sum$Hexp, pair = T, var.equal = FALSE, alternative = "g") # the Hexp and Hobs are significantly different from each other
adults.hwt <- hw.test(nvs_adults, B=10000) # B=number of replicates for MCMC procedure, or B=0 for regular HW test
adults.hwt
pval <- adults.hwt[adults.hwt[,"Pr.exact"] < 0.0100,] # p<0.01, exact test; like in Wolf population structure & canidate genes under selection paper
length(pval[,"Pr.exact"]) #132 SNPs with less than 0.01 probability of being in HWE; this number differs depending on whether Pr(chi^2 >) [138] or Pr.exact is used 

hist(adults.hwt[,4], nclass =100)

sum(is.na(nvs_adults$tab)) #3034
N <- scaleGen(nvs_adults, NA.method = "mean")
dim(N)
class (N)

# make PCA
nvs_pca <- dudi.pca(N,cent=FALSE,scale=FALSE,scannf=FALSE,nf=5) # saved 5 axes
barplot(nvs_pca$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))

nvs_pca
plot(nvs_pca$li[,1], nvs_pca$li[,2], col="blue", xlab = "PC1", ylab = "PC2") # PCA plot of 232 adult summer flounder using 1137 loci
plot(nvs_pca$li, nvs_pca$li)

# Using color to differentiate between 3 geographic locations
# zCol <- function(nCols, Z){
#   cols <- colorRampPalette(c("#000099", "#009E73", "#FF3100"))(nCols)
#   colVec_ind <- cut(Z, breaks=c(32.56, 34.6, 39.6, 41.3)) # manually defining the break points
#   colVec <- cols[colVec_ind]
# } 

zCol <- function(nCols, Z){
  cols <- colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100"))(nCols)
  colVec_ind <- cut(Z, breaks=nCols)
  colVec <- cols[colVec_ind]
} 

Z <- nvs_locs$V2 # using ID number instead of latitude
plot(nvs_pca$li[,1], nvs_pca$li[,2], col=zCol(6, Z), xlab = "PC1", ylab = "PC2") # based roughly on latitute, but colors don't represent populations correctly; see below for correct color breakdown
plot(nvs_pca$li, nvs_pca$li, col=zCol(6,Z)) # plots a 5x5 matrix showing all combinations of axes

# Makes a nice PCA plot using colors to distinguish between populations
png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation/232nvs_pca.png", width=5, height=4.5, res=300, units="in")

par(
  mar=c(5, 4, 4, 2), # panel magin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=10, # point size, which is the font size
  bg=NA
)

plot(nvs_pca$li[1:40,1], nvs_pca$li[1:40,2], pch = 1, col = "violet", xlab = "PC1", ylab = "PC2", xlim = c(-10,30), ylim = c(-40,20))
points(nvs_pca$li[41:94,1], nvs_pca$li[41:94,2], pch = 1, col = "blue")
points(nvs_pca$li[95:135,1], nvs_pca$li[95:135,2], pch=1, col = "green")
points(nvs_pca$li[136:195,1], nvs_pca$li[136:195,2], pch=1, col = "gold")
points(nvs_pca$li[196:232,1], nvs_pca$li[196:232,2], pch=1, col = "tomato")

legend("bottomright",
	legend=c("Population 1", "Population 2", "Population 3", "Population 4", "Population 5"),
	pch=c(1, 1, 1, 1, 1),
	col=c("violet", "blue", "green", "gold", "tomato"))

dev.off()

#### Calculating F statistics ####
fstat(nvs_adults)
#             pop         Ind
# Total 0.002536828(FST)  -0.09950153(FIT)
# pop   0.000000000       -0.10229787(FIS)

# FIT = FIS + FST - (FIS)(FST)
# negative FIS = excess of heterozygotes

pairwise.fst(nvs_adults) # Calculates Nei's pairwise Fst between all pairs of populations, heterozygosities weighted by group size
#        1           2           3           4
# 2 0.006127746                                    
# 3 0.007192645 0.006340812                        
# 4 0.005569441 0.004987458 0.005995512            
# 5 0.007211772 0.005985609 0.007784451 0.005257621

nvs_adults_pegas <- as.loci(nvs_adults) # Weir and Cockerham (1984)
F <- Fst(nvs_adults_pegas) # for each locus in data
F <- as.data.frame(F)
hist(F$Fst, breaks=50, main = "Fst at each locus", xlab = "Fst")
mean(F$Fst) # 0.00278727
mean(F$Fit) # -0.04612725
mean(F$Fis) # -0.04905604

#### I want heterozygozities on a population basis, not together like above ####
#### This will make a heterozygosity plot based on population ####
#### Splitting populations into different data files, since I can't figure out how to subset a genind object ####
pop1 <- read.structure("pop1.str",n.ind = 42, n.loc = 1137, col.lab = 1, 
                        col.pop = 2, row.marknames = 1, onerowperind = FALSE)
pop1.sum <- summary(pop1)
plot(pop1.sum$Hobs ~ pop1.sum$Hexp)
abline(a=0, b=1, col="red")
hist(pop1.sum$Hobs)
pop1.hwt <- hw.test(pop1, B=10000) # B=number of replicates for MCMC procedure (Fisher exact test), or B=0 for regular HW test
# pop1.hwt
# pop1.pval <- pop1.hwt[pop1.hwt[,"Pr.exact"] < 0.0100,] 
# length(pop1.pval[,"Pr.exact"]) 
pop1.avg <- mean(pop1.sum$Hobs) #0.2774038
pop1.sem <- sd(pop1.sum$Hobs)/sqrt(length(pop1.sum$Hobs)) #0.005112956

pop2 <- read.structure("pop2.str",n.ind = 57, n.loc = 1137, col.lab = 1, 
                       col.pop = 2, row.marknames = 1, onerowperind = FALSE)
pop2.sum <- summary(pop2)
plot(pop2.sum$Hobs ~ pop2.sum$Hexp)
abline(a=0, b=1, col="red")
hist(pop2.sum$Hobs)
pop2.avg <- mean(pop2.sum$Hobs) #0.2709195
pop2.sem <- sd(pop2.sum$Hobs)/sqrt(length(pop2.sum$Hobs)) #0.005009899

pop3 <- read.structure("pop3.str",n.ind = 36, n.loc = 1137, col.lab = 1, 
                       col.pop = 2, row.marknames = 1, onerowperind = FALSE)
pop3.sum <- summary(pop3)
plot(pop3.sum$Hobs ~ pop3.sum$Hexp)
abline(a=0, b=1, col="red")
hist(pop3.sum$Hobs)
pop3.avg <- mean(pop3.sum$Hobs) #0.2825785
pop3.sem <- sd(pop3.sum$Hobs)/sqrt(length(pop3.sum$Hobs)) #0.005323603

pop4 <- read.structure("pop4.str",n.ind = 60, n.loc = 1137, col.lab = 1, 
                       col.pop = 2, row.marknames = 1, onerowperind = FALSE)
pop4.sum <- summary(pop4)
plot(pop4.sum$Hobs ~ pop4.sum$Hexp)
abline(a=0, b=1, col="red")
hist(pop4.sum$Hobs)
pop4.avg <- mean(pop4.sum$Hobs) #0.2546314
pop4.sem <- sd(pop4.sum$Hobs)/sqrt(length(pop4.sum$Hobs)) #0.004743794

pop5 <- read.structure("pop5.str",n.ind = 37, n.loc = 1137, col.lab = 1, 
                       col.pop = 2, row.marknames = 1, onerowperind = FALSE)
pop5.sum <- summary(pop5)
plot(pop5.sum$Hobs ~ pop5.sum$Hexp)
abline(a=0, b=1, col="red")
hist(pop5.sum$Hobs)
pop5.avg <- mean(pop5.sum$Hobs) #0.2708982
pop5.sem <- sd(pop5.sum$Hobs)/sqrt(length(pop5.sum$Hobs)) #0.005380584

data.sum <- data.frame(matrix(nrow=5,ncol=3))
names(data.sum) <- c("pop", "mean", "sem")
data.sum$pop <- c(1, 2, 3, 4, 5)
data.sum$mean <- c(pop1.avg, pop2.avg, pop3.avg, pop4.avg, pop5.avg)
data.sum$sem <- c(pop1.sem, pop2.sem, pop3.sem, pop4.sem, pop5.sem)
data.sum


png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation/232_avgh.png", width=3.5, height=4.5, res=300, units="in")

par(
  mar=c(5, 4, 4, 2), # panel magin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=10, # point size, which is the font size
  bg=NA
)

bp <- barplot(data.sum$mean, ylab = "average heterozygosity", names.arg = c("Population 1", "Population 2", "Population 3", "Population 4", "Population 5"), las = 3, col = c("violet", "blue", "green", "gold", "tomato"), border = NA, ylim = c(0,0.30))
arrows(bp, data.sum$mean-data.sum$sem, bp, data.sum$mean+data.sum$sem, length=0.05, angle=90, code=3, lwd =1.5)

dev.off()

#### Bayenv2 implementation to detect Fst outliers ####
# Creating the ENVIROFILE using the average latitude of each population (and longitude for fun)
pop1.avglat <- mean(nvs_locs$lat[1:40])
pop1.avglong <- mean(nvs_locs$long[1:40])
pop1.avgdepth <- mean(nvs_locs$depth[1:40])
pop1.avgb_temp <- mean(nvs_locs$b_temp[1:40])
pop1.avgb_salin <- mean(nvs_locs$b_salin[1:40])

pop2.avglat <- mean(nvs_locs$lat[42:95])
pop2.avglong <- mean(nvs_locs$long[42:95])
pop2.avgdepth <- mean(nvs_locs$depth[42:95])
pop2.avgb_temp <- mean(nvs_locs$b_temp[42:95])
pop2.avgb_salin <- mean(nvs_locs$b_salin[42:95])

pop3.avglat <- mean(nvs_locs$lat[97:137])
pop3.avglong <- mean(nvs_locs$long[97:137])
pop3.avgdepth <- mean(nvs_locs$depth[97:137])
pop3.avgb_temp <- mean(nvs_locs$b_temp[97:137])
pop3.avgb_salin <- mean(nvs_locs$b_salin[97:137])

pop4.avglat <- mean(nvs_locs$lat[139:198])
pop4.avglong <- mean(nvs_locs$long[139:198])
pop4.avgdepth <- mean(nvs_locs$depth[139:198])
pop4.avgb_temp <- mean(nvs_locs$b_temp[139:198])
pop4.avgb_salin <- mean(nvs_locs$b_salin[139:198])

pop5.avglat <- mean(nvs_locs$lat[200:236])
pop5.avglong <- mean(nvs_locs$long[200:236])
pop5.avgdepth <- mean(nvs_locs$depth[200:236])
pop5.avgb_temp <- mean(nvs_locs$b_temp[200:236])
pop5.avgb_salin <- mean(nvs_locs$b_salin[200:236])

ll.sum <- data.frame(matrix(nrow=5,ncol=6))
names(ll.sum) <- c("pop", "avglat", "avglong", "avgdepth", "avgb_temp", "avgb_salin")
ll.sum$pop <- c(1, 2, 3, 4, 5)
ll.sum$avglat <- c(pop1.avglat, pop2.avglat, pop3.avglat, pop4.avglat, pop5.avglat)
ll.sum$avglong <- c(pop1.avglong, pop2.avglong, pop3.avglong, pop4.avglong, pop5.avglong)
ll.sum$avgdepth <- c(pop1.avgdepth, pop2.avgdepth, pop3.avgdepth, pop4.avgdepth, pop5.avgdepth)
ll.sum$avgb_temp <- c(pop1.avgb_temp, pop2.avgb_temp, pop3.avgb_temp, pop4.avgb_temp, pop5.avgb_temp)
ll.sum$avgb_salin <- c(pop1.avgb_salin, pop2.avgb_salin, pop3.avgb_salin, pop4.avgb_salin, pop5.avgb_salin)
ll.sum
plot(ll.sum$avglong ~ ll.sum$avglat, xlab = "Average latitude", ylab = "Average longitude")

# Standardizing the aggregated individuals
stan.ll.sum <- scale(ll.sum[,2:6])
stan.ll.sum
t.stan.ll.sum <- t(stan.ll.sum)
t.stan.ll.sum

# Exporting the envirofile to the Local adaptation folder on my computer
write.table (t.stan.ll.sum, "232stanenvirofile.txt", sep ="\t", col.names = FALSE, row.names = FALSE) # manually redid the tabbing afterwards

# Aggregate the individuals into populations and standardize
by <- list(nvs_locs$bayenv_pop)
pop.avg <- aggregate(nvs_locs[,c(7, 8, 9, 11)], by = by, FUN = mean)
stan.pop.avg <- scale(pop.avg[,2:5])
rownames(stan.pop.avg) <- c("Pop1", "Pop2", "Pop3", "Pop4", "Pop5")
t.stan.pop.avg <- t(stan.pop.avg) # the order is Pop1, Pop2, Pop3, Pop4, Pop5; rows are: depth, b_temp, b_salin, dist

# Exporting the envirofile to my computer
write.table (t.stan.pop.avg, "232stan4envirofile.txt", sep ="\t", col.names = FALSE, row.names = FALSE) # manually redid the tabbing afterwards; order needs to match that of SNP frequencies (Japan, Indonesia, Philippines)


### Below is code for standardizing slightly incorrectly....not totally sure if it acutually makes a differences because the relative differences are the same
g.lat.avg <- mean(nvs_locs$lat, na.rm=TRUE)
g.long.avg <- mean(nvs_locs$long, na.rm=TRUE)
g.depth.avg <- mean(nvs_locs$depth, na.rm=TRUE)
g.b_temp.avg <- mean(nvs_locs$b_temp, na.rm=TRUE)
g.b_salin.avg <- mean(nvs_locs$b_salin, na.rm=TRUE)

g.lat.sd <- sd(nvs_locs$lat, na.rm = TRUE)
g.long.sd <- sd(nvs_locs$long, na.rm = TRUE)
g.depth.sd <- sd(nvs_locs$depth, na.rm = TRUE)
g.b_temp.sd <- sd(nvs_locs$b_temp, na.rm = TRUE)
g.b_salin.sd <- sd(nvs_locs$b_salin, na.rm = TRUE)

# Standardizing the data
std.avglat <- c(((pop1.avglat-g.lat.avg)/g.lat.sd), ((pop2.avglat-g.lat.avg)/g.lat.sd), ((pop3.avglat-g.lat.avg)/g.lat.sd), ((pop4.avglat-g.lat.avg)/g.lat.sd), ((pop5.avglat-g.lat.avg)/g.lat.sd))
std.avglong <- c(((pop1.avglong-g.long.avg)/g.long.sd), ((pop2.avglong-g.long.avg)/g.long.sd), ((pop3.avglong-g.long.avg)/g.long.sd), ((pop4.avglong-g.long.avg)/g.long.sd), ((pop5.avglong-g.long.avg)/g.long.sd))
std.avgdepth <- c(((pop1.avgdepth-g.depth.avg)/g.depth.sd), ((pop2.avgdepth-g.depth.avg)/g.depth.sd), ((pop3.avgdepth-g.depth.avg)/g.depth.sd), ((pop4.avgdepth-g.depth.avg)/g.depth.sd), ((pop5.avgdepth-g.depth.avg)/g.depth.sd))
std.avgb_temp <- c(((pop1.avgb_temp-g.b_temp.avg)/g.b_temp.sd), ((pop2.avgb_temp-g.b_temp.avg)/g.b_temp.sd), ((pop3.avgb_temp-g.b_temp.avg)/g.b_temp.sd), ((pop4.avgb_temp-g.b_temp.avg)/g.b_temp.sd), ((pop5.avgb_temp-g.b_temp.avg)/g.b_temp.sd))
std.avgb_salin <- c(((pop1.avgb_salin-g.b_salin.avg)/g.b_salin.sd), ((pop2.avgb_salin-g.b_salin.avg)/g.b_salin.sd), ((pop3.avgb_salin-g.b_salin.avg)/g.b_salin.sd), ((pop4.avgb_salin-g.b_salin.avg)/g.b_salin.sd), ((pop5.avgb_salin-g.b_salin.avg)/g.b_salin.sd))

envirofile <- rbind(std.avglat, std.avglong, std.avgdepth, std.avgb_temp, std.avgb_salin)
envirofile

# Exporting the envirofile to the Local adaptation folder on my computer
write.table (envirofile, "232envirofile.txt", sep ="\t", col.names = FALSE, row.names = FALSE) # manually redid the tabbing afterwards

sample(-10000:10000, 1) # random number generator to set the random seed

# Visualizing the covariance matrix
library(corrplot)

covar <- matrix (c(3.847203e-03, -5.259565e-06,  1.509697e-05, -5.259565e-06, 3.774831e-03, -2.342748e-05, 1.509697e-05, -2.342748e-05, 3.786794e-03), nrow = 3, ncol = 3) # this is covarmatrix_hwe
covar

cor <- cov2cor(covar)
cor

# Plotting the covariance and correlation matricies
# Correlation matrix should reflect Fst between populations
par(mfrow=c(1,2))
corrplot(covar, method = "color", main="Covariance Matrix", is.corr = FALSE, cl.lim = c(-.00005259570, 3.847203e-03))
corrplot(cor, method = "color", main="Correlation Matrix", is.corr = TRUE, cl.lim = c(-0.01, 1.000000000))


#### Averaging the MCMC output matricies so that I can use this averaged matrix as input for the next step of BayEnv2
# Reading in the matrix.out file from BayEnv2 using the fuction that Ryan wrote [most up to date code is in Ryan_newcode.R]
####
read_jho <- function(file, n_rows, n_matrices, n_comment, n_top_comment=0, ...){
  stopifnot(file.exists(file)) # check that file exists
  
  data_read <- list()
  for(i in 1:n_matrices){
    line_start <- (i-1)*(n_comment+n_rows) + n_comment + n_top_comment
    data_read[[i]] <- scan(file, skip=line_start, nlines=n_rows, ...)
    data_read[[i]] <- as.numeric(data_read[[i]])
    data_read[[i]] <- matrix(data_read[[i]], nrow=n_rows, byrow=TRUE)
  }
  
  data_array <- simplify2array(data_read)
  
  return(data_array)
}

data_array <- read_jho(file="matrix.out", n_rows=5, n_matrices=200, n_comment=2, n_top_comment=13, what="character")

# Take average of the last 40 matricies from the MCMC chain
data_array_mean <- apply(data_array[,,160:200], c(1,2), mean)
data_array_mean # this is the covarmatrix on Amphiprion that I am using to estimate BF for each SNP

########################################################################################
#### This section reads in results from multiple BayEnv runs to look at consistency ####
#### Check the working directory to make sure you're looking at the right results ######
########################################################################################

#### How similar are BF for two independent runs of BayEnv2? ####
setwd("/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation/Archived files from old analyses/results_avgcovarmatrix")

run1 <- read.table("bf_1.txt")
run1 <- run1$V2
run2 <- read.table("bf_2.txt")
run2 <- run2$V2
run3 <- read.table("bf_3.txt")
run3 <- run3$V2
run4 <- read.table("bf_4.txt")
run4 <- run4$V2
run5 <- read.table("bf_5.txt")
run5 <- run5$V2
run6 <- read.table("bf_6.txt")
run6 <- run6$V2
run7 <- read.table("bf_7.txt")
run7 <- run7$V2
run8 <- read.table("bf_8.txt")
run8 <- run8$V2
run9 <- read.table("bf_9.txt")
run9 <- run9$V2
run10 <- read.table("bf_10.txt")
run10 <- run10$V2

bf_mean <- rowMeans(cbind(run1, run2, run3, run4, run5, run6, run7, run8, run9, run10))
bf_mean
write.table(cbind(run1, run2, run3, run4, run5, run6, run7, run8, run9, run10), "~/Desktop/100k.txt", sep="\t")
hist(bf_mean)
which(bf_mean > 3) #16 SNPs: 59   66  102  125  293  316  334  430  511  597  751  868  920 1082 1088 1099

# Does running the MCMC for 500000 iterations make a difference?
setwd("/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation/Archived files from old analyses/results_avgcovarmatrix_500k")

run1_500k <- read.table("bf_1.txt")
run1_500k <- run1_500k$V2
run2_500k <- read.table("bf_2.txt")
run2_500k <- run2_500k$V2
run3_500k <- read.table("bf_3.txt")
run3_500k <- run3_500k$V2
run4_500k <- read.table("bf_4.txt")
run4_500k <- run4_500k$V2
run5_500k <- read.table("bf_5.txt") # SNP 917 is missing...
run5_500k <- run5_500k$V2

bf_500k_mean <- rowMeans(cbind(run1_500k, run2_500k, run3_500k, run4_500k))
bf_500k_mean
write.table(cbind(run1_500k, run2_500k, run3_500k, run4_500k), "~/Desktop/500k.txt", sep="\t")
hist(bf_500k_mean)
which(bf_500k_mean > 3) #13 106  218  282  293  334  430  638  720  831  877  919  944 1133 when only including bf1 and bf2 runs
                        #59  106  218  282  293  334  430  720  885  919  944 1133 when including all 3 runs
                        # log10 > which(log>1), [1]  218  282  430  720 1133

# Does removing loci not in global HWE make a difference/more consistant across independent BayEnv runs?
setwd("/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation/Archived files from old analyses/results_avgcovarmatrix_hwe")

run1_hwe <- read.table("bf1.txt")
run1_hwe <- run1_hwe$V2
run2_hwe <- read.table("bf2.txt")
run2_hwe <- run2_hwe$V2
run3_hwe <- read.table("bf3.txt")
run3_hwe <- run3_hwe$V2
run4_hwe <- read.table("bf4.txt")
run4_hwe <- run4_hwe$V2
run5_hwe <- read.table("bf5.txt")
run5_hwe <- run5_hwe$V2

bf_hwe_mean <- rowMeans(cbind(run1_hwe, run2_hwe, run3_hwe, run4_hwe, run5_hwe))
bf_hwe_mean
hist(bf_hwe_mean)
which(bf_hwe_mean > 3) # 34  125  137  141  244  293  312  334  394  430  631  699  749  907  918  919  995 1082
which(log10(bf_hwe_mean) > 2) # 141  244  293  312  430 1082

# 100k MCMC iterations is still variable across runs, so trying 500k iterations in the hopes for better consistency
setwd("/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation/Archived files from old analyses/results_avgcovarmatrix_hwe_500k")

run1_hwe_500k <- read.table("bf1.txt")
run1_hwe_500k <- run1_hwe_500k$V2
run2_hwe_500k <- read.table("bf2.txt")
run2_hwe_500k <- run2_hwe_500k$V2
run3_hwe_500k <- read.table("bf3.txt")
run3_hwe_500k <- run3_hwe_500k$V2
run4_hwe_500k <- read.table("bf4.txt")
run4_hwe_500k <- run4_hwe_500k$V2
run5_hwe_500k <- read.table("bf5.txt")
run5_hwe_500k <- run5_hwe_500k$V2
run6_hwe_500k <- read.table("bf6.txt")
run6_hwe_500k <- run6_hwe_500k$V2
run7_hwe_500k <- read.table("bf7.txt")
run7_hwe_500k <- run7_hwe_500k$V2
run8_hwe_500k <- read.table("bf8.txt")
run8_hwe_500k <- run8_hwe_500k$V2
run9_hwe_500k <- read.table("bf9.txt")
run9_hwe_500k <- run9_hwe_500k$V2
run10_hwe_500k <- read.table("bf10.txt")
run10_hwe_500k <- run10_hwe_500k$V2


bf_hwe_500k_mean <- rowMeans(cbind(run1_hwe_500k, run2_hwe_500k, run3_hwe_500k, run4_hwe_500k, run5_hwe_500k, run6_hwe_500k, run7_hwe_500k, run8_hwe_500k, run9_hwe_500k, run10_hwe_500k))
bf_hwe_500k_mean
pairs(cbind(run1_hwe_500k, run2_hwe_500k, run3_hwe_500k, run4_hwe_500k, run5_hwe_500k, run6_hwe_500k, run7_hwe_500k, run8_hwe_500k, run9_hwe_500k, run10_hwe_500k)) # plots a scatterplot matrix
hist(bf_hwe_500k_mean)
which(bf_hwe_500k_mean > 3) # Average of bf1-bf5: 125  293  304  334  340  430  486  505  597  755  841  907  919  944  947  950  995 1082
                            # Average of bf1-bf10: 125  155  169  214  228  293  334  340  430  486  718  841  919  944  947  950  995 1082 1133
which(log10(bf_hwe_500k_mean) > 0.5)
# 500k iterations improves consistency between runs, for the most part --> would more MCMC runs lead to even better convergence?

# Plotting boxplots of loci with bf > 3 for either the 100k, 500k or both runs
pdf("canidateloci.pdf")
boxplot(cbind(run1_hwe, run2_hwe, run3_hwe, run4_hwe, run5_hwe)[34,], cbind(run1_hwe_500k, run2_hwe_500k, run3_hwe_500k, run4_hwe_500k, run5_hwe_500k)[34,], main="34")
boxplot(cbind(run1_hwe, run2_hwe, run3_hwe, run4_hwe, run5_hwe)[125,], cbind(run1_hwe_500k, run2_hwe_500k, run3_hwe_500k, run4_hwe_500k, run5_hwe_500k)[125,], main="125")
boxplot(cbind(run1_hwe, run2_hwe, run3_hwe, run4_hwe, run5_hwe)[137,], cbind(run1_hwe_500k, run2_hwe_500k, run3_hwe_500k, run4_hwe_500k, run5_hwe_500k)[137,], main="137")
boxplot(cbind(run1_hwe, run2_hwe, run3_hwe, run4_hwe, run5_hwe)[141,], cbind(run1_hwe_500k, run2_hwe_500k, run3_hwe_500k, run4_hwe_500k, run5_hwe_500k)[141,], main="141")
boxplot(cbind(run1_hwe, run2_hwe, run3_hwe, run4_hwe, run5_hwe)[244,], cbind(run1_hwe_500k, run2_hwe_500k, run3_hwe_500k, run4_hwe_500k, run5_hwe_500k)[244,], main="244")
boxplot(cbind(run1_hwe, run2_hwe, run3_hwe, run4_hwe, run5_hwe)[293,], cbind(run1_hwe_500k, run2_hwe_500k, run3_hwe_500k, run4_hwe_500k, run5_hwe_500k)[293,], main="293")
boxplot(cbind(run1_hwe, run2_hwe, run3_hwe, run4_hwe, run5_hwe)[304,], cbind(run1_hwe_500k, run2_hwe_500k, run3_hwe_500k, run4_hwe_500k, run5_hwe_500k)[304,], main="304")
boxplot(cbind(run1_hwe, run2_hwe, run3_hwe, run4_hwe, run5_hwe)[312,], cbind(run1_hwe_500k, run2_hwe_500k, run3_hwe_500k, run4_hwe_500k, run5_hwe_500k)[312,], main="312")
boxplot(cbind(run1_hwe, run2_hwe, run3_hwe, run4_hwe, run5_hwe)[334,], cbind(run1_hwe_500k, run2_hwe_500k, run3_hwe_500k, run4_hwe_500k, run5_hwe_500k)[334,], main="334")
boxplot(cbind(run1_hwe, run2_hwe, run3_hwe, run4_hwe, run5_hwe)[340,], cbind(run1_hwe_500k, run2_hwe_500k, run3_hwe_500k, run4_hwe_500k, run5_hwe_500k)[340,], main="340")
boxplot(cbind(run1_hwe, run2_hwe, run3_hwe, run4_hwe, run5_hwe)[394,], cbind(run1_hwe_500k, run2_hwe_500k, run3_hwe_500k, run4_hwe_500k, run5_hwe_500k)[394,], main="394")
boxplot(cbind(run1_hwe, run2_hwe, run3_hwe, run4_hwe, run5_hwe)[430,], cbind(run1_hwe_500k, run2_hwe_500k, run3_hwe_500k, run4_hwe_500k, run5_hwe_500k)[430,], main="430")
boxplot(cbind(run1_hwe, run2_hwe, run3_hwe, run4_hwe, run5_hwe)[486,], cbind(run1_hwe_500k, run2_hwe_500k, run3_hwe_500k, run4_hwe_500k, run5_hwe_500k)[486,], main="486")
boxplot(cbind(run1_hwe, run2_hwe, run3_hwe, run4_hwe, run5_hwe)[505,], cbind(run1_hwe_500k, run2_hwe_500k, run3_hwe_500k, run4_hwe_500k, run5_hwe_500k)[505,], main="505")
boxplot(cbind(run1_hwe, run2_hwe, run3_hwe, run4_hwe, run5_hwe)[597,], cbind(run1_hwe_500k, run2_hwe_500k, run3_hwe_500k, run4_hwe_500k, run5_hwe_500k)[597,], main="597")
boxplot(cbind(run1_hwe, run2_hwe, run3_hwe, run4_hwe, run5_hwe)[631,], cbind(run1_hwe_500k, run2_hwe_500k, run3_hwe_500k, run4_hwe_500k, run5_hwe_500k)[631,], main="631")
boxplot(cbind(run1_hwe, run2_hwe, run3_hwe, run4_hwe, run5_hwe)[699,], cbind(run1_hwe_500k, run2_hwe_500k, run3_hwe_500k, run4_hwe_500k, run5_hwe_500k)[699,], main="699")
boxplot(cbind(run1_hwe, run2_hwe, run3_hwe, run4_hwe, run5_hwe)[749,], cbind(run1_hwe_500k, run2_hwe_500k, run3_hwe_500k, run4_hwe_500k, run5_hwe_500k)[749,], main="749")
boxplot(cbind(run1_hwe, run2_hwe, run3_hwe, run4_hwe, run5_hwe)[755,], cbind(run1_hwe_500k, run2_hwe_500k, run3_hwe_500k, run4_hwe_500k, run5_hwe_500k)[755,], main="755")
boxplot(cbind(run1_hwe, run2_hwe, run3_hwe, run4_hwe, run5_hwe)[841,], cbind(run1_hwe_500k, run2_hwe_500k, run3_hwe_500k, run4_hwe_500k, run5_hwe_500k)[841,], main="841")
boxplot(cbind(run1_hwe, run2_hwe, run3_hwe, run4_hwe, run5_hwe)[907,], cbind(run1_hwe_500k, run2_hwe_500k, run3_hwe_500k, run4_hwe_500k, run5_hwe_500k)[907,], main="907")
boxplot(cbind(run1_hwe, run2_hwe, run3_hwe, run4_hwe, run5_hwe)[918,], cbind(run1_hwe_500k, run2_hwe_500k, run3_hwe_500k, run4_hwe_500k, run5_hwe_500k)[918,], main="918")
boxplot(cbind(run1_hwe, run2_hwe, run3_hwe, run4_hwe, run5_hwe)[919,], cbind(run1_hwe_500k, run2_hwe_500k, run3_hwe_500k, run4_hwe_500k, run5_hwe_500k)[919,], main="919")
boxplot(cbind(run1_hwe, run2_hwe, run3_hwe, run4_hwe, run5_hwe)[944,], cbind(run1_hwe_500k, run2_hwe_500k, run3_hwe_500k, run4_hwe_500k, run5_hwe_500k)[944,], main="944")
boxplot(cbind(run1_hwe, run2_hwe, run3_hwe, run4_hwe, run5_hwe)[947,], cbind(run1_hwe_500k, run2_hwe_500k, run3_hwe_500k, run4_hwe_500k, run5_hwe_500k)[947,], main="947")
boxplot(cbind(run1_hwe, run2_hwe, run3_hwe, run4_hwe, run5_hwe)[950,], cbind(run1_hwe_500k, run2_hwe_500k, run3_hwe_500k, run4_hwe_500k, run5_hwe_500k)[950,], main="950")
boxplot(cbind(run1_hwe, run2_hwe, run3_hwe, run4_hwe, run5_hwe)[995,], cbind(run1_hwe_500k, run2_hwe_500k, run3_hwe_500k, run4_hwe_500k, run5_hwe_500k)[995,], main="995")
boxplot(cbind(run1_hwe, run2_hwe, run3_hwe, run4_hwe, run5_hwe)[1082,], cbind(run1_hwe_500k, run2_hwe_500k, run3_hwe_500k, run4_hwe_500k, run5_hwe_500k)[1082,], main="1082")
dev.off()

# Creates a data frame of canidate loci and their averaged BF:
can_loci <- data.frame(matrix(nrow=19,ncol=2))
names(can_loci) <- c("locus", "bf")
can_loci$locus <- c(which(bf_hwe_500k_mean > 3))
can_loci$bf <- c(bf_hwe_500k_mean[125], bf_hwe_500k_mean[155], bf_hwe_500k_mean[169], bf_hwe_500k_mean[214], bf_hwe_500k_mean[228], bf_hwe_500k_mean[293], bf_hwe_500k_mean[334], bf_hwe_500k_mean[340], bf_hwe_500k_mean[430], bf_hwe_500k_mean[486], bf_hwe_500k_mean[718], bf_hwe_500k_mean[841], bf_hwe_500k_mean[919], bf_hwe_500k_mean[944], bf_hwe_500k_mean[947], bf_hwe_500k_mean[950], bf_hwe_500k_mean[995], bf_hwe_500k_mean[1082], bf_hwe_500k_mean[1133])
can_loci

#############################################################################################################
#### Examining consistency for BayEnv runs using a covariance matrix created with 232 fish and 1005 loci ####
#### Testing 1137 loci for associations with latitude                                                    ####
#############################################################################################################
setwd("/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation/results_avg232covarmatrix_hwe")

run1_232hwe_500k <- read.table("bf1.txt")
run1_232hwe_500k <- run1_232hwe_500k$V2
run2_232hwe_500k <- read.table("bf2.txt")
run2_232hwe_500k <- run2_232hwe_500k$V2
run3_232hwe_500k <- read.table("bf3.txt")
run3_232hwe_500k <- run3_232hwe_500k$V2
run4_232hwe_500k <- read.table("bf4.txt")
run4_232hwe_500k <- run4_232hwe_500k$V2
run5_232hwe_500k <- read.table("bf5.txt")
run5_232hwe_500k <- run5_232hwe_500k$V2
run6_232hwe_500k <- read.table("bf6.txt")
run6_232hwe_500k <- run6_232hwe_500k$V2
run7_232hwe_500k <- read.table("bf7.txt")
run7_232hwe_500k <- run7_232hwe_500k$V2
run8_232hwe_500k <- read.table("bf8.txt")
run8_232hwe_500k <- run8_232hwe_500k$V2
run9_232hwe_500k <- read.table("bf9.txt")
run9_232hwe_500k <- run9_232hwe_500k$V2
run10_232hwe_500k <- read.table("bf10.txt")
run10_232hwe_500k <- run10_232hwe_500k$V2


bf_232hwe_500k_mean <- rowMeans(cbind(run1_232hwe_500k, run2_232hwe_500k, run3_232hwe_500k, run4_232hwe_500k, run5_232hwe_500k, run6_232hwe_500k, run7_232hwe_500k, run8_232hwe_500k, run9_232hwe_500k, run10_232hwe_500k))
bf_232hwe_500k_mean
hist(bf_232hwe_500k_mean)
which(bf_232hwe_500k_mean > 3) # Average from runs 1-5: 125 214 442 499 609 615 825 919
                               # Average from runs 1-10: 125 214 396 442 499 609 615 703 825 919
which(log10(bf_232hwe_500k_mean) > 1) # Average from runs 1-5: 125 615
                                      # Average from runs 1-10:  125

########################################################################################################
# Okay, including more environmental variables and rerunning BayEnv. Reading in 10 independent results #
########################################################################################################
setwd("/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation/results_avg232covarmatrix_hwe_5enviro")

library(abind)

fullrun1 <- read.table("bf1.txt")
fullrun1$V1 <- NULL
fullrun2 <- read.table("bf2.txt")
fullrun2$V1 <- NULL
fullrun3 <- read.table("bf3.txt")
fullrun3$V1 <- NULL
fullrun4 <- read.table("bf4.txt")
fullrun4$V1 <- NULL
fullrun5 <- read.table("bf5.txt")
fullrun5$V1 <- NULL
fullrun6 <- read.table("bf6.txt")
fullrun6$V1 <- NULL
fullrun7 <- read.table("bf7.txt")
fullrun7$V1 <- NULL
fullrun8 <- read.table("bf8.txt")
fullrun8$V1 <- NULL
fullrun9 <- read.table("bf9.txt")
fullrun9$V1 <- NULL
fullrun10 <- read.table("bf10.txt")
fullrun10$V1 <- NULL

full_array <- abind(fullrun1, fullrun2, fullrun3, fullrun4, fullrun5, fullrun6, fullrun7, fullrun8, fullrun9, fullrun10, along=3)

full_array_mean <- apply(full_array, c(1,2), mean)
full_array_mean

which(full_array_mean[,1] > 3) #125 214 442 499 579 609 615 703 817 825 919
which(full_array_mean[,2] > 3) #125 214 228 396 442 579 609 647 825
which(full_array_mean[,3] > 3) #125  214  396  609  615  647  708  825  826 1061
which(full_array_mean[,4] > 3) #125  214  396  442  499  524  609  615  626  647  698  703  708  808  825  829  916  919 1050
which(full_array_mean[,5] > 3) #35  125  334  542  609  615  703  743  826  923  990 1031 1090

which(log10(full_array_mean[,1]) > 1) #none
which(log10(full_array_mean[,2]) > 1) #125
which(log10(full_array_mean[,3]) > 1) #none
which(log10(full_array_mean[,4]) > 1) #214
which(log10(full_array_mean[,5]) > 1) #743

quantile(log10(full_array_mean[,1]), .995)

# Maybe taking the median will help get rid of extreme runs that pull the average up
full_array_median <- apply(full_array, c(1,2), median)
full_array_median

which(full_array_median[,1] > 3) # 125 214 609 615 919
which(full_array_median[,2] > 3) # 125 214 396 442 609 825
which(full_array_median[,3] > 3) # 125 214 396 615 825 826
which(full_array_median[,4] > 3) # 125  214  396  442  499  524  609  615  626  703  808  825  829  916  919 1050
which(full_array_median[,5] > 3) # 35 542 609 615 743 826 923 990

#### This is for the correctly standardized environmental variables
setwd("/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation/results_avg232covarmatrix_hwe_5stanenviro")

library(abind)

stanrun1 <- read.table("bf1.txt")
stanrun1$V1 <- NULL
stanrun2 <- read.table("bf2.txt")
stanrun2$V1 <- NULL
stanrun3 <- read.table("bf3.txt")
stanrun3$V1 <- NULL

full_stanarray <- abind(stanrun1, stanrun2, stanrun3, along=3)

full_stanarray_median <- apply(full_stanarray, c(1,2), median)
full_stanarray_median

which(full_stanarray_median[,1] > 3) #125 214 609 615 825 919
which(full_stanarray_median[,2] > 3) #125  214  396  442  579  615  825  829  919 1081
which(full_stanarray_median[,3] > 3) #125  214  334  396  579  615  825  919 1081
which(full_stanarray_median[,4] > 3) #125  214  396  442  524  615  626  808  825  919 1050
which(full_stanarray_median[,5] > 3) #125 542 609 615 703 743 990

# pdf of histograms for each median BF for each environmental variable #
# different x-axis ranges depending on environmental variable
# pdf("histogramplots.pdf")
# 
# par(
#   mfrow = c(3, 2), 
#   mar=c(3.5, 3.5, 2, 1), # panel magin size in "line number" units
#   mgp=c(2, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
#   tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
#   cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
#   ps=10, # point size, which is the font size
#   bg=NA
# )
# 
# hist(full_array_median[,1], xlab="Bayes Factor", main="Latitude")
# rug(line=5, jitter(full_array_median[,1]))
# text(6,800, paste("max BF = 7.16"))
# hist(full_array_median[,2], xlab="Bayes Factor", main="Longitude")
# rug(jitter(full_array_median[,2]))
# text(8,800, paste("max BF = 9.24"))
# hist(full_array_median[,3], xlab="Bayes Factor", main="Depth")
# rug(jitter(full_array_median[,3]))
# text(5.5,700, paste("max BF = 6.79"))
# hist(full_array_median[,4], xlab="Bayes Factor", main="Bottom Temperature")
# rug(jitter(full_array_median[,4]))
# text(12,800, paste("max BF = 15.29"))
# hist(full_array_median[,5], xlab="Bayes Factor", main="Bottom Salinity")
# rug(jitter(full_array_median[,5]))
# text(8,800, paste("max BF = 9.87"))
# 
# dev.off()

# Same x-axis range for each environmental variable
png(file = "/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation/results_avg232covarmatrix_hwe_5enviro/histogramplots_same.png", width=8, height=11, res=300, units="in")

par(
  mfrow = c(3, 2), 
  mar=c(4.5, 3.5, 1, 1), # panel magin size in "line number" units
  mgp=c(2, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=10, # point size, which is the font size
  bg=NA
)

# hist(log(full_array_median[,1]), xlab="Bayes Factor", main="Latitude")

hist(full_array_median[,1], xlab="", ylab = "", main = "", breaks=16, xlim=c(0,16), xaxt="n")
axis(side=1, at=seq(0,16, 1), labels=seq(0,16,1), line=1.3)
rug(jitter(full_array_median[,1]), ticksize = -0.1, line=-0.2)
text(8,800, paste("max BF = 7.16"))
text(8,950, paste("Latitude"), cex = 1.3)
mtext("Bayes Factor", side = 1, line = 3.6)
mtext("Frequency", side = 2, line = 2.7)
hist(full_array_median[,2], xlab="", ylab = "", main = "", breaks=16, xlim=c(0,16), xaxt="n")
axis(side=1, at=seq(0,16, 1), labels=seq(0,16,1), line=1.3)
rug(jitter(full_array_median[,2]), ticksize = -0.1, line=-0.2)
text(8,800, paste("max BF = 9.24"))
text(8,950, paste("Longitude"), cex = 1.3)
mtext("Bayes Factor", side = 1, line = 3.6)
mtext("Frequency", side = 2, line = 2.7)
hist(full_array_median[,3], xlab="", ylab = "", main = "", breaks=16, xlim=c(0,16), xaxt="n")
axis(side=1, at=seq(0,16, 1), labels=seq(0,16,1), line=1.3)
rug(jitter(full_array_median[,3]), ticksize = -0.1, line=-0.2)
text(8,700, paste("max BF = 6.79"))
text(8,850, paste("Depth"), cex = 1.3)
mtext("Bayes Factor", side = 1, line = 3.6)
mtext("Frequency", side = 2, line = 2.7)
hist(full_array_median[,4], xlab="", ylab = "", main = "", breaks=16, xlim=c(0,16), xaxt="n")
axis(side=1, at=seq(0,16, 1), labels=seq(0,16,1), line=1.3)
rug(jitter(full_array_median[,4]), ticksize = -0.1, line=-0.2)
text(8,800, paste("max BF = 15.29"))
text(8,950, paste("Bottom Temperature"), cex = 1.3)
mtext("Bayes Factor", side = 1, line = 3.6)
mtext("Frequency", side = 2, line = 2.7)
hist(full_array_median[,5], xlab="", ylab = "", main = "", breaks=16, xlim=c(0,16), xaxt="n")
axis(side=1, at=seq(0,16, 1), labels=seq(0,16,1), line=1.3)
rug(jitter(full_array_median[,5]), ticksize = -0.1, line=-0.2)
text(8,600, paste("max BF = 9.87"))
text(8,750, paste("Bottom Salinity"), cex = 1.3)
mtext("Bayes Factor", side = 1, line = 3.6)
mtext("Frequency", side = 2, line = 2.7)

dev.off()

#######################################################################################################################################################################################
## Redundancy analysis to see how much variation in allele frequency can be explained by environmental variables ######################################################################
# Redundancy analysis (RDA) is a form of constrained ordination that examines how much of the variation in one set of variables explains the variation in another set of variables. ###
#######################################################################################################################################################################################

library(vegan)

vasc <- read.delim ('http://www.davidzeleny.net/anadat-r/data-download/vasc_plants.txt', row.names = 1)
chem <- read.delim ('http://www.davidzeleny.net/anadat-r/data-download/chemistry.txt', row.names = 1)
vasc.hell <- decostand (vasc, 'hell')  # transform species data by Hellinger transformation
rda.vasc <- rda (vasc.hell ~ ., chem)
plot(rda.vasc)

data(varespec)
data(varechem)

vare.rda <- rda(varespec, varechem)
vare.rda
plot(vare.rda)

# Bringing in adult genotype data
setwd("/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis")

library(ade4)
library(adegenet)
library(devtools)
library("hierfstat")
library(pegas)
library(fields)

# Reading in SNP data file containing only the first SNP at each locus
adults <- read.structure("structure_input_Nov_11_2015.str",
                         n.ind = 241, n.loc = 1137, col.lab = 1, col.pop = 2, row.marknames = 1, 
                         onerowperind = FALSE)

which(adults@loc.n.all > 2) # which snps have more than 2 alelles?

adults.nonas <- scaleGen(adults, center = TRUE, scale = FALSE, NA.method = "mean") # filling in NAs with allele frequencies and centering
# adults.nonas <- scaleGen(adults, center = FALSE, scale = FALSE, NA.method = "mean") # based on Meirimans 2015 and others, it's okay to perform rda on allele frequencies
hist(adults.nonas)
# adults.nonas <- scaleGen(adults, center = TRUE, scale = TRUE, NA.method = "mean") # this fills in NAs w/ allele frequency means, so now everything is in allele frequencies, another way of transforming data using z-scores and allele frequencies
sum(is.na(adults.nonas)) #0 All NAs have been replaced with means

# adults.hell <- decostand(adults.nonas, 'hell')  # transform allele frequency data by Hellinger transformation
# rownames(adults.hell) # lists rownames b/c these need to be in same order for genetic and environmental matrix
# adults.hell <- as.data.frame(adults.hell)

exclu_names <- c("PADE_14230L1439",
                 "PADE_14231L1440",
                 "PADE_14232L1529",
                 "PADE_14233L1588",
                 "PADE_14234L1441",
                 "PADE_14235L1442",
                 "PADE_14236L1530",
                 "PADE_14237L1531",
                 "PADE_14238L1532") # Creating a list of IDs to exclude from the genetic matrix

adults.nonas.232 <- adults.nonas[ ! rownames(adults.nonas) %in% exclu_names, ] # This is genotype data for 232 fish where NAs have been replaced by mean allele frequencies
rownames(adults.nonas.232)
adults.nonas.232.2274 <- adults.nonas.232[, -c(56, 96, 734, 2115)] # SNP 28, 47, 366 and 1056 have 3 alleles. All the rest have 2. Need to remove 3rd allele for SNP 28, 47, 366 and 1056. Removing one with lowest count
dim(adults.nonas.232.2274)
adults.nonas.232.1137 <- adults.nonas.232.2274[, seq(1, ncol(adults.nonas.232.2274),by = 2)] # including every other column
# adults.nonas.232.1137 <- adults.nonas.232.2274[, -seq(1, ncol(adults.nonas.232.2274),by = 2)] # excluding every other column
dim(adults.nonas.232.1137)
# adults.hell.232 <- adults.hell[ ! rownames(adults.hell) %in% exclu_names, ] # This is genotype data for 232 fish where NAs have been replaced by mean allele frequencies and then a Hellinger transformation has been applied
# rownames(adults.hell.232)

# Bringing in environmental data
setwd("/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation")

envi <- read.csv("allpops_combo.csv", header=TRUE)
envi <- envi[-c(41, 96, 138, 199),] # excludes lines with no data
envi.ordered <- envi[with(envi, order(PinskyID)),]
envi.ordered[,2]
envi.ordered
envi.ordered.matrix <- scale(data.matrix(envi.ordered))
envi.ordered.matrix
rownames(envi.ordered.matrix) <- envi.ordered[,"PinskyID"]
envi.ordered.matrix[,1] <- envi.ordered[,1]

# Let's plot environmental variables against each other
pairs(envi.ordered.matrix, main="Bivariate plots")

mod1 <- lm(envi.ordered.matrix[,5] ~ envi.ordered.matrix[,6])
mod2 <- lm(envi.ordered.matrix[,5] ~ envi.ordered.matrix[,7])
mod3 <- lm(envi.ordered.matrix[,5] ~ envi.ordered.matrix[,8])
mod4 <- lm(envi.ordered.matrix[,5] ~ envi.ordered.matrix[,9])
mod5 <- lm(envi.ordered.matrix[,6] ~ envi.ordered.matrix[,7])
mod6 <- lm(envi.ordered.matrix[,6] ~ envi.ordered.matrix[,8])
mod7 <- lm(envi.ordered.matrix[,7] ~ envi.ordered.matrix[,8])
mod8 <- lm(envi.ordered.matrix[,7] ~ envi.ordered.matrix[,9])
mod9 <- lm(envi.ordered.matrix[,8] ~ envi.ordered.matrix[,9])

# Are individuals ordered in the same way in the genetic and environmental data frames? Here are multiple ways to test this.
rownames(adults.nonas.232.1137) == rownames(envi.ordered.matrix)
# rownames(adults.nonas.232) == rownames(envi.ordered.matrix)
# rownames(adults.nonas.232) == envi.ordered[,2]

all.equal(rownames(adults.nonas.232.1137), rownames(envi.ordered.matrix))
# all.equal(rownames(adults.nonas.232), rownames(envi.ordered.matrix))
# all.equal(rownames(adults.nonas.232), envi.ordered[,2])

# rownames(adults.hell.232)[229] # same order? let's test
# envi.ordered[229,2]

# Performing full RDA; pay attention to whether genetic matrix is Hellinger transformed or not
envi.ordered.matrix <- as.data.frame(envi.ordered.matrix)

adults.rda <- rda(adults.nonas.232.1137 ~ lat + long + depth + b_temp + b_salin, envi.ordered.matrix, scale = FALSE)
adults.rda
summary(adults.rda)
RsquareAdj (adults.rda)
R2a.adults.rda <- RsquareAdj(adults.rda)$adj.r.squared
coef(adults.rda)
plot(adults.rda, scaling = 3)
plot(adults.rda, xlim = c(-0.1,0.1), ylim = c(-0.1,0.1))
anova(adults.rda)
anova(adults.rda, by="axis", step=1000)

plot(adults.rda, choices = c(1,2), scaling=2)
plot(adults.rda, choices = c(1,3), xlim = c(-0.6,0.6), ylim = c(-0.6,0.6), scaling=3)
plot(adults.rda, choices = c(1,2), xlim = c(-0.6,0.6), ylim = c(-0.6,0.6), scaling=3)
plot(adults.rda, choices = c(2,3), xlim = c(-0.6,0.6), ylim = c(-0.6,0.6), scaling=3)
plot(adults.rda, choices = c(1,4), xlim = c(-0.6,0.6), ylim = c(-0.6,0.6), scaling=3)
plot(adults.rda, choices = c(2,4), xlim = c(-0.6,0.6), ylim = c(-0.6,0.6), scaling=3)
plot(adults.rda, choices = c(3,4), xlim = c(-0.6,0.6), ylim = c(-0.6,0.6), scaling=3)
plot(adults.rda, choices = c(1,5), xlim = c(-0.6,0.6), ylim = c(-0.6,0.6), scaling=3)
plot(adults.rda, choices = c(2,5), xlim = c(-0.6,0.6), ylim = c(-0.6,0.6), scaling=3)
plot(adults.rda, choices = c(3,5), xlim = c(-0.6,0.6), ylim = c(-0.6,0.6), scaling=3)
plot(adults.rda, choices = c(4,5), xlim = c(-0.6,0.6), ylim = c(-0.6,0.6), scaling=3)

spp.scr <- scores(adults.rda, display = "species", scaling = 2, choices = c(1,2,3,4,5))
arrows(0,0, spp.scr[,1], spp.scr[,2], length = 0, lty=1, col = "blue")

rainbow(22)
points(spp.scr[35,1], spp.scr[35,3], col="#FF0000FF", pch = 16)
points(spp.scr[125,1], spp.scr[125,3], col="#FF4600FF", pch = 16)
points(spp.scr[214,1], spp.scr[214,3], col="#FF8B00FF", pch = 16)
points(spp.scr[396,1], spp.scr[396,3], col="#FFD100FF", pch = 16)
points(spp.scr[442,1], spp.scr[442,3], col="#E8FF00FF", pch = 16)
points(spp.scr[499,1], spp.scr[499,3], col="#A2FF00FF", pch = 16)
points(spp.scr[524,1], spp.scr[524,3], col="#5DFF00FF", pch = 16)
points(spp.scr[542,1], spp.scr[542,3], col="#17FF00FF", pch = 16)
points(spp.scr[609,1], spp.scr[609,3], col="#00FF2EFF", pch = 16)
points(spp.scr[615,1], spp.scr[615,3], col="#00FF74FF", pch = 16)
points(spp.scr[626,1], spp.scr[626,3], col="#00FFB9FF", pch = 16)
points(spp.scr[703,1], spp.scr[703,3], col="#00FFFFFF", pch = 16)
points(spp.scr[743,1], spp.scr[743,3], col="#00B9FFFF", pch = 16)
points(spp.scr[808,1], spp.scr[808,3], col="#0074FFFF", pch = 16)
points(spp.scr[825,1], spp.scr[825,3], col="#002EFFFF", pch = 16)
points(spp.scr[826,1], spp.scr[826,3], col="#1700FFFF", pch = 16)
points(spp.scr[829,1], spp.scr[829,3], col="#5D00FFFF", pch = 16)
points(spp.scr[916,1], spp.scr[916,3], col="#A200FFFF", pch = 16)
points(spp.scr[919,1], spp.scr[919,3], col="#E800FFFF", pch = 16)
points(spp.scr[923,1], spp.scr[923,3], col="#FF00D1FF", pch = 16)
points(spp.scr[990,1], spp.scr[990,3], col="#FF008BFF", pch = 16)
points(spp.scr[1050,1], spp.scr[1050,3], col="#FF0046FF", pch = 16)

spp.scr <- scores(adults.rda, display = "species", scaling = 3, choices = c(1,2,3,4,5))
points(spp.scr[35,1], spp.scr[35,2], col="#FF0000FF", pch = 16)
points(spp.scr[125,1], spp.scr[125,2], col="#FF4600FF", pch = 16)
points(spp.scr[214,1], spp.scr[214,2], col="#FF8B00FF", pch = 16)
points(spp.scr[396,1], spp.scr[396,2], col="#FFD100FF", pch = 16)
points(spp.scr[442,1], spp.scr[442,2], col="#E8FF00FF", pch = 16)
points(spp.scr[499,1], spp.scr[499,2], col="#A2FF00FF", pch = 16)
points(spp.scr[524,1], spp.scr[524,2], col="#5DFF00FF", pch = 16)
points(spp.scr[542,1], spp.scr[542,2], col="#17FF00FF", pch = 16)
points(spp.scr[609,1], spp.scr[609,2], col="#00FF2EFF", pch = 16)
points(spp.scr[615,1], spp.scr[615,2], col="#00FF74FF", pch = 16)
points(spp.scr[626,1], spp.scr[626,2], col="#00FFB9FF", pch = 16)
points(spp.scr[703,1], spp.scr[703,2], col="#00FFFFFF", pch = 16)
points(spp.scr[743,1], spp.scr[743,2], col="#00B9FFFF", pch = 16)
points(spp.scr[808,1], spp.scr[808,2], col="#0074FFFF", pch = 16)
points(spp.scr[825,1], spp.scr[825,2], col="#002EFFFF", pch = 16)
points(spp.scr[826,1], spp.scr[826,2], col="#1700FFFF", pch = 16)
points(spp.scr[829,1], spp.scr[829,2], col="#5D00FFFF", pch = 16)
points(spp.scr[916,1], spp.scr[916,2], col="#A200FFFF", pch = 16)
points(spp.scr[919,1], spp.scr[919,2], col="#E800FFFF", pch = 16)
points(spp.scr[923,1], spp.scr[923,2], col="#FF00D1FF", pch = 16)
points(spp.scr[990,1], spp.scr[990,2], col="#FF008BFF", pch = 16)
points(spp.scr[1050,1], spp.scr[1050,2], col="#FF0046FF", pch = 16)

plot(adults.rda, choices = c(1,2), scaling=3)
spp.scr <- scores(adults.rda, display = "species", scaling = 3, choices = c(1,2,3,4,5))
can_loci <- spp.scr[c(35, 125, 214, 396, 442, 499, 524, 542, 609, 615, 626, 703, 743, 808, 825, 826, 829, 916, 919, 923, 990, 1050),]
for (snp in can_loci){
  points(can_loci, col = "darkblue", pch = 16, cex = 0.8)
  arrows(0,0, can_loci[,1], can_loci[,2], length = 0, lty=1, col = "darkblue", cex = 0.8)
}

points(spp.scr[35,1], spp.scr[35,2], col="darkblue", pch = 16)
points(spp.scr[125,1], spp.scr[125,2], col="darkblue", pch = 16)
points(spp.scr[214,1], spp.scr[214,2], col="darkblue", pch = 16)
points(spp.scr[396,1], spp.scr[396,2], col="darkblue", pch = 16)
points(spp.scr[442,1], spp.scr[442,2], col="darkblue", pch = 16)
points(spp.scr[499,1], spp.scr[499,2], col="darkblue", pch = 16)
points(spp.scr[524,1], spp.scr[524,2], col="darkblue", pch = 16)
points(spp.scr[542,1], spp.scr[542,2], col="darkblue", pch = 16)
points(spp.scr[609,1], spp.scr[609,2], col="darkblue", pch = 16)
points(spp.scr[615,1], spp.scr[615,2], col="darkblue", pch = 16)
points(spp.scr[626,1], spp.scr[626,2], col="darkblue", pch = 16)
points(spp.scr[703,1], spp.scr[703,2], col="darkblue", pch = 16)
points(spp.scr[743,1], spp.scr[743,2], col="darkblue", pch = 16)
points(spp.scr[808,1], spp.scr[808,2], col="darkblue", pch = 16)
points(spp.scr[825,1], spp.scr[825,2], col="darkblue", pch = 16)
points(spp.scr[826,1], spp.scr[826,2], col="darkblue", pch = 16)
points(spp.scr[829,1], spp.scr[829,2], col="darkblue", pch = 16)
points(spp.scr[916,1], spp.scr[916,2], col="darkblue", pch = 16)
points(spp.scr[919,1], spp.scr[919,2], col="darkblue", pch = 16)
points(spp.scr[923,1], spp.scr[923,2], col="darkblue", pch = 16)
points(spp.scr[990,1], spp.scr[990,2], col="darkblue", pch = 16)
points(spp.scr[1050,1], spp.scr[1050,2], col="darkblue", pch = 16)

plot(adults.rda, choices = c(1,3), scaling=3)
spp.scr <- scores(adults.rda, display = "species", scaling = 3, choices = c(1,2,3,4,5))
can_loci <- spp.scr[c(35, 125, 214, 396, 442, 499, 524, 542, 609, 615, 626, 703, 743, 808, 825, 826, 829, 916, 919, 923, 990, 1050), c(1,3)]
for (snp in can_loci){
  points(can_loci, col = "darkblue", pch = 16)
  arrows(0,0, can_loci[,1], can_loci[,2], length = 0, lty=1, col = "darkblue", cex = 0.8)
}

points(spp.scr[35,1], spp.scr[35,3], col="darkblue", pch = 16)
points(spp.scr[125,1], spp.scr[125,3], col="darkblue", pch = 16)
points(spp.scr[214,1], spp.scr[214,3], col="darkblue", pch = 16)
points(spp.scr[396,1], spp.scr[396,3], col="darkblue", pch = 16)
points(spp.scr[442,1], spp.scr[442,3], col="darkblue", pch = 16)
points(spp.scr[499,1], spp.scr[499,3], col="darkblue", pch = 16)
points(spp.scr[524,1], spp.scr[524,3], col="darkblue", pch = 16)
points(spp.scr[542,1], spp.scr[542,3], col="darkblue", pch = 16)
points(spp.scr[609,1], spp.scr[609,3], col="darkblue", pch = 16)
points(spp.scr[615,1], spp.scr[615,3], col="darkblue", pch = 16)
points(spp.scr[626,1], spp.scr[626,3], col="darkblue", pch = 16)
points(spp.scr[703,1], spp.scr[703,3], col="darkblue", pch = 16)
points(spp.scr[743,1], spp.scr[743,3], col="darkblue", pch = 16)
points(spp.scr[808,1], spp.scr[808,3], col="darkblue", pch = 16)
points(spp.scr[825,1], spp.scr[825,3], col="darkblue", pch = 16)
points(spp.scr[826,1], spp.scr[826,3], col="darkblue", pch = 16)
points(spp.scr[829,1], spp.scr[829,3], col="darkblue", pch = 16)
points(spp.scr[916,1], spp.scr[916,3], col="darkblue", pch = 16)
points(spp.scr[919,1], spp.scr[919,3], col="darkblue", pch = 16)
points(spp.scr[923,1], spp.scr[923,3], col="darkblue", pch = 16)
points(spp.scr[990,1], spp.scr[990,3], col="darkblue", pch = 16)
points(spp.scr[1050,1], spp.scr[1050,3], col="darkblue", pch = 16)

# Plotting species and biplot scores
plot(adults.rda, display=c("sp", "bp"))
ordispider(adults.rda)

# Let's look at the loadings for each allele
spp.scr <- scores(adults.rda, display = "species", scaling = 0, choices = c(1,2,3,4,5))
head(spp.scr)
allele_loadings <- spp.scr[c(35, 125, 214, 396, 442, 499, 524, 542, 609, 615, 626, 703, 743, 808, 825, 826, 829, 916, 919, 923, 990, 1050),]
allele_loadings # unscaled species scores

# Forester et al. (2016) used loci with scores +/- 3 SD from mean score for that axis to ID outlier loci. Only used first 3 axes.
mean.rda <- colMeans(spp.scr)
sd.rda1 <- sd(spp.scr[,"RDA1"])
sd.rda2 <- sd(spp.scr[,"RDA2"])
sd.rda3 <- sd(spp.scr[,"RDA3"])
sd.rda4 <- sd(spp.scr[,"RDA4"])
sd.rda5 <- sd(spp.scr[,"RDA5"])

rda1.hi <- mean.rda[1] + 3*sd.rda1
rda1.lo <- mean.rda[1] - 3*sd.rda1
which(spp.scr[,"RDA1"] > rda1.hi)
which(spp.scr[,"RDA1"] < rda1.lo)

rda2.hi <- mean.rda[2] + 3*sd.rda2
rda2.lo <- mean.rda[2] - 3*sd.rda2
which(spp.scr[,"RDA2"] > rda2.hi)
which(spp.scr[,"RDA2"] < rda2.lo)

rda3.hi <- mean.rda[3] + 3*sd.rda3
rda3.lo <- mean.rda[3] - 3*sd.rda3
which(spp.scr[,"RDA3"] > rda3.hi)
which(spp.scr[,"RDA3"] < rda3.lo)

# rda4.hi <- mean.rda[4] + 3*sd.rda4
# rda4.lo <- mean.rda[4] - 3*sd.rda4
# which(spp.scr[,"RDA4"] > rda4.hi)
# which(spp.scr[,"RDA4"] < rda4.lo)
# 
# rda5.hi <- mean.rda[5] + 3*sd.rda5
# rda5.lo <- mean.rda[5] - 3*sd.rda5
# which(spp.scr[,"RDA5"] > rda5.hi)
# which(spp.scr[,"RDA5"] < rda5.lo)

rda.cans <- c((which(spp.scr[,"RDA1"] > rda1.hi)),
              (which(spp.scr[,"RDA1"] < rda1.lo)),
              (which(spp.scr[,"RDA2"] > rda2.hi)),
              (which(spp.scr[,"RDA2"] < rda2.lo)),
              (which(spp.scr[,"RDA3"] > rda3.hi)),
              (which(spp.scr[,"RDA3"] < rda3.lo))
#               (which(spp.scr[,"RDA4"] > rda4.hi)),
#               (which(spp.scr[,"RDA4"] < rda4.lo)),
#               (which(spp.scr[,"RDA5"] > rda5.hi)),
#               (which(spp.scr[,"RDA5"] < rda5.lo))
              )

rda.cans <- unique(sort(rda.cans)) # 125  128  291  380  402  432  466  577  579  639  691  694  695  704  720  743  826  842  937  947  997 1035 1061 1090

# Linear regression between allele frequencies and environmental variables for the RDA outliers
lat125 <- lm(adults.nonas.232.1137[,125] ~ envi.ordered.matrix[,"lat"])
long125 <- lm(adults.nonas.232.1137[,125] ~ envi.ordered.matrix[,"long"])
depth125 <- lm(adults.nonas.232.1137[,125] ~ envi.ordered.matrix[,"depth"])
btemp125 <- lm(adults.nonas.232.1137[,125] ~ envi.ordered.matrix[,"b_temp"])
bsalin125 <- lm(adults.nonas.232.1137[,125] ~ envi.ordered.matrix[,"b_salin"]) # 0.000149

lat128 <- lm(adults.nonas.232.1137[,128] ~ envi.ordered.matrix[,"lat"])
long128 <- lm(adults.nonas.232.1137[,128] ~ envi.ordered.matrix[,"long"])
depth128 <- lm(adults.nonas.232.1137[,128] ~ envi.ordered.matrix[,"depth"])
btemp128 <- lm(adults.nonas.232.1137[,128] ~ envi.ordered.matrix[,"b_temp"])
bsalin128 <- lm(adults.nonas.232.1137[,128] ~ envi.ordered.matrix[,"b_salin"])

lat291 <- lm(adults.nonas.232.1137[,291] ~ envi.ordered.matrix[,"lat"])
long291 <- lm(adults.nonas.232.1137[,291] ~ envi.ordered.matrix[,"long"])
depth291 <- lm(adults.nonas.232.1137[,291] ~ envi.ordered.matrix[,"depth"]) # 2.765e-05
btemp291 <- lm(adults.nonas.232.1137[,291] ~ envi.ordered.matrix[,"b_temp"]) # 2.62e-05
bsalin291 <- lm(adults.nonas.232.1137[,291] ~ envi.ordered.matrix[,"b_salin"])

lat380 <- lm(adults.nonas.232.1137[,380] ~ envi.ordered.matrix[,"lat"])
long380 <- lm(adults.nonas.232.1137[,380] ~ envi.ordered.matrix[,"long"])
depth380 <- lm(adults.nonas.232.1137[,380] ~ envi.ordered.matrix[,"depth"]) 
btemp380 <- lm(adults.nonas.232.1137[,380] ~ envi.ordered.matrix[,"b_temp"]) 
bsalin380 <- lm(adults.nonas.232.1137[,380] ~ envi.ordered.matrix[,"b_salin"])

lat402 <- lm(adults.nonas.232.1137[,402] ~ envi.ordered.matrix[,"lat"])
long402 <- lm(adults.nonas.232.1137[,402] ~ envi.ordered.matrix[,"long"])
depth402 <- lm(adults.nonas.232.1137[,402] ~ envi.ordered.matrix[,"depth"]) 
btemp402 <- lm(adults.nonas.232.1137[,402] ~ envi.ordered.matrix[,"b_temp"]) 
bsalin402 <- lm(adults.nonas.232.1137[,402] ~ envi.ordered.matrix[,"b_salin"])

lat432 <- lm(adults.nonas.232.1137[,432] ~ envi.ordered.matrix[,"lat"])
long432 <- lm(adults.nonas.232.1137[,432] ~ envi.ordered.matrix[,"long"])
depth432 <- lm(adults.nonas.232.1137[,432] ~ envi.ordered.matrix[,"depth"]) 
btemp432 <- lm(adults.nonas.232.1137[,432] ~ envi.ordered.matrix[,"b_temp"]) 
bsalin432 <- lm(adults.nonas.232.1137[,432] ~ envi.ordered.matrix[,"b_salin"])

lat466 <- lm(adults.nonas.232.1137[,466] ~ envi.ordered.matrix[,"lat"])
long466 <- lm(adults.nonas.232.1137[,466] ~ envi.ordered.matrix[,"long"])
depth466 <- lm(adults.nonas.232.1137[,466] ~ envi.ordered.matrix[,"depth"]) 
btemp466 <- lm(adults.nonas.232.1137[,466] ~ envi.ordered.matrix[,"b_temp"]) 
bsalin466 <- lm(adults.nonas.232.1137[,466] ~ envi.ordered.matrix[,"b_salin"])

lat577 <- lm(adults.nonas.232.1137[,577] ~ envi.ordered.matrix[,"lat"])
long577 <- lm(adults.nonas.232.1137[,577] ~ envi.ordered.matrix[,"long"])
depth577 <- lm(adults.nonas.232.1137[,577] ~ envi.ordered.matrix[,"depth"]) 
btemp577 <- lm(adults.nonas.232.1137[,577] ~ envi.ordered.matrix[,"b_temp"]) 
bsalin577 <- lm(adults.nonas.232.1137[,577] ~ envi.ordered.matrix[,"b_salin"])

lat579 <- lm(adults.nonas.232.1137[,579] ~ envi.ordered.matrix[,"lat"])
long579 <- lm(adults.nonas.232.1137[,579] ~ envi.ordered.matrix[,"long"])
depth579 <- lm(adults.nonas.232.1137[,579] ~ envi.ordered.matrix[,"depth"]) 
btemp579 <- lm(adults.nonas.232.1137[,579] ~ envi.ordered.matrix[,"b_temp"]) 
bsalin579 <- lm(adults.nonas.232.1137[,579] ~ envi.ordered.matrix[,"b_salin"])

lat639 <- lm(adults.nonas.232.1137[,639] ~ envi.ordered.matrix[,"lat"])
long639 <- lm(adults.nonas.232.1137[,639] ~ envi.ordered.matrix[,"long"])
depth639 <- lm(adults.nonas.232.1137[,639] ~ envi.ordered.matrix[,"depth"]) 
btemp639 <- lm(adults.nonas.232.1137[,639] ~ envi.ordered.matrix[,"b_temp"]) 
bsalin639 <- lm(adults.nonas.232.1137[,639] ~ envi.ordered.matrix[,"b_salin"])

lat691 <- lm(adults.nonas.232.1137[,691] ~ envi.ordered.matrix[,"lat"])
long691 <- lm(adults.nonas.232.1137[,691] ~ envi.ordered.matrix[,"long"])
depth691 <- lm(adults.nonas.232.1137[,691] ~ envi.ordered.matrix[,"depth"]) 
btemp691 <- lm(adults.nonas.232.1137[,691] ~ envi.ordered.matrix[,"b_temp"]) 
bsalin691 <- lm(adults.nonas.232.1137[,691] ~ envi.ordered.matrix[,"b_salin"])

lat694 <- lm(adults.nonas.232.1137[,694] ~ envi.ordered.matrix[,"lat"])
long694 <- lm(adults.nonas.232.1137[,694] ~ envi.ordered.matrix[,"long"])
depth694 <- lm(adults.nonas.232.1137[,694] ~ envi.ordered.matrix[,"depth"]) 
btemp694 <- lm(adults.nonas.232.1137[,694] ~ envi.ordered.matrix[,"b_temp"]) 
bsalin694 <- lm(adults.nonas.232.1137[,694] ~ envi.ordered.matrix[,"b_salin"])

lat695 <- lm(adults.nonas.232.1137[,695] ~ envi.ordered.matrix[,"lat"])
long695 <- lm(adults.nonas.232.1137[,695] ~ envi.ordered.matrix[,"long"])
depth695 <- lm(adults.nonas.232.1137[,695] ~ envi.ordered.matrix[,"depth"]) 
btemp695 <- lm(adults.nonas.232.1137[,695] ~ envi.ordered.matrix[,"b_temp"]) 
bsalin695 <- lm(adults.nonas.232.1137[,695] ~ envi.ordered.matrix[,"b_salin"])

lat704 <- lm(adults.nonas.232.1137[,704] ~ envi.ordered.matrix[,"lat"])
long704 <- lm(adults.nonas.232.1137[,704] ~ envi.ordered.matrix[,"long"])
depth704 <- lm(adults.nonas.232.1137[,704] ~ envi.ordered.matrix[,"depth"]) 
btemp704 <- lm(adults.nonas.232.1137[,704] ~ envi.ordered.matrix[,"b_temp"]) 
bsalin704 <- lm(adults.nonas.232.1137[,704] ~ envi.ordered.matrix[,"b_salin"])

lat720 <- lm(adults.nonas.232.1137[,720] ~ envi.ordered.matrix[,"lat"])
long720 <- lm(adults.nonas.232.1137[,720] ~ envi.ordered.matrix[,"long"])
depth720 <- lm(adults.nonas.232.1137[,720] ~ envi.ordered.matrix[,"depth"]) 
btemp720 <- lm(adults.nonas.232.1137[,720] ~ envi.ordered.matrix[,"b_temp"]) 
bsalin720 <- lm(adults.nonas.232.1137[,720] ~ envi.ordered.matrix[,"b_salin"])

lat743 <- lm(adults.nonas.232.1137[,743] ~ envi.ordered.matrix[,"lat"])
long743 <- lm(adults.nonas.232.1137[,743] ~ envi.ordered.matrix[,"long"])
depth743 <- lm(adults.nonas.232.1137[,743] ~ envi.ordered.matrix[,"depth"]) 
btemp743 <- lm(adults.nonas.232.1137[,743] ~ envi.ordered.matrix[,"b_temp"]) 
bsalin743 <- lm(adults.nonas.232.1137[,743] ~ envi.ordered.matrix[,"b_salin"]) #0.0006089858

lat826 <- lm(adults.nonas.232.1137[,826] ~ envi.ordered.matrix[,"lat"])
long826 <- lm(adults.nonas.232.1137[,826] ~ envi.ordered.matrix[,"long"])
depth826 <- lm(adults.nonas.232.1137[,826] ~ envi.ordered.matrix[,"depth"]) 
btemp826 <- lm(adults.nonas.232.1137[,826] ~ envi.ordered.matrix[,"b_temp"]) 
bsalin826 <- lm(adults.nonas.232.1137[,826] ~ envi.ordered.matrix[,"b_salin"])

lat842 <- lm(adults.nonas.232.1137[,842] ~ envi.ordered.matrix[,"lat"])
long842 <- lm(adults.nonas.232.1137[,842] ~ envi.ordered.matrix[,"long"])
depth842 <- lm(adults.nonas.232.1137[,842] ~ envi.ordered.matrix[,"depth"]) 
btemp842 <- lm(adults.nonas.232.1137[,842] ~ envi.ordered.matrix[,"b_temp"]) 
bsalin842 <- lm(adults.nonas.232.1137[,842] ~ envi.ordered.matrix[,"b_salin"])

lat937 <- lm(adults.nonas.232.1137[,937] ~ envi.ordered.matrix[,"lat"])
long937 <- lm(adults.nonas.232.1137[,937] ~ envi.ordered.matrix[,"long"])
depth937 <- lm(adults.nonas.232.1137[,937] ~ envi.ordered.matrix[,"depth"]) 
btemp937 <- lm(adults.nonas.232.1137[,937] ~ envi.ordered.matrix[,"b_temp"]) 
bsalin937 <- lm(adults.nonas.232.1137[,937] ~ envi.ordered.matrix[,"b_salin"])

lat947 <- lm(adults.nonas.232.1137[,947] ~ envi.ordered.matrix[,"lat"])
long947 <- lm(adults.nonas.232.1137[,947] ~ envi.ordered.matrix[,"long"])
depth947 <- lm(adults.nonas.232.1137[,947] ~ envi.ordered.matrix[,"depth"]) 
btemp947 <- lm(adults.nonas.232.1137[,947] ~ envi.ordered.matrix[,"b_temp"]) 
bsalin947 <- lm(adults.nonas.232.1137[,947] ~ envi.ordered.matrix[,"b_salin"])

lat997 <- lm(adults.nonas.232.1137[,997] ~ envi.ordered.matrix[,"lat"])
long997 <- lm(adults.nonas.232.1137[,997] ~ envi.ordered.matrix[,"long"])
depth997 <- lm(adults.nonas.232.1137[,997] ~ envi.ordered.matrix[,"depth"]) 
btemp997 <- lm(adults.nonas.232.1137[,997] ~ envi.ordered.matrix[,"b_temp"]) 
bsalin997 <- lm(adults.nonas.232.1137[,997] ~ envi.ordered.matrix[,"b_salin"])

lat1035 <- lm(adults.nonas.232.1137[,1035] ~ envi.ordered.matrix[,"lat"])
long1035 <- lm(adults.nonas.232.1137[,1035] ~ envi.ordered.matrix[,"long"])
depth1035 <- lm(adults.nonas.232.1137[,1035] ~ envi.ordered.matrix[,"depth"]) 
btemp1035 <- lm(adults.nonas.232.1137[,1035] ~ envi.ordered.matrix[,"b_temp"]) 
bsalin1035 <- lm(adults.nonas.232.1137[,1035] ~ envi.ordered.matrix[,"b_salin"])

lat1061 <- lm(adults.nonas.232.1137[,1061] ~ envi.ordered.matrix[,"lat"])
long1061 <- lm(adults.nonas.232.1137[,1061] ~ envi.ordered.matrix[,"long"])
depth1061 <- lm(adults.nonas.232.1137[,1061] ~ envi.ordered.matrix[,"depth"]) 
btemp1061 <- lm(adults.nonas.232.1137[,1061] ~ envi.ordered.matrix[,"b_temp"]) 
bsalin1061 <- lm(adults.nonas.232.1137[,1061] ~ envi.ordered.matrix[,"b_salin"])

lat1090 <- lm(adults.nonas.232.1137[,1090] ~ envi.ordered.matrix[,"lat"])
long1090 <- lm(adults.nonas.232.1137[,1090] ~ envi.ordered.matrix[,"long"])
depth1090 <- lm(adults.nonas.232.1137[,1090] ~ envi.ordered.matrix[,"depth"]) 
btemp1090 <- lm(adults.nonas.232.1137[,1090] ~ envi.ordered.matrix[,"b_temp"]) 
bsalin1090 <- lm(adults.nonas.232.1137[,1090] ~ envi.ordered.matrix[,"b_salin"])

# Function to calculate p-values
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

pvals <- c(lmp(lat125), lmp(long125), lmp(depth125), lmp(btemp125),  lmp(bsalin125), lmp(lat128), lmp(long128), lmp(depth128), lmp(btemp128),  lmp(bsalin128), lmp(lat291), lmp(long291), lmp(depth291), lmp(btemp291), lmp(bsalin291), 
             lmp(lat380), lmp(long380), lmp(depth380),  lmp(btemp380),  lmp(bsalin380), lmp(lat402), lmp(long402), lmp(depth402), lmp(btemp402), lmp(bsalin402), lmp(lat432), lmp(long432), lmp(depth432), lmp(btemp432), lmp(bsalin432), lmp(lat466), 
             lmp(long466), lmp(depth466), lmp(btemp466), lmp(bsalin466), lmp(lat577), lmp(long577), lmp(depth577), lmp(btemp577), lmp(bsalin577), lmp(lat579), lmp(long579), lmp(depth579), lmp(btemp579), lmp(bsalin579), lmp(lat639), lmp(long639), 
             lmp(depth639), lmp(btemp639), lmp(bsalin639), lmp(lat691), lmp(long691), lmp(depth691), lmp(btemp691), lmp(bsalin691), lmp(lat694), lmp(long694), lmp(depth694), lmp(btemp694), lmp(bsalin694), lmp(lat695), lmp(long695), lmp(depth695),  
             lmp(btemp695), lmp(bsalin695), lmp(lat704), lmp(long704), lmp(depth704), lmp(btemp704), lmp(bsalin704), lmp(lat720), lmp(long720), lmp(depth720), lmp(btemp720), lmp(bsalin720), lmp(lat743), lmp(long743), lmp(depth743), lmp(btemp743), 
             lmp(bsalin743), lmp(lat826), lmp(long826), lmp(depth826), lmp(btemp826), lmp(bsalin826), lmp(lat842), lmp(long842), lmp(depth842), lmp(btemp842), lmp(bsalin842), lmp(lat937), lmp(long937),  lmp(depth937),  lmp(btemp937), lmp(bsalin937), 
             lmp(lat947), lmp(long947), lmp(depth947), lmp(btemp947), lmp(bsalin947), lmp(lat997), lmp(long997), lmp(depth997), lmp(btemp997), lmp(bsalin997), lmp(lat1035), lmp(long1035), lmp(depth1035), lmp(btemp1035), lmp(bsalin1035), lmp(lat1061), 
             lmp(long1061), lmp(depth1061), lmp(btemp1061), lmp(bsalin1061), lmp(lat1090), lmp(long1090), lmp(depth1090), lmp(btemp1090), lmp(bsalin1090))

which(pvals < 0.001)

# Plotting the RDA plot with the three loci indicated by redundancy analysis to have strong locus-environmental associations
png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Redundancy analysis/RDA plot.png", width=11, height=6, res=300, units="in")

par(
  mfrow = c(1, 2), 
  mar=c(5, 5, 4, 2), # panel magin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  bg=NA
)

plot(adults.rda, choices = c(1,2), scaling=3, type = "none")
points(adults.rda, choices = c(1,2), display = "sites", col = "gray80", scaling=3, cex = 0.6)
points(adults.rda, choices = c(1,2), display = "sp", col = "gray50", scaling=3, pch = 3, cex = 0.6)
points(adults.rda, choices = c(1,2), display = "bp", scaling=3)
text(adults.rda, choices = c(1,2), display = "bp", scaling=3)

spp.scr <- scores(adults.rda, display = "species", scaling = 3, choices = c(1,2,3,4,5))
# site.scr <- scores(adults.rda, display = "sites", scaling = 3, choices = c(1,2,3,4,5))
# bp.scr <- scores(adults.rda, display = "bp", scaling = 3, choices = c(1,2,3,4,5))
points(spp.scr[125,1], spp.scr[125,2], col = "black", pch = 16, cex = 1.1)
points(spp.scr[291,1], spp.scr[291,2], col = "black", pch = 17, cex = 1.1)
points(spp.scr[743,1], spp.scr[743,2], col = "black", pch = 15, cex = 1.1)
legend("topleft",
	legend=c("contig 8420", "contig 19728", "contig 54288"),
	pch=c(16, 17, 15),
	col=c("black", "black", "black"))
legend("bottomleft",
       legend = expression("A"),
       bty = "n", 
       cex = 1.3)

plot(adults.rda, choices = c(1,3), scaling=3, type = "none")
points(adults.rda, choices = c(1,3), display = "sites", col = "gray80", scaling=3, cex = 0.6)
points(adults.rda, choices = c(1,3), display = "sp", col = "gray50", scaling=3, pch = 3, cex = 0.6)
points(adults.rda, choices = c(1,3), display = "bp", scaling=3)
text(adults.rda, choices = c(1,3), display = "bp", scaling=3)
points(spp.scr[125,1], spp.scr[125,3], col = "black", pch = 16, cex = 1.1)
points(spp.scr[291,1], spp.scr[291,3], col = "black", pch = 17, cex = 1.1)
points(spp.scr[743,1], spp.scr[743,3], col = "black", pch = 15, cex = 1.1)
legend("topleft",
       legend=c("contig 8420", "contig 19728", "contig 54288"),
       pch=c(16, 17, 15),
       col=c("black", "black", "black"))
legend("bottomleft",
       legend = expression("B"),
       bty = "n",
       cex = 1.3)

dev.off()

#### Partial RDAs ####
# Partial RDA to control for confouding effects of geography (lat and long)
adults.prda1 <- rda(adults.nonas.232.1137 ~ depth + b_temp + b_salin + Condition(lat + long), envi.ordered.matrix)
adults.prda1
RsquareAdj(adults.prda1)
anova(adults.prda1) # p-value 0.001
plot(adults.prda1, choices = c(1,2), scaling=3)

# Partial RDA to control for confouding effects of climate (depth, b_temp and b_salin)
adults.prda2 <- rda(adults.nonas.232.1137 ~ lat + long + Condition(depth + b_temp + b_salin), envi.ordered.matrix)
adults.prda2
RsquareAdj(adults.prda2)
anova (adults.prda2) # p-value 0.001
plot(adults.prda2, choices = c(1,2), scaling=3)

# Partitioning the variance componenets
head(summary(adults.rda)) # total explainable variance = inertia of constrained matrix = 31.2
head(summary(adults.prda1)) # variance explained by purely climate = inertia value for constrained matrix = 18.67
head(summary(adults.prda2)) # variance explained by purely geography = inertia value for constrained matrix = 12.62
# Joint effect of climate and geography = proportion of variance in which climate and geography cannot be separated due to collinearity = 31.2-18.67-12.62 = -0.09

# Using adjusted R2 to partition variance
# geography, climate and joint geography & climate = 0.617%
# geography = 0.273%
# climate = 0.372%
# joint geography & climate = 0.617-0.273-0.372 = -0.028%

# Variance inflation factors (VIF) in the full and partial RDAs
vif.cca(adults.rda)
vif.cca(adults.prda1)
vif.cca(adults.prda2)

# Forward selection using vegan's ordistep()
step.forward <- ordistep(rda(adults.nonas.232.1137 ~ 1, data=envi.ordered.matrix), scope=formula(adults.rda), direction="forward", pstep=1000)

# Variance partitioning with 2 sets of explanatory variables
adults.rda.part <- varpart(adults.rda, envi.ordered.matrix)

####################################################################
## Performing full RDA using Hellinger transformed genetic matrix ##
####################################################################
adults.hell.rda <- rda(adults.hell.232 ~ lat + long + depth + b_temp + b_salin, envi)
adults.hell.rda
summary(adults.hell.rda) # inertia is pretty much the same, but the decimal has been moved over 3 places
RsquareAdj (adults.hell.rda)
plot(adults.hell.rda)
plot(adults.hell.rda, xlim=c(-0.015,0.015), ylim=c(-0.015,0.015))
anova(adults.hell.rda) # 0.001

# Partial RDA to control for confouding effects of geography (lat and long)
adults.hell.prda1 <- rda(adults.hell.232 ~ depth + b_temp + b_salin + Condition(lat + long), envi)
adults.hell.prda1
anova(adults.hell.prda1) # p-value 0.001
plot(adults.hell.prda1)

# Partial RDA to control for confouding effects of climate (depth, b_temp and b_salin)
adults.hell.prda2 <- rda(adults.hell.232 ~ lat + long + Condition(depth + b_temp + b_salin), envi)
adults.hell.prda2
anova (adults.hell.prda2) # p-value 0.005 
plot(adults.hell.prda2)

# Partitioning the variance componenets
head(summary(adults.hell.rda)) # total explainable variance = inertia of constrained matrix = 0.003216
head(summary(adults.hell.prda1)) # variance explained by purely climate = inertia value for constrained matrix = 0.001934
head(summary(adults.hell.prda2)) # variance explained by purely geography = inertia value for constrained matrix = 0.001200
# Joint effect of climate and geography = proportion of variance in which climate and geography cannot be separated due to collinearity = 0.003216-0.001934-0.001200 = 0.000082

# Plotting each BayEnv loci separately on RDA1 vs. RDA2
setwd("/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Redundancy analysis")
pdf("can_loci_separately.pdf")

spp.scr <- scores(adults.rda, display = "species", scaling = 3, choices = c(1,2,3,4,5))

plot(adults.rda, choices = c(1,2), scaling=3, main="Locus 35")
points(spp.scr[35,1], spp.scr[35,2], col="darkblue", pch = 16)

plot(adults.rda, choices = c(1,2), scaling=3, main="Locus 125")
points(spp.scr[125,1], spp.scr[125,2], col="darkblue", pch = 16)

plot(adults.rda, choices = c(1,2), scaling=3, main="Locus 214")
points(spp.scr[214,1], spp.scr[214,2], col="darkblue", pch = 16)

plot(adults.rda, choices = c(1,2), scaling=3, main="Locus 396")
points(spp.scr[396,1], spp.scr[396,2], col="darkblue", pch = 16)

plot(adults.rda, choices = c(1,2), scaling=3, main="Locus 442")
points(spp.scr[442,1], spp.scr[442,2], col="darkblue", pch = 16)

plot(adults.rda, choices = c(1,2), scaling=3, main="Locus 499")
points(spp.scr[499,1], spp.scr[499,2], col="darkblue", pch = 16)

plot(adults.rda, choices = c(1,2), scaling=3, main="Locus 524")
points(spp.scr[524,1], spp.scr[524,2], col="darkblue", pch = 16)

plot(adults.rda, choices = c(1,2), scaling=3, main="Locus 542")
points(spp.scr[542,1], spp.scr[542,2], col="darkblue", pch = 16)

plot(adults.rda, choices = c(1,2), scaling=3, main="Locus 609")
points(spp.scr[609,1], spp.scr[609,2], col="darkblue", pch = 16)

plot(adults.rda, choices = c(1,2), scaling=3, main="Locus 615")
points(spp.scr[615,1], spp.scr[615,2], col="darkblue", pch = 16)

plot(adults.rda, choices = c(1,2), scaling=3, main="Locus 626")
points(spp.scr[626,1], spp.scr[626,2], col="darkblue", pch = 16)

plot(adults.rda, choices = c(1,2), scaling=3, main="Locus 703")
points(spp.scr[703,1], spp.scr[703,2], col="darkblue", pch = 16)

plot(adults.rda, choices = c(1,2), scaling=3, main="Locus 743")
points(spp.scr[743,1], spp.scr[743,2], col="darkblue", pch = 16)

plot(adults.rda, choices = c(1,2), scaling=3, main="Locus 808")
points(spp.scr[808,1], spp.scr[808,2], col="darkblue", pch = 16)

plot(adults.rda, choices = c(1,2), scaling=3, main="Locus 825")
points(spp.scr[825,1], spp.scr[825,2], col="darkblue", pch = 16)

plot(adults.rda, choices = c(1,2), scaling=3, main="Locus 826")
points(spp.scr[826,1], spp.scr[826,2], col="darkblue", pch = 16)

plot(adults.rda, choices = c(1,2), scaling=3, main="Locus 829")
points(spp.scr[829,1], spp.scr[829,2], col="darkblue", pch = 16)

plot(adults.rda, choices = c(1,2), scaling=3, main="Locus 916")
points(spp.scr[916,1], spp.scr[916,2], col="darkblue", pch = 16)

plot(adults.rda, choices = c(1,2), scaling=3, main="Locus 919")
points(spp.scr[919,1], spp.scr[919,2], col="darkblue", pch = 16)

plot(adults.rda, choices = c(1,2), scaling=3, main="Locus 923")
points(spp.scr[923,1], spp.scr[923,2], col="darkblue", pch = 16)

plot(adults.rda, choices = c(1,2), scaling=3, main="Locus 990")
points(spp.scr[990,1], spp.scr[990,2], col="darkblue", pch = 16)

plot(adults.rda, choices = c(1,2), scaling=3, main="Locus 1050")
points(spp.scr[1050,1], spp.scr[1050,2], col="darkblue", pch = 16)

dev.off()

# Attempting to write the plotting above as a for loop - not perfected yet
pdf("can_loci_separately_test.pdf")
loci <- c(35, 125, 214, 396, 442, 499, 524, 542, 609, 615, 626, 703, 743, 808, 825, 826, 829, 916, 919, 923, 990, 1050)
can_loci <- spp.scr[c(35, 125, 214, 396, 442, 499, 524, 542, 609, 615, 626, 703, 743, 808, 825, 826, 829, 916, 919, 923, 990, 1050), c(1,2)]

  for (locus in loci){
  plot(adults.rda, choices = c(1,2), scaling=3, main=paste("Locus", loci))
    for (snp in can_loci){
      points(can_loci, col = "darkblue", pch = 16)
      # arrows(0,0, can_loci[snp], can_loci[snp], length = 0, lty=1, col = "darkblue", cex = 0.8)
}}

dev.off()

# What if I did RDA using only 22 loci that BayEnv indicated as important? Would b_temp come through as more important environmental variable? And 16 loci fall on that axis?


