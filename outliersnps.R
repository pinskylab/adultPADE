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

# Write to my computer
write.table(nvs_locs, "232envirowithdist.txt", sep = "\t", col.names = TRUE, row.names = FALSE)

#### Plotting covariance between five environmental variables ####
envirowdist <- read.table("232envirowithdist.txt", header = TRUE)
envs <- envirowdist[, c("dist", "depth", "b_temp", "b_salin")]
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

png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation/results_avg232covarmatrix_hwe_4stanenviro/envs4_covariance.png", width=7.5, height=7.5, res=300, units="in")
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
# Creating the ENVIROFILE using the average depth, bottom temperature, bottom salinity and distance from a southern point
# Aggregate the individuals into populations and standardize
by <- list(nvs_locs$bayenv_pop)
pop.avg <- aggregate(nvs_locs[,c(7, 8, 9, 11)], by = by, FUN = mean)
stan.pop.avg <- scale(pop.avg[,2:5])
rownames(stan.pop.avg) <- c("Pop1", "Pop2", "Pop3", "Pop4", "Pop5")
t.stan.pop.avg <- t(stan.pop.avg) # the order is Pop1, Pop2, Pop3, Pop4, Pop5; rows are: depth, b_temp, b_salin, dist

# Exporting the envirofile to my computer
write.table (t.stan.pop.avg, "232stan4envirofile.txt", sep ="\t", col.names = FALSE, row.names = FALSE)

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

#### This is for the correctly standardized environmental variables ####
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

###############################################################################################################################################
#### BayEnv analysis using 4 standardized environmental variables (depth, bottom temp, bottom salinity and distance from a southern point) ####
###############################################################################################################################################
setwd("/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation/results_avg232covarmatrix_hwe_4stanenviro")

library(abind)

run1 <- read.table("bf1.txt")
run1 <- run1[-c(1:2), ]# first two snps were analyzed twice, I think because I accidently killed the first run. Deleting the first two rows.
run1$V1 <- NULL
rownames(run1) <- 1:length(run1[,1]) # fix row names so that it's less confusing later
run2 <- read.table("bf2.txt")
run2$V1 <- NULL
run3 <- read.table("bf3.txt")
run3$V1 <- NULL
run4 <- read.table("bf4.txt")
run4$V1 <- NULL
run5 <- read.table("bf5.txt")
run5$V1 <- NULL
run6 <- read.table("bf6.txt")
run6$V1 <- NULL
run7 <- read.table("bf7.txt")
run7$V1 <- NULL
run8 <- read.table("bf8.txt")
run8$V1 <- NULL
run9 <- read.table("bf9.txt")
run9$V1 <- NULL
run10 <- read.table("bf10.txt")
run10$V1 <- NULL

full_array <- abind(run1, run2, run3, run4, run5, run6, run7, run8, run9, run10, along=3)

full_array_median <- apply(full_array, c(1,2), median)
colnames(full_array_median) <- c("depth", "b_temp", "b_salin", "dist")
full_array_median

which(full_array_median[,1] > 3) #125 214 396 615 743 825 
which(full_array_median[,2] > 3) #125 214 396 442 499 524 609 615 825 919 
which(full_array_median[,3] > 3) #35 542 609 615 743 990  
which(full_array_median[,4] > 3) #125 214 396 442 499 609 615 825 919 

# Get the actual BF
full_array_median[which(full_array_median[,1] > 3),1]
full_array_median[which(full_array_median[,2] > 3),2]
full_array_median[which(full_array_median[,3] > 3),3]
full_array_median[which(full_array_median[,4] > 3),4]

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
png(file = "/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation/results_avg232covarmatrix_hwe_4stanenviro/histogramplots_same.png", width=8, height=8, res=300, units="in")
# png(file = "/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation/results_avg232covarmatrix_hwe_5enviro/histogramplots_same.png", width=8, height=11, res=300, units="in")

par(
  mfrow = c(2, 2), 
  mar=c(4.5, 3.5, 1, 1), # panel magin size in "line number" units
  mgp=c(2, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=10, # point size, which is the font size
  bg=NA
)

hist(full_array_median[,4], xlab="", ylab = "", main = "", breaks=11, xlim=c(0,11), xaxt="n")
axis(side=1, at=seq(0,11, 1), labels=seq(0,11,1), line=1.3)
rug(jitter(full_array_median[,4]), ticksize = -0.1, line=-0.2)
text(5.5,800, paste("max BF = 8.85"))
text(5.5,950, paste("Distance"), cex = 1.3)
mtext("Bayes Factor", side = 1, line = 3.6)
mtext("Frequency", side = 2, line = 2.7)
hist(full_array_median[,1], xlab="", ylab = "", main = "", breaks=11, xlim=c(0,11), xaxt="n")
axis(side=1, at=seq(0,11, 1), labels=seq(0,11,1), line=1.3)
rug(jitter(full_array_median[,1]), ticksize = -0.1, line=-0.2)
text(5.5,650, paste("max BF = 6.85"))
text(5.5,800, paste("Depth"), cex = 1.3)
mtext("Bayes Factor", side = 1, line = 3.6)
mtext("Frequency", side = 2, line = 2.7)
hist(full_array_median[,2], xlab="", ylab = "", main = "", breaks=11, xlim=c(0,11), xaxt="n")
axis(side=1, at=seq(0,11, 1), labels=seq(0,11,1), line=1.3)
rug(jitter(full_array_median[,2]), ticksize = -0.1, line=-0.2)
text(5.5,800, paste("max BF = 10.59"))
text(5.5,950, paste("Bottom temperature"), cex = 1.3)
mtext("Bayes Factor", side = 1, line = 3.6)
mtext("Frequency", side = 2, line = 2.7)
hist(full_array_median[,3], xlab="", ylab = "", main = "", breaks=12, xlim=c(0,11), xaxt="n")
axis(side=1, at=seq(0,11, 1), labels=seq(0,11,1), line=1.3)
rug(jitter(full_array_median[,3]), ticksize = -0.1, line=-0.2)
text(5.5,800, paste("max BF = 9.86"))
text(5.5,950, paste("Bottom salinity"), cex = 1.3)
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
# write.table(adults.nonas.232.1137, "~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/GAM/232fish1137alleles_notcentered.txt", sep = "\t", col.names = TRUE, row.names = TRUE)
# adults.hell.232 <- adults.hell[ ! rownames(adults.hell) %in% exclu_names, ] # This is genotype data for 232 fish where NAs have been replaced by mean allele frequencies and then a Hellinger transformation has been applied
# rownames(adults.hell.232)

# Bringing in environmental data
setwd("/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation")

# This does not contain distance from a southern point
# envi <- read.csv("allpops_combo.csv", header=TRUE)
# envi.ordered <- envi[with(envi, order(PinskyID)),]
# envi.ordered[,2]
# as.character(envi.ordered[,2]) == rownames(adults.nonas.232.1137)
# envi.ordered
# envi.ordered.matrix <- scale(data.matrix(envi.ordered))
# envi.ordered.matrix
# rownames(envi.ordered.matrix) <- envi.ordered[,"PinskyID"]
# envi.ordered.matrix[,1] <- envi.ordered[,1]

# This contains distance from a southern point
envi <- read.table("232envirowithdist.txt", header = TRUE)
envi.ordered <- envi[with(envi, order(PinskyID)),]
as.character(envi.ordered[,2]) == rownames(adults.nonas.232.1137)
envi.ordered.matrix <- scale(data.matrix(envi.ordered))
envi.ordered.matrix
rownames(envi.ordered.matrix) <- envi.ordered[,"PinskyID"]
envi.ordered.matrix[,1] <- envi.ordered[,1]

# Let's plot environmental variables against each other
pairs(envi.ordered.matrix[,c(11,7,8,9)], main="Bivariate plots")

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

adults.rda <- rda(adults.nonas.232.1137 ~ dist + depth + b_temp + b_salin, envi.ordered.matrix, scale = FALSE)
adults.rda
summary(adults.rda)
RsquareAdj (adults.rda)
R2a.adults.rda <- RsquareAdj(adults.rda)$adj.r.squared
coef(adults.rda)
plot(adults.rda, scaling = 3)
plot(adults.rda, xlim = c(-0.1,0.1), ylim = c(-0.1,0.1))
anova(adults.rda)
anova(adults.rda, by="axis", step=1000)

# To get values proportional to eigenvalues
ev <- sfs_pca$sdev^2
ev.varprop <- ev/sum(ev)
ev.varprop[1:25]

plot(adults.rda, choices = c(1,2), scaling=2)
plot(adults.rda, choices = c(1,3), xlim = c(-0.85,0.8), ylim = c(-0.8,0.7), scaling=3)
plot(adults.rda, choices = c(1,2), xlim = c(-0.85,0.8), ylim = c(-0.8,0.7), scaling=3)
plot(adults.rda, choices = c(2,3), xlim = c(-0.85,0.8), ylim = c(-0.8,0.7), scaling=3)
plot(adults.rda, choices = c(1,4), xlim = c(-0.85,0.8), ylim = c(-0.8,0.7), scaling=3)
plot(adults.rda, choices = c(2,4), xlim = c(-0.85,0.8), ylim = c(-0.8,0.7), scaling=3)
plot(adults.rda, choices = c(3,4), xlim = c(-0.85,0.8), ylim = c(-0.8,0.7), scaling=3)

spp.scr <- scores(adults.rda, display = "species", scaling = 2, choices = c(1,2,3,4))
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

spp.scr <- scores(adults.rda, display = "species", scaling = 3, choices = c(1,2,3,4))
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
spp.scr <- scores(adults.rda, display = "species", scaling = 0, choices = c(1,2,3,4))
head(spp.scr)
allele_loadings <- spp.scr[c(35, 125, 214, 396, 442, 499, 524, 542, 609, 615, 626, 703, 743, 808, 825, 826, 829, 916, 919, 923, 990, 1050),]
allele_loadings # unscaled species scores

# Forester et al. (2016) used loci with scores +/- 3 SD from mean score for that axis to ID outlier loci. Only used first 3 axes.
mean.rda <- colMeans(spp.scr)
sd.rda1 <- sd(spp.scr[,"RDA1"])
sd.rda2 <- sd(spp.scr[,"RDA2"])
sd.rda3 <- sd(spp.scr[,"RDA3"])
# sd.rda4 <- sd(spp.scr[,"RDA4"])

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
              # (which(spp.scr[,"RDA4"] > rda4.hi)),
              # (which(spp.scr[,"RDA4"] < rda4.lo))
#               (which(spp.scr[,"RDA5"] > rda5.hi)),
#               (which(spp.scr[,"RDA5"] < rda5.lo))
              )

rda.cans <- unique(sort(rda.cans)) # 4 enviro variables, including distance from a southern point: 43   47  125  291  380  402  432  466  483  499  577  579  635  695  704  743  826  842  937  947  965  997 1090
# when 5 enviro variables, including lat & lon: 125  128  291  380  402  432  466  577  579  639  691  694  695  704  720  743  826  842  937  947  997 1035 1061 1090

# Linear regression between allele frequencies and environmental variables for the RDA outliers
dist43 <- lm(adults.nonas.232.1137[,43] ~ envi.ordered.matrix[,"dist"])
depth43 <- lm(adults.nonas.232.1137[,43] ~ envi.ordered.matrix[,"depth"])
btemp43 <- lm(adults.nonas.232.1137[,43] ~ envi.ordered.matrix[,"b_temp"])
bsalin43 <- lm(adults.nonas.232.1137[,43] ~ envi.ordered.matrix[,"b_salin"])

dist47 <- lm(adults.nonas.232.1137[,47] ~ envi.ordered.matrix[,"dist"])
depth47 <- lm(adults.nonas.232.1137[,47] ~ envi.ordered.matrix[,"depth"])
btemp47 <- lm(adults.nonas.232.1137[,47] ~ envi.ordered.matrix[,"b_temp"])
bsalin47 <- lm(adults.nonas.232.1137[,47] ~ envi.ordered.matrix[,"b_salin"])

dist125 <- lm(adults.nonas.232.1137[,125] ~ envi.ordered.matrix[,"dist"])
depth125 <- lm(adults.nonas.232.1137[,125] ~ envi.ordered.matrix[,"depth"])
btemp125 <- lm(adults.nonas.232.1137[,125] ~ envi.ordered.matrix[,"b_temp"])
bsalin125 <- lm(adults.nonas.232.1137[,125] ~ envi.ordered.matrix[,"b_salin"]) # 1.490386e-04

dist291 <- lm(adults.nonas.232.1137[,291] ~ envi.ordered.matrix[,"dist"])
depth291 <- lm(adults.nonas.232.1137[,291] ~ envi.ordered.matrix[,"depth"]) # 2.764823e-05
btemp291 <- lm(adults.nonas.232.1137[,291] ~ envi.ordered.matrix[,"b_temp"]) # 2.620017e-05
bsalin291 <- lm(adults.nonas.232.1137[,291] ~ envi.ordered.matrix[,"b_salin"])

dist380 <- lm(adults.nonas.232.1137[,380] ~ envi.ordered.matrix[,"dist"])
depth380 <- lm(adults.nonas.232.1137[,380] ~ envi.ordered.matrix[,"depth"]) 
btemp380 <- lm(adults.nonas.232.1137[,380] ~ envi.ordered.matrix[,"b_temp"]) 
bsalin380 <- lm(adults.nonas.232.1137[,380] ~ envi.ordered.matrix[,"b_salin"])

dist402 <- lm(adults.nonas.232.1137[,402] ~ envi.ordered.matrix[,"dist"])
depth402 <- lm(adults.nonas.232.1137[,402] ~ envi.ordered.matrix[,"depth"]) 
btemp402 <- lm(adults.nonas.232.1137[,402] ~ envi.ordered.matrix[,"b_temp"]) 
bsalin402 <- lm(adults.nonas.232.1137[,402] ~ envi.ordered.matrix[,"b_salin"])

dist432 <- lm(adults.nonas.232.1137[,432] ~ envi.ordered.matrix[,"dist"])
depth432 <- lm(adults.nonas.232.1137[,432] ~ envi.ordered.matrix[,"depth"]) 
btemp432 <- lm(adults.nonas.232.1137[,432] ~ envi.ordered.matrix[,"b_temp"]) 
bsalin432 <- lm(adults.nonas.232.1137[,432] ~ envi.ordered.matrix[,"b_salin"])

dist466 <- lm(adults.nonas.232.1137[,466] ~ envi.ordered.matrix[,"dist"])
depth466 <- lm(adults.nonas.232.1137[,466] ~ envi.ordered.matrix[,"depth"]) 
btemp466 <- lm(adults.nonas.232.1137[,466] ~ envi.ordered.matrix[,"b_temp"]) 
bsalin466 <- lm(adults.nonas.232.1137[,466] ~ envi.ordered.matrix[,"b_salin"])

dist483 <- lm(adults.nonas.232.1137[,483] ~ envi.ordered.matrix[,"dist"])
depth483 <- lm(adults.nonas.232.1137[,483] ~ envi.ordered.matrix[,"depth"]) 
btemp483 <- lm(adults.nonas.232.1137[,483] ~ envi.ordered.matrix[,"b_temp"]) 
bsalin483 <- lm(adults.nonas.232.1137[,483] ~ envi.ordered.matrix[,"b_salin"])

dist499 <- lm(adults.nonas.232.1137[,499] ~ envi.ordered.matrix[,"dist"]) # 5.981221e-04
depth499 <- lm(adults.nonas.232.1137[,499] ~ envi.ordered.matrix[,"depth"]) 
btemp499 <- lm(adults.nonas.232.1137[,499] ~ envi.ordered.matrix[,"b_temp"]) 
bsalin499 <- lm(adults.nonas.232.1137[,499] ~ envi.ordered.matrix[,"b_salin"])

dist577 <- lm(adults.nonas.232.1137[,577] ~ envi.ordered.matrix[,"dist"])
depth577 <- lm(adults.nonas.232.1137[,577] ~ envi.ordered.matrix[,"depth"]) 
btemp577 <- lm(adults.nonas.232.1137[,577] ~ envi.ordered.matrix[,"b_temp"]) 
bsalin577 <- lm(adults.nonas.232.1137[,577] ~ envi.ordered.matrix[,"b_salin"])

dist579 <- lm(adults.nonas.232.1137[,579] ~ envi.ordered.matrix[,"dist"])
depth579 <- lm(adults.nonas.232.1137[,579] ~ envi.ordered.matrix[,"depth"]) 
btemp579 <- lm(adults.nonas.232.1137[,579] ~ envi.ordered.matrix[,"b_temp"]) 
bsalin579 <- lm(adults.nonas.232.1137[,579] ~ envi.ordered.matrix[,"b_salin"])

dist635 <- lm(adults.nonas.232.1137[,635] ~ envi.ordered.matrix[,"dist"])
depth635 <- lm(adults.nonas.232.1137[,635] ~ envi.ordered.matrix[,"depth"]) 
btemp635 <- lm(adults.nonas.232.1137[,635] ~ envi.ordered.matrix[,"b_temp"]) 
bsalin635 <- lm(adults.nonas.232.1137[,635] ~ envi.ordered.matrix[,"b_salin"])

dist695 <- lm(adults.nonas.232.1137[,695] ~ envi.ordered.matrix[,"dist"])
depth695 <- lm(adults.nonas.232.1137[,695] ~ envi.ordered.matrix[,"depth"]) 
btemp695 <- lm(adults.nonas.232.1137[,695] ~ envi.ordered.matrix[,"b_temp"]) 
bsalin695 <- lm(adults.nonas.232.1137[,695] ~ envi.ordered.matrix[,"b_salin"])

dist704 <- lm(adults.nonas.232.1137[,704] ~ envi.ordered.matrix[,"dist"])
depth704 <- lm(adults.nonas.232.1137[,704] ~ envi.ordered.matrix[,"depth"]) 
btemp704 <- lm(adults.nonas.232.1137[,704] ~ envi.ordered.matrix[,"b_temp"]) 
bsalin704 <- lm(adults.nonas.232.1137[,704] ~ envi.ordered.matrix[,"b_salin"])

dist743 <- lm(adults.nonas.232.1137[,743] ~ envi.ordered.matrix[,"dist"])
depth743 <- lm(adults.nonas.232.1137[,743] ~ envi.ordered.matrix[,"depth"]) 
btemp743 <- lm(adults.nonas.232.1137[,743] ~ envi.ordered.matrix[,"b_temp"]) 
bsalin743 <- lm(adults.nonas.232.1137[,743] ~ envi.ordered.matrix[,"b_salin"]) #6.089858e-04

dist826 <- lm(adults.nonas.232.1137[,826] ~ envi.ordered.matrix[,"dist"])
depth826 <- lm(adults.nonas.232.1137[,826] ~ envi.ordered.matrix[,"depth"]) 
btemp826 <- lm(adults.nonas.232.1137[,826] ~ envi.ordered.matrix[,"b_temp"]) 
bsalin826 <- lm(adults.nonas.232.1137[,826] ~ envi.ordered.matrix[,"b_salin"])

dist842 <- lm(adults.nonas.232.1137[,842] ~ envi.ordered.matrix[,"dist"])
depth842 <- lm(adults.nonas.232.1137[,842] ~ envi.ordered.matrix[,"depth"]) 
btemp842 <- lm(adults.nonas.232.1137[,842] ~ envi.ordered.matrix[,"b_temp"]) 
bsalin842 <- lm(adults.nonas.232.1137[,842] ~ envi.ordered.matrix[,"b_salin"])

dist937 <- lm(adults.nonas.232.1137[,937] ~ envi.ordered.matrix[,"dist"])
depth937 <- lm(adults.nonas.232.1137[,937] ~ envi.ordered.matrix[,"depth"]) 
btemp937 <- lm(adults.nonas.232.1137[,937] ~ envi.ordered.matrix[,"b_temp"]) 
bsalin937 <- lm(adults.nonas.232.1137[,937] ~ envi.ordered.matrix[,"b_salin"])

dist947 <- lm(adults.nonas.232.1137[,947] ~ envi.ordered.matrix[,"dist"])
depth947 <- lm(adults.nonas.232.1137[,947] ~ envi.ordered.matrix[,"depth"]) 
btemp947 <- lm(adults.nonas.232.1137[,947] ~ envi.ordered.matrix[,"b_temp"]) 
bsalin947 <- lm(adults.nonas.232.1137[,947] ~ envi.ordered.matrix[,"b_salin"])

dist965 <- lm(adults.nonas.232.1137[,965] ~ envi.ordered.matrix[,"dist"])
depth965 <- lm(adults.nonas.232.1137[,965] ~ envi.ordered.matrix[,"depth"]) 
btemp965 <- lm(adults.nonas.232.1137[,965] ~ envi.ordered.matrix[,"b_temp"]) 
bsalin965 <- lm(adults.nonas.232.1137[,965] ~ envi.ordered.matrix[,"b_salin"])

dist997 <- lm(adults.nonas.232.1137[,997] ~ envi.ordered.matrix[,"dist"])
depth997 <- lm(adults.nonas.232.1137[,997] ~ envi.ordered.matrix[,"depth"]) 
btemp997 <- lm(adults.nonas.232.1137[,997] ~ envi.ordered.matrix[,"b_temp"]) 
bsalin997 <- lm(adults.nonas.232.1137[,997] ~ envi.ordered.matrix[,"b_salin"])

dist1090 <- lm(adults.nonas.232.1137[,1090] ~ envi.ordered.matrix[,"dist"])
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

# cor.test(adults.nonas.232.1137[,43], envi.ordered.matrix[,'depth']) # results in same pvalues and using lmp()

pvals <- c(lmp(dist43),lmp(depth43),lmp(btemp43),lmp(bsalin43),lmp(dist47),lmp(depth47),lmp(btemp47),lmp(bsalin47), lmp(dist125), lmp(depth125),lmp(btemp125),lmp(bsalin125), 
           lmp(dist291),lmp(depth291),lmp(btemp291),lmp(bsalin291), lmp(dist380),lmp(depth380),lmp(btemp380),lmp(bsalin380), lmp(dist402),lmp(depth402),lmp(btemp402),lmp(bsalin402), 
           lmp(dist432),lmp(depth432),lmp(btemp432),lmp(bsalin432), lmp(dist466),lmp(depth466),lmp(btemp466),lmp(bsalin466), lmp(dist483),lmp(depth483),lmp(btemp483),lmp(bsalin483), 
           lmp(dist499),lmp(depth499),lmp(btemp499),lmp(bsalin499), lmp(dist577), lmp(depth577),lmp(btemp577),lmp(bsalin577), lmp(dist579),lmp(depth579),lmp(btemp579), 
           lmp(bsalin579), lmp(dist635),lmp(depth635),lmp(btemp635),lmp(bsalin635), lmp(dist695),lmp(depth695),lmp(btemp695),lmp(bsalin695), 
           lmp(dist704),lmp(depth704),lmp(btemp704),lmp(bsalin704), lmp(dist743),lmp(depth743),lmp(btemp743),lmp(bsalin743), lmp(dist826),lmp(depth826),lmp(btemp826),lmp(bsalin826), lmp(dist842),
           lmp(depth842),lmp(btemp842),lmp(bsalin842), lmp(dist937),lmp(depth937),lmp(btemp937), lmp(bsalin937), lmp(dist947),lmp(depth947),lmp(btemp947),lmp(bsalin947), lmp(dist965), 
           lmp(depth965), lmp(btemp965), lmp(bsalin965), lmp(dist997), lmp(depth997), lmp(btemp997), lmp(bsalin997), lmp(dist1090), lmp(depth1090), lmp(btemp1090), lmp(bsalin1090))

which(pvals < 0.001)

# Same thing as above, but in for loop form. Good to double check #
# For loop through lm for each candidate locus & 4 enviromental variables (23*4 = 92)
# Create a vector of all lm names
names_lms <- vector()
for(i in rda.cans){
  names_lms <- append(names_lms, paste0('dist',i))
  names_lms <- append(names_lms, paste0('depth',i))
  names_lms <- append(names_lms, paste0('btemp',i))
  names_lms <- append(names_lms, paste0('bsalin',i))
}

# For loop to regress all RDA candidates against 4 enviromental variables
for(i in rda.cans){
  assign(paste0('dist',i), lm(adults.nonas.232.1137[,i] ~ envi.ordered.matrix[,"dist"]))
  assign(paste0('depth',i), lm(adults.nonas.232.1137[,i] ~ envi.ordered.matrix[,"depth"]))
  assign(paste0('btemp',i), lm(adults.nonas.232.1137[,i] ~ envi.ordered.matrix[,"b_temp"]))
  assign(paste0('bsalin',i), lm(adults.nonas.232.1137[,i] ~ envi.ordered.matrix[,"b_salin"]))
}

pvalues <- c(lmp(dist43),  lmp(depth43),   lmp(btemp43),  lmp(bsalin43), lmp(dist47),    lmp(depth47),   lmp(btemp47),  lmp(bsalin47),   lmp(dist125),  lmp(depth125),   lmp(btemp125),  lmp(bsalin125),  lmp(dist291)  ,  lmp(depth291), 
           lmp(btemp291),   lmp(bsalin291),  lmp(dist380),  lmp(depth380), lmp(btemp380),   lmp(bsalin380),  lmp(dist402),  lmp(depth402),   lmp(btemp402),  lmp(bsalin402),  lmp(dist432),   lmp(depth432),   lmp(btemp432),   lmp(bsalin432), 
           lmp(dist466),  lmp(depth466),   lmp(btemp466),  lmp(bsalin466),  lmp(dist483),    lmp(depth483),   lmp(btemp483),  lmp(bsalin483),  lmp(dist499),  lmp(depth499),   lmp(btemp499),  lmp(bsalin499),  lmp(dist577)  ,  lmp(depth577), 
           lmp(btemp577),   lmp(bsalin577),  lmp(dist579),   lmp(depth579), lmp(btemp579),   lmp(bsalin579),  lmp(dist635),  lmp(depth635),   lmp(btemp635),  lmp(bsalin635),  lmp(dist695),   lmp(depth695),   lmp(btemp695),   lmp(bsalin695), 
           lmp(dist704),  lmp(depth704),   lmp(btemp704),  lmp(bsalin704),  lmp(dist743),    lmp(depth743),   lmp(btemp743),  lmp(bsalin743),  lmp(dist826),  lmp(depth826),   lmp(btemp826),  lmp(bsalin826),  lmp(dist842) ,  lmp(depth842), 
           lmp(btemp842),   lmp(bsalin842),  lmp(dist937),   lmp(depth937), lmp(btemp937),   lmp(bsalin937),  lmp(dist947),  lmp(depth947),   lmp(btemp947),  lmp(bsalin947),  lmp(dist965),   lmp(depth965),   lmp(btemp965),   lmp(bsalin965), 
           lmp(dist997),  lmp(depth997),   lmp(btemp997),  lmp(bsalin997),  lmp(dist1090),   lmp(depth1090),  lmp(btemp1090),  lmp(bsalin1090))

rda.can.loci <- data.frame(names_lms,pvalues)
rda.can.loci2 <- rda.can.loci[which(rda.can.loci[,2] < 0.001),] # dim is 5 x 2


# Plotting the RDA plot with the three loci indicated by redundancy analysis to have strong locus-environmental associations
png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Redundancy analysis/RDA plotwith4enviros.png", width=11, height=6, res=300, units="in")

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
points(spp.scr[499,1], spp.scr[499,2], col = "black", pch = 18, cex = 1.1)
points(spp.scr[743,1], spp.scr[743,2], col = "black", pch = 15, cex = 1.1)
legend("topleft",
	legend=c("contig 8420", "contig 19728", "contig 35399", "contig 54288"),
	pch=c(16, 17, 18, 15),
	col=c("black", "black", "black", "black"))
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
points(spp.scr[499,1], spp.scr[499,3], col = "black", pch = 18, cex = 1.1)
points(spp.scr[743,1], spp.scr[743,3], col = "black", pch = 15, cex = 1.1)
legend("topleft",
       legend=c("contig 8420", "contig 19728", "contig 35399", "contig 54288"),
       pch=c(16, 17, 18, 15),
       col=c("black", "black", "black", "black"))
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

#####################################################################################################
# Plots of allele frequency vs environment are really hard to visualize --> aggregate by environment
# locus 125 and depth
d<- cbind(adults.nonas.232.1137[,125], envi.ordered.matrix[,"depth"])
d2 <- aggregate(d[, 1], list(d[,2]), mean)
length(unique(d[,2]))
plot(d2$x ~ d2$Group.1, xlab = 'Standardized depth', ylab = 'Aggregated allele freq', ylim = c(-0.8,0.2))
d.lm <- lm(d2$x ~ d2$Group.1)
abline(d.lm, col = 'red')

plot(adults.nonas.232.1137[,125] ~ envi.ordered.matrix[,"depth"]) # in comparision
abline(lm(adults.nonas.232.1137[,125] ~ envi.ordered.matrix[,"depth"]), col = 'red')


# locus 125 and salinity
d<- cbind(adults.nonas.232.1137[,125], envi.ordered.matrix[,"b_salin"])
d2 <- aggregate(d[, 1], list(d[,2]), mean)
length(unique(d[,2]))
plot(d2$x ~ d2$Group.1, xlab = 'Standardized b_salin', ylab = 'Aggregated allele freq', ylim = c(-0.8,0.2))
d.lm <- lm(d2$x ~ d2$Group.1)
abline(d.lm, col = 'red')

plot(adults.nonas.232.1137[,125] ~ envi.ordered.matrix[,"b_salin"]) # in comparision
abline(lm(adults.nonas.232.1137[,125] ~ envi.ordered.matrix[,"b_salin"]), col = 'red')

# locus 291 and depth
d<- cbind(adults.nonas.232.1137[,291], envi.ordered.matrix[,"depth"])
d2 <- aggregate(d[, 1], list(d[,2]), mean)
length(unique(d[,2]))
plot(d2$x ~ d2$Group.1, xlab = 'Standardized depth', ylab = 'Aggregated allele freq')
d.lm <- lm(d2$x ~ d2$Group.1)
abline(d.lm, col = 'red')

plot(adults.nonas.232.1137[,291] ~ envi.ordered.matrix[,"depth"]) # in comparision
abline(lm(adults.nonas.232.1137[,291] ~ envi.ordered.matrix[,"depth"]), col = 'red')

# locus 291 and temperature
d<- cbind(adults.nonas.232.1137[,291], envi.ordered.matrix[,"b_temp"])
d2 <- aggregate(d[, 1], list(d[,2]), mean)
length(unique(d[,2]))
plot(d2$x ~ d2$Group.1, xlab = 'Standardized b_temp', ylab = 'Aggregated allele freq')
d.lm <- lm(d2$x ~ d2$Group.1)
abline(d.lm, col = 'red')

plot(adults.nonas.232.1137[,291] ~ envi.ordered.matrix[,"b_temp"]) # in comparision
abline(lm(adults.nonas.232.1137[,291] ~ envi.ordered.matrix[,"b_temp"]), col = 'red')

# locus 499 and distance
d<- cbind(adults.nonas.232.1137[,499], envi.ordered.matrix[,"dist"])
d2 <- aggregate(d[, 1], list(d[,2]), mean)
length(unique(d[,2]))
plot(d2$x ~ d2$Group.1, xlab = 'Standardized distance', ylab = 'Aggregated allele freq')
d.lm <- lm(d2$x ~ d2$Group.1)
abline(d.lm, col = 'red')

plot(adults.nonas.232.1137[,499] ~ envi.ordered.matrix[,"dist"]) # in comparision
abline(lm(adults.nonas.232.1137[,499] ~ envi.ordered.matrix[,"dist"]), col = 'red')

# locus 499 and depth
d<- cbind(adults.nonas.232.1137[,499], envi.ordered.matrix[,"depth"])
d2 <- aggregate(d[, 1], list(d[,2]), mean)
length(unique(d[,2]))
plot(d2$x ~ d2$Group.1, xlab = 'Standardized depth', ylab = 'Aggregated allele freq')
d.lm <- lm(d2$x ~ d2$Group.1)
abline(d.lm, col = 'red')

plot(adults.nonas.232.1137[,499] ~ envi.ordered.matrix[,"depth"]) # in comparision
abline(lm(adults.nonas.232.1137[,499] ~ envi.ordered.matrix[,"depth"]), col = 'red')

# locus 743 and salinity
d<- cbind(adults.nonas.232.1137[,743], envi.ordered.matrix[,"b_salin"])
d2 <- aggregate(d[, 1], list(d[,2]), mean)
length(unique(d[,2]))
plot(d2$x ~ d2$Group.1, xlab = 'Standardized salinity', ylab = 'Aggregated allele freq')
d.lm <- lm(d2$x ~ d2$Group.1)
abline(d.lm, col = 'red')

plot(adults.nonas.232.1137[,743] ~ envi.ordered.matrix[,"b_salin"]) # in comparision
abline(lm(adults.nonas.232.1137[,743] ~ envi.ordered.matrix[,"b_salin"]), col = 'red')


##################################################################################################
#### Part of RDA is multivarite linear regression. How to bypass this? Rank regression or GAM ####
##################################################################################################
# gam and mcgv packages do not play well
library(gam)

setwd("/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/GAM")
allele.table <- read.table("232fish1137alleles_notcentered.txt", header = TRUE) # allele frequencies
# allele.table <- read.table("232fish1137alleles.txt", header = TRUE)
setwd("/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation")
envi <- read.table("232envirowithdist.txt", header = TRUE)

# Using gam::gam
par(mfrow = c(2,2))
gam1 <- gam(allele.table[,125] ~ s(dist, df = 4) + s(depth, df = 4) + s(b_temp, df = 4) + s(b_salin, df = 4), data = envi)
plot(gam1, se = TRUE)
summary(gam1)

gam2 <- gam(allele.table[,291] ~ s(dist, df = 4) + s(depth, df = 4) + s(b_temp, df = 4) + s(b_salin, df = 4), data = envi)
plot(gam2, se = TRUE)
summary(gam2)

gam3 <- gam(allele.table[,499] ~ s(dist, df = 4) + s(depth, df = 4) + s(b_temp, df = 4) + s(b_salin, df = 4), data = envi)
plot(gam3, se = TRUE)
summary(gam3)

gam4 <- gam(allele.table[,743] ~ s(dist, df = 4) + s(depth, df = 4) + s(b_temp, df = 4) + s(b_salin, df = 4), data = envi)
plot(gam4, se = TRUE)
summary(gam4)

# Automatic smoothness selection
library(mgcv)
mod1 <- mgcv::gam(allele.table[,125] ~ s(dist) + s(depth) + s(b_temp) + s(b_salin), data = envi, method = 'REML', family = quasibinomial) #b_salin
plot(mod1, pages = 1, scheme = 1, all.terms = TRUE, seWithMean = TRUE)
summary(mod1)

mod2 <- mgcv::gam(allele.table[,291] ~ s(dist) + s(depth) + s(b_temp) + s(b_salin), data = envi, method = 'REML') #depth, b_temp
plot(mod2, pages = 1, scheme = 1, all.terms = TRUE, seWithMean = TRUE)
summary(mod2)

mod3 <- mgcv::gam(allele.table[,499] ~ s(dist) + s(depth) + s(b_temp) + s(b_salin), data = envi, method = 'REML') #dist
plot(mod3, pages = 1, scheme = 1, all.terms = TRUE, seWithMean = TRUE)
summary(mod3)

mod4 <- mgcv::gam(allele.table[,743] ~ s(dist) + s(depth) + s(b_temp) + s(b_salin), data = envi, method = 'REML') #b_salin
plot(mod4, pages = 1, scheme = 1, all.terms = TRUE, seWithMean = TRUE)
summary(mod4)


###############################################################################################################################################
#### Single linear regressions are above. Now I need to fit a GAM for each candidate locus-environmental variable. Then calculate p-value. ####
###############################################################################################################################################
envi.ordered <- envi[with(envi, order(PinskyID)),]
as.character(envi.ordered[,2]) == rownames(adults.nonas.232.1137)
rownames(envi.ordered) <- envi.ordered[,"PinskyID"] # so that environmental variable units in GAM plots won't be standardized

library(mgcv)

names_gams <- vector()
for(i in rda.cans){
  names_gams <- append(names_gams, paste0('dist',i,'gam'))
  names_gams <- append(names_gams, paste0('depth',i,'gam'))
  names_gams <- append(names_gams, paste0('btemp',i,'gam'))
  names_gams <- append(names_gams, paste0('bsalin',i,'gam'))
}

# For loop to fit GAM between all RDA candidates and each 4 enviromental variables
for(i in rda.cans){
  assign(paste0('dist',i,'gam'), mgcv::gam(adults.nonas.232.1137[,i] ~ s(dist), data = envi.ordered, method = 'REML', select = TRUE))
  assign(paste0('depth',i,'gam'), mgcv::gam(adults.nonas.232.1137[,i] ~ s(depth), data = envi.ordered, method = 'REML', select = TRUE))
  assign(paste0('btemp',i,'gam'), mgcv::gam(adults.nonas.232.1137[,i] ~ s(b_temp), data = envi.ordered, method = 'REML', select = TRUE))
  assign(paste0('bsalin',i,'gam'), mgcv::gam(adults.nonas.232.1137[,i] ~ s(b_salin), data = envi.ordered, method = 'REML', select = TRUE))
}

# Function to calculate p-values
gamp <- function (modelobject) {
  if (class(modelobject)[1] != "gam") stop("Not an object of class 'gam' ")
  f <- summary(modelobject)$s.pv
  # p <- pf(f[1],f[2],f[3],lower.tail=F)
  # attributes(p) <- NULL
  return(f)
}

gampvalues <- c(gamp(dist43gam),  gamp(depth43gam),   gamp(btemp43gam),  gamp(bsalin43gam), gamp(dist47gam),    gamp(depth47gam),   gamp(btemp47gam),  gamp(bsalin47gam),   gamp(dist125gam),  gamp(depth125gam),   gamp(btemp125gam),  gamp(bsalin125gam),  gamp(dist291gam)  ,  gamp(depth291gam), 
             gamp(btemp291gam),   gamp(bsalin291gam),  gamp(dist380gam),  gamp(depth380gam), gamp(btemp380gam),   gamp(bsalin380gam),  gamp(dist402gam),  gamp(depth402gam),   gamp(btemp402gam),  gamp(bsalin402gam),  gamp(dist432gam),   gamp(depth432gam),   gamp(btemp432gam),   gamp(bsalin432gam), 
             gamp(dist466gam),  gamp(depth466gam),   gamp(btemp466gam),  gamp(bsalin466gam),  gamp(dist483gam),    gamp(depth483gam),   gamp(btemp483gam),  gamp(bsalin483gam),  gamp(dist499gam),  gamp(depth499gam),   gamp(btemp499gam),  gamp(bsalin499gam),  gamp(dist577gam)  ,  gamp(depth577gam), 
             gamp(btemp577gam),   gamp(bsalin577gam),  gamp(dist579gam),   gamp(depth579gam), gamp(btemp579gam),   gamp(bsalin579gam),  gamp(dist635gam),  gamp(depth635gam),   gamp(btemp635gam),  gamp(bsalin635gam),  gamp(dist695gam),   gamp(depth695gam),   gamp(btemp695gam),   gamp(bsalin695gam), 
             gamp(dist704gam),  gamp(depth704gam),   gamp(btemp704gam),  gamp(bsalin704gam),  gamp(dist743gam),    gamp(depth743gam),   gamp(btemp743gam),  gamp(bsalin743gam),  gamp(dist826gam),  gamp(depth826gam),   gamp(btemp826gam),  gamp(bsalin826gam),  gamp(dist842gam) ,  gamp(depth842gam), 
             gamp(btemp842gam),   gamp(bsalin842gam),  gamp(dist937gam),   gamp(depth937gam), gamp(btemp937gam),   gamp(bsalin937gam),  gamp(dist947gam),  gamp(depth947gam),   gamp(btemp947gam),  gamp(bsalin947gam),  gamp(dist965gam),   gamp(depth965gam),   gamp(btemp965gam),   gamp(bsalin965gam), 
             gamp(dist997gam),  gamp(depth997gam),   gamp(btemp997gam),  gamp(bsalin997gam),  gamp(dist1090gam),   gamp(depth1090gam),  gamp(btemp1090gam),  gamp(bsalin1090gam))

rda.can.loci.gam <- data.frame(names_gams,gampvalues)
rda.can.loci2.gam <- rda.can.loci.gam[which(rda.can.loci.gam[,2] < 0.001),]

#### Nice plot of the significant locus-environmental associations based on GAMs ####
png(file = "/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/GAM/GAM_plots.png", width=8, height=9, res=300, units="in")

par(
  mfrow = c(3, 3), 
  mar=c(3, 3, 2, 1), # panel magin size in "line number" units
  mgp=c(2, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=12, # point size, which is the font size
  bg=NA
)

plot(depth125gam, scheme = 1, seWithMean = TRUE, xlab = 'Depth')
mtext('Contig 8420', side = 3, line = 0.5, adj = 0)
plot(bsalin125gam, scheme = 1, seWithMean = TRUE, xlab = 'Bottom salinity')
mtext('Contig 8420', side = 3, line = 0.5, adj = 0)
plot(depth291gam, scheme = 1, seWithMean = TRUE, xlab = 'Depth')
mtext('Contig 19728', side = 3, line = 0.5, adj = 0)
plot(btemp291gam, scheme = 1, seWithMean = TRUE, xlab = 'Bottom temperature')
mtext('Contig 19728', side = 3, line = 0.5, adj = 0)
plot(dist499gam, scheme = 1, seWithMean = TRUE, xlab = 'Distance')
mtext('Contig 35399', side = 3, line = 0.5, adj = 0)
plot(depth499gam, scheme = 1, seWithMean = TRUE, xlab = 'Depth')
mtext('Contig 35399', side = 3, line = 0.5, adj = 0)
plot(bsalin743gam, scheme = 1, seWithMean = TRUE, xlab = 'Bottom salinity')
mtext('Contig 54288', side = 3, line = 0.5, adj = 0)

dev.off()

anova(dist43, test="Chisq")


dist43 <- mgcv::gam(adults.nonas.232.1137[,43] ~ s(dist), data = envi.ordered.matrix, method = 'REML', select = TRUE)
depth43 <- mgcv::gam(adults.nonas.232.1137[,43] ~ s(depth), data = envi.ordered.matrix, method = 'REML', select = TRUE)
btemp43 <- mgcv::gam(adults.nonas.232.1137[,43] ~ s(b_temp), data = envi.ordered.matrix, method = 'REML', select = TRUE)
bsalin43 <- mgcv::gam(adults.nonas.232.1137[,43] ~ s(b_salin), data = envi.ordered.matrix, method = 'REML', select = TRUE)

plot(dist43, pages = 1, scheme = 1, all.terms = TRUE, seWithMean = TRUE)

test125 <- lm(adults.nonas.232.1137[,125] ~ envi.ordered.matrix[,"b_salin"]+envi.ordered.matrix[,"b_temp"]+envi.ordered.matrix[,"depth"]+envi.ordered.matrix[,"dist"])
test125.1 <- lm(adults.nonas.232.1137[,125] ~ envi.ordered.matrix[,"b_salin"]+envi.ordered.matrix[,"b_temp"]+envi.ordered.matrix[,"dist"])
test125.2 <- lm(adults.nonas.232.1137[,125] ~ envi.ordered.matrix[,"b_salin"])
mod1 <- mgcv::gam(adults.nonas.232.1137[,125] ~ s(dist) + s(depth) + s(b_temp) + s(b_salin), data = envi.ordered.matrix, method = 'REML', select = TRUE)
mod1.1 <- mgcv::gam(adults.nonas.232.1137[,125] ~ s(depth) + s(b_salin), data = envi.ordered.matrix, method = 'REML')
mod1.2 <- mgcv::gam(adults.nonas.232.1137[,125] ~ s(b_salin), data = envi.ordered.matrix, method = 'REML')
AIC(test125) #36.47104
AIC(mod1) #29.30607
AIC(test125.1) #34.92524
AIC(mod1.1)
AIC(test125.2) #35.74117
AIC(mod1.2) # 34.87272

test291 <- lm(adults.nonas.232.1137[,291] ~ envi.ordered.matrix[,"b_salin"]+envi.ordered.matrix[,"b_temp"]+envi.ordered.matrix[,"depth"]+envi.ordered.matrix[,"dist"])
test291.1 <- lm(adults.nonas.232.1137[,291] ~ envi.ordered.matrix[,"b_temp"]+envi.ordered.matrix[,"depth"]+envi.ordered.matrix[,"dist"])
mod2 <- mgcv::gam(adults.nonas.232.1137[,291] ~ s(dist) + s(depth) + s(b_temp) + s(b_salin), data = envi.ordered.matrix, method = 'REML')
AIC(test291) #125.0026
AIC(test291.1) #124.7255
AIC(mod2) #125.0028

test499 <- lm(adults.nonas.232.1137[,499] ~ envi.ordered.matrix[,"b_salin"]+envi.ordered.matrix[,"b_temp"]+envi.ordered.matrix[,"depth"]+envi.ordered.matrix[,"dist"])
test499.1 <- lm(adults.nonas.232.1137[,499] ~ envi.ordered.matrix[,"depth"]+envi.ordered.matrix[,"dist"])
mod3 <- mgcv::gam(adults.nonas.232.1137[,499] ~ s(dist) + s(depth) + s(b_temp) + s(b_salin), data = envi.ordered.matrix, method = 'REML')
AIC(test499) #-170.446
AIC(test499.1) #-172.5401
AIC(mod3) #-169.8731

test743 <- lm(adults.nonas.232.1137[,743] ~ envi.ordered.matrix[,"b_salin"]+envi.ordered.matrix[,"b_temp"]+envi.ordered.matrix[,"depth"]+envi.ordered.matrix[,"dist"])
test743.1 <- lm(adults.nonas.232.1137[,743] ~ envi.ordered.matrix[,"b_salin"])
mod4 <- mgcv::gam(adults.nonas.232.1137[,743] ~ s(dist) + s(depth) + s(b_temp) + s(b_salin), data = envi.ordered.matrix, method = 'REML')
AIC(test743) #-8.783139
AIC(test743.1) #-10.24018
AIC(mod4) #-8.424157

