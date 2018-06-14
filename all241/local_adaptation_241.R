setwd("/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation/all241")

library(maps)
library(mapdata)
library(geosphere)
library(ade4)
library(adegenet)
library(devtools)
library(hierfstat)
library(pegas)
library(vegan)

data <- read.csv('allpops_241.csv')

#### Instead of lat & long, add a column for greater-circle distance from a southern point ####
longlat <- data[,c("long", "lat")] # subset data to only longitude and latitude
longlat <- rbind(longlat, c(-80.546297, 28.384993)) # add a southern-most point
geodist <- distm(longlat, fun=distCosine) #matrix is in meters
hist(geodist[lower.tri(geodist)], nclass = 20)

# Distance matrix in km
geodistkm <- geodist * 0.001

# Add greater circle distance to data, minus the south reference point
data$dist <- geodistkm[242,1:241]

write.table(data, "241envirowithdist.txt", sep = "\t", col.names = TRUE, row.names = FALSE)

#### Bayenv2 implementation to detect Fst outliers ####
# Creating the ENVIROFILE using the average depth, bottom temperature, bottom salinity and distance from a southern point
# Aggregate the individuals into populations and standardize
by <- list(data$Population)
pop.avg <- aggregate(data[,6:9], by = by, FUN = mean)
stan.pop.avg <- scale(pop.avg[,2:5])
rownames(stan.pop.avg) <- c("Pop1", "Pop2", "Pop3", "Pop4", "Pop5")
t.stan.pop.avg <- t(stan.pop.avg) # the order is Pop1, Pop2, Pop3, Pop4, Pop5; rows are: depth, b_temp, b_salin, dist
write.table (t.stan.pop.avg, "241stan4envirofile.txt", sep ="\t", col.names = FALSE, row.names = FALSE)

#### Reading in SNP data file containing only the first SNP at each locus ####
adults241 <- read.structure("structure_input_241forbayenv.str",
                             n.ind = 241, n.loc = 1137, col.lab = 1, col.pop = 2, row.marknames = 1, 
                             onerowperind = FALSE)

adults241.df <- as.data.frame(adults241)

#### Which loci are not in HWP? ####
adults.hwt <- hw.test(adults241, B=10000) # B=number of replicates for MCMC procedure, or B=0 for regular HW test; needs to be genind object
adults.hwt
pval <- adults.hwt[adults.hwt[,"Pr.exact"] < 0.0100,] # p<0.01, exact test; like in Wolf population structure & canidate genes under selection paper
length(pval[,"Pr.exact"]) #132 SNPs with less than 0.01 probability of being in HWE; this number differs depending on whether Pr(chi^2 >) [138] or Pr.exact is used 

# Isolate SNP names that I want to exclude
# For genind object
# snp.names <- colnames(adults241.df) # isolate names
# snp.names.split <- do.call(rbind, strsplit(as.character(snp.names), '.', fixed = TRUE)) # split locus from allele
# snp.names.split.sub <- snp.names.split[snp.names.split[,1] %in% rownames(pval),] # subset to only loci not in HWE
# snp.names.exclude <- paste(snp.names.split.sub[,1], snp.names.split.sub[,2], sep = '.') # put locus and allele back together into a single name
# For .str file type
adults242.str <- read.table('structure_input_241forbayenv.txt', header = TRUE)

# Remove these from the the .str file
# adults241_968loci <- adults241.df[, !colnames(adults241.df) %in% snp.names.exclude] # 241 x 1939, for genind object
adults241_968loci <- adults242.str[, ! colnames(adults242.str) %in% rownames(pval)]

# Write to file and save as a .str file afterwards
write.table(adults241_968loci, "structure_input_241forbayenv_968loci.txt", col.names = TRUE, row.names = FALSE, sep = '\t')
# Find and replace all " with nothing

# Test to make sure it's being read in properly. Yes!
iguana <- read.structure("structure_input_241forbayenv_968loci.str",
                            n.ind = 241, n.loc = 968, col.lab = 1, col.pop = 2, row.marknames = 1, 
                            onerowperind = FALSE)

# Average last 40 covariance matrices
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

data_array <- read_jho(file="241matrix_hwe.out", n_rows=5, n_matrices=200, n_comment=2, n_top_comment=13, what="character")

# Take average
data_array_mean <- apply(data_array[,,160:200], c(1,2), mean)
data_array_mean

# Save averaged covariance matrix
write.table(data_array_mean, "241covarmatrix_hwe.txt", row.names = FALSE, col.names = FALSE, sep = '\t')

#### Reading in 10 BayEnv runs (depth, b_temp, b_salin, distance from a southern point) ####
library(abind)
# Randomization #1
setwd("/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation/all241/results_avg241covarmatrix_hwe_4stanenviro")

temp <- list.files(pattern = "*.txt")
list1 <- lapply(temp, read.table)
randarray1 <- array(unlist(list1), dim = c(nrow(list1[[1]]), ncol(list1[[1]]), length(list1)))
randarray1 <- randarray1[,-1,] # gets rid of snpname
randarray1.median <- apply(randarray1, c(1,2), median)
colnames(randarray1.median) <- c("depth", "b_temp", "b_salin", "dist")

which(randarray1.median[,1] > 3) # 35  125  276  432  609  615  703  937 1020 1081 1130
which(randarray1.median[,2] > 3) # 125  214  432  442  524  609  615  626  652  703  919 1050 1096
which(randarray1.median[,3] > 3) # 35 125 432 609 615 703 743 937
which(randarray1.median[,4] > 3) # 35  125  214  432  442  609  615  626  703  937 1020 1081 1096

sort(unique(c(which(randarray1.median[,1] > 3), which(randarray1.median[,2] > 3), which(randarray1.median[,3] > 3), which(randarray1.median[,4] > 3))))

# Get actual BF values
randarray1.median[which(randarray1.median[,1] > 3),1] # depth
randarray1.median[which(randarray1.median[,2] > 3),2] #b_temp
randarray1.median[which(randarray1.median[,3] > 3),3] #b_salin
randarray1.median[which(randarray1.median[,4] > 3),4] #b_salin

###############################################
### Redundancy analysis using all 241 fish ####
###############################################
library(ade4)
library(adegenet)
library(devtools)
library("hierfstat")
library(pegas)
library(fields)
library(vegan)

setwd("/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/")

# Reading in SNP data file containing only the first SNP at each locus
adults <- read.structure("structure_input_Nov_11_2015.str",
                          n.ind = 241, n.loc = 1137, col.lab = 1, col.pop = 2, row.marknames = 1, 
                          onerowperind = FALSE)

which(adults@loc.n.all > 2) # which snps have more than 2 alelles?
adults.nonas <- scaleGen(adults, center = TRUE, scale = FALSE, NA.method = "mean") # filling in NAs with allele frequencies and centering
hist(adults.nonas)
sum(is.na(adults.nonas)) #0 All NAs have been replaced with means
adults.nonas.241.2274 <- adults.nonas[, -c(56, 96, 734, 2115)]
dim(adults.nonas.241.2274)
adults.nonas.241.1137 <- adults.nonas.241.2274[, seq(1, ncol(adults.nonas.241.2274),by = 2)] # including every other column
dim(adults.nonas.241.1137)

# Bringing in environmental data
# This contains distance from a southern point
envi <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation/all241/241envirowithdist.txt", header = TRUE)
envi.ordered <- envi[with(envi, order(PinskyID)),]
as.character(envi.ordered[,3]) == rownames(adults.nonas.241.1137)
envi.ordered.matrix <- scale(data.matrix(envi.ordered))
envi.ordered.matrix
rownames(envi.ordered.matrix) <- envi.ordered[,"PinskyID"]
envi.ordered.matrix[,1] <- envi.ordered[,1]

# Perform RDA
envi.ordered.matrix <- as.data.frame(envi.ordered.matrix)
adults.rda <- rda(adults.nonas.241.1137 ~ dist + depth + b_temp + b_salin, envi.ordered.matrix, scale = FALSE)
adults.rda

spp.scr <- scores(adults.rda, display = "species", scaling = 0, choices = c(1,2,3,4))
head(spp.scr)

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

rda.cans <- unique(sort(rda.cans)) # 4 enviro variables, including distance from a southern point: 35  43 221 291 363 466 568 704 716 842 859 947 965 997 999

# Create a vector of all lm names - needs to reflect order that lms are in the list
names_lms <- vector()
for(i in rda.cans){
  names_lms <- append(names_lms, paste0('dist',i))
}

for(i in rda.cans){
  names_lms <- append(names_lms, paste0('depth',i))
}

for(i in rda.cans){
  names_lms <- append(names_lms, paste0('btemp',i))
}

for(i in rda.cans){
  names_lms <- append(names_lms, paste0('bsalin',i))
}

# For loop to regress all RDA candidates against 4 enviromental variables
lms.dist <- lapply(rda.cans, function(x) lm(adults.nonas.241.1137[,x] ~ envi.ordered.matrix[,"dist"]))
lms.depth <- lapply(rda.cans, function(x) lm(adults.nonas.241.1137[,x] ~ envi.ordered.matrix[,"depth"]))
lms.bt <- lapply(rda.cans, function(x) lm(adults.nonas.241.1137[,x] ~ envi.ordered.matrix[,"b_temp"]))
lms.bs <- lapply(rda.cans, function(x) lm(adults.nonas.241.1137[,x] ~ envi.ordered.matrix[,"b_salin"]))

all.lms <- c(lms.dist, lms.depth, lms.bt, lms.bs)

# Using the list of lm objects to put each into the lmp function to calculate p-values
# Function to calculate p-values that's needed in the for loop below
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

pvalues <- vector()
for (i in 1:length(all.lms)){
  pvalues[i] <- lmp(all.lms[[i]])
}

rda.can.loci <- data.frame(names_lms,pvalues)
rda.can.loci2 <- rda.can.loci[which(rda.can.loci[,2] < 0.001),]


