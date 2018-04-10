setwd("/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis")

library(ade4)
library(adegenet)
library(devtools)
library("hierfstat")
library(pegas)
library(fields)

# Reading in SNP data file containing only the first SNP at each locus for 1137 loci across 232 fish
adults232 <- read.structure("structure_input_Nov_11_2015_minus9.str",
                            n.ind = 232, n.loc = 1137, col.lab = 1, col.pop = 2, row.marknames = 1, 
                            onerowperind = FALSE)

# Reading in adult locations
adults_locs <- read.table("adults_locations.txt", header=TRUE)
adults_locs <- adults_locs[-c(217:225),] # remove 9 of 10 GB fish

is.genind(adults232)
head(indNames(adults232),10)

sum(is.na(adults232$tab)) #3034
X <- scaleGen(adults232, NA.method = "mean")
dim(X)
class (X)

# make PCA
pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3)
barplot(pca1$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))

pca1

# Calculate variation explained by PC axes
eig_percent <- round((pca1$eig/(sum(pca1$eig)))*100,2)
eig_percent [1:3]

# Plotting PC1 and PC2
s.class(pca1$li, pop(adults232))
title("PCA of summer flounder dataset\naxes 1-2")
add.scatter.eig(pca1$eig[1:20], 3,1,2)

# add heat map based on latitude to plot using grey scale (for publication)
# grey.colors(6, start = 0, end = 0.8)
# zCol <- function(nCols, Z){
#   cols <- colorRampPalette(c("#000000", "#626262", "#878787", "#A2A2A2", "#B8B8B8", "#CCCCCC"))(nCols)
#   colVec_ind <- cut(Z, breaks=nCols)
#   colVec <- cols[colVec_ind]
# }

# add heat map colors based on latitude to plot using color blind friendly colors (reviewers want color)
zCol <- function(nCols, Z){
  cols <- colorRampPalette(c("#0072B2", "#56B4E9", "#009E73", "#F0E442", "#D55E00", "#CC79A7"))(nCols)
  colVec_ind <- cut(Z, breaks=nCols)
  colVec <- cols[colVec_ind]
}

### To make a B&W plot of the PCA using 241 fish, 1137 loci and fish with color based on latitude###
# png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Popstructure_and_spatialoutliers_results/232fish_1137loci_pca_b&w.png", width=5, height=4.55, res=300, units="in")
png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Popstructure_and_spatialoutliers_results/232fish_1137loci_pca_color.png", width=5, height=4.55, res=300, units="in")

par(
  mar=c(5, 4, 4, 2), # panel magin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  bg=NA
)

# zCol <- function(nCols, Z){
#   cols <- colorRampPalette(c("#000000", "#626262", "#878787", "#A2A2A2", "#B8B8B8", "#CCCCCC"))(nCols)
#   colVec_ind <- cut(Z, breaks=nCols)
#   colVec <- cols[colVec_ind]
# }

zCol <- function(nCols, Z){
  cols <- colorRampPalette(c("#0072B2", "#56B4E9", "#009E73", "#F0E442", "#D55E00", "#CC79A7"))(nCols)
  colVec_ind <- cut(Z, breaks=nCols)
  colVec <- cols[colVec_ind]
}

Z <- adults_locs[,2]
plot(pca1$li[,1], pca1$li[,2], col=zCol(6, Z), xlab = "PC1 (1.42%)", ylab = "PC2 (1.18%)")

dev.off()

# Some code for making a color bar for the PCA, not perfected yet. Makes plot and then I used PPT to paste it on top of PCA, would be better if I could do PCA and color bar all in one
# png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/colorbar_b&w.png", width=4, height=4.5, res=300, units="in")
png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/colorbar_color.png", width=4, height=4.5, res=300, units="in")
par(
  mar=c(5, 4, 4, 2), # panel magin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  bg=NA,
  pin=c(0.5,3)
)

# Code to make the color bar in black & white
# my.colors = colorRampPalette(c("#000000", "#626262", "#878787", "#A2A2A2", "#B8B8B8", "#CCCCCC")) #b&w
my.colors = colorRampPalette(c("#0072B2", "#56B4E9", "#009E73", "#F0E442", "#D55E00", "#CC79A7")) #color-blind colors
z=matrix(1:6,nrow=1)
x=0.25
y=seq(28.96,41.55,len=7) # range of your data
image(x,y,z,col=my.colors(6),axes=FALSE,xlab='', ylab='',main='')
mtext("Latitude", side=4, line=2.5)
axis(4)
dev.off()

#### Calculate HWE for each locus so that these can be removed for STRUCTURE analysis ####
adults.hwt <- hw.test(adults232, B=10000)
adults.hwt
hist(adults.hwt, nclass =30)
pval <- adults.hwt[adults.hwt[,"Pr.exact"] < 0.0100,] # p<0.01, exact test; like in Wolf population structure & canidate genes under selection paper
length(pval[,"Pr.exact"]) 

# Reading in SNP data file containing only the first SNP at each locus for 1137 loci across 232 fish as a text file. 
# This is only to make the final str file creation easier because it preserves the format of having two rows per individual and one column per locus
adults232_2rows <- read.table("structure_input_Nov_11_2015_minus9.str", skip = 1)

# For loop to create locus names
loci <- 1:1137
for (i in 1:length(loci)){
  loci.names <- as.vector(paste('SNP_', loci, sep = ''))
}

# Add these names to the SNP file
colnames(adults232_2rows) <- c('name', 'pop', loci.names)

# # Fix names in genind object
# allele.names <- data.frame(colnames(adults232@tab))
# snp.names <- do.call(rbind, strsplit(as.character(allele.names[,1]), '[.]'))
# snp.names.hwe <- snp.names[!snp.names[,1] %in% rownames(pval),]
# 
# colnames(adults232@tab) <- as.vector(snp.names[,1])# Stick these fixed names back in the genind@tab object

# Remove loci not in HWE
adults232hwe <- adults232_2rows[,!as.character(colnames(adults232_2rows)) %in% rownames(pval)] # Exclude the loci not in HWE from the SNP data text file
dim(adults232hwe) # 464 x 1007

#### Rebuild the str file ####
# Write str file for 232 fish and 1005 loci in HWE
write.table(adults232hwe, "~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/structure_input_232fish_1005loci.str", sep="\t", col.names = TRUE, row.names = FALSE)
# Look at this file in a text editor and delete 'names' and 'pop' column headers and all the "'s

#### Isolation by distance for 1137 loci across 232 fish ####
setwd("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/")

library(geosphere)
rousset <- read.csv('Rousset_indivdist.csv', header = TRUE, check.names = FALSE, row.names = 1)

GB9 <- c("14230L1439", "14231L1440", "14232L1529", "14233L1588", "14234L1441", "14235L1442", "14236L1530", "14237L1531", "14238L1532")
rousset232 <- rousset[!rownames(rousset) %in% GB9,!colnames(rousset) %in% GB9] # remove 9 of 10 GB fish
dim(rousset232)

locations <- read.table("adults_locations3.txt", header = TRUE, row.names = 1)
locations.names <- do.call(rbind, strsplit(as.character(locations[,1]), '_'))
locations [,1] <- locations.names[,2]
locations232 <- locations[!locations[,1] %in% GB9,] # remove 9 of 10 GB fish
geodist <- distm(locations232[,2:3], fun = distCosine)
geodistkm <- geodist * 0.001

rownames(rousset232) == locations232[,1] # make sure order of genetic and geographic distances are the same: they are

# Makes a nice plot of Rousset's distance vs geographic distance 
png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Popstructure_and_spatialoutliers_results/IBD232.png", width=5, height=4.5, res=300, units="in")
par(
  mar=c(5, 4, 4, 2), # panel magin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  bg=NA
)

plot(rousset232[lower.tri(rousset232)] ~ geodistkm[lower.tri(geodistkm)], xlab = 'Geographic distance (km)', ylab = 'Genetic distance (Rousset)')
linreg <- lm(rousset232[lower.tri(rousset232)] ~ geodistkm[lower.tri(geodistkm)])
summary(linreg)
abline(linreg, col = "gray65")

dev.off()

# Mantel test for IBD using Rousset's individual distances and greater circle distances (works)
class(geodistkm)
class(rousset232)
geo <- as.dist(geodistkm) # b/c mantel.randtest requires objects of class dist
gen <- as.dist(rousset232)
man <- mantel.randtest(geo, gen, nrepet=10000)
man
# plot(man, nclass=30, main='Mantel Test', sub='p=0.8029')
plot(man, nclass=30, main='Mantel Test', sub='p=0.8595')

#### Incorporating MARMAP into IBD ####
setwd("/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis")

library(marmap)
library(geosphere)
library(ade4)

# Read in genetic distances
rousset <- read.csv('Rousset_indivdist.csv', header = TRUE, check.names = FALSE, row.names = 1)

GB9 <- c("14230L1439", "14231L1440", "14232L1529", "14233L1588", "14234L1441", "14235L1442", "14236L1530", "14237L1531", "14238L1532")
rousset232 <- rousset[!rownames(rousset) %in% GB9,!colnames(rousset) %in% GB9] # remove 9 of 10 GB fish
dim(rousset232)

# Reading in adult locations
locations <- read.table("adults_locations3.txt", header = TRUE, row.names = 1)
locations.names <- do.call(rbind, strsplit(as.character(locations[,1]), '_'))
locations [,1] <- locations.names[,2]
locations232 <- locations[!locations[,1] %in% GB9,] # remove 9 of 10 GB fish

# Get bathymetry data
useast <- getNOAA.bathy(lon1 = -85, lon2 = -60, lat1 = 23, lat2 = 48, resolution = 1)

plot(useast)
plot(useast, image = TRUE)
# blues <- colorRampPalette(c("red","purple","blue",
                            # "cadetblue1","white"))
blues <- c("lightsteelblue4", "lightsteelblue3", "lightsteelblue2", "lightsteelblue1")
greys <- c(grey(0.6), grey(0.93), grey(0.99))

png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Maps/marmap_plot.png", width=7, height=7, res=300, units="in")
par(
  mar=c(5, 4, 4, 2), # panel magin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  bg=NA
)

plot(useast, image = TRUE, land = TRUE, lwd = 0.1, bpal = list(c(0, max(useast), greys), c(min(useast), 0, blues)), drawlabels = FALSE)
plot(useast, lwd = 0.8, deep = 0, shallow = 0, step = 0, add = TRUE) #highlight coastline
scaleBathy(useast, deg = 2, x = 'bottomright', inset = 5)
points(locations232$long, locations232$lat, pch = 21, col = 'black', bg = 'white')

# adults_locs_withdepth <- get.depth(useast, locations232$long, locations232$lat, locator = FALSE)
plot(useast, deep = -200, shallow = -200, step = 0, lwd = 0.5, add = TRUE, drawlabels = TRUE) # draw 200 meter isobath (delimiting shelf)

dev.off()

# Doesn't seem there is an easy way to calculate greater circle distance between two points within an isobath, but I can do least cost path within an isobath
# Compute least cost paths
lonlat <- locations232[,-1] # lon column then lat column for dist cost analysis
lonlat2 <- unique(lonlat) # only unique lon/lat for path cost analysis

tr <- trans.mat(useast, min.depth = -0.5, max.depth = -200)
cost <- lc.dist(tr, lonlat2[1:10,], res = 'path') # this is just for plotting on a map; probably want the 'dist' option for IBD; non-unique lat/lon seems to give errors, but maybe also a resolution thing, def a resolution thing

# Add least cost paths to map (if I can get it to work, seems to be computational)
paths <- lapply(cost, lines, col = 'black', lwd = 0.8, lty = 1)

cost2 <- lc.dist(tr, as.data.frame(lonlat), res = 'dist') # for IBD calculations
cost2 <- as.matrix(cost2)

# Makes a nice plot of Rousset's distance vs least cost distance 
png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Popstructure_and_spatialoutliers_results/IBD232_marmap.png", width=5, height=4.5, res=300, units="in")
par(
  mar=c(5, 4, 4, 2), # panel magin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  bg=NA
)

plot(rousset232[lower.tri(rousset232)] ~ cost2[lower.tri(cost2)], xlab = 'Geographic distance (km)', ylab = 'Genetic distance (Rousset)')
linreg <- lm(rousset232[lower.tri(rousset232)] ~ cost2[lower.tri(cost2)])
summary(linreg)
abline(linreg, col = "gray65")

dev.off()

# Mantel test for IBD using Rousset's individual distances and greater circle distances (works)
class(cost2)
class(rousset232)
geo <- as.dist(cost2) # b/c mantel.randtest requires objects of class dist
gen <- as.dist(rousset232)
man <- mantel.randtest(geo, gen, nrepet=10000)
man
plot(man, nclass=30, main='Mantel Test', sub='p=0.8618')
