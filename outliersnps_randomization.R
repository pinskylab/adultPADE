#### Code to randomize allele frequencies for use in BayEnv and RDA sensitivity analyses ####

setwd("/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation")

library(ade4)
library(adegenet)
library(devtools)
library(hierfstat)
library(pegas)

# Reading in SNP data file containing only the first SNP at each locus
adults <- read.table("structure_input_232forbayenv.txt", header = TRUE)
adults[adults == -9] <- NA # missing data is currently coded as a -9, replace all -9's with NA
adults <- adults[,-c(1:2)]
popnames <- adults[,c(1:2)] # put population and ID columns into an object called names

#### Randomize each column ####
# test code
df1 <- data.frame(A=c(1,1,0,0), B=c(0,0,0,1), C=c(1,0,0,0), D = c(1,1,1,0), E= c(1,0,1,1)) # randomization test
df2 <- apply(df1,2, sample)

# Code for real PADE data
adults2 <- apply(adults, 2, function(x){sample(x[!is.na(x)])}) # randomizes each column, except NAs

# Column sums stay the same? yes.
sums1 <- colSums(adults[,-c(1:2)],na.rm = TRUE) # sum of all snps
sum(adults2[[1]]) # sum of first snp
sums1[1000] == sum(adults2[[1000]])

##### Now I need to recreate the df with NAs and the randomized alleles. Basically, randomize everything but NAs ####
# test code
x <- data.frame(a = c(1,1,NA,NA,1,2,1,1,1,2), b = c(NA,NA,2,3,2,3,2,3,NA,NA), c = c(3,4,3,4,3,4,3,4,3,4)) 
x2 <- apply(x, 2, function(x){sample(x[!is.na(x)])})
x$a <- replace(x$a, !is.na(x$a), x2[,"a"])
x$b <- replace(x$b, !is.na(x$b), x2[,"b"])
x$c <- replace(x$c, !is.na(x$c), x2[,"c"])

for (i in 1:ncol(x)){
  x[,i] <- replace(x[,i], !is.na(x[,i]), x2[[i]])
}

# Code for real PADE data  
for (i in 1:ncol(adults)){
  adults[,i] <- replace(adults[,i], !is.na(adults[,i]), adults2[[i]])
}

# For SNPs with no missing data, these should all be TRUE
adults[,1] == adults2[[1]]

#### Final data manipulation steps to get randomized .str file ####
# Cbind Fish ID and populaltion number to randomized alleles
adults_rand <- cbind(popnames, adults)
dim(adults_rand) # 464 x 1139

# Replace -9's with NA's
adults_rand[is.na(adults_rand)] <- -9

# Write randomized adult allele frequencies to the Local adaptation folder on my computer. Modify on my computer into .str file so that it can be read into PGDspider.
write.table (adults_rand, "structure_input_232forbayenv_randomized.txt", sep ="\t", col.names = TRUE, row.names = FALSE) # manually deleted ""s and column names that don't denote SNP names

# Test to see whether .str file is correctly formatted by reading into R as a genind object
# Reading in SNP data file containing only the first SNP at each locus
adults_test <- read.structure("structure_input_232forbayenv_randomized.str",
                             n.ind = 232, n.loc = 1137, col.lab = 1, col.pop = 2, row.marknames = 1, 
                             onerowperind = FALSE)

#####################################################
#### Randomizing environmental variables instead ####
#####################################################
setwd("/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation")

# Read in environmental data
nvs_locs <- read.table("232envirowithdist.txt", header=TRUE) #all 241-9=232 fish divided into 5 populations

# Subset to 4 environmental variables
envir <- nvs_locs[,c('bayenv_pop', 'dist', 'depth', 'b_temp', 'b_salin')]

#### For loop to randomize, standardize and write 10 environmental matrices ####
for (i in 1:10){
  
# Randomize environmental variables
envir2 <- apply(envir[,-1],2, sample)

# Aggregate the individuals into populations and standardize
by <- list(nvs_locs$bayenv_pop)
pop.avg <- aggregate(envir2, by = by, FUN = mean)
stan.pop.avg <- scale(pop.avg[,2:5])
rownames(stan.pop.avg) <- c("Pop1", "Pop2", "Pop3", "Pop4", "Pop5")
t.stan.pop.avg <- t(stan.pop.avg) # the order is Pop1, Pop2, Pop3, Pop4, Pop5; rows are: dist, depth, b_temp, b_salin

# Exporting the randomized envirofiles to my computer in a randomization folder
write.table (t.stan.pop.avg, paste("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation/randomization/232stan4envirofile", i,".txt", sep = ""), sep ="\t", col.names = FALSE, row.names = FALSE)

}

#### Reading in 10 runs for each of 10 BayEnv analyses (distance from a southern point, depth, b_temp, b_salin) ####
library(abind)
# Randomization #1
setwd("/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation/randomization/rand1")

temp <- list.files(pattern = "*.txt")
list1 <- lapply(temp, read.table)
randarray1 <- array(unlist(list1), dim = c(nrow(list1[[1]]), ncol(list1[[1]]), length(list1)))
randarray1 <- randarray1[,-1,]
randarray1.median <- apply(randarray1, c(1,2), median)
colnames(randarray1.median) <- c("dist", "depth", "b_temp", "b_salin")

which(randarray1.median[,1] > 3) #10 125 228 442 914
which(randarray1.median[,2] > 3) #79 118 706 755 785 951
which(randarray1.median[,3] > 3) #118 785 794 911
which(randarray1.median[,4] > 3) #153  228  808  934 1080

# Randomization #2
setwd("/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation/randomization/rand2")

temp <- list.files(pattern = "*.txt")
list2 <- lapply(temp, read.table)
randarray2 <- array(unlist(list2), dim = c(nrow(list2[[1]]), ncol(list2[[1]]), length(list2)))
randarray2 <- randarray2[,-1,]
randarray2.median <- apply(randarray2, c(1,2), median)
colnames(randarray2.median) <- c("dist", "depth", "b_temp", "b_salin")

which(randarray2.median[,1] > 3) #263 581 725 743 882
which(randarray2.median[,2] > 3) #228 263 396 825 889
which(randarray2.median[,3] > 3) #22   71  153  228 1003 1080
which(randarray2.median[,4] > 3) #78  104  125  442  499  609  626  703  829  845  911  914  916 1074 1096

# Randomization #3
setwd("/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation/randomization/rand3")

temp <- list.files(pattern = "*.txt")
list3 <- lapply(temp, read.table)
randarray3 <- array(unlist(list3), dim = c(nrow(list3[[1]]), ncol(list3[[1]]), length(list3)))
randarray3 <- randarray3[,-1,]
randarray3.median <- apply(randarray3, c(1,2), median)
colnames(randarray3.median) <- c("dist", "depth", "b_temp", "b_salin")

which(randarray3.median[,1] > 3) #10  228  826  968  981 1125
which(randarray3.median[,2] > 3) #609 615 703 743 990
which(randarray3.median[,3] > 3) #263 480 725 743 785 809 889
which(randarray3.median[,4] > 3) #103 414 524 635 706 757 881

# Randomization #4
setwd("/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation/randomization/rand4")

temp <- list.files(pattern = "*.txt")
list4 <- lapply(temp, read.table)
randarray4 <- array(unlist(list4), dim = c(nrow(list4[[1]]), ncol(list4[[1]]), length(list4)))
randarray4 <- randarray4[,-1,]
randarray4.median <- apply(randarray4, c(1,2), median)
colnames(randarray4.median) <- c("dist", "depth", "b_temp", "b_salin")

which(randarray4.median[,1] > 3) #125  228  396  825  826  981 1090 1125
which(randarray4.median[,2] > 3) #125  153  214  442  499  524  609  626  703  829  916  919 1050 1096
which(randarray4.median[,3] > 3) #104 120 435 911
which(randarray4.median[,4] > 3) #153  214  496  524  626  790  916  919 1050

# Randomization #5
setwd("/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation/randomization/rand5")

temp <- list.files(pattern = "*.txt")
list5 <- lapply(temp, read.table)
randarray5 <- array(unlist(list5), dim = c(nrow(list5[[1]]), ncol(list5[[1]]), length(list5)))
randarray5 <- randarray5[,-1,]
randarray5.median <- apply(randarray5, c(1,2), median)
colnames(randarray5.median) <- c("dist", "depth", "b_temp", "b_salin")

which(randarray5.median[,1] > 3) #609 615 743 865 886 990
which(randarray5.median[,2] > 3) #114  132  496  751  794  826  968 1125
which(randarray5.median[,3] > 3) #10  20  71 228 396
which(randarray5.median[,4] > 3) #71  245  743  865  990 1003

# Randomization #6
setwd("/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation/randomization/rand6")

temp <- list.files(pattern = "*.txt")
list6 <- lapply(temp, read.table)
randarray6 <- array(unlist(list6), dim = c(nrow(list6[[1]]), ncol(list6[[1]]), length(list6)))
randarray6 <- randarray6[,-1,]
randarray6.median <- apply(randarray6, c(1,2), median)
colnames(randarray6.median) <- c("dist", "depth", "b_temp", "b_salin")

which(randarray6.median[,1] > 3) #71  79 609 706 757 990
which(randarray6.median[,2] > 3) #10  20  71 228 462 706 755 820
which(randarray6.median[,3] > 3) #543  664  755  809  826  981 1051 1090
which(randarray6.median[,4] > 3) #125  153  214  442  499  524  609  615  626  703  916  919 1050 1096

# Randomization #7
setwd("/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation/randomization/rand7")

temp <- list.files(pattern = "*.txt")
list7 <- lapply(temp, read.table)
randarray7 <- array(unlist(list7), dim = c(nrow(list7[[1]]), ncol(list7[[1]]), length(list7)))
randarray7 <- randarray7[,-1,]
randarray7.median <- apply(randarray7, c(1,2), median)
colnames(randarray7.median) <- c("dist", "depth", "b_temp", "b_salin")

which(randarray7.median[,1] > 3) #153 524 543 826 916
which(randarray7.median[,2] > 3) #135 263 751 794 882
which(randarray7.median[,3] > 3) #153 570 790 916
which(randarray7.median[,4] > 3) #132  826  882  968 1125

# Randomization #8
setwd("/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation/randomization/rand8")

temp <- list.files(pattern = "*.txt")
list8 <- lapply(temp, read.table)
randarray8 <- array(unlist(list8), dim = c(nrow(list8[[1]]), ncol(list8[[1]]), length(list8)))
randarray8 <- randarray8[,-1,]
randarray8.median <- apply(randarray8, c(1,2), median)
colnames(randarray8.median) <- c("dist", "depth", "b_temp", "b_salin")

which(randarray8.median[,1] > 3) #120  435  599  609  794  911 1074
which(randarray8.median[,2] > 3) #125 153 214 228 396 442 570 808 916 919
which(randarray8.median[,3] > 3) #71  615  743  886  990 1003 1031
which(randarray8.median[,4] > 3) #114  120  132  435  609  703  751  794  911 1074 1118

# Randomization #9
setwd("/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation/randomization/rand9")

temp <- list.files(pattern = "*.txt")
list9 <- lapply(temp, read.table)
randarray9 <- array(unlist(list9), dim = c(nrow(list9[[1]]), ncol(list9[[1]]), length(list9)))
randarray9 <- randarray9[,-1,]
randarray9.median <- apply(randarray9, c(1,2), median)
colnames(randarray9.median) <- c("dist", "depth", "b_temp", "b_salin")

which(randarray9.median[,1] > 3) #581 882
which(randarray9.median[,2] > 3) #79  114  125  226  299  609  703  751  794 1009 1074
which(randarray9.median[,3] > 3) #132  794  826  882  968  981 1125
which(randarray9.median[,4] > 3) #125  334  609  615  647  703  826 1074 1090 1125

# Randomization #10
setwd("/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation/randomization/rand10")

temp <- list.files(pattern = "*.txt")
list10 <- lapply(temp, read.table)
randarray10 <- array(unlist(list10), dim = c(nrow(list10[[1]]), ncol(list10[[1]]), length(list10)))
randarray10 <- randarray10[,-1,]
randarray10.median <- apply(randarray10, c(1,2), median)
colnames(randarray10.median) <- c("dist", "depth", "b_temp", "b_salin")

which(randarray10.median[,1] > 3) #263 581 725 743 882
which(randarray10.median[,2] > 3) #609 615 703 881 990
which(randarray10.median[,3] > 3) #153 496 524 751 794
which(randarray10.median[,4] > 3) #22  71 706 990

#### Cumulative density function of all median BF's for each environmental variable ####
par(mfrow=c(2,2))
plot.ecdf(log10(randarray1.median[,1]), xlab = expression(paste('log'[10], '(Bayes Factor)')), ylab = 'Fraction of BFs', col = 'gray60', main = 'Distance from southern point')
plot.ecdf(log10(randarray2.median[,1]), add = TRUE, col = 'gray60')
plot.ecdf(log10(randarray3.median[,1]), add = TRUE, col = 'gray60')
plot.ecdf(log10(randarray4.median[,1]), add = TRUE, col = 'gray60')
plot.ecdf(log10(randarray5.median[,1]), add = TRUE, col = 'gray60')
plot.ecdf(log10(randarray6.median[,1]), add = TRUE, col = 'gray60')
plot.ecdf(log10(randarray7.median[,1]), add = TRUE, col = 'gray60')
plot.ecdf(log10(randarray8.median[,1]), add = TRUE, col = 'gray60')
plot.ecdf(log10(randarray9.median[,1]), add = TRUE, col = 'gray60')
plot.ecdf(log10(randarray10.median[,1]), add = TRUE, col = 'gray60')
plot.ecdf(log10(full_array_median[,4]), add = TRUE, col = 'tomato') # plot cdf of real data

plot.ecdf(log10(randarray1.median[,2]), xlab = expression(paste('log'[10], '(Bayes Factor)')), ylab = 'Fraction of BFs', col = 'gray60', main = 'Depth')
plot.ecdf(log10(randarray2.median[,2]), add = TRUE, col = 'gray60')
plot.ecdf(log10(randarray3.median[,2]), add = TRUE, col = 'gray60')
plot.ecdf(log10(randarray4.median[,2]), add = TRUE, col = 'gray60')
plot.ecdf(log10(randarray5.median[,2]), add = TRUE, col = 'gray60')
plot.ecdf(log10(randarray6.median[,2]), add = TRUE, col = 'gray60')
plot.ecdf(log10(randarray7.median[,2]), add = TRUE, col = 'gray60')
plot.ecdf(log10(randarray8.median[,2]), add = TRUE, col = 'gray60')
plot.ecdf(log10(randarray9.median[,2]), add = TRUE, col = 'gray60')
plot.ecdf(log10(randarray10.median[,2]), add = TRUE, col = 'gray60')
plot.ecdf(log10(full_array_median[,1]), add = TRUE, col = 'tomato') # plot cdf of real data

plot.ecdf(log10(randarray1.median[,3]), xlab = expression(paste('log'[10], '(Bayes Factor)')), ylab = 'Fraction of BFs', col = 'gray60', main = 'Bottom Temp')
plot.ecdf(log10(randarray2.median[,3]), add = TRUE, col = 'gray60')
plot.ecdf(log10(randarray3.median[,3]), add = TRUE, col = 'gray60')
plot.ecdf(log10(randarray4.median[,3]), add = TRUE, col = 'gray60')
plot.ecdf(log10(randarray5.median[,3]), add = TRUE, col = 'gray60')
plot.ecdf(log10(randarray6.median[,3]), add = TRUE, col = 'gray60')
plot.ecdf(log10(randarray7.median[,3]), add = TRUE, col = 'gray60')
plot.ecdf(log10(randarray8.median[,3]), add = TRUE, col = 'gray60')
plot.ecdf(log10(randarray9.median[,3]), add = TRUE, col = 'gray60')
plot.ecdf(log10(randarray10.median[,3]), add = TRUE, col = 'gray60')
plot.ecdf(log10(full_array_median[,2]), add = TRUE, col = 'tomato') # plot cdf of real data

plot.ecdf(log10(randarray1.median[,4]), xlab = expression(paste('log'[10], '(Bayes Factor)')), ylab = 'Fraction of BFs', col = 'gray60', main = 'Bottom Salinity')
plot.ecdf(log10(randarray2.median[,4]), add = TRUE, col = 'gray60')
plot.ecdf(log10(randarray3.median[,4]), add = TRUE, col = 'gray60')
plot.ecdf(log10(randarray4.median[,4]), add = TRUE, col = 'gray60')
plot.ecdf(log10(randarray5.median[,4]), add = TRUE, col = 'gray60')
plot.ecdf(log10(randarray6.median[,4]), add = TRUE, col = 'gray60')
plot.ecdf(log10(randarray7.median[,4]), add = TRUE, col = 'gray60')
plot.ecdf(log10(randarray8.median[,4]), add = TRUE, col = 'gray60')
plot.ecdf(log10(randarray9.median[,4]), add = TRUE, col = 'gray60')
plot.ecdf(log10(randarray10.median[,4]), add = TRUE, col = 'gray60')
plot.ecdf(log10(full_array_median[,3]), add = TRUE, col = 'tomato') # plot cdf of real data

#### Two sample K-S test comparing real data to each randomized data set ####
list <- list(randarray1.median[,1],
          randarray2.median[,1],
          randarray3.median[,1],
          randarray4.median[,1],
          randarray5.median[,1],
          randarray6.median[,1],
          randarray7.median[,1],
          randarray8.median[,1],
          randarray9.median[,1],
          randarray10.median[,1])

distance <- vector()

for (i in 1:length(list)){
  kstest <- ks.test(full_array_median[,4], list[[i]])
  distance[i] <- kstest$p.value
}

distance.adj <- round(p.adjust(distance, method = "BH"), 5)

#### Comparing all randomized BFs for a given environmental variable to that of real data ####
all.randomized.medians <- abind(randarray1.median,randarray2.median, randarray3.median, randarray4.median, 
                              randarray5.median, randarray6.median, randarray7.median, randarray8.median, 
                              randarray9.median, randarray10.median, along = 3) # dist, depth, b_temp, b_salin

ks1 <- ks.test(all.randomized.medians[,"dist",], full_array_median[,"dist"])
ks2 <- ks.test(all.randomized.medians[,"depth",], full_array_median[,"depth"])
ks3 <- ks.test(all.randomized.medians[,"b_temp",], full_array_median[,"b_temp"])
ks4 <- ks.test(all.randomized.medians[,"b_salin",], full_array_median[,"b_salin"])

par(mfrow = c(2,2))
plot.ecdf(log10(all.randomized.medians[,"dist",]), main = 'Distance from southern point', ylab = 'Proportion of BFs', xlab = expression(paste('log'[10], '(Bayes Factor)')))
plot.ecdf(log10(full_array_median[,"dist"]), col = 'tomato', add = TRUE)
text(0.5, 0.2,  paste('D = ', round(ks1$statistic, 3), '\n p-value =', round(ks1$p.value, 3)))
plot.ecdf(log10(all.randomized.medians[,"depth",]), main = 'Depth', ylab = 'Proportion of BFs', xlab = expression(paste('log'[10], '(Bayes Factor)')))
plot.ecdf(log10(full_array_median[,"depth"]), col = 'tomato', add = TRUE)
text(0.5, 0.2,  paste('D = ', round(ks2$statistic, 4), '\n p-value =', round(ks2$p.value, 3)))
plot.ecdf(log10(all.randomized.medians[,"b_temp",]), main = 'Bottom Temperature', ylab = 'Proportion of BFs', xlab = expression(paste('log'[10], '(Bayes Factor)')))
plot.ecdf(log10(full_array_median[,"b_temp"]), col = 'tomato', add = TRUE)
text(0.5, 0.2,  paste('D = ', round(ks3$statistic, 4), '\n p-value =', round(ks3$p.value, 3)))
plot.ecdf(log10(all.randomized.medians[,"b_salin",]), main = 'Bottom Salinity', ylab = 'Proportion of BFs', xlab = expression(paste('log'[10], '(Bayes Factor)')))
plot.ecdf(log10(full_array_median[,"b_salin"]), col = 'tomato', add = TRUE)
text(0.5, 0.2,  paste('D = ', round(ks4$statistic, 3), '\n p-value =', round(ks4$p.value, 3)))

#### What about exact binomial test for each environmental variable? ####
# First need to figure out the 'probability of success' from the randomized data

# Distance from southern point
list <- list(randarray1.median[,1],
             randarray2.median[,1],
             randarray3.median[,1],
             randarray4.median[,1],
             randarray5.median[,1],
             randarray6.median[,1],
             randarray7.median[,1],
             randarray8.median[,1],
             randarray9.median[,1],
             randarray10.median[,1])

distance <- vector()

for (i in 1:length(list)){
  distance[i] <- length(which(list[[i]] > 3))
}

# Depth
list <- list(randarray1.median[,2],
             randarray2.median[,2],
             randarray3.median[,2],
             randarray4.median[,2],
             randarray5.median[,2],
             randarray6.median[,2],
             randarray7.median[,2],
             randarray8.median[,2],
             randarray9.median[,2],
             randarray10.median[,2])

depth <- vector()

for (i in 1:length(list)){
  depth[i] <- length(which(list[[i]] > 3))
}

# Bottom Temperature
list <- list(randarray1.median[,3],
             randarray2.median[,3],
             randarray3.median[,3],
             randarray4.median[,3],
             randarray5.median[,3],
             randarray6.median[,3],
             randarray7.median[,3],
             randarray8.median[,3],
             randarray9.median[,3],
             randarray10.median[,3])

bottomtemp <- vector()

for (i in 1:length(list)){
  bottomtemp[i] <- length(which(list[[i]] > 3))
}

# Bottom Salinity
list <- list(randarray1.median[,4],
             randarray2.median[,4],
             randarray3.median[,4],
             randarray4.median[,4],
             randarray5.median[,4],
             randarray6.median[,4],
             randarray7.median[,4],
             randarray8.median[,4],
             randarray9.median[,4],
             randarray10.median[,4])

bottomsalin <- vector()

for (i in 1:length(list)){
  bottomsalin[i] <- length(which(list[[i]] > 3))
}

mean(distance)
mean(depth)
mean(bottomtemp)
mean(bottomsalin)

# of 'successes', # of trials, prob of success from randomized data
binom.test(9, 1137, mean(distance)/1137) # Distace from southern point
binom.test(6, 1137, mean(depth)/1137) # Depth
binom.test(10, 1137, mean(bottomtemp)/1137) # Bottom Temp
binom.test(6, 1137, mean(bottomsalin)/1137) # Bottom Salinity

# Exact binomial test treating all locus x permutations as independent, this results in the same pvalues as above
binom.test(9, 1137, sum(distance)/(10*1137)) # Distace from southern point
binom.test(6, 1137, sum(depth)/(10*1137)) # Depth
binom.test(10, 1137, sum(bottomtemp)/(10*1137)) # Bottom Temp
binom.test(6, 1137, sum(bottomsalin)/(10*1137)) # Bottom Salinity

binom.test(14, 1137, mean(c(17, 29, 24, 27, 21, 32, 15, 29, 22, 18))/1137) # all loci together across environmental variables, prob not best method because certain environmental variables (bottom temp, distance from a southern point) may be swamped out

#####################################################
#### Barplots of observed vs permuted BFs counts ####
#####################################################
all.randomized.medians <- abind(randarray1.median,randarray2.median, randarray3.median, randarray4.median, 
                                randarray5.median, randarray6.median, randarray7.median, randarray8.median, 
                                randarray9.median, randarray10.median, along = 3) # dist, depth, b_temp, b_salin

save(all.randomized.medians, file = "/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation/randomization/randomized_array10.RData")

load(file = "/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation/randomization/randomized_array10.RData")
all.randomized.medians.mean <- apply(all.randomized.medians, c(1,2), mean)
# se <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x))) # Function to calculate SE for permuted data

# Plot all the barplots together
png(file = "/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation/randomization/BF_counts.png", width=8, height=8, res=300, units="in")

par(
  mfrow = c(2, 2), 
  mar=c(5, 4, 1.3, 0.7), # panel magin size in "line number" units
  mgp=c(2, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=12, # point size, which is the font size
  bg=NA
)

# Distance
dist.perm <- table(factor(findInterval(all.randomized.medians.mean[,'dist'], seq(0,9, by = 1)), levels = 1:9))
dist.obs <- table(factor(findInterval(full_array_median[,"dist"], seq(0,9, by = 1)), levels = 1:9))
dist.diffs <- dist.obs - dist.perm # needs to contain zeros
dist.perm[5:9] <- NA # so count isn't zero when log10
dist.obs[c(6,8)] <- NA # so count isn't zero when log10
dist.all <- log10(rbind(dist.perm, dist.obs))
barplot(dist.all, beside = TRUE, ylab = '', main = 'Distance', col = c('gray60', 'white'), yaxt = 'n', xaxt = 'n', ylim = c(0,3.2))
axis(1, at=seq(2,30, by = 3), labels = FALSE)
text(c(2, 5,  8, 11, 14, 17, 20, 23, 26), -0.4, c('0.0-0.99','1.0-1.99', '2.0-2.99','3.0-3.99', '4.0-4.99', '5.0-5.99','6.0-6.99', '7.0-7.99', '> 8.0'), srt = 45, xpd = TRUE, cex = 0.89)
mtext("Bayes Factor", side = 1, line = 3.3)
axis(2, at=seq(0,3, by = 1), labels = c('0','10','100', '1000'), las = 2)
mtext("Number of loci", side = 2, line = 2.5)
dist.sd <- aggregate(all.randomized.medians.mean[,'dist'], list(findInterval(all.randomized.medians.mean[,'dist'], seq(0,9, by = 1))), sd) # calculate SD for permuted data by bin
arrows(c(1.5, 4.5, 7.5), log10(dist.perm[c(1:4)]), c(1.5, 4.5, 7.5), log10(dist.perm[c(1:4)] + dist.sd$x), length = 0.05, angle = 90, code = 3)
text(c(2, 5,  8, 11, 14, 17, 20, 23, 26), c(3.13, 1.69, 1, 0.6, 0.6, 0.1, 0.1, 0.1, 0.45), labels = dist.diffs)

legend(16, 2.5,
       legend=c("permuted", "observed"),
       pch=c(22, 22),
       col=c('black', 'black'),
       bty = "n",
       pt.bg = c('gray60','white'),
       y.intersp = 1)

# Depth
depth.perm <- table(factor(findInterval(all.randomized.medians.mean[,'depth'], seq(0,8, by = 1)), levels = 1:8))
depth.obs <- table(factor(findInterval(full_array_median[,"depth"], seq(0,8, by = 1)), levels = 1:8))
depth.diffs <- depth.obs - depth.perm # needs to contain zeros
depth.perm[5:7] <- NA # so count isn't zero when log10
depth.obs[8] <- NA # so count isn't zero when log10
depth.all <- log10(rbind(depth.perm, depth.obs))
barplot(depth.all, beside = TRUE, ylab = '', axisnames = FALSE, main = 'Depth', col = c('gray60', 'white'), yaxt = 'n', xaxt = 'n', ylim = c(0,3.2))
axis(1, at=seq(2,30, by = 3), labels = FALSE)
text(c(2, 5,  8, 11, 14, 17, 20, 23), -0.4, c('0.0-0.99','1.0-1.99', '2.0-2.99','3.0-3.99', '4.0-4.99', '5.0-5.99','6.0-6.99', '> 7.0'), srt = 45, xpd = TRUE, cex = 0.89)
mtext("Bayes Factor", side = 1, line = 3.3)
axis(2, at=seq(0,3, by = 1), labels = c('0','10','100', '1000'), las = 2)
mtext("Number of loci", side = 2, line = 2.5)
depth.sd <- aggregate(all.randomized.medians.mean[,'depth'], list(findInterval(all.randomized.medians.mean[,'depth'], seq(0,8, by = 1))), sd) # calculate SD for permuted data by bin
arrows(c(1.5, 4.5, 7.5, 10.5, 22.5), log10(depth.perm[c(1,2,3,4,8)]), c(1.5, 4.5, 7.5, 10.5, 22.5), log10(depth.perm[c(1,2,3,4,8)] + depth.sd$x), length = 0.05, angle = 90, code = 3)
text(c(2, 5,  8, 11, 14, 17, 20, 23), c(3.13, 1.71, 0.99, 0.76, 0.1, 0.1, 0.41, 0.1), labels = depth.diffs)


legend(16, 2.5,
       legend=c("permuted", "observed"),
       pch=c(22, 22),
       col=c('black', 'black'),
       bty = "n",
       pt.bg = c('gray60','white'),
       y.intersp = 1)

# Bottom temperature
bt.perm <- table(factor(findInterval(all.randomized.medians.mean[,'b_temp'], seq(0,11, by = 1)), levels = 1:11))
bt.obs <- table(factor(findInterval(full_array_median[,"b_temp"], seq(0,11, by = 1)), levels = 1:11))
bt.diffs <- bt.obs - bt.perm # needs to contain zeros
bt.perm[5:11] <- NA # so count isn't zero when log10
bt.obs[c(8, 9, 10)] <- NA # so count isn't zero when log10
bt.all <- log10(rbind(bt.perm, bt.obs))
barplot(bt.all, beside = TRUE, ylab = '', axisnames = FALSE, main = 'Bottom temperature', col = c('gray60', 'white'), yaxt = 'n', xaxt = 'n', ylim = c(0, 3.2))
axis(1, at=seq(2,33, by = 3), labels = FALSE)
text(c(2, 5,  8, 11, 14, 17, 20, 23, 26, 29, 32), -0.4, c('0.0-0.99','1.0-1.99', '2.0-2.99','3.0-3.99', '4.0-4.99', '5.0-5.99','6.0-6.99', '7.0-7.99', '8.00-8.99', '9.00-9.99', '> 10.0'), srt = 45, xpd = TRUE, cex = 0.89)
mtext("Bayes Factor", side = 1, line = 3.3)
axis(2, at=seq(0,3, by = 1), labels = c('0','10','100', '1000'), las = 2)
mtext("Number of loci", side = 2, line = 2.5)
bt.sd <- aggregate(all.randomized.medians.mean[,'b_temp'], list(findInterval(all.randomized.medians.mean[,'b_temp'], seq(0,11, by = 1))), sd) # calculate SD for permuted data by bin
arrows(c(1.5, 4.5, 7.5, 10.5), log10(bt.perm[c(1,2,3,4)]), c(1.5, 4.5, 7.5, 10.5), log10(bt.perm[c(1,2,3,4)] + bt.sd$x), length = 0.05, angle = 90, code = 3)
text(c(2, 5,  8, 11, 14, 17, 20, 23, 26, 29, 32), c(3.14, 1.69, 0.89, 0.8, 0.1, 0.42, 0.1, 0.1, 0.1, 0.1, 0.1), labels = bt.diffs)

legend(16, 2.5,
       legend=c("permuted", "observed"),
       pch=c(22, 22),
       col=c('black', 'black'),
       bty = "n",
       pt.bg = c('gray60', 'white'),
       y.intersp = 1)

# Bottom salinity
bs.perm <- table(factor(findInterval(all.randomized.medians.mean[,'b_salin'], seq(0,10, by = 1)), levels = 1:10))
bs.obs <- table(factor(findInterval(full_array_median[,"b_salin"], seq(0,10, by = 1)), levels = 1:10))
bs.diffs <- bs.obs - bs.perm # needs to contain zeros
bs.perm[c(5,7,9,10)] <- NA # so count isn't zero when log10
bs.obs[5:9] <- NA # so count isn't zero when log10
bs.all <- log10(rbind(bs.perm, bs.obs))
barplot(bs.all, beside = TRUE, ylab = '', axisnames = FALSE, main = 'Bottom salinity', col = c('gray60', 'white'), yaxt = 'n', xaxt = 'n', ylim = c(0, 3.2))
axis(1, at=seq(2,30, by = 3), labels = FALSE)
text(c(2, 5,  8, 11, 14, 17, 20, 23, 26, 29), -0.4, c('0.0-0.99','1.0-1.99', '2.0-2.99','3.0-3.99', '4.0-4.99', '5.0-5.99','6.0-6.99', '7.0-7.99', '8.00-8.99', '> 9.0'), srt = 45, xpd = TRUE, cex = 0.89)
mtext("Bayes Factor", side = 1, line = 3.3)
axis(2, at=seq(0,3, by = 1), labels = c('0','10','100', '1000'), las = 2)
mtext("Number of loci", side = 2, line = 2.5)
bs.sd <- aggregate(all.randomized.medians.mean[,'b_salin'], list(findInterval(all.randomized.medians.mean[,'b_salin'], seq(0,10, by = 1))), sd) # calculate SD for permuted data by bin
arrows(c(1.5, 4.5, 7.5, 10.5), log10(bs.perm[c(1,2,3,4)]), c(1.5, 4.5, 7.5, 10.5), log10(bs.perm[c(1,2,3,4)] + bs.sd$x[1:4]), length = 0.05, angle = 90, code = 3)
text(c(2, 5,  8, 11, 14, 17, 20, 23, 26, 29), c(3.13, 1.69, 1.05, 0.71, 0.1, 0.1, 0.1, 0.1, 0.1, 0.42), labels = bs.diffs)

legend(16, 2.5,
       legend=c("permuted", "observed"),
       pch=c(22, 22),
       col=c('black', 'black'),
       bty = "n",
       pt.bg = c('gray60', 'white'),
       y.intersp = 1)

dev.off()

# KS test using mean of 10 random permutations vs observed data for each environmental variable
ks1 <- ks.test(all.randomized.medians.mean[,'dist'], full_array_median[,"dist"])
ks2 <- ks.test(all.randomized.medians.mean[,'depth'], full_array_median[,"depth"])
ks3 <- ks.test(all.randomized.medians.mean[,'b_temp'], full_array_median[,"b_temp"])
ks4 <- ks.test(all.randomized.medians.mean[,'b_salin'], full_array_median[,"b_salin"])

png(file = "/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation/randomization/CDF_plots.png", width=8, height=8, res=300, units="in")
par(mfrow = c(2,2))
plot.ecdf(log10(all.randomized.medians.mean[,'dist']), main = 'Distance', ylab = 'Proportion of BFs', xlab = expression(paste('log'[10], '(Bayes Factor)')), xlim = c(-1,1))
plot.ecdf(log10(full_array_median[,"dist"]), col = 'tomato', add = TRUE)
# text(0.1, 0.2, paste('D = ', round(ks1$statistic, 3), '\n p-value <',round(ks1$p.value, 9)))
text(0.32, 0.6, paste('D = ', round(ks1$statistic, 3)))
text(0.32, 0.52, expression('p-value = 5.57x10'^-07))

plot.ecdf(log10(all.randomized.medians.mean[,'depth']), main = 'Depth', ylab = 'Proportion of BFs', xlab = expression(paste('log'[10], '(Bayes Factor)')), xlim = c(-1,1))
plot.ecdf(log10(full_array_median[,"depth"]), col = 'tomato', add = TRUE)
# text(0.1, 0.2,  paste('D = ', round(ks2$statistic, 4), '\n p-value <', round(ks2$p.value, 9)))
text(0.32, 0.6, paste('D = ', round(ks2$statistic, 3)))
text(0.32, 0.52, expression('p-value = 2.40x10'^-08))

plot.ecdf(log10(all.randomized.medians.mean[,'b_temp']), main = 'Bottom temperature', ylab = 'Proportion of BFs', xlab = expression(paste('log'[10], '(Bayes Factor)')), xlim = c(-1,1))
plot.ecdf(log10(full_array_median[,"b_temp"]), col = 'tomato', add = TRUE)
# text(0.1, 0.2,  paste('D = ', round(ks3$statistic, 4), '\n p-value <', round(ks3$p.value, 9)))
text(0.32, 0.6, paste('D = ', round(ks3$statistic, 3)))
text(0.32, 0.52, expression('p-value = 2.29x10'^-09))

plot.ecdf(log10(all.randomized.medians.mean[,'b_salin']), main = 'Bottom salinity', ylab = 'Proportion of BFs', xlab = expression(paste('log'[10], '(Bayes Factor)')), xlim = c(-1,1))
plot.ecdf(log10(full_array_median[,"b_salin"]), col = 'tomato', add = TRUE)
# text(0.1, 0.2,  paste('D = ', round(ks4$statistic, 3), '\n p-value <', round(ks4$p.value, 14)))
text(0.32, 0.6, paste('D = ', round(ks4$statistic, 3)))
text(0.32, 0.52, expression('p-value = 6.32x10'^-14))

legend("bottomright", legend = c('permuted', 'observed'), lty = c(1,1), col = c('black', 'tomato'))

dev.off()

# Exact binomial test using mean of 10 random permutations vs observed data for each environmental variable
# of 'successes', # of trials, prob of success from randomized data
binom.test(9, 1137, length(which(all.randomized.medians.mean[,'dist'] > 3))/1137, alternative = 'greater') # Distace from southern point
binom.test(6, 1137, length(which(all.randomized.medians.mean[,'depth'] > 3))/1137, alternative = 'greater') # Depth
binom.test(10, 1137, length(which(all.randomized.medians.mean[,'b_temp'] > 3))/1137, alternative = 'greater') # Bottom Temp
binom.test(6, 1137, length(which(all.randomized.medians.mean[,'b_salin'] > 3))/1137, alternative = 'greater') # Bottom Salinity

binom.test(9, 1137, length(which(all.randomized.medians.mean[,'dist'] > 3))/1137) # Distace from southern point
binom.test(6, 1137, length(which(all.randomized.medians.mean[,'depth'] > 3))/1137) # Depth
binom.test(10, 1137, length(which(all.randomized.medians.mean[,'b_temp'] > 3))/1137) # Bottom Temp
binom.test(6, 1137, length(which(all.randomized.medians.mean[,'b_salin'] > 3))/1137) # Bottom Salinity

binom.test(34, 1825, 9/1825, alternative = 'g')
x <- all.randomized.medians.mean[,'dist']
y <- dbinom(all.randomized.medians.mean[,'dist'], 1137, length(which(all.randomized.medians.mean[,'dist'] > 3))/1137)
plot(x, y)

#############################################
####           Randomizing RDA           ####
#############################################
# Bringing in adult genotype data
setwd("/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis")

library(ade4)
library(adegenet)
library(devtools)
library("hierfstat")
library(pegas)
library(fields)
library(vegan)
library(plotrix)

# Reading in SNP data file containing only the first SNP at each locus
adults <- read.structure("structure_input_Nov_11_2015.str",
                         n.ind = 241, n.loc = 1137, col.lab = 1, col.pop = 2, row.marknames = 1, 
                         onerowperind = FALSE)

which(adults@loc.n.all > 2) # which snps have more than 2 alelles?

adults.nonas <- scaleGen(adults, center = TRUE, scale = FALSE, NA.method = "mean") # filling in NAs with allele frequencies and centering
hist(adults.nonas)
sum(is.na(adults.nonas)) #0 All NAs have been replaced with means

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
dim(adults.nonas.232.1137)

# On Amphiprion
setwd("~/24_RDA_rand")
# write.table(adults.nonas.232.1137, "/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Redundancy analysis/adults.nonas.232.1137.txt", sep = '\t')
adults.nonas.232.1137 <- read.table("adults.nonas.232.1137.txt", header = TRUE)

# Bringing in environmental data
# setwd("/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation")
envi <- read.table("232envirowithdist.txt", header = TRUE)
envi.ordered <- envi[with(envi, order(PinskyID)), c('PinskyID', 'dist', 'depth', 'b_temp', 'b_salin')] # Subset to 4 environmental variables
as.character(envi.ordered[,1]) == rownames(adults.nonas.232.1137)

#### For loop to do emperical p-values ####
# Function to calculate p-values that's needed in the for loop below
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

rda.cans.dist <- vector()
locusenvi.cans <- vector()
for (j in 1:10000){
envi.ordered2 <- apply(envi.ordered[,-1],2, sample) # randomize environmental variables
envi.ordered.matrix <- scale(data.matrix(envi.ordered2))
rownames(envi.ordered.matrix) <- envi.ordered[,"PinskyID"]
# envi.ordered.matrix[,1] <- envi.ordered[,1]

envi.ordered.matrix <- as.data.frame(envi.ordered.matrix)
adults.rda2 <- rda(adults.nonas.232.1137 ~ dist + depth + b_temp + b_salin, envi.ordered.matrix, scale = FALSE)

spp.scr2 <- scores(adults.rda2, display = "species", scaling = 0, choices = c(1,2,3,4))

# Forester et al. (2016) used loci with scores +/- 3 SD from mean score for that axis to ID outlier loci. Only used first 3 axes.
mean.rda <- colMeans(spp.scr2)
sd.rda1 <- sd(spp.scr2[,"RDA1"])
sd.rda2 <- sd(spp.scr2[,"RDA2"])
sd.rda3 <- sd(spp.scr2[,"RDA3"])

rda1.hi <- mean.rda[1] + 3*sd.rda1
rda1.lo <- mean.rda[1] - 3*sd.rda1
which(spp.scr2[,"RDA1"] > rda1.hi)
which(spp.scr2[,"RDA1"] < rda1.lo)

rda2.hi <- mean.rda[2] + 3*sd.rda2
rda2.lo <- mean.rda[2] - 3*sd.rda2
which(spp.scr2[,"RDA2"] > rda2.hi)
which(spp.scr2[,"RDA2"] < rda2.lo)

rda3.hi <- mean.rda[3] + 3*sd.rda3
rda3.lo <- mean.rda[3] - 3*sd.rda3
which(spp.scr2[,"RDA3"] > rda3.hi)
which(spp.scr2[,"RDA3"] < rda3.lo)

rda.cans <- c((which(spp.scr2[,"RDA1"] > rda1.hi)),
              (which(spp.scr2[,"RDA1"] < rda1.lo)),
              (which(spp.scr2[,"RDA2"] > rda2.hi)),
              (which(spp.scr2[,"RDA2"] < rda2.lo)),
              (which(spp.scr2[,"RDA3"] > rda3.hi)),
              (which(spp.scr2[,"RDA3"] < rda3.lo)))
              
rda.cans.dist[j] <- length(unique(sort(rda.cans)))

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
lms.dist <- lapply(rda.cans, function(x) lm(adults.nonas.232.1137[,x] ~ envi.ordered.matrix[,"dist"]))
lms.depth <- lapply(rda.cans, function(x) lm(adults.nonas.232.1137[,x] ~ envi.ordered.matrix[,"depth"]))
lms.bt <- lapply(rda.cans, function(x) lm(adults.nonas.232.1137[,x] ~ envi.ordered.matrix[,"b_temp"]))
lms.bs <- lapply(rda.cans, function(x) lm(adults.nonas.232.1137[,x] ~ envi.ordered.matrix[,"b_salin"]))

all.lms <- c(lms.dist, lms.depth, lms.bt, lms.bs)

# for(i in rda.cans){
#   assign(paste0('dist',i), lm(adults.nonas.232.1137[,i] ~ envi.ordered.matrix[,"dist"]))
#   assign(paste0('depth',i), lm(adults.nonas.232.1137[,i] ~ envi.ordered.matrix[,"depth"]))
#   assign(paste0('btemp',i), lm(adults.nonas.232.1137[,i] ~ envi.ordered.matrix[,"b_temp"]))
#   assign(paste0('bsalin',i), lm(adults.nonas.232.1137[,i] ~ envi.ordered.matrix[,"b_salin"]))
# }


# Using the list of lm objects to put each into the lmp function to calculate p-values
pvalues <- vector()
for (i in 1:length(all.lms)){
  pvalues[i] <- lmp(all.lms[[i]])
}

rda.can.loci <- data.frame(names_lms,pvalues)
rda.can.loci2 <- rda.can.loci[which(rda.can.loci[,2] < 0.001),]

locusenvi.cans[j] <- nrow(rda.can.loci2)
}

#### Calculate emperical p-values ####
# Get R to tell me the number of simulations (# of loci > 3 & # of loci > 3 and p-value < 0.001) that are larger than observed values & calculate proportion (p value = (r+1)/((n+1)))
# Number of loci > 3
rda.cans.dist # should be of length 10,000
locusenvi.cans # should be of length 10,000

# save(rda.cans.dist, file = "/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Redundancy analysis/randomization/rda.cans.dist.RData")
# save(locusenvi.cans, file = "/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Redundancy analysis/randomization/locusenvi.cans.RData")

rda.cans.dist.pvalue <- (length(which(rda.cans.dist >= 23))+1)/(length(rda.cans.dist)+1) # p-value = 0.3617638

# Number of loci > 3 & p-value < 0.001
locusenvi.cans.pvalue <- (length(which(locusenvi.cans >= 5))+1)/(length(locusenvi.cans)+1) # p-value = 0.04429557

#### Plot emperical p-value distributions ####
png(file = "~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Redundancy analysis/randomization/loci3SDfrommean.png", width=6, height=10, res=300, units="in")

par(
  mfrow = c(2, 1), 
  mar=c(3.5, 4, 2, 0.7), # panel magin size in "line number" units
  mgp=c(2, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=13.5, # point size, which is the font size
  bg=NA
)

hist(rda.cans.dist, main = '', breaks = 25, xlab = 'Number of RDA outliers', ylab = 'Count', ylim = c(0,1200))
ablineclip(v = 23, col = 'tomato', y1=0, y2=1100) # observed number of loci
text(23,1100, pos = 4, paste('Observed # of loci = 23 \n p-value =', round(rda.cans.dist.pvalue,3)))
mtext(expression(bold('A')), side = 2, line = 3, padj = -13, las = 1)

hist(locusenvi.cans, breaks = 10, main = '', xlab = 'Number of RDA outliers with a strong environmental association', ylab = 'Count')
ablineclip(v = 5, col = 'tomato', y1=0, y2=3500) # oberved number of loci w/ strong enviornmental associations
text(5,3500, pos = 4, paste('Observed # of loci = 5 \n p-value =', round(locusenvi.cans.pvalue,3)))
mtext(expression(bold('B')), side = 2, line = 3, padj = -13, las = 1)

dev.off()

################################

# Long way of calculating p-values for RDA outliers
pvalues <- c(lmp(dist7), lmp(depth7), lmp(btemp7), lmp(bsalin7), lmp(dist168), lmp(depth168),  lmp(btemp168),  lmp(bsalin168),  lmp(dist248),  lmp(depth248),   lmp(btemp248),  lmp(bsalin248),
             lmp(dist290), lmp(depth290), lmp(btemp290), lmp(bsalin290), lmp(dist300), lmp(depth300),  lmp(btemp300),  lmp(bsalin300),  lmp(dist311),  lmp(depth311),   lmp(btemp311),  lmp(bsalin311),
             lmp(dist334), lmp(depth334), lmp(btemp334), lmp(bsalin334), lmp(dist353), lmp(depth353),  lmp(btemp353),  lmp(bsalin353),  lmp(dist360),  lmp(depth360),   lmp(btemp360),  lmp(bsalin360),
             lmp(dist392), lmp(depth392), lmp(btemp392), lmp(bsalin392), lmp(dist395), lmp(depth395),  lmp(btemp395),  lmp(bsalin395),  lmp(dist524),  lmp(depth524),   lmp(btemp524),  lmp(bsalin524),
             lmp(dist551), lmp(depth551), lmp(btemp551), lmp(bsalin551), lmp(dist565), lmp(depth565),  lmp(btemp565),  lmp(bsalin565),  lmp(dist578),  lmp(depth578),   lmp(btemp578),  lmp(bsalin578),
             lmp(dist696), lmp(depth696), lmp(btemp696), lmp(bsalin696), lmp(dist705), lmp(depth705),  lmp(btemp705),  lmp(bsalin705),  lmp(dist810),  lmp(depth810),   lmp(btemp810),  lmp(bsalin810),
             lmp(dist849), lmp(depth849), lmp(btemp849), lmp(bsalin849), lmp(dist939), lmp(depth939),  lmp(btemp939),  lmp(bsalin939),  lmp(dist985),  lmp(depth985),   lmp(btemp985),  lmp(bsalin985),
             lmp(dist997), lmp(depth997), lmp(btemp997), lmp(bsalin997), lmp(dist1070), lmp(depth1070),  lmp(btemp1070),  lmp(bsalin1070)
)

rda.can.loci <- data.frame(names_lms,pvalues)
rda.can.loci2 <- rda.can.loci[which(rda.can.loci[,2] < 0.001),]


#### Why do some environmental variables in the BayEnv null model have more spurious associations than others? ####
setwd("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation/randomization")
list <- list.files(pattern = "*.txt")
enviro_files <- lapply(list, read.table)
enviro_files_array <- array(unlist(enviro_files), dim = c(nrow(enviro_files[[1]]), ncol(enviro_files[[1]]), length(enviro_files))) # the order is Pop1, Pop2, Pop3, Pop4, Pop5; rows are: dist, depth, b_temp, b_salin

plot(enviro_files_array[2,,1]) # For example, plots distance of 1st permutation

par(mfrow = c(2,2))
hist(enviro_files_array[1,,], main = "Distance")
hist(enviro_files_array[2,,], main = "Depth")
hist(enviro_files_array[3,,], main = "Bottom Temp")
hist(enviro_files_array[4,,], main = "Bottom Salinity")

#### Let's look at the BF files ####
load(file = "/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation/randomization/randomized_array10.RData")
all.randomized.medians.mean <- apply(all.randomized.medians, c(1,2), mean)

# Check out skewness and kurtosis to see if there is a pattern
library(moments)

skewness(enviro_files_array[1,,])
kurtosis(enviro_files_array[1,,])
all.randomized.medians[which(all.randomized.medians.mean[,1]>3),1,]

skewness(enviro_files_array[2,,])
kurtosis(enviro_files_array[2,,])
all.randomized.medians[which(all.randomized.medians.mean[,2]>3),2,]

skewness(enviro_files_array[3,,])
kurtosis(enviro_files_array[3,,])
all.randomized.medians[which(all.randomized.medians.mean[,3]>3),3,]

skewness(enviro_files_array[4,,])
kurtosis(enviro_files_array[4,,])
all.randomized.medians[which(all.randomized.medians.mean[,4]>3),4,]

skewness(as.vector(enviro_files_array[1,,]))
skewness(as.vector(enviro_files_array[2,,]))
skewness(as.vector(enviro_files_array[3,,]))
skewness(as.vector(enviro_files_array[4,,]))

kurtosis(as.vector(enviro_files_array[1,,])) # pretty similar
kurtosis(as.vector(enviro_files_array[2,,]))
kurtosis(as.vector(enviro_files_array[3,,]))
kurtosis(as.vector(enviro_files_array[4,,]))

#### Let's look at allele frequencies ####
library(abind)

genes <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation/232forbayenv.txt')

# Split data into odd and even dataframes
even_indexes<-seq(2,2274,2)
odd_indexes<-seq(1,2273,2)

odds <- data.frame(genes[odd_indexes,]) # 1137 x 5
evens <- data.frame(genes[even_indexes,]) # 1137 x 5

totals <- odds + evens

odds.freqs <- as.matrix(odds/totals)
rownames(odds.freqs) <- 1:1137
evens.freqs <- as.matrix(evens/totals)
rownames(evens.freqs) <- 1:1137

# I need to fit a line to standardized enviromental variables & each allele freq, for each environmental variables and simulation
# Here is how to do it for a single simulation
# n <- 1137
# sim.dist10 <- lapply(1:n, function(x) lm(odds.freqs[x,] ~ enviro_files_array[1,,10]))
# sim.dist10.r2 <- sapply(sim.dist10, function(x) summary(x)$r.squared) # extract r2 
# sim.depth1 <- lapply(1:n, function(x) lm(odds.freqs[x,] ~ enviro_files_array[2,,1]))
# sim.depth1.r2 <- sapply(sim.depth1, function(x) summary(x)$r.squared) # extract r2 
# sim.btemp1 <- lapply(1:n, function(x) lm(odds.freqs[x,] ~ enviro_files_array[3,,1]))
# sim.btemp1.r2 <- sapply(sim.btemp1, function(x) summary(x)$r.squared) # extract r2 
# sim.bsalin1 <- lapply(1:n, function(x) lm(odds.freqs[x,] ~ enviro_files_array[4,,1]))
# sim.bsalin1.r2 <- sapply(sim.bsalin1, function(x) summary(x)$r.squared) # extract r2 

# for loop through each simulation
n <- 1137
for (m in 1:10){
  assign(paste0('dist',m,'r2'), sapply(lapply(1:n, function(x) lm(odds.freqs[x,] ~ enviro_files_array[1,,m])), function(x) summary(x)$r.squared)) # extract r2 and fit lm in one step
  assign(paste0('depth',m,'r2'), sapply(lapply(1:n, function(x) lm(odds.freqs[x,] ~ enviro_files_array[2,,m])), function(x) summary(x)$r.squared))
  assign(paste0('btemp',m,'r2'), sapply(lapply(1:n, function(x) lm(odds.freqs[x,] ~ enviro_files_array[3,,m])), function(x) summary(x)$r.squared))
  assign(paste0('bsalin',m,'r2'), sapply(lapply(1:n, function(x) lm(odds.freqs[x,] ~ enviro_files_array[4,,m])), function(x) summary(x)$r.squared))
}
  
sim1r2 <- cbind(dist1r2, depth1r2, btemp1r2, bsalin1r2)
sim2r2 <- cbind(dist2r2, depth2r2, btemp2r2, bsalin2r2)
sim3r2 <- cbind(dist3r2, depth3r2, btemp3r2, bsalin3r2)
sim4r2 <- cbind(dist4r2, depth4r2, btemp4r2, bsalin4r2)
sim5r2 <- cbind(dist5r2, depth5r2, btemp5r2, bsalin5r2)
sim6r2 <- cbind(dist6r2, depth6r2, btemp6r2, bsalin6r2)
sim7r2 <- cbind(dist7r2, depth7r2, btemp7r2, bsalin7r2)
sim8r2 <- cbind(dist8r2, depth8r2, btemp8r2, bsalin8r2)
sim9r2 <- cbind(dist9r2, depth9r2, btemp9r2, bsalin9r2)
sim10r2 <- cbind(dist10r2, depth10r2, btemp10r2, bsalin10r2)

simr2.array <- abind(sim1r2, sim2r2, sim3r2, sim4r2, sim5r2, sim6r2, sim7r2, sim8r2, sim9r2, sim10r2, along = 3)
save(simr2.array, file = "/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation/randomization/locus.enviro.rsquared.RData")
load(file = "/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation/randomization/locus.enviro.rsquared.RData")

# R2 mean for each environmental variable & simulation
simr2.means <- matrix(ncol = 4, nrow = 10)

for (t in 1:10){
  simr2.means[t,] <-colMeans(simr2.array[,,t])
}

par(mfrow = c(1,4))
boxplot(simr2.means[,1], main = 'Distance', ylim = c(0.2, 0.3), ylab = 'Mean R2')
points(simr2.means.envi[1])
points(mean(dist.obs.r2), col = 'tomato', pch = 19) # observed mean R2
boxplot(simr2.means[,2], main = 'Depth', ylim = c(0.2, 0.3))
points(simr2.means.envi[2])
points(mean(depth.obs.r2), col = 'tomato', pch = 19) # observed mean R2
boxplot(simr2.means[,3], main = 'Bottom Temp', ylim = c(0.2, 0.3))
points(simr2.means.envi[3])
points(mean(btemp.obs.r2), col = 'tomato', pch = 19) # observed mean R2
boxplot(simr2.means[,4], main = 'Bottom Salin', ylim = c(0.2, 0.3))
points(simr2.means.envi[4])
points(mean(bsalin.obs.r2), col = 'tomato', pch = 19) # observed mean R2

sd(simr2.means[,1])
sd(simr2.means[,2])
sd(simr2.means[,3])
sd(simr2.means[,4])

# R2 averaged through/across simulations for each locus/environmental variable
simr2.sim.means <- apply(simr2.array, c(1,2), mean)

par(mfrow = c(1,4))
boxplot(simr2.sim.means[,1], main = 'Distance', ylim = c(0,0.5), ylab = expression(paste('mean simulated R'^2, ' (locus basis)')))
points(mean(dist.obs.r2), col = 'tomato', pch = 19) # observed mean R2
boxplot(simr2.sim.means[,2], main = 'Depth', ylim = c(0,0.5))
points(mean(depth.obs.r2), col = 'tomato', pch = 19) # observed mean R2
boxplot(simr2.sim.means[,3], main = 'Bottom Temp', ylim = c(0,0.5))
points(mean(btemp.obs.r2), col = 'tomato', pch = 19) # observed mean R2
boxplot(simr2.sim.means[,4], main = 'Bottom Salin', ylim = c(0,0.5))
points(mean(bsalin.obs.r2), col = 'tomato', pch = 19) # observed mean R2

sd(simr2.sim.means[,1])
sd(simr2.sim.means[,2])
sd(simr2.sim.means[,3])
sd(simr2.sim.means[,4])

# Plotting SD of 1137 R2 values that have been averaged through the 10 simulations against # of unique values
plot(c(sd(simr2.sim.means[,1]), sd(simr2.sim.means[,2]), sd(simr2.sim.means[,3]), sd(simr2.sim.means[,4])) ~ c(85, 47, 69, 60), xlab = 'Number of unique values', ylab = expression(paste('SD of simulated R'^2)))

# Mean R2 within an environmental variable & across simulations
simr2.means.envi <- apply(simr2.array, 2, mean)
simr2.sd <- apply(simr2.array, 2, sd) # this must be sd across 1137*10 values, instead of 10

# Plot/R2 BFs against R2
plot(log10(all.randomized.medians[,4,]) ~ log10(simr2.array[,4,]), xlab = 'log10(randomized R2)', ylab = 'log10(randomized BF)')
plot(log10(all.randomized.medians[,1,1]) ~ log10(simr2.array[,1,1]), xlab = 'randomized R2', ylab = 'randomized BF')

# Fit linear regression to untransformed BFs and R2
for (m in 1:10){
  assign(paste0('dist',m,'bf'), lm(log10(simr2.array[,1,m]) ~ log10(all.randomized.medians[,'dist',m]))) # fit lm for each bf/r2 pair for distance
  assign(paste0('depth',m,'bf'), lm(log10(simr2.array[,1,m]) ~ log10(all.randomized.medians[,'depth',m])))
  assign(paste0('btemp',m,'bf'), lm(log10(simr2.array[,1,m]) ~ log10(all.randomized.medians[,'b_temp',m])))
  assign(paste0('bsalin',m,'bf'), lm(log10(simr2.array[,1,m]) ~ log10(all.randomized.medians[,'b_salin',m])))
}

# Plot dist BFs vs R2
par(mfrow = c(2,2))
plot(log10(simr2.array[,1,]) ~ log10(all.randomized.medians[,'dist',]), xlab = 'log10(randomized BF)', ylab = 'log10(randomized R2)', main = 'Distance')
abline(dist1bf, col = 'gray50')
abline(dist2bf, col = 'gray50')
abline(dist3bf, col = 'gray50')
abline(dist4bf, col = 'gray50')
abline(dist5bf, col = 'gray50')
abline(dist6bf, col = 'gray50')
abline(dist7bf, col = 'gray50')
abline(dist8bf, col = 'gray50')
abline(dist9bf, col = 'gray50')
abline(dist10bf, col = 'gray50')

# Plot depth BFs vs R2
plot(log10(simr2.array[,2,]) ~ log10(all.randomized.medians[,'depth',]), xlab = 'log10(randomized BF)', ylab = 'lgo10(randomized R2)', main = 'Depth')
abline(depth1bf, col = 'gray50')
abline(depth2bf, col = 'gray50')
abline(depth3bf, col = 'gray50')
abline(depth4bf, col = 'gray50')
abline(depth5bf, col = 'gray50')
abline(depth6bf, col = 'gray50')
abline(depth7bf, col = 'gray50')
abline(depth8bf, col = 'gray50')
abline(depth9bf, col = 'gray50')
abline(depth10bf, col = 'gray50')

# Plot btemp BFs vs R2
plot(log10(simr2.array[,3,]) ~ log10(all.randomized.medians[,'b_temp',]), xlab = 'log10(randomized BF)', ylab = 'log10(randomized R2)', main = 'Bottom Temp')
abline(btemp1bf, col = 'gray50')
abline(btemp2bf, col = 'gray50')
abline(btemp3bf, col = 'gray50')
abline(btemp4bf, col = 'gray50')
abline(btemp5bf, col = 'gray50')
abline(btemp6bf, col = 'gray50')
abline(btemp7bf, col = 'gray50')
abline(btemp8bf, col = 'gray50')
abline(btemp9bf, col = 'gray50')
abline(btemp10bf, col = 'gray50')

# Plot btemp BFs vs R2
plot(log10(simr2.array[,4,]) ~ log10(all.randomized.medians[,'b_salin',]), xlab = 'log10(randomized BF)', ylab = 'log10(randomized R2)', main = 'Bottom Salin')
abline(bsalin1bf, col = 'gray50')
abline(bsalin2bf, col = 'gray50')
abline(bsalin3bf, col = 'gray50')
abline(bsalin4bf, col = 'gray50')
abline(bsalin5bf, col = 'gray50')
abline(bsalin6bf, col = 'gray50')
abline(bsalin7bf, col = 'gray50')
abline(bsalin8bf, col = 'gray50')
abline(bsalin9bf, col = 'gray50')
abline(bsalin10bf, col = 'gray50')

# Plot mean R2 through simulations against mean BF through simulations
plot(simr2.sim.means[,1] ~ all.randomized.medians.mean[,1], xlab = 'Mean BF', ylab = 'Mean R2', main = 'Distance')
lm1 <- lm(simr2.sim.means[,1] ~ all.randomized.medians.mean[,1])
abline(lm1, col = 'tomato')
plot(simr2.sim.means[,2] ~ all.randomized.medians.mean[,2], xlab = 'Mean BF', ylab = 'Mean R2', main = 'Depth')
lm2 <- lm(simr2.sim.means[,2] ~ all.randomized.medians.mean[,2])
abline(lm2, col = 'tomato')
plot(simr2.sim.means[,3] ~ all.randomized.medians.mean[,3], xlab = 'Mean BF', ylab = 'Mean R2', main = 'Bottom Temp')
lm3 <- lm(simr2.sim.means[,3] ~ all.randomized.medians.mean[,3])
abline(lm3, col = 'tomato')
plot(simr2.sim.means[,4] ~ all.randomized.medians.mean[,4], xlab = 'Mean BF', ylab = 'Mean R2', main = 'Bottom Salinity')
lm4 <- lm(simr2.sim.means[,4] ~ all.randomized.medians.mean[,4])
abline(lm4, col = 'tomato')

# for loop through each simulation for each environmental variable (40 comparisons) to get R2 between BF and R2
distbf <- vector()
depthbf <- vector()
btempbf <- vector()
bsalinbf <- vector()

n <- 1137
for (m in 1:10){
  distbf[m] <-  summary(lm(simr2.array[,1,m] ~ all.randomized.medians[,"dist",m]))$r.squared # extract r2 and fit lm in one step
  depthbf[m] <-  summary(lm(simr2.array[,2,m] ~ all.randomized.medians[,"depth",m]))$r.squared
  btempbf[m] <-  summary(lm(simr2.array[,3,m] ~ all.randomized.medians[,"b_temp",m]))$r.squared
  bsalinbf[m] <-  summary(lm(simr2.array[,4,m] ~ all.randomized.medians[,"b_salin",m]))$r.squared
}

# Observed R2 for locus-environmental associations? 
# Read in observed standardized envirofile
envir.obs <- read.table("/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation/232stan4envirofile.txt", header = FALSE)

# Calculate r2 for each observed environmental-allele freq pair (1137 x 4)
n <- 1137

assign(paste0('dist.obs.r2'), sapply(lapply(1:n, function(x) lm(odds.freqs[x,] ~ unlist(envir.obs[1,]))), function(x) summary(x)$r.squared)) # extract r2 and fit lm in one step
assign(paste0('depth.obs.r2'), sapply(lapply(1:n, function(x) lm(odds.freqs[x,] ~ unlist(envir.obs[2,]))), function(x) summary(x)$r.squared))
assign(paste0('btemp.obs.r2'), sapply(lapply(1:n, function(x) lm(odds.freqs[x,] ~ unlist(envir.obs[3,]))), function(x) summary(x)$r.squared))
assign(paste0('bsalin.obs.r2'), sapply(lapply(1:n, function(x) lm(odds.freqs[x,] ~ unlist(envir.obs[4,]))), function(x) summary(x)$r.squared))

mean(dist.obs.r2) # 0.2616155
mean(depth.obs.r2) # 0.2632938
mean(btemp.obs.r2) # 0.2419458
mean(bsalin.obs.r2) # 0.2624925

#### Plot difference between mean observed R2 and mean simulated R2 across all simulations ####
r2.diffs <- simr2.means.envi - (c(mean(dist.obs.r2), mean(depth.obs.r2), mean(btemp.obs.r2), mean(bsalin.obs.r2)))

plot(c(85, 47, 69, 60) ~ c(1, 5, 2, 4), xlab = 'Number of expected \nlocus-environmental associations', ylab = 'Number of unique values', main = 'Observed')
plot(abs(r2.diffs) ~ c(85, 47, 69, 60), xlab = 'Number of unique values', ylab = expression(paste('Abs diff between obs and mean simul. R'^2)))
plot(abs(r2.diffs) ~ c(1, 5, 2, 4), xlab = 'Number of expected \nlocus-environmental associations', ylab = expression(paste('Abs diff between obs and mean simul. R'^2)))
plot(simr2.sd ~ c(1, 5, 2, 4), xlab = 'Number of expected \nlocus-environmental associations', ylab = expression(paste('SD of simulated R'^2)))
plot(c(sd(simr2.means[,1]), sd(simr2.means[,2]), sd(simr2.means[,3]),sd(simr2.means[,4])) ~ c(85, 47, 69, 60), xlab = 'Number of unique values', ylab = expression(paste('SD of simulated R'^2)))

#####################################################################
#### Generate fake data with different numbers of unique numbers ####
#####################################################################
fake1 <- round(rnorm(120), 2)
length(unique(fake1)) # 110
fake1 <- append(fake1, sample(fake1, 112, replace = TRUE))
length(fake1)

fake2 <- round(rnorm(100),2)
length(unique(fake2)) #85
fake2 <- append(fake2, sample(fake2, 132, replace = TRUE))
length(fake2)

fake3 <- round(rnorm(80),2)
length(unique(fake3)) #74
fake3 <- append(fake3, sample(fake3, 152, replace = TRUE))
length(fake3)

fake4 <- round(rnorm(60),2)
length(unique(fake4)) #55
fake4 <- append(fake4, sample(fake4, 172, replace = TRUE))
length(fake4)

fake5 <- round(rnorm(40),2)
length(unique(fake5)) #40
fake5 <- append(fake5, sample(fake5, 192, replace = TRUE))
length(fake5)

fake6 <- round(rnorm(30),2)
length(unique(fake6)) #27
fake6 <- append(fake6, sample(fake6, 202, replace = TRUE))
length(fake6)

fakes <- cbind(fake1, fake2, fake3, fake4, fake5, fake6) #232 x 6
write.table (fakes, "~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation/randomization/Fake Data/232fakedata.txt", sep ="\t", col.names = FALSE, row.names = FALSE)

# For loop to randomize, standardize and write 10 fake matrices #
for (i in 1:10){
  
  # Randomize environmental variables
  fakes.rand <- apply(fakes,2, sample)
  
  # Aggregate the individuals into populations and standardize
  by <- list(c(rep(1,40), rep(2,54), rep(3,41), rep(4,60), rep(5,37))) # generate vector of bayenv pop numbers
  pop.avg <- aggregate(fakes.rand, by = by, FUN = mean)
  stan.pop.avg <- scale(pop.avg[,2:7])
  rownames(stan.pop.avg) <- c("Pop1", "Pop2", "Pop3", "Pop4", "Pop5")
  t.stan.pop.avg <- t(stan.pop.avg) # the order is Pop1, Pop2, Pop3, Pop4, Pop5; rows are: 110 unique values, 85 unique values, 74 unique values, 55 unique values, 40 unique values, 27 unique values
  
  # Exporting the randomized envirofiles to my computer in a randomization folder
  write.table (t.stan.pop.avg, paste("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation/randomization/Fake Data/232stanfake", i,".txt", sep = ""), sep ="\t", col.names = FALSE, row.names = FALSE)
  
}

# Read in fake data and genetic data
library(abind)

setwd("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation/randomization/Fake Data")
list <- list.files(pattern = "232stan*")
fake_files <- lapply(list, read.table)
fake_files_array <- array(unlist(fake_files), dim = c(nrow(fake_files[[1]]), ncol(fake_files[[1]]), length(fake_files))) # the columns are Pop1, Pop2, Pop3, Pop4, Pop5; rows are different numbers of unique values

genes <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation/232forbayenv.txt')

# Split data into odd and even dataframes
even_indexes<-seq(2,2274,2)
odd_indexes<-seq(1,2273,2)

odds <- data.frame(genes[odd_indexes,]) # 1137 x 5
evens <- data.frame(genes[even_indexes,]) # 1137 x 5

totals <- odds + evens

odds.freqs <- as.matrix(odds/totals) # columns are BayEnv populations, rows are loci
rownames(odds.freqs) <- 1:1137
evens.freqs <- as.matrix(evens/totals)
rownames(evens.freqs) <- 1:1137

# I need to fit a line for each standardized fake dataset & each allele freq, for each environmental variables and simulation (1137*6*10)
# for loop through each simulation
n <- 1137
for (m in 1:10){
  assign(paste0('fake1.',m,'.r2'), sapply(lapply(1:n, function(x) lm(odds.freqs[x,] ~ fake_files_array[1,,m])), function(x) summary(x)$r.squared)) # extract r2 and fit lm in one step
  assign(paste0('fake2.',m,'.r2'), sapply(lapply(1:n, function(x) lm(odds.freqs[x,] ~ fake_files_array[2,,m])), function(x) summary(x)$r.squared))
  assign(paste0('fake3.',m,'.r2'), sapply(lapply(1:n, function(x) lm(odds.freqs[x,] ~ fake_files_array[3,,m])), function(x) summary(x)$r.squared))
  assign(paste0('fake4.',m,'.r2'), sapply(lapply(1:n, function(x) lm(odds.freqs[x,] ~ fake_files_array[4,,m])), function(x) summary(x)$r.squared))
  assign(paste0('fake5.',m,'.r2'), sapply(lapply(1:n, function(x) lm(odds.freqs[x,] ~ fake_files_array[5,,m])), function(x) summary(x)$r.squared))
  assign(paste0('fake6.',m,'.r2'), sapply(lapply(1:n, function(x) lm(odds.freqs[x,] ~ fake_files_array[6,,m])), function(x) summary(x)$r.squared))
}

fake1r2 <- cbind(fake1.1.r2, fake2.1.r2, fake3.1.r2, fake4.1.r2, fake5.1.r2, fake6.1.r2)
fake2r2 <- cbind(fake1.2.r2, fake2.2.r2, fake3.2.r2, fake4.2.r2, fake5.2.r2, fake6.2.r2)
fake3r2 <- cbind(fake1.3.r2, fake2.3.r2, fake3.3.r2, fake4.3.r2, fake5.3.r2, fake6.3.r2)
fake4r2 <- cbind(fake1.4.r2, fake2.4.r2, fake3.4.r2, fake4.4.r2, fake5.4.r2, fake6.4.r2)
fake5r2 <- cbind(fake1.5.r2, fake2.5.r2, fake3.5.r2, fake4.5.r2, fake5.5.r2, fake6.5.r2)
fake6r2 <- cbind(fake1.6.r2, fake2.6.r2, fake3.6.r2, fake4.6.r2, fake5.6.r2, fake6.6.r2)
fake7r2 <- cbind(fake1.7.r2, fake2.7.r2, fake3.7.r2, fake4.7.r2, fake5.7.r2, fake6.7.r2)
fake8r2 <- cbind(fake1.8.r2, fake2.8.r2, fake3.8.r2, fake4.8.r2, fake5.8.r2, fake6.8.r2)
fake9r2 <- cbind(fake1.9.r2, fake2.9.r2, fake3.9.r2, fake4.9.r2, fake5.9.r2, fake6.9.r2)
fake10r2 <- cbind(fake1.10.r2, fake2.10.r2, fake3.10.r2, fake4.10.r2, fake5.10.r2, fake6.10.r2)

fake.array <- abind(fake1r2, fake2r2, fake3r2, fake4r2, fake5r2, fake6r2, fake7r2, fake8r2, fake9r2, fake10r2, along = 3) #1137 x 6 x 10
save(fake.array, file = "/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation/randomization/Fake Data/locus.fake.rsquared.RData")
load(file = "/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation/randomization/Fake Data/locus.fake.rsquared.RData")

# R2 mean for each environmental variable & simulation (6*10)
faker2.means <- matrix(ncol = 6, nrow = 10)

for (t in 1:10){
  faker2.means[t,] <-colMeans(fake.array[,,t])
}

# Plot boxplots of mean r2 for fake data - hard to actually see...I wonder if outliers are messing things up
boxplot(faker2.means[,1], ylim = c(0.21, 0.29), main = '110 uniques', ylab = expression(paste('Mean R'^2)))
boxplot(faker2.means[,2], ylim = c(0.21, 0.29), main = '85 uniques')
boxplot(faker2.means[,3], ylim = c(0.21, 0.29), main = '74 uniques')
boxplot(faker2.means[,4], ylim = c(0.21, 0.29), main = '55 uniques')
boxplot(faker2.means[,5], ylim = c(0.21, 0.29), main = '40 uniques')
boxplot(faker2.means[,6], ylim = c(0.21, 0.29), main = '27 uniques')

# How about plots of SD? Removing 'outliers' based on boxplot
sd(faker2.means[-5,1])
sd(faker2.means[-10,2])
sd(faker2.means[,3])
sd(faker2.means[-5,4])
sd(faker2.means[,5])
sd(faker2.means[,6 ])

plot(c(sd(faker2.means[,1]), sd(faker2.means[,2]), sd(faker2.means[,3]), sd(faker2.means[,4]), sd(faker2.means[,5]), sd(faker2.means[,6])) ~ c(110, 85, 74, 55, 40, 27), xlab = 'Number of unique values', ylab = expression(paste('SD of simulated R'^2)))

faker2.sd <- apply(fake.array, 2, sd) # this must be sd across 1137*10 values, instead of 10

##############################################################################
#### Try again to make fake data using environmental variable means & SDs ####
##############################################################################
setwd("/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation")

# Read in environmental data
nvs_locs <- read.table("232envirowithdist.txt", header=TRUE) #all 241-9=232 fish divided into 5 populations

# Subset to 4 environmental variables
envir <- nvs_locs[,c('bayenv_pop', 'dist', 'depth', 'b_temp', 'b_salin')]
mean(envir$dist)
mean(envir$depth)
mean(envir$b_temp)
mean(envir$b_salin)
sd(envir$dist)
sd(envir$depth)
sd(envir$b_temp)
sd(envir$b_salin)

#### Generate fake data with different numbers of unique numbers ####
fake1 <- round(rnorm(85, mean(envir$dist), sd(envir$dist)), 2)
length(unique(fake1)) # 85
fake1 <- append(fake1, sample(fake1, 147, replace = TRUE))
length(fake1)

fake2 <- round(rnorm(47, mean(envir$depth), sd(envir$depth)),2)
length(unique(fake2)) #47
fake2 <- append(fake2, sample(fake2, 185, replace = TRUE))
length(fake2)

fake3 <- round(rnorm(69, mean(envir$b_temp), sd(envir$b_temp)),2)
length(unique(fake3)) #69
fake3 <- append(fake3, sample(fake3, 163, replace = TRUE))
length(fake3)

fake4 <- round(rnorm(65, mean(envir$b_salin), sd(envir$b_salin)),2)
length(unique(fake4)) #60
fake4 <- append(fake4, sample(fake4, 167, replace = TRUE))
length(fake4)

fakes <- cbind(fake1, fake2, fake3, fake4) #232 x 4
write.table (fakes, "~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation/randomization/fake2/232fakedata.txt", sep ="\t", col.names = FALSE, row.names = FALSE)

# For loop to randomize, standardize and write 10 fake matrices #
for (i in 1:10){
  
  # Randomize environmental variables
  fakes.rand <- apply(fakes,2, sample)
  
  # Aggregate the individuals into populations and standardize
  by <- list(c(rep(1,40), rep(2,54), rep(3,41), rep(4,60), rep(5,37))) # generate vector of bayenv pop numbers
  pop.avg <- aggregate(fakes.rand, by = by, FUN = mean)
  stan.pop.avg <- scale(pop.avg[,2:5])
  rownames(stan.pop.avg) <- c("Pop1", "Pop2", "Pop3", "Pop4", "Pop5")
  t.stan.pop.avg <- t(stan.pop.avg) # the order is Pop1, Pop2, Pop3, Pop4, Pop5; rows are: 85 unique values, 47 unique values, 69 unique values, 60 unique values
  
  # Exporting the randomized envirofiles to my computer in a randomization folder
  write.table (t.stan.pop.avg, paste("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation/randomization/fake2/232stanfake", i,".txt", sep = ""), sep ="\t", col.names = FALSE, row.names = FALSE)
  
}

# Read in fake data and genetic data
library(abind)

setwd("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation/randomization/fake2")
list <- list.files(pattern = "232stan*")
fake_files <- lapply(list, read.table)
fake_files_array <- array(unlist(fake_files), dim = c(nrow(fake_files[[1]]), ncol(fake_files[[1]]), length(fake_files))) # the columns are Pop1, Pop2, Pop3, Pop4, Pop5; rows are different numbers of unique values; 4 x 5 x 10

genes <- read.table('~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation/232forbayenv.txt')

# Split data into odd and even dataframes
even_indexes<-seq(2,2274,2)
odd_indexes<-seq(1,2273,2)

odds <- data.frame(genes[odd_indexes,]) # 1137 x 5
evens <- data.frame(genes[even_indexes,]) # 1137 x 5

totals <- odds + evens

odds.freqs <- as.matrix(odds/totals) # columns are BayEnv populations, rows are loci
rownames(odds.freqs) <- 1:1137
evens.freqs <- as.matrix(evens/totals)
rownames(evens.freqs) <- 1:1137

# I need to fit a line for each standardized fake dataset & each allele freq, for each environmental variables and simulation (1137*6*10)
# for loop through each simulation
n <- 1137
for (m in 1:10){
  assign(paste0('fake1.',m,'.r2'), sapply(lapply(1:n, function(x) lm(odds.freqs[x,] ~ fake_files_array[1,,m])), function(x) summary(x)$r.squared)) # extract r2 and fit lm in one step
  assign(paste0('fake2.',m,'.r2'), sapply(lapply(1:n, function(x) lm(odds.freqs[x,] ~ fake_files_array[2,,m])), function(x) summary(x)$r.squared))
  assign(paste0('fake3.',m,'.r2'), sapply(lapply(1:n, function(x) lm(odds.freqs[x,] ~ fake_files_array[3,,m])), function(x) summary(x)$r.squared))
  assign(paste0('fake4.',m,'.r2'), sapply(lapply(1:n, function(x) lm(odds.freqs[x,] ~ fake_files_array[4,,m])), function(x) summary(x)$r.squared))
}

fake1r2 <- cbind(fake1.1.r2, fake2.1.r2, fake3.1.r2, fake4.1.r2)
fake2r2 <- cbind(fake1.2.r2, fake2.2.r2, fake3.2.r2, fake4.2.r2)
fake3r2 <- cbind(fake1.3.r2, fake2.3.r2, fake3.3.r2, fake4.3.r2)
fake4r2 <- cbind(fake1.4.r2, fake2.4.r2, fake3.4.r2, fake4.4.r2)
fake5r2 <- cbind(fake1.5.r2, fake2.5.r2, fake3.5.r2, fake4.5.r2)
fake6r2 <- cbind(fake1.6.r2, fake2.6.r2, fake3.6.r2, fake4.6.r2)
fake7r2 <- cbind(fake1.7.r2, fake2.7.r2, fake3.7.r2, fake4.7.r2)
fake8r2 <- cbind(fake1.8.r2, fake2.8.r2, fake3.8.r2, fake4.8.r2)
fake9r2 <- cbind(fake1.9.r2, fake2.9.r2, fake3.9.r2, fake4.9.r2)
fake10r2 <- cbind(fake1.10.r2, fake2.10.r2, fake3.10.r2, fake4.10.r2)

fake.array <- abind(fake1r2, fake2r2, fake3r2, fake4r2, fake5r2, fake6r2, fake7r2, fake8r2, fake9r2, fake10r2, along = 3) #1137 x 4 x 10
save(fake.array, file = "/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation/randomization/fake2/locus.fake.rsquared.RData")
load(file = "/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local adaptation/randomization/fake2/locus.fake.rsquared.RData")

# R2 averaged through/across simulations for each locus/environmental variable
simr2.sim.means <- apply(fake.array, c(1,2), mean)

par(mfrow = c(1,4))
boxplot(simr2.sim.means[,1], main = 'Distance', ylim = c(0,0.6), ylab = expression(paste('mean simulated R'^2, ' (locus basis)')))
points(mean(dist.obs.r2), col = 'tomato', pch = 19) # observed mean R2
boxplot(simr2.sim.means[,2], main = 'Depth', ylim = c(0,0.6))
points(mean(depth.obs.r2), col = 'tomato', pch = 19) # observed mean R2
boxplot(simr2.sim.means[,3], main = 'Bottom Temp', ylim = c(0,0.6))
points(mean(btemp.obs.r2), col = 'tomato', pch = 19) # observed mean R2
boxplot(simr2.sim.means[,4], main = 'Bottom Salin', ylim = c(0,0.6))
points(mean(bsalin.obs.r2), col = 'tomato', pch = 19) # observed mean R2

sd(simr2.sim.means[,1])
sd(simr2.sim.means[,2])
sd(simr2.sim.means[,3])
sd(simr2.sim.means[,4])

# Plotting SD of 1137 R2 values that have been averaged through the 10 simulations against # of unique values
plot(c(sd(simr2.sim.means[,1]), sd(simr2.sim.means[,2]), sd(simr2.sim.means[,3]), sd(simr2.sim.means[,4])) ~ c(85, 47, 69, 60), xlab = 'Number of unique values', ylab = expression(paste('SD of simulated R'^2)))



# R2 mean for each environmental variable & simulation (4*10)
faker2.means <- matrix(ncol = 4, nrow = 10)

for (t in 1:10){
  faker2.means[t,] <-colMeans(fake.array[,,t])
}

# Plot boxplots of mean r2 for fake data
par(mfrow = c(1,4))
boxplot(faker2.means[,1], ylim = c(0.20, 0.29), main = '85 uniques', ylab = expression(paste('Mean R'^2)))
boxplot(faker2.means[,2], ylim = c(0.20, 0.29), main = '47 uniques')
boxplot(faker2.means[,3], ylim = c(0.20, 0.29), main = '69 uniques')
boxplot(faker2.means[,4], ylim = c(0.20, 0.29), main = '60 uniques')

# How about plots of SD? Removing 'outliers' based on boxplot
sd(faker2.means[,1])
sd(faker2.means[,2])
sd(faker2.means[,3])
sd(faker2.means[,4])

plot(c(sd(faker2.means[,1]), sd(faker2.means[,2]), sd(faker2.means[,3]), sd(faker2.means[,4])) ~ c(85, 47, 69, 60), xlab = 'Number of unique values', ylab = expression(paste('SD of simulated R'^2)))

