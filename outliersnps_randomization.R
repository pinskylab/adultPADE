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




# Randomize everything but NA's??
df1 <- data.frame(A=c(1,1,2,2), B=c(NA,NA,1,3), C=c(4,4,2,2), D = c(1,1,NA,NA), E= c(1,3,3,3))
df4 <- apply(df1, 2, function(x){sample(x[!is.na(x)])})
df2 <- apply(df1, 2, sample)

df3 <- df1[is.na(df1)]

df1[df1 > 0] <- 0
  
x <- c(NA, 3, NA, 5)
x[!is.na(x)]


