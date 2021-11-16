# Analysis of methylation data

# Data Origin:
# Old Data comes from /datasets/work/af-mlai-bio/work/Heat_Stress_Use-case/processed_data/HS_Methylation_Data3_Msp1_Restricted_Mapping/
# New Data comes from /datasets/work/af-mlai-bio/work/Heat_Stress_Use-case/HS_ML3_data

# Purpose:  to identify regions which have changed in methylation value 

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#
#  BiocManager::install("DSS")
#  BiocManager::install("bsseq")


# The low read counts and filtering is the biggest issue we have with 
# having so few sites in common. 


rm(list=ls())
DIR <- "../Results/"
library(DSS)
library(bsseq)
library(data.table)

source("./functions.R")

# Stability selection of fisher score tests
# stability_selection()



# form two  data sets, one for reps in HE1 and one for reps in HE2 
DataDIR <- "../New_Data/"
animals <- c(1:36, 40:45)
animalsHE1 <- c(seq(1,34,3), 40, 43)
animalsHE2 <- c(seq(2,35,3), 41, 44)
animalsHE3 <- c(seq(3,36,3), 42,45)


for(CHRM in 1:29){
 cat(" Processing chromosome ", CHRM, "\n")


finaldfHE1 <- all_together(animals=animalsHE1, CHRM=CHRM)
finaldfHE2 <- all_together(animals=animalsHE2, CHRM=CHRM)
finaldfHE3 <- all_together(animals=animalsHE3, CHRM=CHRM)


# base  on finaldfHE1 
Ncols <- paste0("N", animalsHE1)
Ncols <- which(names(finaldfHE1) %in% Ncols )
Xcols <- paste0("X", animalsHE1)
Xcols <- which(names(finaldfHE1) %in% Xcols )

df <- data.frame(Pos=finaldfHE1$Pos, Prop=methylationValues <- rowSums(finaldfHE1[, Xcols], na.rm=TRUE)/rowSums(finaldfHE1[, Ncols], na.rm=TRUE) )



#----------------------------------------------------------
# Find sites that are common to animalsHE1 and animalsHE2 and animalsHE3
#----------------------------------------------------------

# form indx of rows in finaldfHE1 that are also in finaldfHE2
indx <- which( finaldfHE2$Pos %in%  finaldfHE1$Pos )
finaldfHE2 <- finaldfHE2[indx,]
indx <- which( finaldfHE3$Pos %in%  finaldfHE1$Pos )
finaldfHE3 <- finaldfHE3[indx,]


indx <- which( finaldfHE1$Pos %in%  finaldfHE2$Pos )
finaldfHE1 <- finaldfHE1[indx,]
indx <- which( finaldfHE3$Pos %in%  finaldfHE2$Pos )
finaldfHE3 <- finaldfHE3[indx,]

indx <- which( finaldfHE1$Pos %in%  finaldfHE3$Pos )
finaldfHE1 <- finaldfHE1[indx,]
indx <- which( finaldfHE2$Pos %in%  finaldfHE3$Pos )
finaldfHE2 <- finaldfHE2[indx,]

print(dim(finaldfHE1))
print(dim(finaldfHE2))
print(dim(finaldfHE3))







#=====> write files to disk

DIR_tai <- "../Results/"
fn1 <-  paste0(DIR_tai, "filtered_meth_HE1_chrm", CHRM, ".dat")
fn2 <-  paste0(DIR_tai, "filtered_meth_HE2_chrm", CHRM, ".dat")
fn3 <-  paste0(DIR_tai, "filtered_meth_HE3_chrm", CHRM, ".dat")

write.table(finaldfHE1, file=fn1, col.names=TRUE, row.names=FALSE, quote=FALSE)
write.table(finaldfHE2, file=fn2, col.names=TRUE, row.names=FALSE, quote=FALSE)
write.table(finaldfHE3, file=fn3, col.names=TRUE, row.names=FALSE, quote=FALSE)

}





