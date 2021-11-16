# Functions

## Functions .... 
#----------------------------

get_pos <- function(q, a1, a2){

  # get pos for chromosome 1
  ans <- c(a1, a2)
  pos <-  q[[ans[1] ]]$pos
  for(ii in 2:length(ans)){
    pos <- c(pos, q[[ ans[ii] ]]$pos)
    pos <- unique(pos)
    pos <- sort(pos)
  }

 return(pos)
}

all_together <- function(animals=NULL, CHRM=NULL){

  df <- list()
 q <- list()
 originalSize <- NULL
 filteredSize <- NULL
 

 for(ii in animals)
    {
      cat(" Processing animal ... ", ii, "\n")
      df[[ii]] <- fread(file=paste0(DataDIR, "HS",ii, ".csv"),  header=TRUE)
      df[[ii]] <- as.data.frame(df[[ii]])
      df[[ii]] <- subset(df[[ii]], Chromosome==CHRM)
      originalSize <- c(nrow(df[[ii]]), originalSize)
      cat("Mean total read count for replicate is ", median(df[[ii]]$"Context coverage", na.rm=TRUE), "\n")
      # QC 
      # Remove loci outside of bounds
      cat(" Removing lower coverage SNPs ... \n")
      indx <- with(df[[ii]], which(`Context coverage` >= 5 & `Context coverage` < 300 ))
      cat(" % CpG sites kept",  length(indx)/nrow(df[[ii]]), "\n")
      df[[ii]] <- df[[ii]][indx,]

      # Remove loci where strand coverage neq context coverage
      #indx <- which(df[[ii]]$`Context coverage` != df[[ii]]$`Strand coverage`)
      #cat(" % CpG sites removed",  length(indx)/nrow(df[[ii]] ), "\n")
      #df[[ii]] <- df[[ii]][-indx,]

      filteredSize <- c(nrow(df[[ii]]), filteredSize)

      # Form data set in preparation for analysis
      q[[ii]] <- data.frame("chr"=df[[ii]]$Chromosome, "pos"=df[[ii]]$Region,
                                 "N"=df[[ii]]$`Context coverage`,
                                 "X" = df[[ii]]$`Methylated coverage`)

      # turn q$pos into a numeric position value
      q[[ii]]$pos <- gsub("complement\\(", "", q[[ii]]$pos)
      q[[ii]]$pos <- as.numeric(gsub("\\)", "", q[[ii]]$pos))
    }

  # create object to hold all data
  #pos <- get_pos(q=q, a1=animals, a2=animalsHE2)
  pos <- get_pos(q=q, a1=animals, a2=animals)  #tricky way of using get_pos just for one set of replicates

  finaldf <- as.data.frame(matrix(data=NA, nrow=length(pos), ncol=(2*length(animals) + 1 ) ) )
  names(finaldf) <- c("Pos", paste(c("N","X"), rep(animals,each=2) , sep=""))
  finaldf[, "Pos"] <- pos
  for (ii in animals){
    indx <- which(pos %in% q[[ii]]$pos)    
    finaldf[indx, paste0("N",ii) ] <- q[[ii]]$N
    finaldf[indx, paste0("X",ii) ] <- q[[ii]]$X

  }


  return (finaldf)

}

