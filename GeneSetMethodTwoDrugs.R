#####################################################
#     Joelle Strom, last updated 02/17/2021         #
# Program to predict a synergy score for two        #
# drugs or toxins, begin by defining functional     #
# gene sets, then sort through expression files to  #
# find commonly expressed genes across sets between #
# the two treatments. Method derived from Hsu et al.#
# (2016). A simple gene set-based method accurately #
# predicts the synergy of drug pairs.               #
# DOI: 10.1186/s12918-016-0310-3                    #
#####################################################

##Set working directory, import library to read gene set datafile
setwd("D:/Documents/Applied Stats MS/Research")
library(qusage)
library(httr)
library(data.table)

#Import one expression dataset to extract list of genes
example <- read.csv("Acetaminophen-compoundsData.csv")
genes <- example$Gene[which(example$Dataset == "Open TG-GATEs Human" &
                            example$Dose == "Middle" &
                            example$Time == 24)]

#Open .gmt file from MSigDB containing hallmark gene sets
sets <- read.gmt("h.all.v7.2.symbols.gmt")
#Arrange genes from dataset into the hallmark gene sets
mysets <- list()
for(i in 1:length(sets)){
  temp <- intersect(sets[[i]],genes)
  mysets[[names(sets[i])]] <- temp
}

#Read in datafile containing list of compounds (see CreateCompoundList.R)
compounds <- read.table("CompoundsList.txt")
compounds$mark <- NULL
tox <- compounds[which(compounds$class_in_vitro=="Genotoxic" | compounds$class_in_vivo=="Genotoxic"),]
carc <- compounds[which(compounds$carcinogenicity=="Carcinogenic"),]

#Create a dataframe with fillers to store final co-gene/GS scores with compound pair labels
allcogene.GS <- data.frame(
  comp1 = rep('char',dim(compounds)[1]**2),
  comp2 = rep('char',dim(compounds)[1]**2),
  score = rep(0,dim(compounds)[1]**2),
  stringsAsFactors = FALSE
)

start.time <- proc.time()

i <- 1
for(i in 1:dim(compounds)[1]){
  #Retrieve expression data for a drug with ToxicoDB id number i via API
  id <- as.character(compounds$id[i])
  path1 <- paste("https://www.toxicodb.ca/api/v1/compounds/",id,"/analysis",sep="")
  drug1json <- GET(url=path1)
  drug1json <- content(drug1json)
  
  #Reformat into a dataframe with proper column attributes
  drug1 <- rbindlist(drug1json$data)
  drug1 <- transform(drug1, fdr = as.numeric(fdr), p_value = as.numeric(p_value))
  
  for(j in i:dim(compounds)[1]){
    id2 <- as.character(compounds$id[j])
    path2 <- paste("https://www.toxicodb.ca/api/v1/compounds/",id2,"/analysis",sep="")
    drug2json <- GET(url=path2)
    drug2json <- content(drug2json)
    drug2 <- rbindlist(drug2json$data)
    drug2 <- transform(drug2, fdr = as.numeric(fdr), p_value = as.numeric(p_value))
  
    #Create new dataframes for all significantly positively and negatively expressed genes for each drug
    # at middle dose and 24 hours
    drug1pos <- drug1[which(drug1$fold_change > 0 &
                            drug1$fdr <= 0.05 & 
                            drug1$dataset_name == "Open TG-GATEs Human" & 
                            drug1$dose == "Middle" & 
                            drug1$time == 24),]
    drug1neg <- drug1[which(drug1$fold_change < 0 &
                            drug1$fdr <= 0.05 & 
                            drug1$dataset_name == "Open TG-GATEs Human" & 
                            drug1$dose == "Middle" & 
                            drug1$time == 24),]
    drug2pos <- drug2[which(drug2$fold_change > 0 &
                            drug2$fdr <= 0.05 &
                            drug2$dataset_name == "Open TG-GATEs Human" &
                            drug2$dose == "Middle" &
                            drug2$time == 24),]
    drug2neg <- drug2[which(drug2$fold_change < 0 &
                            drug2$fdr <= 0.05 &
                            drug2$dataset_name == "Open TG-GATEs Human" &
                            drug2$dose == "Middle" &
                            drug2$time == 24),]
  
    #Arrange positively and negatively expressed genes from drug 1 into hallmark gene sets
    drug1setspos <- list()
    for(k in 1:length(mysets)){
      temp <- intersect(mysets[[k]],drug1pos$gene_name)
      drug1setspos[[names(mysets[k])]] <- temp
    }
  
    drug1setsneg <- list()
    for(k in 1:length(mysets)){
      temp <- intersect(mysets[[k]],drug1neg$gene_name)
      drug1setsneg[[names(mysets[k])]] <- temp
    }
  
    #Find the common positively and negatively expressed genes between drug 1 and 2 for each gene set
    commonpos <- list()
    for(k in 1:length(mysets)){
      temp <- intersect(drug1setspos[[k]],drug2pos$gene_name)
      commonpos[[names(mysets[k])]] <- temp
    }
  
    commonneg <- list()
    for(k in 1:length(mysets)){
      temp <- intersect(drug1setsneg[[k]],drug2neg$gene_name)
      commonneg[[names(mysets[k])]] <- temp
    }
  
    #Compute average number of commonly expressed genes for each gene set: the co-gene score
    cogenescore <- list()
    for(k in 1:length(mysets)){
      a <- length(commonpos[[k]])
      b <- length(commonneg[[k]])
      c <- length(mysets[[k]])
      cogenescore[[names(mysets[k])]] <- as.numeric((a+b)/c)
    }
  
    #Create new dataframes for all expression values at middle dose and 24 hours for each drug
    drug1abbr <- drug1[which(drug1$dataset_name == "Open TG-GATEs Human" & 
                             drug1$dose == "Middle" & 
                             drug1$time == 24),]
    drug2abbr <- drug2[which(drug2$dataset_name == "Open TG-GATEs Human" &
                             drug2$dose == "Middle" &
                             drug2$time == 24),]
  
    #Sum expression values across each gene set for each drug separately, then create TRUE item if the sign of
    # the summed expression value agrees across each drug for a given gene set
    setenrich <- rep(FALSE,50)
    for(k in 1:length(mysets)){
      a <- 0
      b <- 0
      for(l in 1:length(mysets[[k]])){
        a <- a + drug1abbr$fold_change[which(drug1abbr$gene_name==mysets[[k]][l])]
        b <- b + drug2abbr$fold_change[which(drug2abbr$gene_name==mysets[[k]][l])]
      }
      setenrich[k] <- (a>0 && b>0) || (a<0 && b<0)
    }
  
    #Average the number of commonly enriched gene sets to get co-GS score
    coGS <- mean(setenrich)
  
    #Average co-gene scores across all commonly enriched gene sets to get co-gene/GS score
    a <- 0
    for(k in 1:length(setenrich)){
      if(setenrich[k]){
      a <- a + cogenescore[[k]]
      }
    }
    cogene.GS <- a/sum(setenrich)
        
    allcogene.GS$comp1[dim(compounds)[1]*(i-1) + j] <- compounds$name[i]
    allcogene.GS$comp2[dim(compounds)[1]*(i-1) + j] <- compounds$name[j]
    allcogene.GS$score[dim(compounds)[1]*(i-1) + j] <- cogene.GS
  }
}
proc.time()-start.time