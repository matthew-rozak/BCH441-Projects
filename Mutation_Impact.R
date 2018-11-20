# MUTATION_IMPACT.R
#
# Purpose: Explore the impact of mutations on dna in resulting proteins
#
# Version: 0.1
#
# Date:    2018  11  18
# Author:  Matthew Rozak (matthew.rozak@mail.utoronto.ca)
#
# V 0.1    First code
#
#
#

#Use this package for the genetic code and DNA analysis
if (! require(Biostrings, quietly=TRUE)) {
  if (! exists("biocLite")) {
    source("https://bioconductor.org/biocLite.R")
  }
  biocLite("Biostrings")
  library(Biostrings)
}
#use this to load the TSV file
if (! require(readr, quietly = TRUE)) {
  install.packages("readr")
  library(readr)
}

#The following function has been inspired by functions in
#FND-Genetic_code.R and RPR-Genetic_code_optimality.R.
evalMut <- function(FA, N) {
  # Purpose: evaluate the distribution of silent, missense and nonsense
  # codon changes in cDNA read from FA for N random mutation trials.
  # Parameters:
  #     FA   chr      Filename of a FASTA formatted sequence file of cDNA
  #                     beginning with a start codon.
  #     N    integer  The number of point mutation trials to perform
  # Value:   list     List with the following elements:
  #                      FA    chr  the input file
  #                      N     num  same as the input parameter
  #                      nSilent    num  the number of silent mutations
  #                      nMissense  num  the number of missense mutations
  #                      nNonsense  num  the number of nonsense mutations
  if(TRUE){
    #counters for future iff statements
    nSilent <- 0
    nMissense <- 0
    nNonsense <- 0

    #get file name
    fileName <- FA

    #vector for selection of change nucleotide
    change <- c('A','T','G','C')

    #puts input file data into a varaible
    data <- readLines(FA,)

    #removes the header
    data <- data[-1]

    #shorten into a single string
    data <- paste(data,sep = "",collapse = "")

    #split into individula leters
    data <- unlist(strsplit(data,''))

    #copy for laters
    data1 <- data

    #shorten again
    data <- paste(data,sep = "",collapse = "")

    #make vectors for translation
    dataCodons <- as.character(codons(DNAString(data)))
    dataAA <- character(length = length(dataCodons))

    #loop N times for N mutations
    for(j in 1:N){
    #reset previous change
    mutdata <- data1
    #select a rnadom nucleotide to mutate
    muts <- sample(1:length(mutdata),1)
    #swap selected nucleotide for another nucleotide
    mutN <- sample(change[change != mutdata[muts]],1)
    mutdata[muts] <- mutN
    #shorten mutated data into a single string
    mutdata <- paste(mutdata,sep = "",collapse = "")
    #seperate into a sequence fo string with 3 setters each
    mutdataCodons <- as.character(codons(DNAString(mutdata)))
    mutdataAA <- character(length = length(dataCodons))
    for(i in 1:length(dataCodons)){
      #translate
      dataAA[i] <- GENETIC_CODE[dataCodons[i]]
      mutdataAA[i] <- GENETIC_CODE[mutdataCodons[i]]
      if(dataAA[i] == mutdataAA[i] & dataCodons[i] != mutdataCodons[i] ){
        #silent mutation
        nSilent <- nSilent + 1
      } else if(dataAA[i] != mutdataAA[i] & mutdataAA[i] != "*"){
        #missense mutation
        nMissense <- nMissense + 1
      } else if(mutdataAA[i] == "*" & dataAA[i] != "*"){
        #nonsense mutation
        nNonsense <- nNonsense + 1
      } else {
        #no mutation
      }
    }
    }
    #put results into a dataframe
    frame <- data.frame("FA" = FA, "N" = N, "nSilent" = nSilent, "nMissense" = nMissense, "nNonsense" = nNonsense)
    #return results
    return(frame)
  }
}


readIntOGen <- function(IN) {
  # Purpose: read and parse an IntOGen mutation data file. Return only the
  #            number of silent, missense, and nonsense point mutations.
  #            All indels are ignored.
  # Parameters:
  #     IN   chr      Filename of an IntOGen mutation data file.
  # Value:   list     List with the following elements:
  #                      nSilent    num the number of silent mutations
  #                      nMissense  num the number of missense mutations
  #                      nNonsense  num the number of nonsense mutations
  if(TRUE){
  #read data
  data <- read_tsv(IN)
  #remove header
  data <- data[-1]
  #counter variables
  nSilent <- 0
  nMissense <- 0
  nNonsense <- 0
  #loop over data
  for(i in 1:length(data$AA_CHANGE)){
    a <- unlist(strsplit(data$REF[i],""))
    b <- unlist(strsplit(data$ALT[i],""))
    if(data$MOST_SEVERE[i] == 'synonymous_variant'){
      nSilent <- nSilent+1
    } else if (data$MOST_SEVERE[i] == 'missense_variant'){
      nMissense <- nMissense + 1
    } else if(data$MOST_SEVERE[i] == 'stop_gained'){
      nNonsense <- nNonsense + 1
    }
  }
  #put results into a data frame
  frame <- data.frame("IN" = IN,  "nSilent" = nSilent, "nMissense" = nMissense, "nNonsense" = nNonsense)
  #return results
  return(frame)
  }
}

if (FALSE){
IN <- '/Users/Matt/Documents/BCH441/ABC-units/data/intogen-KRAS-distribution-data.tsv'
IN <- '/Users/Matt/Documents/BCH441/ABC-units/data/intogen-OR1A1-distribution-data.tsv'
IN <- '/Users/Matt/Documents/BCH441/ABC-units/data/intogen-PTPN11-distribution-data.tsv'
readIntOGen(IN)

set.seed(2398)
file <- '/Users/Matt/Documents/BCH441/ABC-units/data/PTPN11_HSa_coding.fa'
file <- '/Users/Matt/Documents/BCH441/ABC-units/data/OR1A1_HSa_coding.fa'
file <- '/Users/Matt/Documents/BCH441/ABC-units/data/KRAS_HSa_coding.fa'
evalMut(file,10000)
}
if(FALSE){
  IN <- '/Users/Matt/Documents/BCH441/ABC-units/data/intogen-PTPN11-distribution-data.tsv'
  readIntOGen(IN)
  file <- '/Users/Matt/Documents/BCH441/ABC-units/data/PTPN11_HSa_coding.fa'
  evalMut(file,10000)
}
if(FALSE{
  File <- '/Users/Matt/Documents/BCH441/ABC-units/data/GCsample.fa'
  evalMut(File,10000)
})
