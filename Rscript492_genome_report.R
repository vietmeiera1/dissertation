

#how many genes have been identified in the annotated sequence?
#3900 observations were made indicating 3900 genes

#Dimensions of data frame
dim(genome_report)
#3900    9

#PART 2: Basic Exploration of Data 
#Data classes
class(genome_report$V1) #"factor"
class(genome_report$V2) #"factor"
class(genome_report$V3) #"factor"
class(genome_report$V4) #integer
class(genome_report$V5) #integer
class(genome_report$V6) #factor
class(genome_report$V7) #factor
class(genome_report$V8) #factor


#Q2.b levels of the cdsStartStat and cdsEndStat 
levels(genome_report$V1) #lists all 70 nodes used for assembly
levels(genome_report$V2) # "."
levels(genome_report$V3) # "gene"
levels(genome_report$V4) # NULL
levels(genome_report$V5) #NULL
levels(genome_report$V6) # "."
levels(genome_report$V7) # "-" +"
levels(genome_report$V8) # "."

#Add new column to dataframe that subtracts V4 & V5
#genome_report$DNALength=
#|genome_report$V4 - genome_report$V5|
#C = data$columnV4 - data$columnV5

mycolumnafteradding <- genome_report$V4 - genome_report$V5 #this did work, but didnt add it as a column
abs(mycolumnafteradding)
#now need to add this data as a column
print(mycolumnafteradding)

#Q2.c minimum, maximum, and mean of the number of Protein Length 
min(genome_report$V4) # 1
max(genome_report$V4) #424982
mean(genome_report$V4) #77725.37



############################################



#minimum, maximum, and mean of the number of gene length
min(mv492cytochrome$geneLength) #48
max(mv492cytochrome$geneLength) #2458176
mean(mv492cytochrome$geneLength) #66354.44

#Q2.f How many genes are likely not coding based on the absense of cds start and end site?
summary(mv492cytochrome$Feature.id) #not meaningful
summary(mv492cytochrome$Type)#not meaningful
summary(mv492cytochrome$ProductFunction) #not meaningful
summary(mv492cytochrome$Location) #not meaningful
summary(mv492cytochrome$ProteinLength)#give min, 1st Qu, Median, Mean, 3rd Qu, Max... This could be graphed!
summary(mv492cytochrome$ProtineSeq) #not meaningful
summary(mv492cytochrome$DNALength) #gives compl 41267 none 4789
summary(mv492cytochrome$DNAseq) #gives length 0 class NULL mode NULL
summary(mv492cytochrome$geneLength) #give min, 1st Qu, Median, Mean, 3rd Qu, Max... This could be graphed!


#make object vector of protein length, dna length, gene length
proteinlength <- mv492cytochrome$ProteinLength
dnalength <- mv492cytochrome$DNALength
genelength <- mv492cytochrome$geneLength

#density of above vectors
proteinlengthdistr<-density(proteinlength) #okay
dnalengthdistr<-density(dnalength) #gives error that this is NOT numeric?
genelengthdistr<-density(genelength) #okay




#Q3.a Describe the distribution of exon numbers per genes
plot(proteinlengthdistr, main="Distribution of Protein Length per MV492", xlab = "protein length", ylab = "density")
#^^okay

#Q3.b Describe the distribution of gene length across all genes
plot(genelengthdistr, main="Distribution of Gene Length for MV492", xlab = "gene length", ylab = "density")
#^^okay





#how many genes are coded for on + strand?

#how many genes are coded for on - strand?

#how many genes as associated with cytochrome?


#Q2.e how many genes annoated for your species
#under dataframe in environment click the blue play button that opens a dropdown menu
#under name2 lists the actual number of genes that have been annotated



#PART 3 Plotting the distribution of variables
#make object vector of exon count 
exonnumber<- elephant.genes$exonCount
#make object vector of gene length
genelength<- elephant.genes$geneLength

#density of above vectors
exondistr<-density(exonnumber)
lengthdistr<-density(genelength)

#create a plot of density for both vectors
#label x & y axis
#plot(density(exonnumber), main="Distribution of Exons per Gene.Elephant", xlab = "density", ylab = "exon count")
#plot(density(genelength), main="Distribution of Gene Length.Elephant", xlab = "density", ylab = "exon count")

#Q3.a Describe the distribution of exon numbers per genes
plot(exondistr, main="Distribution of Exons per Gene.Elephant", xlab = "exon count", ylab = "density")

#Q3.b Describe the distribution of gene length across all genes
plot(lengthdistr, main="Distribution of Gene Length.Elephant", xlab = "gene length", ylab = "density")


#PART 4 Plotting the two variables against each other 
plot(x=exonnumber, y=genelength, main="Exon Number vs. Gene Length for Loxodonta africana (Elephant)", xlab ="Number of Exons", ylab="Gene Length")

