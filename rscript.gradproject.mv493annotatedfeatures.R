genome_report <- read.delim("~/Downloads/gradproject/genome_report.tab", stringsAsFactors=TRUE)
   View(genome_report)



   #how many genes have been identified in the annotated sequence?
#5048 observations were made indicating 5048 genes
   
    #Q1.b dimensions of data frame
#dim(mv493.genomereport)
dim(genome_report)

#PART 2: Basic Exploration of Data 
#Q2.a data classes
class(genome_report$Feature.type) #factor
class(genome_report$Contig) #factor
class(genome_report$Location) #factor
class(genome_report$Strand) #factor
class(genome_report$Feature.function) #factor

#Q2.b levels 
levels(genome_report$Feature.ID) #not meaningful
levels(genome_report$Contig) #provides list of all Nodes, might be meaningful
levels(genome_report$Location) #not meaningful
levels(genome_report$Strand) # gives answer of +/-, might be meaningful
levels(genome_report$Feature.function) #not meaningful


#how many genes are coded for on + strand?

#how many genes are coded for on - strand?

#how many genes as associated with cytochrome?
   
   

   
   #Q2.c minimum, maximum, and mean of the number of exons per gene 
   min(elephant.genes$exonCount)
   max(elephant.genes$exonCount)
   mean(elephant.genes$exonCount)
   
   #Q2.d shortest, longest, and mean length of all genes in species
   min(elephant.genes$geneLength)
   max(elephant.genes$geneLength)
   mean(elephant.genes$geneLength)
   
   #Q2.e how many genes annoated for your species
   #under dataframe in environment click the blue play button that opens a dropdown menu
   #under name2 lists the actual number of genes that have been annotated
   
   #Q2.f How many genes are likely not coding based on the absense of cds start and end site?
   summary(elephant.genes$cdsStartStat)
   summary(elephant.genes$cdsEndStat)
   
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
   
   