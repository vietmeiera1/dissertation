#Insert target size is ~460bp 
#Load dada2 

install.packages("dada2") #wont work for R version 4.0.2

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.11")

install.packages("devtools")
library("devtools")
devtools::install_github("benjjneb/dada2", ref="v1.16") # change the ref argument to get other versions

install.packages("~/github/dada2",
                 repos = NULL,
                 type = "source",
                 dependencies = c("Depends", "Suggests","Imports"))



library(dada2)

#Load Sequences in to R
Boyce16SPond4_Seq<-"All.Sequences_16Seq2020" 
                            #"DATA FILE"

#Confirm Sequences loaded in ok 
list.files(Boyce16SPond4_Seq)

#Assign names to forward and reverse reads 
fn_Pond4SUFs <-sort(list.files(Boyce16SPond4_Seq, pattern="SU-B4-A_S68_L001_R1_001.fastq", full.names = TRUE))
fn_Pond4SURs <-sort(list.files(Boyce16SPond4_Seq, pattern="SU-B4-A_S68_L001_R2_001.fastq", full.names = TRUE))

length(fn_Pond4SUFs) #2
length(fn_Pond4SURs) #2 
# ^^ These should match 

#Forward read quality 
plotQualityProfile(fn_Pond4SUFs[1:2])

#Reverse read quality 
plotQualityProfile(fn_Pond4SURs[1:2])

#Extract names from the reads for downstream analysis 
sample.namesF <-sapply(strsplit(basename(fn_Pond4SUFs), "_"), `[`, 1)
sample.namesR <-sapply(strsplit(basename(fn_Pond4SURs), "_"), `[`, 1)


#Place Filtered Files into a new folder 
filtFs <-file.path(~Download, "filtered", paste0(sample.namesF, "_F_filt.fastq"))
filtRs <-file.path(~Download, "filtered", paste0(sample.namesR, "_R_filt.fastq"))
names(filtFs) <-sample.namesF
names(filtRs) <-sample.namesR

#Perform necessary trimming based on quality profiles
# trimLeft = based on quality profile generated above ^^
# trimRight = based on quality profile generated above ^^
# maxN = 
# maxEE = 
# truncQ = 
out20 <-filterAndTrim(fn_Pond4SUFs, filtFs, fn_Pond4SURs, filtRs, trimLeft=c(20,20),
                      trimRight=c(100,120),
                      maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                      compress=TRUE, multithread=TRUE, matchIDs=TRUE)


#FilteredSequences Result –tells you the “reads.in” and “reads.out”
list(out20)

#parametric error model (err)
#Depending on number of sequences this step can take a long time to run 
errF <-learnErrors(filtFs, multithread=TRUE)
errR <-learnErrors(filtRs, multithread=TRUE)

dadaFs <-dada(filtFs, err=errF, multithread=TRUE) 
dadaRs <-dada(filtRs, err=errR, multithread=TRUE) 

#Merge the forward and reverse reads 
mergers <-mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

#amplicon sequence variant table (ASV) table
seqtab <-makeSequenceTable(mergers)
dim(seqtab)#tells you how many samples you have and how many unique ASV you have

#Distribution of fragment size
table(nchar(getSequences(seqtab)))

#Chimeras 
seqtab.nochim <-removeBimeraDenovo(seqtab,method="consensus", multithread=TRUE, 
                                   verbose=TRUE)

#number of ASVs after chimeras removed
dim(seqtab.nochim)

#percent of sequences remaining
sum(seqtab.nochim)/sum(seqtab)

#summarize the quality filtering 
getN <-function(x) sum(getUniques(x))
track <-cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
colnames(track) <-c("input", "filtered", "denoisedF", "nonchim")
rownames(track) <-sample.names

#lists the results of the quality filtering (provides input, filtered, denoised, and nonchim)
list(track) 

#assign taxa (This takes a very long time, set this up to run over night)
taxa <-assignTaxonomy(seqtab.nochim, "~/Documents/silva_nr_v132_train_set.fa.gz", multithread=TRUE)

#turn ASV taxa table into a text file that is saved on your computer 
write.table(taxa, "~/Documents/SulfurCyclingTaxa.txt", sep="\t")#example

#to get species, you run a second taxa file called “taxaSpecies” that only labels species identity if it is 100% match 
taxaSpecies <-addSpecies(taxa, "~/Documents/silva_species_assignment_v132.fa.gz")

SavetaxSpecies into a text file
write.table(taxaSpecies, "~/Documents/SulfurCyclingTaxa_wSpecies.txt", sep="\t")#example

###This ends the dada2 portion of analysis, now phyloseq R package will be utilized
####PHYLOSEQ
#activating necessary packages library(phyloseq); packageVersion("phyloseq")library(Biostrings); packageVersion("Biostrings")library(ggplot2); packageVersion("ggplot2")theme_set(theme_bw())
#Read in the metadata file that characterizes all of your samples meta =read.csv(“FILE LOCATION”) 
#identify the sample names in the meta table row.names(meta) <-meta$Sample_ID
#construct phyloseq object directly from outputs created in dada2> ps <-phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), sample_data(meta), tax_table(taxa))
#Rename sequences to ASVs anf store the raw sequence as a refseq slot for later reference by entering refseq(ps)dna <-Biostrings::DNAStringSet(taxa_names(ps))names(dna) <-taxa_names(ps)ps <-merge_phyloseq(ps, dna)taxa_names(ps) <-paste0("ASV", seq(ntaxa(ps)))ps
#you can remove a sample by using the following script (i.e. if it does not meet the quality standards or belongs to another experiment)LessSample<-prune_samples(sample_names(ps) != "Sample", ps)
#Save OTU table, taxa table, and refseq (raw sequences) as text files and send to computer write.table(otu_table(ps), "~/Documents/SulfurCycling_OTU.txt", sep="\t")
#examplewrite.table(tax_table(ps), "~/Documents/SulfurCycling_taxa_OTU.txt", sep="\t")
#examplewrite.table(refseq(ps), "~/Documents/SulfurCycling_seq_OTU.txt", sep="\t")
#example#Prevalance of Phyla 
#Compute prevalence of each feature, store as data.frameprevdf = apply(X = otu_table(ps), MARGIN= ifelse(taxa_are_rows(ps), yes = 1, no = 2),FUN = function(x){sum(x > 0)})
#Add taxonomy and total read counts to this data.frameprevdf = data.frame(Prevalence = prevdf, TotalAbundance = taxa_sums(ps), tax_table(ps))
530plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
#Normalized PCoA Plot by transforming into proportions ps2Sulfur.prop <-transform_sample_counts(ps, function(otu) otu/sum(otu))Norm.pcoa <-ordinate(physeq = ps2Sulfur.prop, method = “PCoA”, distance = “bray”)plot_ordination(physeq = ps2Sulfur.prop, ordination = Norm.pcoa, shape = “Type”, color = “Sample_ID”, title = “PCoA of Bacterial Communities”)
#NMDS Plot#Transform data to proportions as appropriate for Bray Curtis psSulfur.prop <-transform_sample_counts(ps, function(otu) otu/sum(otu))ord.nmds.bray <-ordinate(psSulfur.prop, method=”NMDS”, distance=”bray”)
#Diversity between samplesplot_richness(ps, x = "Type", title = "Diversity by Type") + geom_boxplot()
#Diversity between samples with just observed metricplot_richness(ps, x = "Sample_ID", title = "Sample Diversity", measures=c("Observed"))
#With color coding with boxplotplot_richness(ps, x = "Type", color = "Sample_Location", title = "Diversity by Type", measures=c("Observed")) + geom_boxplot()
#With color coding without boxplotplot_richness(ps, x = "Type", color = "Sample_Location", title = "Diversity by Type", measures=c("Observed")) 
#Top 20 Bar Plot top20 <-names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]ps.top20 <-transform_sample_counts(ps, function(OTU) OTU/sum(OTU))ps.top20 <-prune_taxa(top20, ps.top20)plot_bar(ps.top20, x="Type", fill="Family")
#relative abundance table psSulfur_ra = transform_sample_counts(ps, function(x){x / sum(x)})#turn it into text file write.table(otu_table(psSulfur_ra), "~/Documents/SulfurCycling_taxa_RelativeAbundance.txt", sep="\t")
#example
#Take Subset of Samples –example is “raw” but can be any metadata category rawSamples <-subset_samples(ps, Type==”raw”)
#Filter specific bacterial families –example is SRB families SRB <-subset_taxa(ps, Family=="Campylobacteraceae" | Family=="Desulfarculaceae" | Family=="Desulfobacteraceae" | Family=="Desulfobulbaceae" |
531Family=="Desulfomicrobiaceae" | Family=="Desulfovibrionaceae" | Family=="Desulfuromonadaceae" | Family=="Syntrophaceae" |Family=="Syntrophobacteraceae")Filter for SRB families from full relative abundance sheet based on list compiled SRB.ra <-subset_taxa(psSulfur_ra, Family=="Campylobacteraceae" | Family=="Desulfarculaceae" | Family=="Desulfobacteraceae" | Family=="Desulfobulbaceae" | Family=="Desulfomicrobiaceae" | Family=="Desulfovibrionaceae" | Family=="Desulfuromonadaceae" | Family=="Syntrophaceae" | Family=="Syntrophobacteraceae")
#Turn both into text files for reading write.table(otu_table(SRB.ra), "~/Documents/SRB_RelativeAbundance.txt", sep="\t")
#examplewrite.table(otu_table(SRB), "~/Documents/SRB_Raw.txt", sep="\t")
#example
#Filter for SOB families based on list compiledSOB <-subset_taxa(ps, Family=="Acidithiobacillaceae" | Family=="Beggiatoaceae" | Family=="Bradyrhizobiaceae" | Family=="Burkholderiaceae" | Family=="Chlorobiaceae" | Family=="Chromatiaceae" | Family=="Comamonadaceae" | Family=="Ectothiorhodospiraceae" | Family=="Helicobacteraceae" | Family=="Hydrogenophilaceae" | Family=="Hyphomicrobiaceae"| Family=="Piscirickettsiaceae" | Family=="Rhodobacteraceae" | Family=="Rhodocyclaceae" | Family=="Thiotrichaceae")
#Filter for SOB families from full relative abundance sheet based on list compiled SOB.ra <-subset_taxa(psSulfur_ra, Family=="Acidithiobacillaceae" | Family=="Beggiatoaceae" | Family=="Bradyrhizobiaceae" | Family=="Burkholderiaceae" | Family=="Chlorobiaceae" | Family=="Chromatiaceae" | Family=="Comamonadaceae" | Family=="Ectothiorhodospiraceae" | Family=="Helicobacteraceae" | Family=="Hydrogenophilaceae" | Family=="Hyphomicrobiaceae" | Family=="Piscirickettsiaceae" | Family=="Rhodobacteraceae" | Family=="Rhodocyclaceae" | Family=="Thiotrichaceae")
#Plot tables as bar graphs of SRB and SOBplot_bar(SRB, x="Sample_Location", fill="Family")plot_bar(SRB.ra, x="Sample_ID", fill="Family")plot_bar(SOB, x="Sample_Location", fill="Family")plot_bar(SOB.ra, x="Sample_ID", fill="Family")
#Turn both into text files for reading write.table(otu_table(SOB.ra), "~/Documents/SOB_RelativeAbundance.txt", sep="\t")
#examplewrite.table(otu_table(SOB), "~/Documents/SOB_Raw.txt", sep="\t")

#example
#Increases size on all samples + geom_point(size = 4)
#Increasing font size+ theme(text = element_text(size = 15))

532

