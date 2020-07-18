
#Select workplace
#loading used packages
library(phyloseq)
library(ggplot2)
library(vegan)
library(Tax4Fun)

#Functional prediction
QIIMESingleData <- importQIIMEData("Input_Tax4fun.txt")
folderReferenceData <- "SILVA115/"
Tax4FunOutput <- Tax4Fun(QIIMESingleData, folderReferenceData, 
                         fctProfiling = T, refProfile = "UProC", shortReadMode = TRUE, 
                         normCopyNo = TRUE)
Tax4FunProfile<-Tax4FunOutput$Tax4FunProfile
Tax4FunProfile <- data.frame(t(Tax4FunOutput$Tax4FunProfile))
#View(Tax4FunProfile)

write.table(Tax4FunProfile, "IAC_L0_L3.txt",  dec = ".")



##Load biom file manually
otumat<-read.csv2("otu_mat.csv", header=T, row.names=1, sep=";", dec=".")
taxmat<-read.csv2("tax_mat.csv", header=T, row.names = 1, sep=";", dec=".")

#converting to phyloseq class
otumat<-as.matrix(otumat)
taxmat<-as.matrix(taxmat)
class(otumat) #matrix
class(taxmat)#matrix
OTU<-otu_table(otumat, taxa_are_rows = T)
TAX<-tax_table(taxmat)
physeq<-phyloseq(OTU,TAX)

#loading metadata
meta<-read.table("mapping_file.txt",row.names = 1, header=T, sep="")

# Check that the rownames match the sample names
all(rownames(meta) %in% sample_names(physeq))
#View(meta)

#if TRUE
meta<-sample_data(meta)
physeq<-merge_phyloseq(physeq,meta)

#check the complete file
class(physeq) #phyloseq file
rank_names(physeq)
physeq

#inserting tree
tree<-read_tree("rep_set.tre")
physeq<-merge_phyloseq(physeq, tree)


#Rarefazendo os dados
phyloseq_raref <- rarefy_even_depth(physeq, rngseed = T)
phyloseq_raref

###Ordination Total
all.ord<-ordinate(phyloseq_raref , distance="bray", method="CAP", ~ Cultivar)
anova(all.ord)

#Ordination plot of samples
plot<-plot_ordination(phyloseq_raref, all.ord, color="Dose")
###Separing in plant x phosphorus treatment
physeq_P1F1 = subset_samples(phyloseq_raref, Subset_P1 == "P1F1")
physeq_P1F1
head(sample_data(physeq_P1F1)$Description, 20) 

physeq_P1F2 = subset_samples(phyloseq_raref, Subset_P1 == "P1F2")
physeq_P1F2
head(sample_data(physeq_P1F2)$Description, 20) 

physeq_P2F1 = subset_samples(phyloseq_raref, Subset_P2 == "P2F1")
physeq_P2F1
head(sample_data(physeq_P2F1)$Description, 20) 

physeq_P2F2 = subset_samples(phyloseq_raref, Subset_P2 == "P2F2")
physeq_P2F2
head(sample_data(physeq_P2F2)$Description, 20) 


##ordination for each treatment changing "physeq_P1F1/P1F2/P2F1/P2F2"
all.ord<-ordinate(physeq_P1F1, distance="bray", method="CAP", ~ Dose)
anova(all.ord)

#Ordination plot of samples
plot<-plot_ordination(physeq_P1F1, all.ord, color="Dose")

######################################
###################################### P source ordination
######################################

#Subset physeq by cultivar 
physeq_IAC = subset_samples(phyloseq_raref, Cultivar == "IAC")
physeq_IAC
head(sample_data(physeq_IAC)$Cultivar, 20) 

physeq_Dor = subset_samples(phyloseq_raref, Cultivar == "Dor")
physeq_Dor
head(sample_data(physeq_Dor)$Cultivar, 20) 

#class(physeq)
physeq_IAC #10638 OTUs
physeq_Dor #10638 OTUs

#Removing zero sum OTUS's 

bac_IAC = prune_taxa (taxa_sums(physeq_IAC) > 0, physeq_IAC)
bac_IAC #9645

bac_DOR = prune_taxa (taxa_sums(physeq_Dor) > 0, physeq_Dor)
bac_DOR #9665

#Removing L0

bac_IAC_P = subset_samples(bac_IAC, TP != "P0")
bac_IAC_P
head(sample_data(bac_IAC_P) $ TP, 3) #Levels PR e SPT

bac_DOR_P = subset_samples(bac_DOR, TP != "P0")
bac_DOR_P
head(sample_data(bac_DOR_P) $ TP, 3) #Levels PR e SPT

#Transform the count in relative abundance

prop_IAC <- transform_sample_counts(bac_IAC_P,  function(x) 100 * x/sum(x))
prop_IAC #9645 OTUs 30 Samples 11 variables
prop_DOR <- transform_sample_counts(bac_DOR_P,  function(x) 100 * x/sum(x))
prop_DOR #9665 OTUs 30 Samples 11 variables
#Ordination
ord_IAC<-ordinate(prop_IAC, distance="bray", method="PCoA")#, ~TP)

plot_IAC<-plot_ordination(prop_IAC, ord_IAC, color="TP", shape="Dose")

ord_DOR<-ordinate(prop_DOR, distance="bray", method="PCoA")#, ~TP)

plot_DOR<-plot_ordination(prop_DOR, ord_DOR, color="TP", shape="Dose")


#For further analysis if needed 
#summary(ord_IAC)
#pcoa_IAC_Scores<-as.data.frame(ord_IAC$vectors)
#summary(pcoa_IAC_Scores)
#pcoa_IAC_values<-as.data.frame(ord_IAC$values) #brings broken stick values
#write.table(pcoa_IAC_Scores, file = "PCOA_IAC.txt", sep = " ", dec = ",")



###Differentially enriched OTUs were acessed with QIIME 1.9.1
###Ploted with iTol
###Available at https://itol.embl.de/login.cgi?logout=1

###############################################
##### Metagenome Prediction using Tax4Fun #####
###############################################

QIIMESingleData <- importQIIMEData("Input_Tax4fun.txt")
folderReferenceData <- "SILVA115/"
Tax4FunOutput <- Tax4Fun(QIIMESingleData, folderReferenceData, 
                         fctProfiling = T, refProfile = "UProC", shortReadMode = TRUE, 
                         normCopyNo = TRUE)
Tax4FunProfile<-Tax4FunOutput$Tax4FunProfile
Tax4FunProfile <- data.frame(t(Tax4FunOutput$Tax4FunProfile))
#View(Tax4FunProfile)


write.table(Tax4FunProfile, "Predicted_Kegg_Functions.txt",  dec = ".")

###Differentially enriched functions were acessed with QIIME 1.9.1
