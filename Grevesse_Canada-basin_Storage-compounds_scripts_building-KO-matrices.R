########## Code to build a table containing metagenome information: genes, scaffolds, annotations, etc##########
### Notes: For each metagenome, all the files are retrieved from JGI - A loop is run to create one table per samples

### Directory and packages
dir<-"/Volumes/Extreme SSD/Research"

library("dplyr")
library("stringr")
library("tidyr")

### Load needed file that contains the metagenome sample names and information
metadata<-read.csv(paste(dir,"/Analysis-Results/Arctic2017/Metagenomics/Gene_abundance_matrices/Names_Table_noAG5.csv",sep=""))

### Build a table with the genes, genes info, genes annotation and their scaffold for each sample
for (i in 1:length(metadata$Station)){
  
## Download all the necessary files with different annotations for genes
# coverage files and scaffold names map comes from JGI right away
cov<-read.delim(header=TRUE,quote="",paste(dir,"/Analysis-Results/Arctic2017/Metagenomics/Gene_abundance_matrices/",metadata$Station[i],"/","reads.v.final.contigs.sorted.bam.cov",sep=""))
names.map<-read.delim(header=FALSE,quote="",paste(dir,"/Analysis-Results/Arctic2017/Metagenomics/Gene_abundance_matrices/",metadata$Name_ID[i],".assembled.names_map",sep=""))

# all the annotation files
gff<-read.delim(header=FALSE,quote="",paste(dir,"/Analysis-Results/Arctic2017/Metagenomics/data/Annotation-NMDC/",metadata$Station[i],"/",metadata$Station[i],"_structural_annotation.gff",sep=""))
scaffold.phyl<-read.delim(header = FALSE,quote="",paste(dir,"/Analysis-Results/Arctic2017/Metagenomics/data/Annotation-NMDC/",metadata$Station[i],"/",metadata$Station[i],"_scaffold_lineage.tsv",sep=""))
gene.phyl<-read.delim(header = FALSE, quote="",paste(dir,"/Analysis-Results/Arctic2017/Metagenomics/data/Annotation-NMDC/",metadata$Station[i],"/",metadata$Station[i],"_gene_phylogeny.tsv",sep=""))
ko<-read.delim(header = FALSE,quote="",paste(dir,"/Analysis-Results/Arctic2017/Metagenomics/data/Annotation-NMDC/",metadata$Station[i],"/",metadata$Station[i],"_ko.tsv",sep=""))
ec<-read.delim(header = FALSE,quote="",paste(dir,"/Analysis-Results/Arctic2017/Metagenomics/data/Annotation-NMDC/",metadata$Station[i],"/",metadata$Station[i],"_ec.tsv",sep=""))
cog<-read.delim(header = FALSE,quote="",paste(dir,"/Analysis-Results/Arctic2017/Metagenomics/data/Annotation-NMDC/",metadata$Station[i],"/",metadata$Station[i],"_cog.gff",sep=""))
prod.name<-read.delim(header = FALSE,quote="",paste(dir,"/Analysis-Results/Arctic2017/Metagenomics/data/Annotation-NMDC/",metadata$Station[i],"/",metadata$Station[i],"_product_names.tsv",sep=""))
tigr<-read.delim(header = FALSE,quote="",paste(dir,"/Analysis-Results/Arctic2017/Metagenomics/data/Annotation-NMDC/",metadata$Station[i],"/",metadata$Station[i],"_tigrfam.gff",sep=""))
pfam<-read.delim(header = FALSE,quote="",paste(dir,"/Analysis-Results/Arctic2017/Metagenomics/data/Annotation-NMDC/",metadata$Station[i],"/",metadata$Station[i],"_pfam.gff",sep=""))

## Prepare the files to modify and keep what we need
cov.temp<-cov %>%
  rename(Scaf_name=X.ID,Scaf_length=Length) %>%
  select(Scaf_name,Avg_fold,Scaf_length)

names.map.temp<-names.map %>%
  rename(Scaf_name=V1,Scaf_id=V2)

gff.temp<-gff %>%
  rename(Scaf_id=V1,Seq_type=V3,Start=V4,Stop=V5,Gene_direction=V7) %>%
  mutate(Gene_id=paste(Scaf_id,"_",Start,"_",Stop,sep="")) %>%
  select(Scaf_id,Gene_id,Seq_type,Start,Stop,Gene_direction)

scafold.phyl.temp<-scaffold.phyl %>%
  rename(Scaf_id=V1,Scaf_phyl=V2) %>%
  select(Scaf_id,Scaf_phyl)

gene.phyl.temp<-gene.phyl %>%
  rename(Gene_id=V1,Gene_phyl=V5,Gene_phyl_perc=V4) %>%
  select(Gene_id,Gene_phyl,Gene_phyl_perc)

ko.temp<-ko %>%
  rename(Gene_id=V1,KO=V3,KO_SimPerc=V4) %>%
  mutate(KO=str_replace(KO, "KO:", "")) %>%
  select(Gene_id,KO,KO_SimPerc)

ec.temp<-ec %>%
  rename(Gene_id=V1,EC=V3,EC_SimPerc=V4) %>%
  select(Gene_id,EC,EC_SimPerc)

cog.temp<-cog %>%
  rename(Gene_id=V1,COG=V3) %>%
  select(Gene_id,COG)

prod.name.temp<-prod.name %>%
  rename(Gene_id=V1,Prod_name=V2,Prod_name_source=V3) %>%
  select(Gene_id,Prod_name,Prod_name_source)

tigr.temp<-tigr %>%
  rename(Gene_id=V1,Tigr=V3) %>%
  select(Gene_id,Tigr)

pfam.temp<-pfam %>%
  rename(Gene_id=V1,Pfam=V3) %>%
  select(Gene_id,Pfam)

## Add all the info for genes, scaffolds to the scaffold abundance table
final.table.temp<-cov.temp %>%
  left_join(names.map.temp,by="Scaf_name") %>%
  left_join(gff.temp,by="Scaf_id") %>%
  left_join(scafold.phyl.temp,by="Scaf_id") %>%
  left_join(gene.phyl.temp,by="Gene_id") %>%
  left_join(ko.temp,by="Gene_id") %>%
  left_join(ec.temp,by="Gene_id") %>%
  left_join(cog.temp,by="Gene_id") %>%
  left_join(prod.name.temp,by="Gene_id") %>%
  left_join(tigr.temp,by="Gene_id") %>%
  left_join(pfam.temp,by="Gene_id")

## save one table for each metagenome sample
write.csv(final.table.temp,paste(dir,"/Analysis-Results/Arctic2017/Metagenomics/results/",metadata$Station[i],"_fully_annotated_gene_table.csv",sep=""),row.names = FALSE)
} 


##################################################################################################


########## Build a table with total KO abundance per metagenome sample with normalization step##########
### Directory and packages
dir<-"/Volumes/ExtremeSSD/Research"

library("dplyr")
library("stringr")
library("tidyr")
library("reshape2")

### Load needed file with sample information
metadata<-read.csv(paste(dir,"/Analysis-Results/Arctic2017/Metagenomics/Gene_abundance_matrices/Names_Table_noAG5.csv",sep=""))

### Loop to download all the tables, sum up abundances by KO and store in a list

ko.aggregation<-list()

for (i in 1:length(metadata$Station)){
  temp<-read.csv(paste(dir,"/Analysis-Results/Arctic2017/Metagenomics/results/",metadata$Station[i],"_fully_annotated_gene_table.csv",sep=""))
  ko.abundance<-temp %>%
    filter(Scaf_length>=1500) %>%
    select(Gene_id,KO,Avg_fold) %>%
    distinct(Gene_id,KO,.keep_all = TRUE) %>%
    aggregate(Avg_fold~KO, FUN=sum) %>%
    mutate(Sample=metadata$Station[i])
  
  ko.aggregation[[i]]<-ko.abundance
}


### Put all samples with aggregated KO abundance value together and normalize to get a "per cell abundance"
# the per cell abudnance is obtained by dividing the abundance of KO by the median abundance of 10 single-copy conserved genes (identified by their KO)
ko.abundance.table<-do.call(rbind , ko.aggregation) %>%
  dcast(KO~Sample, value.var="Avg_fold") %>%
  replace(is.na(.),0)

ko.universal<-c("K06942", "K01889", "K01887", "K01875", "K01883", "K01869", "K01873", "K01409", "K03106", "K03110")

sample.universal<-ko.abundance.table[ko.abundance.table$KO %in% ko.universal,2:ncol(ko.abundance.table)]
medians.universal<-apply(sample.universal,2,median)

ko.abundance.table.nmd<-ko.abundance.table[,2:ncol(ko.abundance.table)]
rownames(ko.abundance.table.nmd)<-ko.abundance.table$KO

for (i in 1:ncol(ko.abundance.table.nmd)){
  ko.abundance.table.nmd[,i]<-ko.abundance.table.nmd[,i]/medians.universal[i]
} 

write.csv(ko.abundance.table.nmd,paste(dir,"/Analysis-Results/Arctic2017/Metagenomics/results/Canada-basin-mg_ko-abundance_normalization-per-cell.csv",sep=""),row.names = TRUE)
##################################################################################################


########## Make the global ocean abundance table - Use this step if separate KO abundance matrices are created for separate groups of metagenomic samples ##########
### Directory and packages
dir<-"/Volumes/ExtremeSSD/Research"

library("dplyr")
library("stringr")
library("tidyr")
library("reshape2")

### Load needed files
cb<-read.csv(paste(dir,"/Analysis-Results/Arctic2017/Metagenomics/results/Canada-basin-mg_ko-abundance_normalization-per-cell.csv",sep=""))
tara<-read.csv(paste(dir,"/Analysis-Results/Commonfiles/Metagenomes/TARA/Results/TARA-mg_ko-abundance_normalization-per-cell.csv",sep=""))
malaspina<-read.csv(paste(dir,"/Analysis-Results/Commonfiles/Metagenomes/MALASPINA/Results/","Malaspina-mg_ko-abundance_normalization-per-cell.csv",sep=""))

colnames(tara)[1]<-"KO"
colnames(cb)[1]<-"KO"
### make the global ocean file
global<-cb %>%
  left_join(tara,by="KO") %>%
  left_join(malaspina,by="KO")

write.csv(global,paste(dir,"/Analysis-Results/Arctic2017/Microbiome-DOM/MT-MG/Results/Global-ocean-mg_ko-abundance_normalization-per-cell.csv",sep=""),row.names = FALSE)
#######################################################################################################