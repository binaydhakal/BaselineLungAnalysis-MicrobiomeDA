library(mia)
library(vegan)
library(ggpubr)
library(tidyverse)
library(scater)
library(lubridate)
library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(ggsignif)
library(ggpubr)
library(vegan)
library(picante)
require(cowplot)
require(ggrepel)
library(ecodist)
library(pheatmap)
library(dendextend)
#library(ANCOMBC)
#library(ALDEx2)
library(knitr)
library(reshape2)
library(palmerpenguins)
library(ggforce)
#library(DESeq2)
library(scran)
library(ggordiplots)
library(bluster)
library(DirichletMultinomial)
library(parallel)
library(tableone)
library(edgeR)
library(survival)
library(survminer)
library(data.table)
library(GGally)
library(networkD3)
library(htmlwidgets)


setwd(getwd())
all_counts <- read.csv('../MbioData/0102_ALIR_2023.taxa.genus.cmF.cln.SampleID.summary_table.tsv',
                       sep = '\t')  

row.names(all_counts)<-all_counts$sample_id

MetaData <- read.csv(file = "../MbioData/alirmbio_6.29.23.csv")

MetaData <- MetaData %>% distinct(sample_id, .keep_all = TRUE)
row.names(MetaData) <- MetaData$sample_id

MetaData_healthy_COVID <- read.csv(file = "../MbioData/Mbiome_SampleMetadata_2.24.23.csv")

MetaData_healthy_COVID <- MetaData_healthy_COVID %>% distinct(sample_id, .keep_all = TRUE)
MetaData_healthy_COVID <- MetaData_healthy_COVID[-which(is.na(MetaData_healthy_COVID$sample_id)),]
row.names(MetaData_healthy_COVID) <- MetaData_healthy_COVID$sample_id
included_studies = which(MetaData_healthy_COVID$StudyName == "alir" | 
                           MetaData_healthy_COVID$StudyName == "alirseeds"|
                           MetaData_healthy_COVID$StudyName == "ExtrNeg"|
                           MetaData_healthy_COVID$StudyName == "fmt"|
                           MetaData_healthy_COVID$StudyName == "lhmp"|
                           MetaData_healthy_COVID$StudyName == "PCRNeg"|
                           MetaData_healthy_COVID$StudyName == "PCRPosZymo")


MetaData_healthy_COVID = MetaData_healthy_COVID[included_studies,]

alirSampleID_3_3_23 <- read_csv("../MbioData/alirSampleID_3.3.23.csv")

### look at covid missingness
length(which(is.na(alirSampleID_3_3_23$COVID)))

alirSampleID_3_3_23$SubjectIDDay[is.na(alirSampleID_3_3_23$COVID)]
### 5130 is covid
alirSampleID_3_3_23$COVID[alirSampleID_3_3_23$SubjectID=="5130"] <- "Yes"
### rest is noncovid
alirSampleID_3_3_23$COVID[is.na(alirSampleID_3_3_23$COVID)] <- "No"


MetaData_healthy_COVID_merged = merge(alirSampleID_3_3_23,MetaData_healthy_COVID, by = "SubjectIDDay", all.y = TRUE)

lhmp_OW_indx = which(MetaData_healthy_COVID_merged$StudyName=="lhmp" & MetaData_healthy_COVID_merged$SampleType == "oral")

lhmp_BAL_indx = which((MetaData_healthy_COVID_merged$StudyName=="lhmp" & MetaData_healthy_COVID_merged$SampleType == "IS") | 
                        (MetaData_healthy_COVID_merged$StudyName=="lhmp" & MetaData_healthy_COVID_merged$SampleType == "BAL"))
fmt_indx = which(MetaData_healthy_COVID_merged$StudyName=="fmt" & MetaData_healthy_COVID_merged$SampleType == "stool")

MetaData_healthy_COVID_merged$compartment = NA

MetaData_healthy_COVID_merged$group[lhmp_OW_indx] = "healthyCtrl"
MetaData_healthy_COVID_merged$group[lhmp_BAL_indx] = "healthyCtrl"
MetaData_healthy_COVID_merged$group[fmt_indx] = "healthyCtrl"

MetaData_healthy_COVID_merged$SampleType[lhmp_OW_indx] = "OW"
MetaData_healthy_COVID_merged$SampleType[lhmp_BAL_indx] = "BAL/IS"
MetaData_healthy_COVID_merged$SampleType[fmt_indx] = "FMT"


MetaData_healthy_COVID_merged$compartment[lhmp_OW_indx] = "oral"
MetaData_healthy_COVID_merged$compartment[lhmp_BAL_indx] = "lung"
MetaData_healthy_COVID_merged$compartment[fmt_indx] = "gut"

MetaData_healthy_COVID_merged$studydayfollowup[lhmp_OW_indx] = "healthyCtrl"
MetaData_healthy_COVID_merged$studydayfollowup[lhmp_BAL_indx] = "healthyCtrl"
MetaData_healthy_COVID_merged$studydayfollowup[fmt_indx] = "healthyCtrl"

MetaData_healthy = MetaData_healthy_COVID_merged[which(MetaData_healthy_COVID_merged$group=="healthyCtrl"),] 
MetaData_COVD = MetaData_healthy_COVID_merged[which(MetaData_healthy_COVID_merged$COVID=="Yes"),] 

MetaData_ExprtCtrl = MetaData_healthy_COVID_merged[which(MetaData_healthy_COVID_merged$SampleType%in%c("BALCon","TCON","OWCon",
                                                                                                       "ExtrNeg","PCRNeg",
                                                                                                       "PCRPosZymo")),] 

table(MetaData_healthy_COVID_merged$SampleType)

names(MetaData)[which(names(MetaData) == "SubjectID.x")] = "SubjectID"
names(MetaData_healthy)[which(names(MetaData_healthy) == "SubjectID.x")] = "SubjectID"
names(MetaData_COVD)[which(names(MetaData_COVD) == "SubjectID.x")] = "SubjectID"
names(MetaData_ExprtCtrl)[which(names(MetaData_ExprtCtrl) == "SubjectID.x")] = "SubjectID"

names(MetaData)[which(names(MetaData) == "SampleType.x")] = "SampleType"


table(MetaData$rectalswab_soiled)

gut.indx1 = which(MetaData$SampleType =="stool")
gut.indx2 = which(MetaData$SampleType == "rectal" & MetaData$rectalswab_soiled == "yes")
gut.indx3 = which(MetaData$SampleType == "rectal" & is.na(MetaData$rectalswab_soiled))

gut.indx = unique(c(gut.indx1,gut.indx2,gut.indx3))

meta_ETA = MetaData[which(MetaData$SampleType=="ETA"),] #| MetaData$SampleType=="BAL"
table(meta_ETA$SampleType)
ETA_subjects = unique(meta_ETA$SubjectID)
indx = list()
sub = "1031"
for(sub in ETA_subjects){
  x = MetaData$RespSample_Depletion[which(MetaData$SubjectID==sub)]
  
  undep_indx = which(x=="Undep")
  dep_indx =  which(x=="Dep")
  if(length(undep_indx)>0){
    indx = c(indx,which(MetaData$SubjectID==sub)[undep_indx])
  }else{
    indx = c(indx,which(MetaData$SubjectID==sub )[dep_indx])
  }
  
}

ETA.indx = unlist(indx)

oral.indx = which( MetaData$SampleType=="oral")
MetaData$compartment = MetaData$SampleType

MetaData$compartment[oral.indx] = "oral"
MetaData$compartment[ETA.indx] = "lung"
MetaData$compartment[gut.indx] = "gut"

MetaData = MetaData[which(MetaData$compartment %in% c("oral","lung","gut")),]

table(MetaData$studydayfollowup, MetaData$compartment)
# table(MetaData$studydayfollowup, MetaData$compartment, is.na(MetaData$ShannonIndex))

x = intersect(names(MetaData),names(MetaData_healthy))
MetaData = MetaData[,x]
MetaData_healthy = MetaData_healthy[,x]
MetaData_COVD = MetaData_COVD[,x]
MetaData_ExprtCtrl = MetaData_ExprtCtrl[,x]

row.names(MetaData_healthy) = MetaData_healthy$sample_id
row.names(MetaData_COVD) = MetaData_COVD$sample_id
row.names(MetaData_ExprtCtrl) = MetaData_ExprtCtrl$sample_id

MetaData = rbind(MetaData,MetaData_healthy,MetaData_COVD,MetaData_ExprtCtrl)

identical(row.names(MetaData), MetaData$sample_id)
# matched_ids = intersect(all_counts$sample_id,MetaData$sample_id)

# remove samples without clinical information
counts = all_counts[MetaData$sample_id,-c(1,2)]
identical(row.names(MetaData), row.names(counts))


# MetaData = MetaData[-which(is.na(MetaData$compartment)),]
table(MetaData$compartment)
length(which(MetaData$total16sreads<300))

MetaData1 = MetaData[-which(MetaData$total16sreads<300),]
MetaData2 = MetaData[-which(MetaData$total16sreads<500),]
MetaData3 = MetaData[-which(MetaData$total16sreads<1000),]
MetaData4 = MetaData[-which(MetaData$total16sreads<2000),]


counts1 = counts[MetaData1$sample_id,]
counts2 = counts[MetaData2$sample_id,]
counts3 = counts[MetaData3$sample_id,]
counts4 = counts[MetaData4$sample_id,]


counts_rare = counts
metadata_rare = MetaData


metadata_rare$compartment[metadata_rare$SampleType=="BAL/IS"] = "BAL/IS"
metadata_rare$compartment[metadata_rare$SampleType=="BALCon"] = "BALCon"
metadata_rare$compartment[metadata_rare$SampleType=="ExtrNeg"] = "ExtrNeg"
metadata_rare$compartment[metadata_rare$SampleType=="FMT"] = "FMT"
metadata_rare$compartment[metadata_rare$SampleType=="OW"] = "OW"
metadata_rare$compartment[metadata_rare$SampleType=="OWCon"] = "OWCon"
metadata_rare$compartment[metadata_rare$SampleType=="PCRNeg"] = "PCRNeg"
metadata_rare$compartment[metadata_rare$SampleType=="PCRPosZymo"] = "PCRPosZymo"
metadata_rare$compartment[metadata_rare$SampleType=="TCON"] = "TCON"
metadata_rare = metadata_rare[-which(is.na(metadata_rare$compartment)),]

metadata_rare$compartment =  as.factor(metadata_rare$compartment)
metadata_rare$compartment = ordered(metadata_rare$compartment,levels = c("oral","lung","gut",
                                                                         "OW","BAL/IS","FMT",
                                                                         "TCON","OWCon","BALCon",
                                                                         "ExtrNeg","PCRNeg","PCRPosZymo"))
table(metadata_rare$compartment)

metadata_rare = metadata_rare[which(metadata_rare$compartment=="oral" | metadata_rare$compartment=="lung" |metadata_rare$compartment=="gut"),]
counts_rare = counts[metadata_rare$sample_id,]

table(metadata_rare$compartment)


metadata_rare$compartment = factor(metadata_rare$compartment, levels = c("oral","lung","gut"))


baseline_lung_metadata <- metadata_rare[which(metadata_rare$compartment == 'lung' & metadata_rare$studydayfollowup == "baseline"), ]
table(baseline_lung_metadata$studydayfollowup, baseline_lung_metadata$compartment)

# Subset counts_rare to include only rows (sample_ids) present in baseline_lung_metadata
baseline_lung_counts <- counts_rare[rownames(counts_rare) %in% rownames(baseline_lung_metadata), ]


df <- data.frame(x = names(baseline_lung_counts))
taxa <- df %>% tidyr::separate(x, c("Kingdom", "Phylum", "Class", "Order", 
                                    "Family", "Genus"), sep = "\\.")
# use Genus as column names
names(counts_rare) <- taxa$Genus


# Creating a TreeSummarizedExperiment Object

#baseline_lung_counts <- as.matrix(baseline_lung_counts)

tse <- TreeSummarizedExperiment(
  assays =  list(counts = t(baseline_lung_counts)),
  colData = baseline_lung_metadata,
  rowData = taxa)


# Subsetting

# Filtering out zero variance features

rowData(tse)[["sd"]] <- rowSds(assay(tse, "counts"))
# Plot
hist(log(rowData(tse)[["sd"]]))

thresold <- 1

selected <- rowData(tse)[['sd']] > thresold

tse_selected <- tse[selected, ]
tse_selected

# Found that subsetting reduces the dimension of the data from 1459*418 to 237*418
#----------------------------------------------------------------------------------------


# Transformations & Subsetting

tse_transformed <- transformAssay(tse, method = 'relabundance')

tse_transformed <- subsetByPrevalentFeatures(tse_transformed, rank="Genus", detection = 1/10000, abund_values = "relabundance",
                                         prevalence = 5/100)

tse_transformed <- transformAssay(x=tse_transformed, method = "clr", pseudocount = 1, name = "clr")

# Transformation & Subsetting reduces the dimension of the data from 1459*418 to 127*418
#----------------------------------------------------------------------------------------

# EDA
# Abundance Analysis.

tse_abundance_transform <- transformAssay(tse, method = "relabundance")

plotAbundanceDensity(
  tse_abundance_transform, layout = "jitter",
  assay.type = "relabundance",
  n = 40, point_size=1, point_shape=19,
  point_alpha=0.1) +
  scale_x_log10(label=scales::percent)

plotAbundanceDensity(
  tse_abundance_transform, layout = "density",
  assay.type = "relabundance",
  n = 5, colour_by="Mortality30Day",
  point_alpha=1/10) +
  scale_x_log10()

getPrevalence(
  tse, detection = 1, sort = TRUE, assay.type = "counts",
  as.relative = FALSE) |>
  head()


# Diversity Analysis

tse <- estimateDiversity(
  tse, abund_values="relabundance", index = "shannon", name="ShannonIndex")

tse <- estimateDiversity(
  tse, abund_values="relabundance", index = "faith", name="Faith")

###  Beta Diversity Analysis using dimensional reduction techniques like PCA, PCoA, tSNE

set.seed(5252)

tse <- transformAssay(tse, method = "relabundance")

counts <- assay(tse, "counts")
libsizes <- colSums(baseline_lung_counts)
size.factors <- libsizes/mean(libsizes)
logcounts(tse) <- log2(t(t(baseline_lung_counts)/size.factors) + 1)

pca_data <- prcomp(t(logcounts(tse)), rank=50)
