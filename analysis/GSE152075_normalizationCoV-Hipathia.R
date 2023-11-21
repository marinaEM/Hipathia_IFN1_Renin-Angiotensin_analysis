
###########################################################################################
###   Normalization of GSE152075 SARS-CoV-2 dataset for CoV-Hipathia web tool         #####
###   at  http://hipathia.babelomics.org/covid19/                                     #####
###   Author: Marina Esteban-Medina marina.esteban@juntadeandalucia.es                #####
###                                                                                   #####
###########################################################################################

# setwd() directory containing the GEO dataset


pacman::p_load("here", "dplyr", "openxlsx", "preprocessCore", "edgeR")

#### Normalize Data  #####

# Load  RNA-seq expression set from raw counts file.
download.file(url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE152075&format=file&file=GSE152075%5Fraw%5Fcounts%5FGEO%2Etxt%2Egz",
              destfile = "GSE152075_raw_counts_GEO.txt.gz")

expreset_raw <- read.delim(file = gzfile("GSE152075_raw_counts_GEO.txt.gz"), header = T, sep = " ")

identifiers <- rownames(expreset_raw)

identifiers_df <- data.frame(id = identifiers,
                             entrez = mapIds(org.Hs.eg.db, keys = identifiers, column = "ENTREZID", keytype = "SYMBOL"), stringsAsFactors = F) %>% .[!is.na(.$entrez),] ## identify weird codes which do not have entrez ID

expreset_raw <- expreset_raw[rownames(expreset_raw) %in% identifiers_df$id, ] ## here we clean the dataset from pseudogenes and rows with none entrez ids whichwill not be used by Hipathia

print("read and clean from non entrez codes ...done")

# Explore how our data is organized
hist(as.numeric(expreset_raw[2,]),breaks=100)
var(as.numeric(expreset_raw[2,]))

getVari <- apply(expreset_raw, 1, var)
hist(getVari,100)

# Normalization by TMM with "edgeR" package
dge <- DGEList(counts=expreset_raw)
print("dge...done")
tmm <-  calcNormFactors(dge, method="TMM")
print("tmm...done")
logcpm <- cpm(tmm, prior.count=3, log=TRUE)
print("logcpm...done")

#### Do Hipathia rescale to speed up the web process ####
gExp = logcpm
gExp = normalize.quantiles(gExp)
rownames(gExp)<- rownames(logcpm)
colnames(gExp) <- colnames(logcpm)

## Construct the metadata file with the info from GEO GSE152075
metadata <- data.frame( fileName = colnames(gExp), type = c(rep("V",430), rep("C",54))) 

trans_data <- translate_data(gExp, "hsa")
trans_data_export <- tibble::rownames_to_column(as.data.frame(trans_data), var ="X")


#### Export results table for the Cov-Hipathia Disease Maps differential signaling example.

write.table(trans_data_export, file = here("GSE152075_SARS-COV2_nasopharingeal", "GSE152075trans_data_corrected.tsv"), col.names = T, row.names = F, quote = F, sep = "\t")
write.table(metadata, file = here("GSE152075_SARS-COV2_nasopharingeal", "GSE152075samples.tsv"), row.names = F, col.names = T, quote = F, sep = "\t" )


