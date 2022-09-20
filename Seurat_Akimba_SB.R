setwd("~/R/msAkimba_scRNA")

library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)

# Using data from https://www.ebi.ac.uk/gxa/sc/experiments/E-MTAB-9061/downloads?plotType=tsne&plotOption=20

# akimba <- Seurat::ReadMtx("E-MTAB-9061.aggregated_filtered_counts.mtx", )

# Akimba.seurat <- Seurat::ReadMtx(
#   mtx = akimba,
#   cells = "E-MTAB-9061.aggregated_filtered_counts.mtx_cols",
#   features = "E-MTAB-9061.aggregated_filtered_counts.mtx_rows"
# )

# Error in parse_url(url = uri) : length(url) == 1 is not TRUE
# I think I need to remove one of the columns in mtx_rows since they're identical
rows.tmp <- read.delim("E-MTAB-9061.aggregated_filtered_counts.mtx_rows", header = FALSE)
head(rows.tmp)

rows.tmp <- rows.tmp[-c(2)]

# Akimba.seurat <- Seurat::ReadMtx(
#   mtx = akimba,
#   cells = "E-MTAB-9061.aggregated_filtered_counts.mtx_cols",
#   features = rows.tmp,
#   feature.column = 1
# )
# Same error as above
# Let's throw mtx_cols into a variable and see what happens

cols.tmp <- read.delim("E-MTAB-9061.aggregated_filtered_counts.mtx_cols", header = FALSE)

# Akimba.seurat <- Seurat::ReadMtx(
#   mtx = akimba,
#   cells = cols.tmp,
#   features = rows.tmp,
#   feature.column = 1
# )
# Same error. :<
# Since the two *.tmp variables are one column, the problem must be in akimba
# Let's try some different Read options

#akimba <- Seurat::Read10X("C:/Users/Sam/Documents/R/msAkimba_scRNA")
# Looks like it's expecting specific file names. Let's see if renaming my stuff
# will make a difference.

akimba <- Seurat::Read10X("C:/Users/Sam/Documents/R/msAkimba_scRNA")
# I think it worked!

Akimba.seurat <- CreateSeuratObject(counts = akimba, project = "Akimba",
                                    assay = "scRNA",
                                    min.cells = 3, min.features = 200)

# Works, but instead of gene names I've got ensembl numbers, which are basically unreadable.
# Going to try https://www.reddit.com/r/bioinformatics/comments/cc5db6/help_with_converting_ensembl_id_to_gene_name/

BiocManager::install('org.Hs.eg.db')
library(AnnotationDbi)
library(org.Hs.eg.db)
ensembl2entrez <- as.list(org.Hs.egENSEMBL)
# This maps things from ENSEMBL numbers to entrez IDs

entrez2genesymbol <- as.list(org.Hs.egSYMBOL)
# This maps from entrez IDs to gene symbols, not the full gene names (too long)

ensembl_list <- Akimba.seurat@assays[['scRNA']]
# that gets me all the data, now I just need to figure out how to subset to
# get only Ensembl numbers

ensembl_rownames <- rownames(ensembl_list@data)
# this gets me all the ensembl numbers from the dataset

# now to set up a loop to iterate over the list and swap out the numbers for gene names
for (x in 1:length(ensembl_rownames)) {
  iter_variable <- ensembl_rownames[x]
  iter_variable <- sub("ENSMUSG", "", iter_variable) # strip leading "ENSMUSG" characters
  iter_variable <- sub("^0+", "", iter_variable) # strip leading zeroes
  if(is.null(entrez2genesymbol[[iter_variable]])) {next} # simple error correction
  gene_name <- entrez2genesymbol[[iter_variable]] # put gene name into gene_name
  ensembl_rownames[x] <- gene_name # switch ensembl row name to gene name
}

rownames(ensembl_list@data) <- ensembl_rownames


