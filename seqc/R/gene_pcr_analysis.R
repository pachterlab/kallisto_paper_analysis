library(data.table)
library(dplyr)
library(mamabear)

mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
dataset = "hsapiens_gene_ensembl", host="sep2015.archive.ensembl.org")
transcripts_to_genes <- biomaRt::getBM(
  attributes = c("ensembl_transcript_id", "ensembl_gene_id",
    "external_gene_name"), mart = mart)
transcripts_to_genes <- dplyr::rename(transcripts_to_genes,
  target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id)

summarized_by_genes <- function(abundance, mapping, which_column) {
  abundance <- as.data.frame(abundance, stringsAsFactors = FALSE)
  abundance <- inner_join(abundance, mapping, by = 'target_id')
  abundance <- group_by_(abundance, which_column)

  summarize(abundance, tpm = sum(tpm))
}

#
# load data
#
pre <- file.path('../results/A/BGI', 1:4)

cuff <- lapply(file.path(pre, 'cufflinks', 'genes.fpkm_tracking'), read_cufflinks, 200)
cuff <- lapply(cuff, rename, ens_gene = target_id)

cuff_bias <- lapply(file.path(pre, 'cufflinks_bias', 'genes.fpkm_tracking'), read_cufflinks, 200)
cuff_bias <- lapply(cuff_bias, rename, ens_gene = target_id)

sailfish <- lapply(file.path(pre, 'sailfish', 'quant.sf'), read_sailfish)
sailfish <- lapply(sailfish, summarized_by_genes, transcripts_to_genes, 'ens_gene')

sailfish_bias <- lapply(file.path(pre, 'sailfish_bias', 'quant_bias_corrected.sf'), read_sailfish)
sailfish_bias <- lapply(sailfish_bias, summarized_by_genes, transcripts_to_genes, 'ens_gene')

xprs <- lapply(file.path(pre, 'express', 'results.xprs'), read_xprs)
for(i in 1:length(xprs)){
      xprs[[i]]$tpm = 10^6*xprs[[i]]$fpkm/sum(xprs[[i]]$fpkm)
      xprs[[i]] <- data.frame(xprs[[i]])
}
xprs <- lapply(xprs, summarized_by_genes, transcripts_to_genes, 'ens_gene')

rsem <- lapply(file.path(pre, 'rsem', 'out.isoforms.results'), read_rsem)
rsem <- lapply(rsem, summarized_by_genes, transcripts_to_genes, 'ens_gene')

kal <- lapply(file.path(pre, 'kallisto', 'abundance.tsv'), read_kallisto_rename)
kal <- lapply(kal, summarized_by_genes, transcripts_to_genes, 'ens_gene')

kal_bias <- lapply(file.path(pre, 'kallisto_bias', 'abundance.tsv'), read_kallisto_rename)
kal_bias <- lapply(kal_bias, summarized_by_genes, transcripts_to_genes, 'ens_gene')

read_emsar <- function(fname) {
  df <- read.table(fname, stringsAsFactors = FALSE, header = TRUE)
  df <- dplyr::select(df, target_id = transcriptID, est_counts = iReadcount,
    tpm = TPM)

  df
}

emsar <- lapply(file.path(pre, 'emsar', 'emsar.0.fpkm'), read_emsar)
emsar <- lapply(emsar, summarized_by_genes, transcripts_to_genes, 'ens_gene')

#
# qPCR
#

# qpcr = readRDS("qpcra_processed.rds")
qpcr <- readRDS('../data/gene_pcr.rds')
names(qpcr) <- c('ens_gene', 'value')
# names(qpcr) = c("MIQE_location", "target_id", "value")

compare_qpcr <- function(pred, qpcr, id){
  df <- inner_join(pred, qpcr, by='ens_gene')
  data.frame(method=id,
    spearman = with(df, cor(value, tpm, meth="spear")),
    pearson = with(df, cor(value, log2(tpm + 1), method = "pearson")),
    stringsAsFactors=FALSE)
}

all_qpcr <- c(
  lapply(cuff, compare_qpcr, qpcr, "cufflinks"),
  lapply(cuff_bias, compare_qpcr, qpcr, "cufflinks (bias)"),
  lapply(sailfish, compare_qpcr, qpcr, "sailfish"),
  lapply(sailfish_bias, compare_qpcr, qpcr, "sailfish (bias)"),
  lapply(xprs, compare_qpcr, qpcr, "express"),
  lapply(rsem, compare_qpcr, qpcr, "rsem"),
  lapply(emsar, compare_qpcr, qpcr, "emsar"),
  lapply(kal, compare_qpcr, qpcr, "kallisto"),
  lapply(kal_bias, compare_qpcr, qpcr, "kallisto (bias)")
  ) %>% bind_rows()
all_qpcr$rep <- 1:4
all_qpcr <- dplyr::arrange(
  dplyr::group_by(all_qpcr, rep),
  method, rep)
#all_qpcr %>% print(n = nrow(.))

write.table(all_qpcr, file = '../results/gene_pcr_analysis.tsv',
  quote = FALSE, sep = '\t', col.names = TRUE)
