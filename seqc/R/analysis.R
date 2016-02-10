library(data.table)
library(dplyr)
library(mamabear)


#
# load data
#

pre <- file.path('../results/A/BGI', 1:4)

cuff <- lapply(file.path(pre, 'cufflinks', 'isoforms.fpkm_tracking'), read_cufflinks, 200)
cuff_bias <- lapply(file.path(pre, 'cufflinks_bias', 'isoforms.fpkm_tracking'), read_cufflinks, 200)
sailfish <- lapply(file.path(pre, 'sailfish', 'quant.sf'), read_sailfish)
sailfish_bias <- lapply(file.path(pre, 'sailfish_bias', 'quant_bias_corrected.sf'), read_sailfish)
xprs <- lapply(file.path(pre, 'express', 'results.xprs'), read_xprs)
for(i in 1:length(xprs)){
      xprs[[i]]$tpm = 10^6*xprs[[i]]$fpkm/sum(xprs[[i]]$fpkm)
      xprs[[i]] <- data.frame(xprs[[i]])
}
rsem <- lapply(file.path(pre, 'rsem', 'out.isoforms.results'), read_rsem)
kal <- lapply(file.path(pre, 'kallisto', 'abundance.tsv'), read_kallisto_rename)
kal_bias <- lapply(file.path(pre, 'kallisto_bias', 'abundance.tsv'), read_kallisto_rename)

read_emsar <- function(fname) {
  df <- read.table(fname, stringsAsFactors = FALSE, header = TRUE)
  df <- dplyr::select(df, target_id = transcriptID, est_counts = iReadcount,
    tpm = TPM)

  df
}
emsar <- lapply(file.path(pre, 'emsar', 'emsar.0.fpkm'), read_emsar)


#
# qPCR
#

qpcr <- readRDS("../data/qpcra_processed.rds")
names(qpcr) = c("MIQE_location", "target_id", "value")

compare_qpcr <- function(pred, qpcr, id){
 df <- inner_join(pred, qpcr, by='target_id')
 data.frame(method=id, spearman=with(df, cor(-value, tpm, meth="spear")), stringsAsFactors=FALSE)
}


all_qpcr <- c(lapply(cuff, compare_qpcr, qpcr, "cufflinks"),
              lapply(cuff_bias, compare_qpcr, qpcr, "cufflinks (bias)"),
              lapply(sailfish, compare_qpcr, qpcr, "sailfish"),
              lapply(sailfish_bias, compare_qpcr, qpcr, "sailfish (bias)"),
              lapply(xprs, compare_qpcr, qpcr, "express"),
              lapply(rsem, compare_qpcr, qpcr, "rsem"),
              lapply(emsar, compare_qpcr, qpcr, "emsar"),
              lapply(kal, compare_qpcr, qpcr, "kallisto"),
              lapply(kal_bias, compare_qpcr, qpcr, "kallisto (bias)")) %>% bind_rows()
all_qpcr$rep = 1:4

all_qpcr %>% arrange(spearman) %>% print(n = nrow(.))

write.table(all_qpcr, file = '../results/isoform_pcr_analysis.tsv',
  quote = FALSE, sep = '\t', col.names = TRUE)
