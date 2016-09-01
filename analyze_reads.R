#library calls to set up my R terminal environment
library(sleuth)

#modify directory to Zika data
base_dir <- "/Users/lynnyi/workspace/Zika Data/"
sample_id <- dir(file.path(base_dir, "kallisto_output"))
sample_id
#delete some sample_ids

#modify directory to Zika data
kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, "kallisto_output", id))
kal_dirs

#modify directory of design matrix
s2c <- read.table("/Users/lynnyi/workspace/sleuth/zikareadsbasedmatrix.txt", header = TRUE, stringsAsFactors = FALSE)
s2c <- dplyr::mutate(s2c, path = kal_dirs)
s2c

mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "hsapiens_gene_ensembl",
  host = 'ensembl.org')

t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
  ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

so <- sleuth_prep(s2c, ~read, target_mapping = t2g)
so <- sleuth_fit(so)
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')

results_table <- sleuth_results(so, 'reduced:full', test_type = 'lrt')

#modify sleuth_live call?
sleuth_live(so)


