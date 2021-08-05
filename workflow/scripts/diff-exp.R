library(tidyverse)
library(DESeq2)

## Args

#FILES = Sys.glob("results/star/*/ReadsPerGene.out.tab")
FILES = snakemake@input
FORMULA = as.formula("~condition")

message(FILES)

## Diffexp

# set names
FILES <- FILES %>% set_names(.,str_extract(.,"(?<=star\\/).+(?=\\/Reads)"))

# make coldata
cd <- enframe(FILES,"sample","file")

# yb is treated as control for germline - this is also done in the original manuscript
cd <- cd %>%
  mutate(condition = map_lgl(sample,~any(str_detect(.,c("ovary","w_","yb"))))) %>%
  #mutate(condition = map_lgl(sample,~any(str_detect(.,c("w_","yb"))))) %>%
  #mutate(condition = map_lgl(sample,~any(str_detect(.,c("w_"))))) %>%
  mutate(condition = ifelse(condition,"control",sample)) %>%
  column_to_rownames("sample")

dat <- FILES %>%
  map_df(~read_tsv(.,skip = 4,col_names = c("gene_id","unstranded","first_strand","second_strand")),.id="sample")

strand_counts <- list("unstranded","first_strand","second_strand") %>%
  set_names(.,.) %>%
  map(~sum(dat[,.]))

which_col <- ifelse(between(strand_counts$first_strand/strand_counts$unstranded,0.3,0.7),"unstranded",
                    names(which.max(strand_counts[c("first_strand","second_strand")])))

# take only the counts from the appropriate strand
dat <- dat %>% dplyr::select(sample,gene_id,all_of(which_col))

# spread to wide shape for matrix conv
mat <- pivot_wider(dat,names_from = "sample",values_from = which_col,id_cols = "gene_id") %>%
  column_to_rownames("gene_id") %>%
  as.matrix()

# make deseq2 obj
dds <- DESeqDataSetFromMatrix(mat,colData = cd, design = FORMULA)

dds$condition <- relevel(dds$condition,"control")

dds <- dds[rowSums(counts(dds)) > 1,]

dds <- DESeq(dds)

#rdds <- rlog(dds)
#plotPCA(rdds,intgroup="condition")

res <- resultsNames(dds) %>%
  set_names(.,.) %>%
  .[names(.)!="Intercept"] %>%
  map(~results(dds,name=.,alpha=0.05)) %>%
  imap(~as_tibble(lfcShrink(dds,res=.x,coef = .y,type = "normal",),rownames = "gene_id")) %>%
  #map(as_tibble, rownames="gene_id") %>%
  bind_rows(.id="comparison")

## export

saveRDS(dds, snakemake@output[['dds']])
write_tsv(res,snakemake@output[['tsv']])
