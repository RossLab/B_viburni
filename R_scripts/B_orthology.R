viburni_B_genes_file <- 'viburni/B_genes.tsv'
B_genes_tab <- read.table(viburni_B_genes_file, col.names = c('gene', 'chr'))

rename_genes <- function(gene_names){
  gene_name_len <- 6
  zeros_counts <- gene_name_len - nchar(gene_names)
  zeros_to_add <- sapply(zeros_counts, function(x){ paste0(rep('0', x), collapse = '') })
  paste0('g', zeros_to_add, substr(gene_names, 2, nchar(gene_names)))
}

row.names(B_genes_tab) <- rename_genes(B_genes_tab$gene)

####
# missing gene length potentially?

##################
# WITHIN VIBURNI #
##################

within_viburni_file <- 'blastout/Within_viburni_reciprocal_OG_pairs.tsv'
viburni_orthologs <- read.table(within_viburni_file, header = T)

viburni_orthologs$gene1 <- sapply(strsplit(viburni_orthologs$gene1, '.t'), function(x){ x[1] })
viburni_orthologs$gene2 <- sapply(strsplit(viburni_orthologs$gene2, '.t'), function(x){ x[1] })

# filtering all self-hits of alternative transcripts
viburni_orthologs <- viburni_orthologs[viburni_orthologs$gene1 != viburni_orthologs$gene2, ]

# getting chromosomal assignments to the table
viburni_orthologs$ch1 <- B_genes_tab[rename_genes(viburni_orthologs$gene1), 'chr']
viburni_orthologs$ch2 <- B_genes_tab[rename_genes(viburni_orthologs$gene2), 'chr']

viburni_orthologs[is.na(viburni_orthologs$ch1), 'ch1'] <- 'A'
viburni_orthologs[is.na(viburni_orthologs$ch2), 'ch2'] <- 'A'

# removing within autosome
B_viburni <- viburni_orthologs[!c(viburni_orthologs$ch1 == 'A' & viburni_orthologs$ch2 == 'A'), ]
# removing alternative transcripts
B_viburni <- B_viburni[!duplicated(paste0(B_viburni$gene1, B_viburni$gene2)), ]

# retaining only Autosome Bchromosome paralogs
B_A_viburni <- B_viburni[c(B_viburni$ch1 == 'A' | B_viburni$ch2 == 'A'), ]

Bs_with_orthologs <- as.character(apply(B_A_viburni, 1, function(x){x[c(2,3)][which(c(x[8], x[9]) != 'A' )]}))
ortholog_names <- as.character(apply(B_A_viburni, 1, function(x){x[c(2,3)][which(c(x[8], x[9]) == 'A' )]}))
B_genes_tab[rename_genes(Bs_with_orthologs), 'core_paralog'] <- ortholog_names
B_genes_tab[rename_genes(Bs_with_orthologs), 'core_similarity'] <- round(B_A_viburni$identity, 2)

##############
# LNGISPINUS #
##############

longispinus_file <- 'blastout/longispinus_and_B_viburni_reciprocal_OG_pairs.tsv'
longispinus_orthologs <- read.table(longispinus_file, header = T)
longispinus_orthologs <- longispinus_orthologs[grepl('Pvib', longispinus_orthologs$gene1) | grepl('Pvib', longispinus_orthologs$gene2), ]
longispinus_orthologs <- longispinus_orthologs[!(grepl('Pvib', longispinus_orthologs$gene1) & grepl('Pvib', longispinus_orthologs$gene2)), ]
long_orthologs <- table(longispinus_orthologs$gene1)
# table(substr(names(long_orthologs), nchar(names(long_orthologs)) - 1, nchar(names(long_orthologs))))

long_ortholog_names <- as.character(sapply(names(long_orthologs), function(x){ glo = longispinus_orthologs[longispinus_orthologs$gene1 == x, ]; glo[which.max(glo[, 'aln_len']), 'gene2'] }))
long_ortholog_names <- sapply(strsplit(long_ortholog_names, '[.]'), function(x){ x[1] } )
long_sim <- round(sapply(names(long_orthologs), function(x){ glo = longispinus_orthologs[longispinus_orthologs$gene1 == x, ]; glo[which.max(glo[, 'aln_len']), 'identity'] }), 2)
Bgenes_with_longi_orth <- sapply(strsplit(names(long_orthologs), "_"), function(x){ x[2] })
B_genes_tab[rename_genes(Bgenes_with_longi_orth), 'longispinus_ortholog'] <- long_ortholog_names
B_genes_tab[rename_genes(Bgenes_with_longi_orth), 'longispinus_similarity'] <- long_sim

##############
# SOLENOPSIS #
##############

solenopsis_file <- 'blastout/solenopsis_and_B_viburni_reciprocal_OG_pairs.tsv'
solenopsis_orthologs <- read.table(solenopsis_file, header = T)
solenopsis_orthologs <- solenopsis_orthologs[grepl('Pvib', solenopsis_orthologs$gene1) | grepl('Pvib', solenopsis_orthologs$gene2), ]
solenopsis_orthologs <- solenopsis_orthologs[!(grepl('Pvib', solenopsis_orthologs$gene1) & grepl('Pvib', solenopsis_orthologs$gene2)), ]
sol_orthologs <- table(solenopsis_orthologs$gene2)
# table(substr(names(sol_orthologs), nchar(names(sol_orthologs)) - 1, nchar(names(sol_orthologs))))

Bgenes_with_sol_orth <- sapply(strsplit(solenopsis_orthologs[, 'gene2'], "_"), function(x){ x[2] })
B_genes_tab[rename_genes(Bgenes_with_sol_orth), 'solenopsis_ortholog'] <- solenopsis_orthologs$gene1
B_genes_tab[rename_genes(Bgenes_with_sol_orth), 'solenopsis_similarity'] <- solenopsis_orthologs$identity

############
# HOMOLOGY #
############

ncbi <- read.table('blast_taxon_annot/B_genes.blast_taxa_overview.tsv', sep = '\t')
Bgenes_with_ncbi_hit <- sapply(strsplit(ncbi[, 1], "_"), function(x){ x[2] })
B_genes_tab[rename_genes(Bgenes_with_ncbi_hit), 'NCBI_hits'] <- sapply(strsplit(ncbi[, 2], ','), function(x){ paste(names(table(x)), collapse = ', ') } )

uniprot <- read.table('blast_taxon_annot/B_genes.diamond_taxa_overview.tsv', sep = '\t')
Bgenes_with_uniprot_hit <- sapply(strsplit(uniprot[, 1], "_"), function(x){ x[2] })
B_genes_tab[rename_genes(Bgenes_with_uniprot_hit), 'uniprot_hits'] <- sapply(strsplit(uniprot[, 2], ','), function(x){ paste(names(table(x)), collapse = ', ') } )

B_genes_tab <- B_genes_tab[order(B_genes_tab$gene), ]
B_genes_tab <- B_genes_tab[order(B_genes_tab$chr, B_genes_tab$core_similarity, B_genes_tab$longispinus_similarity, B_genes_tab$solenopsis_similarity, B_genes_tab$NCBI, B_genes_tab$uniprot), ]

write.table(B_genes_tab, 'B_orthology_final.tsv', row.names = F, sep = '\t', quote = F)
