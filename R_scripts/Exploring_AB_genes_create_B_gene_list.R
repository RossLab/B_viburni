

# master anno
freeze.v0.genes.anno <- read.csv("output/freeze.v0.genes.anno.complete.csv")

# B chromosome assignment
asn_table <- read.table('output/scaffolds.final.assignment.tsv', header = T, sep = '\t')
asn_table[asn_table$seq %in% c("scaffold_360", "scaffold_957"), 'b.status.final'] <- 'B-A'

table(asn_table$b.status.final)

rownames(asn_table) <- asn_table$seq
freeze.v0.genes.anno$b.status.final <- asn_table[freeze.v0.genes.anno$seq, 'b.status.final']

# Total genes on A/B/Bc
table(freeze.v0.genes.anno$b.status)
#     A     B   B-A    Bc
# 23305    55   149   120

# Total genes on A/B/Bc with functional annotation
table(freeze.v0.genes.anno$b.status[freeze.v0.genes.anno$anno == 'Y'])
#     A     B   B-A    Bc
# 13402    20    40    53


# Fraction genes on A/B/Bc with functional annotation
round(100 * (table(freeze.v0.genes.anno$b.status[freeze.v0.genes.anno$anno == 'Y']) / table(freeze.v0.genes.anno$b.status)), 2)
#     A     B    Bc
# 57.51 29.41 44.17

write.table(freeze.v0.genes.anno[, c('gene', 'seq', 'b.status.final')], 'output/genes.by.scaffolds.tsv', sep = '\t', quote = F, row.names = F, col.names = T)
