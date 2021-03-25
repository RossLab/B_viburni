from Bio import SeqIO
from collections import defaultdict

gene2asn = defaultdict(lambda: "A")
with open('viburni/B_genes.tsv', 'r') as gene_file:
    for line in gene_file:
        line = line.split('\t')
        gene = line[0]
        asn = line[1].rstrip('\n')
        gene2asn[gene] = asn

for seq in SeqIO.parse('viburni/transcripts.fa', "fasta"):
    gene,iso = seq.name.split('.')
    if gene2asn[gene] != 'A' and iso == 't1':
        print(">Pvib_" + gene + '_' + gene2asn[gene])
        print(seq.seq)

with open("viburni/B_proteins.faa", 'w') as B_prot_file:
    for seq in SeqIO.parse('viburni/proteins.faa', "fasta"):
        gene,iso = seq.name.split('.')
        if gene2asn[gene] != 'A' and iso == 't1':
            B_prot_file.write(">Pvib_" + gene + '_' + gene2asn[gene] + '\n')
            B_prot_file.write(str(seq.seq) + "\n")
