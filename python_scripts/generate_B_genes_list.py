#!/usr/bin/env python3

from collections import defaultdict
from sys import stdout

scf2asn = defaultdict(lambda :"A")

with open('output/scaffolds.final.assignment.csv') as scf_asn_file:
    for line in scf_asn_file:
        line_tab = line.rstrip('\n').split(',')
        asn = line_tab[-1][1:-1]
        if asn.startswith("B"):
            scf2asn[line_tab[0][1:-1]] = asn



with open('1_pacbio_assembly/7_braker/PVIB.BRAKER1/braker.gff3') as vcf_file:
    for line in vcf_file:
        if line.startswith('#'):
            continue
        annot_tab = line.rstrip('\n').split('\t')
        if annot_tab[2] != "gene":
            continue
        asn = scf2asn[annot_tab[0]]
        if asn.startswith("B"):
            gene = annot_tab[-1].split(';')[0].split('=')[1]
            stdout.write("{}\t{}\n".format(gene, asn))
