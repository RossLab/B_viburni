import argparse
from operator import itemgetter

parser = argparse.ArgumentParser(description='Annotate blast/diamond')
parser.add_argument('blastfile', help='blast file (formated as 6 outfmt qseqid staxids)')
args = parser.parse_args()

id2sp = dict()

# blastout/
with open('data/id2names.tsv') as tx_file:
  for line in tx_file:
    id, sp = itemgetter(0,2)(line.split('\t'))
    id2sp[id] = sp

# blastout/B_genes.blast.out
with open(args.blastfile, 'r') as blast_file:
  prev_gene = ''
  hits = []
  for line in blast_file:
    gene, txid = itemgetter(0,1)(line.split('\t'))
    if prev_gene != gene and prev_gene != '':
      hit_string = ','.join(hits)
      hits = []
      print('{}\t{}'.format(gene, hit_string))
    prev_gene = gene
    try:
        hits.append(id2sp[txid])
    except KeyError:
        hits.append('Unknown')
  print('{}\t{}'.format(gene, hit_string))
