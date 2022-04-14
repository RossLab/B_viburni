#!/usr/bin/env python

import sys
from collections import defaultdict

scf2asn = defaultdict(lambda: 'A')

with open('B1_scaffolds.list', 'r') as f:
    for line in f:
        scf2asn[line.rstrip('\n')] = 'B1'

with open('B2_scaffolds.list', 'r') as f:
    for line in f:
        scf2asn[line.rstrip('\n')] = 'B2'

with open('B3_scaffolds.list', 'r') as f:
    for line in f:
        scf2asn[line.rstrip('\n')] = 'B3'

for line in sys.stdin:
    scf, pos, cov1, cov2 = line.rstrip('\n').split('\t')
    asn = scf2asn[scf]
    if asn != 'A':
        sys.stdout.write('{}\t{}\t{}\t{}\t{}\n'.format(scf, pos, cov1, cov2, asn))
