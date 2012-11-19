# -*- coding: utf-8 -*-
# <nbformat>2</nbformat>

# <codecell>

file='test.vcf'
fh=open(file, 'r')

for line in fh:
    if '#' in line:
        print line.strip()
        continue
    fields=line.strip().split('\t')
    affygstring=fields[9]
    if affygstring != './.':
        affygstring+=':.'
    fields[9]=affygstring
    outstring="\t".join( fields )
    print outstring

# <codecell>


