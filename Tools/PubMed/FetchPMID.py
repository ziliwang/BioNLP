#!/usr/bin/env python3
import urllib.request
import re
import time


def efetch(pmids, num):
    url = (
      'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi' +
      '?db=pubmed&id=' + ','.join(pmids) + '&retmode=xml')

    html = urllib.request.urlopen(url).read().decode(
      'unicode_escape', errors='ignore')
    with open(str(num) + '.xml', 'w') as f:
        f.write(html)

pmidlist = 'esc_pmids.txt'
pmids = []
with open(pmidlist, 'r') as f:
    for line in f.readlines():
        if re.match(r'\d+', line):
            pmid = line.strip('\n')
            if pmid not in pmids:
                pmids.append(pmid)
iter_num = 0
print('start')
batch_bum = 100
start = 0
while len(pmids) > batch_bum:
    post_pmids = [pmids.pop() for i in range(batch_bum)]  # 400 items
    if iter_num >= start:
        print('Downloading {0}~{1}'.format(iter_num*batch_bum,(iter_num +1)*batch_bum))
        efetch(post_pmids, iter_num)
    iter_num += 1
print('Downloading last {0}'.format(batch_bum))
efetch(pmids, iter_num)
