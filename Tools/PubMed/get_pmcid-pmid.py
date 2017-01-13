#!/usr/bin/env python3
import urllib.request
import argparse
import json

FUNCTION = '''
query pmid by pmcid, output pmid-pmcid pairs 
'''
INFO = '''Copyright wzlnot@gmail.com All Rights Reserved. \
Licensed under the MIT License'''

TRYTIMES = 10


def parse_args():
    parser = argparse.ArgumentParser(
                 description=FUNCTION, epilog=INFO,
                 formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-l', '--list', type=str,
                              help='pmcid list file', required=True)
    args = parser.parse_args()
    return args


def query_pmcid(query_pmcid):
    prefix = 'https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?ids='
    pmcid = ','.join(query_pmcid)
    suffix = '&format=json'
    url = prefix + pmcid + suffix
    html = None
    i = 0
    while not html and i < TRYTIMES:
        i += 1
        try:
            html = urllib.request.urlopen(url).read().decode(
                                         'utf-8', errors='ignore')
        except BaseException as e:
            print('Err:{0}.try again({1})'.format(e, i))
    result = json.loads(html)
    record_pairs = []
    for record in result['records']:
        if 'pmid' in record.keys():
            record_pairs.append((record['pmcid'], record['pmid']))
    return record_pairs


if __name__ == "__main__":
    arg = parse_args()
    with open(arg.list, 'r') as f:
        pmid = [i for i in f.read().split('\n') if i]
    post_pmid = []
    id_pairs = []
    while len(pmid) > 200:
        post_pmid = [pmid.pop() for i in range(200)]  # 200 items
        id_pairs += query_pmcid(post_pmid)
    id_pairs += query_pmcid(pmid)  # last
    with open('pmcid-pmid.pairs', 'w') as f:
        for pmc, pm in id_pairs:
            f.write('{0}\t{1}\n'.format(pmc, pm))
