#!/usr/bin/env python3
import urllib.request
import re
import time
import argparse
import os

FUNCTION = '''
Get pubmed xml by pmid based on NCBI API.
'''
INFO = '''Copyright wzlnot@gmail.com All Rights Reserved. \
Licensed under the MIT License'''


def parse_args():
    parser = argparse.ArgumentParser(
                 description=FUNCTION, epilog=INFO,
                 formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('pmid_file', type=str,
                              help='A file contained pmid in one colum.')
    parser.add_argument('-n', '--number', metavar='N', type=int,
                              default=100,
                              help='batch number per request.[100]')
    parser.add_argument('-s', '--start', metavar='N', type=int, default=0,
                              help='Continue to download, start batch number\
                                   , start=item%%batch.[0]')
    parser.add_argument('-p', '--path', metavar='path', type=str,
                        default='./download',
                              help='the path to store XML files.[./download]')
    args = parser.parse_args()
    return args


def efetch(pmids, num, path):
    '''efetch(pmids[list], num_tag[int], path[str]). num_tag is the btach \
    number for indentify;
    '''
    url = (
      'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi' +
      '?db=pubmed&id=' + ','.join(pmids) + '&retmode=xml')
    times = 0
    while True:
        try:
            html = urllib.request.urlopen(url).read().decode(
                'unicode_escape', errors='ignore')
        except:
            times += 1
        # give up when get the response or try more then 5 times
        if html or times > 5:
            break
    with open(os.path.join(path, str(num) + '.xml'), 'w') as f:
        f.write(html)


def main():
    args = parse_args()
    os.makedirs(args.path, exist_ok=True)  # build path if not exist
    with open(args.pmid_file, 'r') as f:
        pmids = [i.strip('\n') for i in f.readlines() if re.match(r'\d+', i)]
    pmids = sorted(set(pmids))  # remove repeat pmids
    iter_num = 0  # for download control, item = iter_num * batch_num
    print('start')
    batch_bum = args.number
    start = args.start
    while len(pmids) > batch_bum:
        post_pmids = [pmids.pop() for i in range(batch_bum)]  # post pmids
        if iter_num >= start:  # control
            print('Downloading {0}~{1}'.format(
                        iter_num*batch_bum, (iter_num + 1)*batch_bum))
            efetch(post_pmids, iter_num, args.path)
        iter_num += 1
    # last
    print('Downloading last {0}'.format(batch_bum))
    efetch(pmids, iter_num, args.path)


if __name__ == '__main__':
    main()
