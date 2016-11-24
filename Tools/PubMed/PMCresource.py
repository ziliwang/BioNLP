#!/usr/bin/env python3
import queue
import threading
import urllib.request as urlr
import xml.etree.ElementTree as etree
import PubMedDB
import argparse
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
import os

FUNCTION = '''
Download free full text from PMC.
'''
INFO = '''Copyright wzlnot@gmail.com All Rights Reserved. \
Licensed under the MIT License'''
TRYTIMES = 10


def cbk(a, b, c):
    per = 100.0 * a * b / c
    if per > 100:
        per = 100
    print('%.2f%%' % per)


def query_pmcid(DBSession):
    query_session = DBSession()
    pmcids = [str(i[0]) for i in query_session.query(
                                              PubMedDB.PM_to_PMC.pmcid).all()]
    query_session.close
    return pmcids


class downloader(object):
    def __init__(self, pmcid, store_path):
        self.pmcid = pmcid
        self.store_path = store_path

    def work(self):
        # url: https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?id=PMC13901
        baseurl = 'https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?id='
        resp = urlr.urlopen(baseurl + self.pmcid).read().decode(
                                             'unicode_escape', errors='ignore')
        tree = etree.fromstring(resp)
        try:
            durl = tree[2][0][0].attrib['href']
            urlr.urlretrieve(
                 durl, self.store_path + '/{0}.tar.gz'.format(self.pmcid))
            print("Download: {0} finished".format(self.pmcid))
        except:
            print('{0} not in NCBI oa_package'.format(self.pmcid))


class downloadTASK(threading.Thread):
    def __init__(self, myqueue):
        self.queue = myqueue
        super(downloadTASK, self).__init__()

    def run(self):
        while True:
            if self.queue.qsize() > 0:
                self.queue.get().work()
            else:
                break


def parse_args():
    parser = argparse.ArgumentParser(
                 description=FUNCTION, epilog=INFO,
                 formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-d', '--database', type=str,
                              help='postgrel database', required=True)
    parser.add_argument('-P', '--PATH', type=str, required=True,
                              help='PATH to store the full text')
    parser.add_argument('-p', '--process', type=int, default=4,
                              help='thread number')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()
    if not os.path.exists(args.PATH):
        os.mkdir(args.PATH)
    db = args.database
    engine = create_engine('postgresql://parser:parser@localhost/' + db)
    DBSession = sessionmaker(bind=engine)
    myqueue = queue.Queue(0)
    for i in query_pmcid(DBSession):
        s_download = downloader(i, args.PATH)
        myqueue.put_nowait(s_download)
    threads = []
    for i in range(args.process):
        thread = downloadTASK(myqueue)
#        thread.start()
        threads.append(thread)
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
