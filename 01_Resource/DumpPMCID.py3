#!/usr/bin/env python3
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
import PubMedDB
import urllib.request
import argparse
import json

FUNCTION = '''
Dump PMCID into database.
'''
INFO = '''Copyright wzlnot@gmail.com All Rights Reserved. \
Licensed under the MIT License'''

TRYTIMES = 10


def parse_args():
    parser = argparse.ArgumentParser(
                 description=FUNCTION, epilog=INFO,
                 formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-d', '--database', type=str,
                              help='postgrel database', required=True)
    args = parser.parse_args()
    return args


def query_pmcid(query_pmid):
    prefix = 'https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?ids='
    pmid = ','.join(query_pmid)
    suffix = '&format=json'
    url = prefix + pmid + suffix
    html = None
    i = 0
    while not html and i < TRYTIMES:
        i += 1
        try:
            html = urllib.request.urlopen(url).read().decode(
                                         'unicode_escape', errors='ignore')
        except BaseException as e:
            print('Err:{0}.try again({1})'.format(e, i))
    result = json.loads(html)
    record_pairs = []
    for record in result['records']:
        if 'pmcid' in record.keys():
            record_pairs.append((record['pmcid'], record['pmid']))
    return record_pairs


def query_pmid(DBSession):
    query_session = DBSession()
    pmid = [str(i[0]) for i in query_session.query(PubMedDB.Citation.pmid).all()]
    query_session.close
    return pmid


def pmcid_dump(DBSession, pairs):
    qury_session = DBSession()
    pmcids = [i[0] for i in qury_session.query(PubMedDB.PM_to_PMC.pmcid).all()]
    qury_session.close
    for pmid, pmcid in pairs:
        if pmcid in pmcids:
            print('{0} exist: pass')
        else:
            session = DBSession()
            record = PubMedDB.PM_to_PMC(pmcid=pmcid, fk_pmid=pmid,
                                        text_path='')
            session.add(record)
            session.commit()
            session.close


if __name__ == "__main__":
    arg = parse_args()
    db = arg.database
    engine = create_engine('postgresql://parser:parser@localhost/' + db)
    DBSession = sessionmaker(bind=engine)
    pmid = query_pmid(DBSession)
    post_pmid = []
    id_pairs = []
    while len(pmid) > 200:
        post_pmid = [pmid.pop() for i in range(200)]  # 200 items
        id_pairs += query_pmcid(post_pmid)
    id_pairs += query_pmcid(pmid)  # last
    # dump into database
    pmcid_dump(DBSession, id_pairs)
