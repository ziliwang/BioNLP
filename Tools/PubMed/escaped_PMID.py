#!/usr/bin/env python3
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
import PubMedDB
import argparse
import re

FUNCTION = '''
find escaped pmid.
'''
INFO = '''Copyright wzlnot@gmail.com All Rights Reserved. \
Licensed under the MIT License'''


def parse_args():
    parser = argparse.ArgumentParser(
                 description=FUNCTION, epilog=INFO,
                 formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-d', '--database', type=str,
                              help='postgrel database', required=True)
    parser.add_argument('-f', '--file', type=str, help='total pmids',
                              required=True)
    args = parser.parse_args()
    return args


def query_pmid(DBSession):
    query_session = DBSession()
    pmid = [str(i[0]) for i in query_session.query(PubMedDB.Citation.pmid).all()]
    query_session.close
    return pmid

if __name__ == "__main__":
    arg = parse_args()
    db = arg.database
    engine = create_engine('postgresql://parser:parser@localhost/' + db)
    DBSession = sessionmaker(bind=engine)
    pmids_in = query_pmid(DBSession)
    with open(arg.file, 'r') as f:
        for line in f.readlines():
            if re.match(r'\d+', line):
                pmid = line.strip('\r\n')
                if pmid not in pmids_in:
                    print(pmid)
