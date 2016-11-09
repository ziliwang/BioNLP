#!/usr/bin/env python3
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
import PubMedDB
import argparse
import json

FUNCTION = '''
Database to table.
'''
INFO = '''Copyright wzlnot@gmail.com All Rights Reserved. \
Licensed under the MIT License'''


def parse_args():
    parser = argparse.ArgumentParser(
                 description=FUNCTION, epilog=INFO,
                 formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-d', '--database', type=str,
                              help='postgrel database', required=True)
    args = parser.parse_args()
    return args


def query_pmid(DBSession):
    query_session = DBSession()
    pmids = [str(i[0]) for i in query_session.query(
                                          PubMedDB.Citation.pmid).all()]
    output = {}
    for pmid in pmids:
        s = query_session.query(PubMedDB.Journal).filter(
                       PubMedDB.Journal.fk_pmid == pmid).one()
        try:
            ab = query_session.query(
                        PubMedDB.Abstract).filter(
                                  PubMedDB.Abstract.fk_pmid == pmid).one()
            ab_text = ab.abstract_text
        except:
            ab_text = None
        tit = query_session.query(PubMedDB.Citation.article_title).filter(
                                  PubMedDB.Citation.pmid == pmid).one()
        try:
            auth = query_session.query(PubMedDB.Author).filter(
                                  PubMedDB.Author.fk_pmid == pmid).all()
            authlist = [' '.join([i.last_name, i.fore_name]) for i in auth]
        except:
            authlist = None
        output[pmid] = {'title': tit[0],
                        'j_issn': s.issn,
                        'j_volume': s.volume,
                        'j_issue': s.issue,
                        'j_pud_date_year': s.pub_date_year,
                        'j_pub_date_month': s.pub_date_month,
                        'j_pub_date_day': s.pub_date_day,
                        'j_title': s.title,
                        'j_iso_abbr': s.iso_abbreviation,
                        'abstract': ab_text,
                        'author_list': authlist}
    query_session.close
    return output


if __name__ == "__main__":
    arg = parse_args()
    db = arg.database
    engine = create_engine('postgresql://parser:parser@localhost/' + db)
    DBSession = sessionmaker(bind=engine)
    data = query_pmid(DBSession)
    with open('pmids.json', 'w') as f:
        json.dump(data, f)
