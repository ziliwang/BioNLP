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
        tit = query_session.query(PubMedDB.Citation).filter(
                                  PubMedDB.Citation.pmid == pmid).one()
        try:
            auth = query_session.query(PubMedDB.Author).filter(
                                  PubMedDB.Author.fk_pmid == pmid).all()
            authlist = [' '.join([i.last_name, i.fore_name]) for i in auth]
        except:
            authlist = None
        output[pmid] = {'title': tit.article_title,  # 文章.标题(isNull=False)
                        'j_start_page': tit.start_page,  # 杂志.起始页码
                        'j_end_page': tit.end_page,  # 杂志.终止页码
                        'j_medline_pgn': tit.medline_pgn,  # 杂志.缩略页码
                        'j_issn': s.issn,  # 杂志.issn号
                        'j_volume': s.volume,  # 杂志.卷号
                        'j_issue': s.issue,  # 杂志.期号
                        'j_pud_date_year': s.pub_date_year,  # 杂志.出版年(isNull=False)
                        'j_pub_date_month': s.pub_date_month,  # 杂志.出版月
                        'j_pub_date_day': s.pub_date_day,  # 杂志.出版日
                        'j_title': s.title,  # 杂志.全名
                        'j_iso_abbr': s.iso_abbreviation,  # 杂志.缩略名
                        'abstract': ab_text,  # 文章.摘要
                        'author_list': authlist}  # 文章.作者（列表）
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
