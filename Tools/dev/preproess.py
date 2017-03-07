from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
import re
import PubMedDB
import os
import sys
import pickle
__all__ = ["db", "pdf_ls_to_input"]


def NCBI_input_text_clean(string):
    # replace '|' to ' '
    # replace '\t' to ' '
    # repalce '\n' to ' '
    string = "".join([i for i in string if ord(i) < 128])
    return re.sub(r'\||\t|\n', ' ', string)


class db():
    '''operate db. extract info from database'''

    def __init__(self, db):
        engine = create_engine('postgresql://parser:parser@localhost/' + db)
        self.DBSession = sessionmaker(bind=engine)

    def get_text(self, output_prefix):
        input_for_tmvar = output_prefix + '.tmvarI'
        input_for_GNP = output_prefix + '.GNormPlusI'
        input_for_DN = output_prefix + '.DNormI'
        input_for_EMU = output_prefix + '.EMU'
        qurey_session = self.DBSession()
        articles = qurey_session.query(PubMedDB.Citation).all()
        # title and keywords
        with open(input_for_tmvar + '.tk', 'w') as tmvar:
            with open(input_for_DN + '.tk', 'w') as dn:
                with open(input_for_EMU + '.tk', 'w') as emu:
                    with open(input_for_GNP + '.tk', 'w') as gnp:
                        for article in articles:
                            keywords = qurey_session.query(PubMedDB.Keyword
                                                           ).filter_by(
                                                    fk_pmid=article.pmid)
                            keyword_ls = [i.keyword for i in keywords if i]
                            tmp_str = article.article_title
                            if len(keyword_ls) > 0:
                                tmp_str += ', '.join(keyword_ls)
                            pmid = article.pmid
                            tmp_str = text_clean(tmp_str, 'pubmed')
                            tmvar.write('{0}|a|{1}\n\n'.format(pmid, tmp_str))
                            dn.write('{0}\t{1}\n'.format(pmid, tmp_str))
                            emu.write('''{0}\t{0}-{1}\t"{2}"\n'''.format(
                                        pmid, 'tk', tmp_str))
                            gnp.write('{0}|a|{1}\n\n'.format(pmid, tmp_str))
        with open(input_for_tmvar + '.abstract', 'w') as tmvar:
            with open(input_for_DN + '.abstract', 'w') as dn:
                with open(input_for_EMU + '.abstract', 'w') as emu:
                    with open(input_for_GNP + '.abstract', 'w') as gnp:
                        for article in articles:
                            pmid = article.pmid
                            try:
                                abstract = qurey_session.query(
                                    PubMedDB.Abstract).filter_by(fk_pmid=pmid)
                                tmp_str = abstract[0].abstract_text
                            except BaseException as e:
                                '''some article have no abstract'''
                                print('{0} have no abstract: {1}'.format(
                                      pmid, e))
                                continue
                            tmp_str = text_clean(tmp_str, 'pubmed')
                            tmvar.write('{0}|a|{1}\n\n'.format(pmid, tmp_str))
                            dn.write('{0}\t{1}\n'.format(pmid, tmp_str))
                            emu.write('''{0}\t{0}-{1}\t"{2}"\n'''.format(
                                        pmid, 'abstract', tmp_str))
                            gnp.write('{0}|a|{1}\n\n'.format(pmid, tmp_str))
        qurey_session.close()


def pdf_ls_to_input(ls, output_prefix):
    '''read text and write into input format'''
    input_for_tmvar = output_prefix + '.tmvarI'
    input_for_GNP = output_prefix + '.GNormPlusI'  # dir
    input_for_DN = output_prefix + '.DNormI'
    input_for_EMU = output_prefix + '.EMU'
    if not os.path.exists(input_for_GNP):
        os.makedirs(input_for_GNP)
    with open(input_for_tmvar + '.results', 'w') as tmvar:
        with open(input_for_DN + '.results', 'w') as dn:
            with open(input_for_EMU + '.results', 'w') as emu:
                with open(input_for_GNP + '.results', 'w') as gnp:
                    for pmid, text in ls:
                        if len(text) == 0:
                            continue
                        tmp_str = text_clean(text, 'pdf')
                        tmvar.write('{0}|a|{1}\n\n'.format(pmid, tmp_str))
                        dn.write('{0}\t{1}\n'.format(pmid, tmp_str))
                        emu.write('''{0}\t{0}-{1}\t"{2}"\n'''.format(
                                                    pmid, 'results', tmp_str))
                        gnp.write('{0}|a|{1}\n\n'.format(pmid, tmp_str))


def __load_aa_table():
    try:
        with open('../data/AA.table.pkl', 'rb') as f:
            return pickle.load(f)
    except FileNotFoundError as e:
            print('aa.table missing: {0}'.format(e))
            sys.exit(0)


def text_clean(text, at):
    '''
    # replace '|' to ' '
    # replace '\t' to ' '
    # repalce '\n' to ' '
    '''
    text = "".join([i for i in text if ord(i) < 128])
    text = re.sub(r'\||\t|\n', ' ', text)
    aa_tuple = __load_aa_table()
    if at == 'pdf':
        '''fix '> to . ' issue'''
        pattern = re.compile('[ATCG]\s?(\.)\s?[ATCG]')
        text = pattern.sub('\g<1>>\g<2>', text)
    else:
        return text
