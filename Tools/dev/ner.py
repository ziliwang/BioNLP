import re
import os
import sys
import Base
import pickle
import normalization
from collections import defaultdict
__all__ = ['output_parser']


def output_parser(outdir, prefix):
    '''output: article obj list'''
    tm_tk = outdir + '/' + prefix + '.tmvarI.tk.PubTator'
    tm_ab = outdir + '/' + prefix + '.tmvarI.abstract.PubTator'
    tm_re = outdir + '/' + prefix + '.tmvarI.results.PubTator'
    dn_tk = outdir + '/' + prefix + '.DNormO.tk'
    dn_ab = outdir + '/' + prefix + '.DNormO.abstract'
    dn_re = outdir + '/' + prefix + '.DNormO.results'
    gn_dir = outdir + '/' + prefix + '.GNormPlusO'
    emu_tk = outdir + '/EMU_1.19_HUGO_' + prefix + '.EMU.tk'
    emu_ab = outdir + '/EMU_1.19_HUGO_' + prefix + '.EMU.abstract'
    emu_re = outdir + '/EMU_1.19_HUGO_' + prefix + '.EMU.results'
    articles = []
    for pmid, string, mutation_entry_ls in var_parser(tm_tk, emu_tk):
        article_obj = Base.Article(pmid)  # create article_obj
        tk_obj = Base.TK(string)  # create tk_obj
        for e in mutation_entry_ls:
            tk_obj.add_e(e)  # add entry
        article_obj.add_tk(tk_obj)
        articles.append(article_obj)
    # abstract
    for pmid, string, mutation_entry_ls in var_parser(tm_ab, emu_ab):
        article_obj = found_article(pmid, articles)  # found article_obj
        abstract_obj = Base.Abstract(string)  # create abstract_obj
        for e in mutation_entry_ls:
            abstract_obj.add_e(e)  # add entry
        article_obj.add_abstract(abstract_obj)
    # results
    has_results = True
    try:
        for pmid, string, mutation_entry_ls in var_parser(tm_re, emu_ab):
            article_obj = found_article(pmid, articles)  # found article_obj
            results_obj = Base.Results(string)  # create tk_obj
            for e in mutation_entry_ls:
                results_obj.add_e(e)  # add entry
            article_obj.add_results(results_obj)
    except FileNotFoundError as e:
        has_results = False
        print('no results: {0}'.format(e))
    # DNorm
    for pmid, entry in DNorm_parser(dn_tk):
        article_obj = found_article(pmid, articles)
        article_obj.tk.add_e(entry)
    for pmid, entry in DNorm_parser(dn_ab):
        article_obj = found_article(pmid, articles)
        article_obj.abstract.add_e(entry)
    if has_results:
        for pmid, entry in DNorm_parser(dn_re):
            article_obj = found_article(pmid, articles)
            article_obj.results.add_e(entry)
    # GNormPlus
    for pmid, entry, fn in GNormPlus_parser(gn_dir):
        article_obj = found_article(pmid, articles)
        if re.search('tk', fn):
            article_obj.tk.add_e(entry)
        elif re.search('abstract', fn):
            article_obj.abstract.add_e(entry)
        else:
            article_obj.results.add_e(entry)

    '''sort entry'''
    for i in articles:
        if i.abstract:  # possible no abstract
            i.abstract.sort_entry()
        if i.results:  # possible no results
            i.results.sort_entry()
            pass
    return articles


def found_article(pmid, articles):
    for article in articles:
        if pmid == article.pmid:
            return article
    raise ValueError('{0} have no title'.format(pmid))


def var_parser(outputfile, emu_out):
    '''output: pmid, string, mutation_entry'''
    '''mutation integration'''
    with open(outputfile, 'r') as f:
        mutation_results = [i.split('\n') for i in f.read().split('\n\n') if i]
    emu = defaultdict(lambda: defaultdict(str))
    for pmid, string, norm in emu_to_tmvar(emu_out):
        emu[int(pmid)][string] = norm
    for article_result in mutation_results:
        text = ''
        pmid = False
        mutation_entry_ls = []
        for record in article_result:
            is_text = re.match(r'(\d+)\|a\|(.+)$', record)
            if is_text:
                # id, string
                text = is_text.group(2)
                pmid = int(is_text.group(1))
            else:
                rpmid, start, end, string, Mtype, norm = record.split('\t')
                start = int(start)
                end = int(end)
                if re.search('[A-T]\d+[A-T]', string.upper()):
                    if re.match('c', norm):
                        norm = norm.replace('c', 'p', 1)
                try:
                    my_nm = normalization.mutation_string_normalizition(
                                string, norm)
                except ValueError as e:
                    print('{0}\t{1:10}'.format(e, 'skip'))
                    continue
                mutation_entry_ls.append(
                    Base.BioEntry(start, end, 'mutation', string, norm, my_nm))
        for string in emu[pmid].keys():
            emu_nm = emu[pmid][string]
            is_exist = 0
            for tmvar_mut in mutation_entry_ls:
                tm_nm = tmvar_mut.id
                if emu_nm in tm_nm or tm_nm in emu_nm:  # 重复
                    is_exist = 1
                    break
            if not is_exist:
                '''add mutation entry'''
                for start, end in __emu_get_pos(string, text):
                    my_nm = normalization.mutation_string_normalizition(
                                    string, emu_nm)
                    mutation_entry_ls.append(
                        Base.BioEntry(start, end, 'mutation', string, emu_nm,
                                      my_nm)
                        )
        yield(pmid, text, mutation_entry_ls)


def DNorm_parser(outputfile):
    '''output: pmid, disease_entry'''
    dnor = normalization.disease_normor('../data/CTD_diseases.tsv')
    with open(outputfile, 'r') as f:
        dnorm_results = [i.split('\t') for i in f.read().split('\n') if i]
    for i in dnorm_results:
        if len(i) == 5:
            pmid, start, end, string, id = i
        else:
            continue
            # pmid, start, end, string = i
            # id = 'null'
        pmid = int(pmid)
        start = int(start)
        end = int(end)
        yield(pmid, Base.BioEntry(start, end, 'disease', string, id,
                                  dnor.norm(id)))


def GNormPlus_parser(outputdir):
    '''output: pmid, gene_entry, article_part'''
    gnor = normalization.gene_normor('../data/hgnc_complete_set.txt')
    for root, dirs, files in os.walk(outputdir):
        for fn in files:
            with open(os.path.join(root, fn), 'r') as f:
                records = [i for i in f.read().split('\n') if i]
                for r in records:
                    is_text = re.match(r'\d+\|a\|(.+)$', r)
                    if is_text:
                        continue
                    else:
                        pmid, start, end, string, type, ids = r.split('\t')
                        if not re.match('gene|protein', type, re.IGNORECASE):
                            continue
                        pmid = int(pmid)
                        start = int(start)
                        end = int(end)
                        for did in ids.split(';'):
                            '''one tring map to mutil-gene'''
                            try:
                                norm = gnor.norm(did)
                            except KeyError as e:
                                print('gene {0} skip: no gene in hgnc'.format(
                                        string))
                                continue
                            yield(pmid,
                                  Base.BioEntry(start, end, type, string, id,
                                                norm),
                                  fn)


def __load_aa_tuple(path):
    with open(path, 'rb') as f:
        aa_tuple = pickle.load(f)
    return aa_tuple


def get_aa_abbr(aa):
    aa_tuple = __load_aa_tuple('../data/AA.table.pkl')
    for abbr, tripe in aa_tuple:
        if aa.upper() == tripe.upper() or aa.upper() == abbr.upper():
            return(abbr)
    raise ValueError('illegal amino acid symble: {0}'.format(aa))


def __get_mut_type(mut_string, emu_mu_type):
    '''INS DEL INDEL classifier'''
    if emu_mu_type:
        return emu_mu_type
    for i in ['INS', 'DEL', 'INDEL']:
        if re.search(i, mut_string.upper()):
            return i
    raise ValueError('can not recognize mutation type: {0}\t{1}'.format(
                        mut_string, emu_mu_type))
    return False


def __emu_level(level_set, mut_string):
    '''judge level: DNA or AAs or MiTHO'''
    if level_set == 'DNA' or level_set == 'DNA;RNA':
        return('c')
    elif level_set == 'PROTEIN':
        return('p')
    elif level_set == 'MITHO':
        return('m')
    else:
        if not re.search('[A-T]\d+[A-T]', mut_string):
            return('c')
        else:
            return('p')


def __emu_ref_mut(string, mut_type, level, wtaa, mtaa):
    '''ref ale check'''
    if mut_type == 'MISSENSE':
        if level in ['c', 'm']:
            if wtaa.upper() in 'ATCG' and mtaa.upper() in 'ATCG':
                return(wtaa, mtaa)
            else:
                raise ValueError('DNA level with err symbol'
                                 ': {0}, {1}'.format(wtaa, mtaa))
                return(False, False)
        else:
            return(get_aa_abbr(wtaa), get_aa_abbr(mtaa))
    else:
        ref = __get_mut_type(string, wtaa)
        if level in ['c', 'm']:
            if set([i.upper() for i in mtaa]).issubset(
                set([i.upper() for i in 'atcg'])) or (
                            re.match('\d+$', mtaa)
                            ):
                return (ref, mtaa)
            else:
                raise ValueError('alle detect failed: {0}'.format(mtaa))
        else:
            if re.match('\d+$', mtaa):
                return (ref, mtaa)  # number
            else:
                return (ref, get_aa_abbr)  # abbr


def __emu_pos(emu_pos, string):
    '''pos check'''
    if re.match('^\d+$', emu_pos):
        return(emu_pos)
    elif re.match(r'^[\d\-\+_\s\t]+$', emu_pos):  # legal chr
        if re.match(r'\d+\-\-\d+$', emu_pos):
            return(emu_pos.replace('--', '_'))
        tmp = re.sub('[^_\-\+\d]', '', string)  # remove space
        if re.search(re.escape(emu_pos), tmp):
            return emu_pos
        else:
            emu_legal_pos = re.sub('\-+', '_', emu_pos, 1)
            if re.search(re.escape(emu_legal_pos), tmp):
                return emu_legal_pos
            else:
                raise ValueError('position detected failed: {0}'.format(
                                                                    emu_pos))
    else:
        raise ValueError('Illegal position: {0}'.format(emu_pos))


def __emu_tm_string(mut_type, level, ref, pos, ale):
    if mut_type == 'MISSENSE':
        return('{0}|SUB|{1}|{2}|{3}'.format(level, ref, pos, ale))
    else:
        return('{0}|{1}|{2}|{3}'.format(level, ref, pos, ale))


def __emu_porcessor(mut_type, level_set, string, wtaa, mtaa, raw_pos):
    '''emu misssense processor'''
    level = __emu_level(level_set, string)
    ref, ale = __emu_ref_mut(string, mut_type, level, wtaa, mtaa)
    pos = __emu_pos(raw_pos, string)
    return __emu_tm_string(mut_type, level, ref, pos, ale)


def emu_to_tmvar(emu_out_put_file):
    with open(emu_out_put_file, 'r') as f:
        f.readline()  # header
        raw = [i.split('\t') for i in f.read().split('\n') if i]
    # process
    for (pmid, organism, mut_pat1, pos_patt,
         mutation_type, wtaa, mtaa, pos, genes, gtype) in raw:
        try:
            yield(pmid, mut_pat1, __emu_porcessor(mutation_type, gtype,
                                                  mut_pat1, wtaa, mtaa, pos))
        except ValueError as e:
            print(pmid, e)
            continue


def __emu_get_pos(string, text):
    reg = re.escape(string)
    for i in re.finditer(reg, text):
        yield(i.span())
