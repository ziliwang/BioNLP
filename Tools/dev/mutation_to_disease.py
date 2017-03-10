import Base
import re
from collections import defaultdict
import numpy as np
import dependency_tree
import normalization
__all__ = ['topic_anno']


def mutation_disease(artilce_ls):
    dnormor = normalization.disease_normor('../data/CTD_diseases.tsv')
    qdic = defaultdict(
            lambda: defaultdict(
                lambda: defaultdict(dict)))
    for artilce in artilce_ls:
        pmid = artilce.pmid
        relation_in_tk(pmid, artilce.tk, 'mutation', 'disease', qdic)
        relation_in_part(pmid, artilce.abstract, 'mutation', 'disease', qdic)
        relation_in_part(pmid, artilce.results, 'mutation', 'disease', qdic)
    qdic = default_to_regular(qdic)
    qdic = dependency_tree.dependency_tree_shortest_path(qdic)
    with open('tmp', 'w', encoding='utf-8') as f:
        for pmid in qdic.keys():
            for A in qdic[pmid].keys():
                for B in pairs(qdic[pmid][A]):
                # tk = [qdic[pmid][A][B]['tk'][2]] if 'tk' in qdic[pmid][A][B].keys() else []
                # for p in [i[2] for i in qdic[pmid][A][B]['evidences']] + tk:
                #     print('{0}\t{1}\t{2}'.format(A, dnormor.norm(B)['name'], p))
                    try:
                        for d, s in [i[1:] for i in qdic[pmid][A][B]['evidences']]:
                            f.write('{4}\t{0}\t{1}\t{2}\t{3}\n'.format(A, dnormor.norm(B)['name'], d, s, pmid))
                    except:
                        pass


def relation_in_tk(pmid, tk, A_type, B_type, qdic):
    for A in tk.bioe:
        if A.type == A_type:
            for B in tk.bioe:
                if B.type == B_type:
                    qdic[pmid][A.norm][B.id]['tk'] = (A.string,
                                                      B.string,
                                                      tk.text)


def relation_in_part(pmid, part, A_type, B_type, qdic):
    for A in part.bioe:
        if A.type == A_type:
            for a, b, s in same_stentenc(part, A, B_type):
                if 'evidences' not in qdic[pmid][a.norm][b.id].keys():
                    qdic[pmid][a.norm][b.id]['evidences'] = []
                if s not in [tripe[2]
                             for tripe in
                             qdic[pmid][a.norm][b.id]['evidences']]:
                    qdic[pmid][a.norm][b.id]['evidences'].append(
                            (a.string, b.string, s))


def pairs(A):
    shortest = 5
    tk = False
    for i in A.keys():
        if 'tk' in A[i].keys():
            yield(i)
        if 'dt_path_lens' in A[i] and min(A[i]['dt_path_lens']) < shortest:
            shortest = min(A[i]['dt_path_lens'])
    for i in A.keys():
        if 'dt_path_lens' in A[i] and min(A[i]['dt_path_lens']) == shortest:
            yield(i)


def same_stentenc(part, target, o_type):
    o_ls = [i for i in part.bioe if i.type == o_type]
    for i in sent_boundary(part.text):
        if i[0] <= target.start and i[1] >= target.end:
            # print('get sent: {0} - {1}'.format(i[0], i[1]))
            for j in o_ls:
                # print('q: {0}, {1}'.format(j.start, j.end))
                if i[0] <= j.start and i[1] >= j.end:
                    yield (target, j, part.text[i[0]:i[1]])


def sent_boundary(sent):
    reg = '[^\s][\w,\),\]][!.?]\s+[A-Z]'
    start = 0
    for i in re.finditer(reg, sent):
        end = i.span()[0] + 3
        yield(start, end)
        start = i.span()[1] - 1
    '''fixed bug: Name abbr, i. e.: A is provided by M. Korc.'''


def default_to_regular(d):
    if isinstance(d, defaultdict):
        d = {k: default_to_regular(v) for k, v in d.items()}
    return d


def topic_model(artilce_ls):
    for a in artilce_ls:
        tk_d, ab_d, re_d = False, False, False
        tk_d = [i for i in a.tk.bioe if i.type == 'disease']
        if a.results:
            re_d = [i for i in a.results.bioe if i.type == 'disease']
        if a.abstract:
            ab_d = [i for i in a.abstract.bioe if i.type == 'disease']
        tmp = tk_d if tk_d else ab_d if ab_d else re_d
        if not tmp:
            continue
        for i in get_topic(tmp, 2):
            yield((a.pmid, i.id, i.norm['name']))


def get_topic(obs, fold):
    '''return topic ob, remove synonymy'''
    tmp = dict()
    for ob in obs:
        '''synonymy?'''
        is_assign = False
        if tmp:
            for i in tmp.keys():
                if __nodes_cont(tmp[i]['ob'].norm['nodes'], ob.norm['nodes']):
                    tmp[i]['times'] += 1
                    is_assign = True
                    break
                elif __nodes_cont(ob.norm['nodes'], tmp[i]['ob'].norm['nodes']):
                    c_t = tmp[i]['times'] + 1
                    del tmp[i]
                    tmp[ob.id] = {'times': c_t, 'ob': ob}
                    is_assign = True
                    break
            if not is_assign:
                tmp[ob.id] = {'times': 1, 'ob': ob}
        else:
            tmp[ob.id] = {'times': 1, 'ob': ob}
    ap_max = 1
    for i in tmp.keys():
        if tmp[i]['times'] > ap_max:
            ap_max = tmp[i]['times']
    for i in tmp.keys():
        if tmp[i]['times'] >= float(ap_max/fold):
            yield tmp[i]['ob']


def __nodes_cont(A, B):
    '''A contains B'''
    for i in A.split('|'):
        for ii in B.split('|'):
            if re.match(re.escape(i.split('/')[0]), ii):
                return True
    return False


def topic_anno(artilce_ls, ann_file):
    d = {}
    for i in topic_model(artilce_ls):
        if i[0] in d.keys():
            d[i[0]].append(i[1:])
        else:
            d[i[0]] = [i[1:]]
    with open(ann_file, 'r') as f:
        with open('disease_' + ann_file, 'w') as ff:
            for i in f.readlines():
                items = i.split('\t')
                if int(items[0]) in d.keys():
                    for disease in d[int(items[0])]:
                        ff.write('\t'.join([items[0]] + list(disease) + items[1:]))
                else:
                    ff.write('\t'.join([items[0]] + ['None']*2 + items[1:]))
