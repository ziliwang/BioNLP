import link
import re
from collections import defaultdict
import csv
import subprocess
__all__ = ['rank_with_dt']


def rank_with_dt(artilce_ls, transVar_out_path):
    '''document components rank with dependency tree'''
    local_linkor = link.linker(transVar_out_path)
    dic = defaultdict(
            lambda: defaultdict(
                lambda: defaultdict(lambda: defaultdict(list))))
    fuzzy_dic = defaultdict(
                    lambda: defaultdict(lambda: defaultdict(list)))
    """dic[pmid][mutation][rank][gene] = evidence.list"""
    with open('gene-mutation.tsv', 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        for i in artilce_ls:
            if i.get_abstract_relations():
                for ii in i.get_abstract_relations():
                    mut_bioe = ii[0]
                    B_bioe = ii[1][1]
                    rank = ii[1][0]
                    evidence = ii[1][2]
                    """ii[0] is mutation, ii[1][1] is gene, ii[1][0] is rank,
                    ii[1][2] is sentence"""
                    if B_bioe.type == 'Gene':
                        """ok: rs, and in link,[futurn add mt g support]"""
                        links = local_linkor.parser(B_bioe.norm + ':' +
                                                    mut_bioe.norm)
                        if re.match('rs', mut_bioe.norm.lower()):
                            out = [i.pmid, mut_bioe.string, mut_bioe.norm, 'r0',
                                   mut_bioe.string, mut_bioe.norm,
                                   mut_bioe.norm, evidence]
                            writer.writerow(out)
                        elif links:
                            '''rank'''
                            dic[i.pmid][mut_bioe.norm][rank][B_bioe.norm].append(
                                (mut_bioe, B_bioe, evidence))
                        else:
                            pass
            if i.get_results_relations():
                for ii in i.get_results_relations():
                    mut_bioe = ii[0]
                    B_bioe = ii[1][1]
                    rank = ii[1][0]
                    evidence = ii[1][2]
                    """ii[0] is mutation, ii[1][1] is gene, ii[1][0] is rank,
                    ii[1][2] is sentence"""
                    if B_bioe.type == 'Gene':
                        """ok: rs, and in link,[futurn add mt g support]"""
                        links = local_linkor.parser(B_bioe.norm + ':' +
                                                    mut_bioe.norm)
                        if re.match('rs', mut_bioe.norm.lower()):
                            out = [i.pmid, mut_bioe.string, mut_bioe.norm, 'r0',
                                   mut_bioe.string, mut_bioe.norm,
                                   mut_bioe.norm, evidence]
                            writer.writerow(out)
                        elif links:
                            '''rank'''
                            dic[i.pmid][mut_bioe.norm][rank][B_bioe.norm].append(
                                (mut_bioe, B_bioe, evidence))
                        else:
                            pass
        for i in dic.keys():
            for j in dic[i].keys():
                for o in ['r0', 'r1', 'r2', 'r3', 'r4', 'r5', 'r6']:
                    if o in dic[i][j].keys():
                        if len(dic[i][j][o].keys()) > 1 and o == 'r0':
                            for p in dic[i][j][o].keys():
                                for q in dic[i][j][o][p]:  # tripe
                                    fuzzy_dic[i][j][p].append(q)
                        else:
                            for p in dic[i][j][o].keys():
                                for q in dic[i][j][o][p]:
                                    a, b, s = q
                                    '''pmid, rawA, normA, rank, rawB, normB,
                                    hgvs, evidences'''
                                    out = [i, a.string, a.norm, o, b.string,
                                           b.norm, '|'.join(local_linkor.parser(
                                             b.norm + ':' + a.norm
                                           )), s]
                                    writer.writerow(out)
                        break
        for pmid, i in get_DT_step(fuzzy_dic):
            for j in i:
                a, b, s = j
                out = [pmid, a.string, a.norm, o, b.string,
                       b.norm, '|'.join(local_linkor.parser(
                         b.norm + ':' + a.norm
                       )), s]
                writer.writerow(out)


def get_DT_step(qdic):
    ls = []
    # i, pmid; j, n_mut; o, n_gene; p, [a_o, b_o, record]
    for i in qdic.keys():
        for j in qdic[i].keys():
            for o in qdic[i][j].keys():
                for p in qdic[i][j][o]:
                    raw_text = p[2]
                    raw_mut_text = p[0].string
                    raw_gen_text = p[1].string
                    # set_trace()
                    raw_text = raw_text.replace(raw_mut_text, 'MMMoneMMM')
                    raw_text = raw_text.replace(raw_gen_text, 'GGGoneGGG')
                    raw_text = raw_text.strip('\n')  # remove new line
                    raw_text = raw_text.strip(' ')  # remove space in end
                    raw_text = raw_text.strip('.')  # remove .
                    raw_text = raw_text.replace('.', '')
                    raw_text = raw_text.replace('>', '')
                    raw_text = raw_text.replace('*', '')
                    raw_text = raw_text + '.'
                    ls.append([i, j, o, raw_text])
    with open('nlpin.tmp.txt', 'w') as f:
        for i in ls:
            f.write(i[-1]+'\n')
    sp = subprocess.Popen('java -cp "../soft/stanford-corenlp-full-2016-10-31/*" '
                          '-Xmx2g edu.stanford.nlp.pipeline.StanfordCoreNLP'
                          ' -annotators tokenize,ssplit,pos,depparse '
                          '-ssplit.eolonly true -file '
                          'nlpin.tmp.txt -outputDirectory ./ '
                          '-outputFormat json',
                          shell=True,
                          stdout=subprocess.DEVNULL,
                          stderr=subprocess.STDOUT)
    sp.wait()
    import json
    with open('nlpin.tmp.txt.json', 'r') as f:
        nlp = json.load(f)
    out = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict())))
    for i, v in enumerate(ls):
        m_ind_set = set([i for i in
                         nlp_get_word_index(nlp['sentences'][i]['tokens'],
                                            'MMMoneMMM')])
        g_ind_set = set([i for i in
                         nlp_get_word_index(nlp['sentences'][i]['tokens'],
                                            'GGGoneGGG')])
        if not m_ind_set or not g_ind_set:
            continue
        DT = nlp['sentences'][i]['basicDependencies']
        p_seq = get_short_path(DT, m_ind_set, g_ind_set)
        min_step = 9999
        for i in p_seq:
            if len(i) < min_step:
                min_step = len(i)
        pmid, n_mut, n_gene, record = v
        out[pmid][n_mut][min_step][n_gene] = qdic[pmid][n_mut][n_gene]
    for pmid in out.keys():
        for n_mut in out[pmid].keys():
            min_step = min(list(out[pmid][n_mut].keys()))
            for n_gene in out[pmid][n_mut][min_step].keys():
                yield (pmid, out[pmid][n_mut][min_step][n_gene])


def nlp_get_word_index(ls_dic, word):
    for i in ls_dic:
        if re.search(word, i['word']):
            yield i['index']


def get_pairs(dt, q):
    for i in dt:
        if q == i['dependent']:
            yield(i['governor'])
        if q == i['governor']:
            yield(i['dependent'])


def get_short_path(dt, start_set, end_set):
    out = []
    if not start_set:
        raise ValueError('start set is empty:{0}, {1}\n{2}'.format(start_set, end_set, dt))
    for i in start_set:
        seq_list = [[i]]  # rote
        end = end_set.copy()  # end set
        while end:
            new_seq_list = []
            for j in seq_list:  # j, seq
                if len(j) == 1:
                    if j[0] in end:
                        out.append(j)
                        end.remove(j[0])
                        new_seq_list.append(j)
                        continue
                for ans in get_pairs(dt, j[-1]):  # j, seq
                    if ans in j:
                        pass
                    elif ans in end_set:
                        out.append(j + [ans])
                        end.remove(ans)
                        new_seq_list.append(j + [ans])
                    else:
                        new_seq_list.append(j + [ans])
            if not new_seq_list and end:
                raise ValueError('list is empty: {0}, {1}, {2} \n{3}'.format(seq_list, end, out, dt))
            seq_list = new_seq_list
    return out
