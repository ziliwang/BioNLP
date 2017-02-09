#!/usr/bin/env python3
import re
import os
import pickle
import csv
import numpy as np
import subprocess
from ipdb import set_trace
from collections import defaultdict
from pprint import pprint


def test(f):
    with open(f, 'rb') as fh:
        data = pickle.loads(fh.read())
    pairs_ls = []
    for i in data:
        if i[6].lower() == 'gene':
            if i[9]:
                pairs_ls.append((i[4], i[9]))
            else:
                pairs_ls.append((i[4], i[7]))
    with open('pass_to_transvar_canno', 'w') as f:
        with open('pass_to_transvar_panno', 'w') as ff:
            for m, g in set(pairs_ls):
                if re.match(r'c\.', m):
                    f.write(g + ':' + m + '\n')
                else:
                    ff.write(g + ':' + m + '\n')


def transvar_out_anno(canno, panno, gene_norm):
    """anno transvar_out_anno with hgvs format"""
    cdic = defaultdict(list)
    with open(canno, 'r') as f:
        for i in f.readlines():
            items = i.split('\t')
            if items[1] != '.':
                trans = items[1].split(' ')[0]
                cdic[items[0]].append(trans + ':' + items[4])
    pdic = defaultdict(list)
    with open(panno, 'r') as f:
        for i in f.readlines():
            items = i.split('\t')
            if items[1] != '.':
                trans = items[1].split(' ')[0]
                pdic[items[0]].append(trans + ':' + items[4])
    with open(gene_norm, 'rb') as fh:
        data = pickle.loads(fh.read())
    output = []
    for i in data:
        if i[6].lower() == 'gene':
            if i[9]:
                query_key = i[9] + ':' + i[4]
            else:
                query_key = i[7] + ':' + i[4]
            if re.match('c', i[4]):
                output.append([dic_search(cdic, query_key)] + i)
            else:
                output.append([dic_search(pdic, query_key)] + i)
    with open('transvar_anno', 'wb') as f:
        pickle.dump(output, f)
    with open('transvar_anno.txt', 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(output)
        IndexError


def dic_search(dic, key):
    if key in dic.keys():
        refseq = [i for i in dic[key] if re.match('[XN]M', i.upper())]
        print(refseq)
        if refseq:
            return('|'.join(refseq))
        else:
            return('|'.join(dic[key]))
    else:
        return('Failed')


def simple_rule(fl):
    with open(fl, 'rb') as f:
        data = pickle.load(f)
    dic = defaultdict(
            lambda: defaultdict(
                lambda: defaultdict(lambda: defaultdict(list))))
    fuzzy_dic = defaultdict(
                    lambda: defaultdict(lambda: defaultdict(list)))
    for i in data:
        if i[0] != 'Failed':
            dic[i[1]][i[5]][i[6]][i[10]].append(i)
    with open('gene-mutation.csv', 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        for i in dic.keys():
            for j in dic[i].keys():
                for o in ['r0', 'r1', 'r2', 'r3', 'r4', 'r5', 'r6']:
                    if o in dic[i][j].keys():
                        if len(dic[i][j][o].keys()) > 1 and o == 'r0':
                            for p in dic[i][j][o].keys():
                                for q in dic[i][j][o][p]:
                                    fuzzy_dic[i][j][p].append(q)
                        else:
                            for p in dic[i][j][o].keys():
                                for q in dic[i][j][o][p]:
                                    writer.writerow(q)
                        break
        for i in get_DT_step(fuzzy_dic):
            for j in i:
                writer.writerow(j)


def get_DT_step(qdic):
    ls = []
    # i, pmid; j, n_mut; o, n_gene; p, record
    for i in qdic.keys():
        for j in qdic[i].keys():
            for o in qdic[i][j].keys():
                for p in qdic[i][j][o]:
                    raw_text = p[11]
                    raw_mut_text = p[3]
                    raw_gen_text = p[8]
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
    with open('tmp/nlpin.tmp.txt', 'w') as f:
        for i in ls:
            f.write(i[-1]+'\n')
    sp = subprocess.Popen('java -cp "./stanford-corenlp-full-2016-10-31/*" '
                          '-Xmx2g edu.stanford.nlp.pipeline.StanfordCoreNLP'
                          ' -annotators tokenize,ssplit,pos,depparse '
                          '-ssplit.eolonly true -file '
                          'tmp/nlpin.tmp.txt -outputDirectory tmp/ '
                          '-outputFormat json',
                          shell=True,
                          stdout=subprocess.DEVNULL,
                          stderr=subprocess.STDOUT)
    sp.wait()
    import json
    with open('tmp/nlpin.tmp.txt.json', 'r') as f:
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
        print('{0:20} -> {1:20} -> {2:3}: {3}'.format(n_mut, n_gene, min_step, [i[11] for i in qdic[pmid][n_mut][n_gene]]))
        out[pmid][n_mut][min_step][n_gene] = qdic[pmid][n_mut][n_gene]
    for pmid in out.keys():
        for n_mut in out[pmid].keys():
            min_step = np.min(list(out[pmid][n_mut].keys()))
            for n_gene in out[pmid][n_mut][min_step].keys():
                yield out[pmid][n_mut][min_step][n_gene]


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
