import re
import os
import subprocess
from collections import defaultdict


def dependency_tree_shortest_path(qdic):
    '''
    dic[pmid][A][B]['sent_list'] = sentence_list: A, tagA; B, tagB; sentence
    dic[pmid][A][B]['path_list'] = path_list:
    '''
    qls = []
    for pmid in qdic.keys():
        for A in qdic[pmid].keys():
            for B in qdic[pmid][A].keys():
                if 'evidences' not in qdic[pmid][A][B].keys():
                    continue
                    # raise ValueError('data structure err:{0} {1} {2} {3}'.format(pmid, A, B, qdic[pmid][A][B]))
                for tagA, tagB, sentence in qdic[pmid][A][B]['evidences']:
                    raw_text = sentence.replace(tagA, 'AAAoneAAA')
                    raw_text = raw_text.replace(tagB, 'BBBoneBBB')
                    raw_text = raw_text.strip('\n')  # remove new line
                    raw_text = raw_text.strip(' ')  # remove space in end
                    raw_text = raw_text.strip('.')  # remove .
                    raw_text = raw_text.replace('.', '')
                    raw_text = raw_text.replace('>', '')
                    raw_text = raw_text.replace('*', '')
                    raw_text = raw_text + '.'
                    qls.append([pmid, A, B, raw_text])
    tmp_dir = './.tmp/'
    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)
    with open(os.path.join(tmp_dir, 'dt_raw.txt'), 'w') as f:
        for i in qls:
            f.write(i[-1]+'\n')
    sp = subprocess.Popen('java -cp '
                          '"../soft/stanford-corenlp-full-2016-10-31/*" '
                          '-Xmx2g edu.stanford.nlp.pipeline.StanfordCoreNLP'
                          ' -annotators tokenize,ssplit,pos,depparse '
                          '-ssplit.eolonly true -file '
                          '{0}/dt_raw.txt -outputDirectory {0}/ '
                          '-outputFormat json'.format(tmp_dir),
                          shell=True,
                          stdout=subprocess.DEVNULL,
                          stderr=subprocess.STDOUT)
    sp.wait()
    import json
    with open('{0}/dt_raw.txt.json'.format(tmp_dir), 'r') as f:
        nlp = json.load(f)
    for i, v in enumerate(qls):
        m_ind_set = set([i for i in
                         nlp_get_word_index(nlp['sentences'][i]['tokens'],
                                            'AAAoneAAA')])
        g_ind_set = set([i for i in
                         nlp_get_word_index(nlp['sentences'][i]['tokens'],
                                            'BBBoneBBB')])
        if not m_ind_set or not g_ind_set:
            continue
        DT = nlp['sentences'][i]['basicDependencies']
        p_seq = get_short_path(DT, m_ind_set, g_ind_set)
        min_step = 9999
        for i in p_seq:
            if len(i) < min_step:
                min_step = len(i)
        pmid, A, B, record = v
        if 'dt_path_lens' not in qdic[pmid][A][B].keys():
            qdic[pmid][A][B]['dt_path_lens'] = []
        qdic[pmid][A][B]['dt_path_lens'].append(min_step)
    return qdic


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
