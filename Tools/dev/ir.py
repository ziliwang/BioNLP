import Base


def mutation_disease(artilce):
    '''input a article obj, return tuple (mutation - disease), and evidence
    '''
    '''three level: single sentence, content, topic;
    topic disease is global, is desided by titlle key, then the most frequency
    disease in abstract; then the most frequency in full text.
    content: same mutation share content, co-occuption for score
    sentence: DT path, with a -ln function
    '''
    tk_m = [i for i in artilce.tk.bioe if i.type == 'mutation']
    ab_m = [i for i in artilce.results.bioe if i.type == 'mutation']
    re_m = [i for i in artilce.abstract.bioe if i.type == 'mutation']
    m_set = set([i.norm for i in tk_m + ab_m + re_m])
    print(m_set)
    '''topic'''
    tmp = [i for i in artilce.tk.bioe if i.type == 'disease']
    tmp_v = 10
    if not tmp:
        tmp = [i for i in artilce.abstract.bioe if i.type == 'disease']
        tmp_v = 5
    if not tmp:
        full = [i for i in artilce.result.bioe if i.type == 'disease']
        tmp_v = 1
    topic = Base.mesh_tree()
    for i in tmp:
        for j in i.norm['nodes'].split('|'):
            topic.add(j.split('/')[0], tmp_v)
    topic_node, topic_sore = topic.get_weight_path()
    topic_sore = float(topic_sore/len(tmp))
    print(topic_node, topic_sore)
    for i in m_set:
        for j in tk_m:
            if j.norm == i:
                print(i, topic_node, 'm - g in title and key words')
                break
        for j in ab_m:
            if j.norm == i:
                '''sentence'''
                if have_d(j, artilce):
                    pass
                '''else content'''
                elif have_cont(j, artilce, up, down):
                    get nearst up sentence, and its
