#!/usr/bin/env python3
import nltk
from ipdb import set_trace


class BioEntry():
    def __init__(self, start, end, type, string, id, norm):
        self.start = start
        self.end = end
        self.type = type
        self.string = string
        self.id = id
        self.norm = norm


class TK():
    def __init__(self, string):
        self.text = string
        self.bioe = []

    def add_e(self, bioe):
        self.bioe.append(bioe)


class Abstract():
    def __init__(self, string):
        self.text = string
        self.bioe = []
        self.sorted = False

    def add_e(self, bioe):
        self.bioe.append(bioe)

    def sort_entry(self):
        sort_ls = []
        for i in self.bioe:
            if len(sort_ls) < 1:
                sort_ls.append(i)
            else:
                for index in range(len(sort_ls)-1, -1, -1):
                        end = sort_ls[index].end
                        if i.end > end:
                            sort_ls.insert(index + 1, i)
                            break
                        elif index == 0:
                            sort_ls.insert(0, i)
        self.bioe = sort_ls
        self.sorted = True


class Results():
    def __init__(self, string):
        self.text = string
        self.bioe = []
        self.sorted = False

    def add_e(self, bioe):
        self.bioe.append(bioe)

    def sort_entry(self):
        sort_ls = []
        for i in self.bioe:
            if len(sort_ls) < 1:
                sort_ls.append(i)
            else:
                for index in range(len(sort_ls)-1, -1, -1):
                        end = sort_ls[index].end
                        if i.end > end:
                            sort_ls.insert(index + 1, i)
                            break
                        elif index == 0:
                            sort_ls.insert(0, i)
        self.bioe = sort_ls
        self.sorted = True


class Article():
    def __init__(self, pmid):
        self.tk = False
        self.abstract = False
        self.results = False
        self.pmid = pmid

    def add_tk(self, tk_obj):
        self.tk = tk_obj

    def add_abstract(self, abstract_obj):
        self.abstract = abstract_obj

    def add_results(self, results_obj):
        self.results = results_obj

    def get_abstract_relations(self):
        if not self.abstract:
            return False
        len_ls = self.__token_sent(self.abstract.text)
        keys = self.tk.bioe
        event = []
        for up_gene, up_disease, i, down_gene, down_disease in self.__get_hyponymy(self.abstract):
            relat = []
            mode_gene_in_abstract, mode_gene_in_result = self.__mode_entry()
            mode_disease_in_abstract, mode_disease_in_result = self.__mode_entry('disease')
            if mode_gene_in_abstract:
                relat.append(('r4', mode_gene_in_abstract[0], nltk.sent_tokenize(self.abstract.text)[self.__index_s(i, len_ls)] + ' Frequency of occurrence: {0}'.format(len(mode_gene_in_abstract))))
            if mode_disease_in_abstract:
                relat.append(('r4', mode_disease_in_abstract[0], nltk.sent_tokenize(self.abstract.text)[self.__index_s(i, len_ls)] + ' Frequency of occurrence: {0}'.format(len(mode_disease_in_abstract))))
            if mode_gene_in_result:
                relat.append(('r5', mode_gene_in_result[0], nltk.sent_tokenize(self.abstract.text)[self.__index_s(i, len_ls)] + ' Frequency of occurrence: {0}'.format(len(mode_gene_in_result))))
            if mode_disease_in_result:
                relat.append(('r5', mode_disease_in_result[0], nltk.sent_tokenize(self.abstract.text)[self.__index_s(i, len_ls)] + ' Frequency of occurrence: {0}'.format(len(mode_disease_in_result))))
            # 1 level
            for ii in [up_gene, up_disease, down_gene, down_disease]:
                if self.__index_s(i, len_ls) == self.__index_s(ii, len_ls):
                    relat.append(('r0', ii, nltk.sent_tokenize(self.abstract.text)[self.__index_s(i, len_ls)]))
                elif abs(self.__index_s(i, len_ls) - self.__index_s(ii, len_ls)) == 1:
                    if self.__index_s(i, len_ls) < self.__index_s(ii, len_ls):
                        relat.append(('r1', ii, nltk.sent_tokenize(self.abstract.text)[self.__index_s(i, len_ls):self.__index_s(i, len_ls) + 2]))
                    else:
                        relat.append(('r1', ii, nltk.sent_tokenize(self.abstract.text)[self.__index_s(i, len_ls)-1:self.__index_s(i, len_ls)+1]))
                else:
                    pass
            if up_gene:
                relat.append(('r6', up_gene, nltk.sent_tokenize(self.abstract.text)[self.__index_s(i, len_ls)] + 'Distance: {0}'.format(self.__index_s(i, len_ls) - self.__index_s(up_gene, len_ls))))
            if up_disease:
                relat.append(('r6', up_disease, nltk.sent_tokenize(self.abstract.text)[self.__index_s(i, len_ls)] + 'Distance: {0}'.format(self.__index_s(i, len_ls) - self.__index_s(up_disease, len_ls))))
            for iii in keys:
                relat.append(('r3', iii, nltk.sent_tokenize(self.abstract.text)[self.__index_s(i, len_ls)] + ' Title and Keywords: ' + self.tk.text))
            for item in relat:
                event.append((i, item))
        return event

    def get_results_relations(self):
        if not self.results:
            return False
        len_ls = self.__token_sent(self.results.text)
        keys = self.tk.bioe
        event = []
        for up_gene, up_disease, i, down_gene, down_disease in self.__get_hyponymy(self.results):
            '''a mutation try to find the gene and disease, [r4, r5, r6 supported]'''
            relat = []
            mode_gene_in_abstract, mode_gene_in_result = self.__mode_entry()
            mode_disease_in_abstract, mode_disease_in_result = self.__mode_entry('disease')
            if mode_gene_in_abstract:
                relat.append(('r4', mode_gene_in_abstract[0], nltk.sent_tokenize(self.results.text)[self.__index_s(i, len_ls)] + ' Frequency of occurrence: {0}'.format(len(mode_gene_in_abstract))))
            if mode_disease_in_abstract:
                relat.append(('r4', mode_disease_in_abstract[0], nltk.sent_tokenize(self.results.text)[self.__index_s(i, len_ls)] + ' Frequency of occurrence: {0}'.format(len(mode_disease_in_abstract))))
            if mode_gene_in_result:
                relat.append(('r5', mode_gene_in_result[0], nltk.sent_tokenize(self.results.text)[self.__index_s(i, len_ls)] + ' Frequency of occurrence: {0}'.format(len(mode_gene_in_result))))
            if mode_gene_in_result:
                relat.append(('r5', mode_disease_in_result[0], nltk.sent_tokenize(self.results.text)[self.__index_s(i, len_ls)] + ' Frequency of occurrence: {0}'.format(len(mode_disease_in_result))))
            # 1 level
            for ii in [up_gene, up_disease, down_gene, down_disease]:
                if self.__index_s(i, len_ls) == self.__index_s(ii, len_ls):
                    relat.append(('r0', ii, nltk.sent_tokenize(self.results.text)[self.__index_s(i, len_ls)]))
                elif abs(self.__index_s(i, len_ls) - self.__index_s(ii, len_ls)) == 1:
                    if self.__index_s(i, len_ls) < self.__index_s(ii, len_ls):
                        relat.append(('r1', ii, nltk.sent_tokenize(self.results.text)[self.__index_s(i, len_ls):self.__index_s(i, len_ls) + 2]))
                    else:
                        relat.append(('r1', ii, nltk.sent_tokenize(self.results.text)[self.__index_s(i, len_ls)-1:self.__index_s(i, len_ls)+1]))
                else:
                    pass
            if up_gene:
                relat.append(('r6', up_gene, nltk.sent_tokenize(self.results.text)[self.__index_s(i, len_ls)] + 'Distance: {0}'.format(self.__index_s(i, len_ls) - self.__index_s(up_gene, len_ls))))
            if up_disease:
                relat.append(('r6', up_disease, nltk.sent_tokenize(self.results.text)[self.__index_s(i, len_ls)] + 'Distance: {0}'.format(self.__index_s(i, len_ls) - self.__index_s(up_disease, len_ls))))
            for iii in keys:
                relat.append(('r3', iii, nltk.sent_tokenize(self.results.text)[self.__index_s(i, len_ls)] + ' Title and Keywords: ' + self.tk.text))
            for item in relat:
                event.append((i, item))
                # print('{0}\t{1}'.format(art['text'][i_s:i_e], art['text'][s:e]))
        return event

    def __get_hyponymy(self, part):
        if not part.sorted:
            part.sort_entry()
        for i in part.bioe:
            if i.type == 'mutation':
                # nearly gene or disease entry
                n = part.bioe.index(i)
                up_disease = False
                up_gene = False
                down_disease = False
                down_gene = False
                for tmp in reversed(part.bioe[0:n]):
                    if tmp.type == 'disease':
                            up_disease = tmp
                            break
                for tmp in reversed(part.bioe[0:n]):
                    if tmp.type == 'Gene' or tmp.type == 'Protein':
                            up_gene = tmp
                            break

                for tmp in part.bioe[n+1:]:
                    if tmp.type == 'disease':
                        down_disease = tmp
                        break
                for tmp in part.bioe[n+1:]:
                    if tmp.type == 'Gene' or tmp.type == 'Protein':
                        down_gene = tmp
                        break
                yield(up_gene, up_disease, i, down_gene, down_disease)

    def __token_sent(self, string):
        import nltk
        sent_len = []
        sents = []
        for i in nltk.sent_tokenize(string):
            if sent_len:
                sent_len.append(sent_len[-1] + len(i) + 1)
            else:
                sent_len.append(len(i))
            sents.append(i)
        return(sent_len)

    def __index_s(self, e, len_ls):
        if e:
            for i in len_ls:
                if e.end <= i:
                    return(len_ls.index(i))
        else:
            return 1000

    def __mode_entry(self, biotype='gene'):
        '''(ab(num, string),re(num, string))'''
        if self.abstract and self.abstract.bioe:
            abs_mode = self.__get_mode([i.string for i in self.abstract.bioe if i.type == biotype])
            abs_re = [i for i in self.abstract.bioe if i.string == abs_mode[1]]
        else:
            abs_re = False
        if self.results and self.results.bioe:
            res_mode = self.__get_mode([i.string for i in self.results.bioe if i.type == biotype])
            res_re = [i for i in self.results.bioe if i.string == res_mode[1]]
        else:
            res_re = False
        return(abs_re, res_re)

    def __get_mode(self, ls):
        mode_sth = (0, False)
        for i in set(ls):
            if ls.count(i) > mode_sth[0]:
                mode_sth = (ls.count(i), i)
        return(mode_sth)

    def get_topic_disease(self):
        tk = [i for i in self.tk.bioe if i.type == 'disease']
        ab = [i for i in self.abstract.bioe if i.type == 'disease']
        if len(tk):
            return tk
        else:
            pass


class w_tree():

    def __init__(self, v):
        self.value = v
        self.children = []

    def add_child(self, c, w):
        self.children.append([c, float(w)])

    def children_set(self):
        return(i[0].value for i in self.children)

    def get_child(self, name):
        if name in self.children_set():
            for i in self.children:
                if i[0].value == name:
                    return i[0]
        else:
            return False

    def odd_child_weight(self, name, weight):
        if name in self.children_set():
            for i in self.children:
                if i[0].value == name:
                    i[1] += float(weight)
        else:
            raise ValueError('odd child weight failed')


class mesh_tree():

    def __init__(self):
        self.root = w_tree('ROOT')

    def add(self, string, weight):
        tmp = []
        p = self.root
        weight = float(weight / len(string.split('.')))
        for i in string.split('.'):
            tmp.append(i)
            c = p.get_child('.'.join(tmp))
            if c:
                p.odd_child_weight('.'.join(tmp), weight)
                p = c
            else:
                n = w_tree('.'.join(tmp))
                p.add_child(n, weight)
                p = n

    def show(self):
        print('graph TB')
        self.__show_child(self.root)

    def __show_child(self, n):
        for i in n.children:
            print('{0}--{1}-->{2}'.format(n.value, i[1], i[0].value))
            self.__show_child(i[0])

    def get_weight_path(self):
        nodes = list(self.getw(self.root, 0))
        w_max = 0
        for n, w in nodes:
            if w > w_max:
                w_max = w
        for n, w in nodes:
            if w == w_max:
                yield n, w

    def getw(self, node, v):
        if node.children:
            for i in node.children:
                yield from self.getw(i[0], v+i[1])
        else:
            yield((node.value, v))
