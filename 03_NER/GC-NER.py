#!/usr/bin/env python3
import NLTK

'''
sentence operation
duanluo -> fenju -> tagger ->
'''


class article():

    def __int__(self):
        self.data = {}

    def add(self, pmid, text):
        if pmid in self.data.keys():
            print('pmid: {0} existed. pass'.format(pmid))
            return
        self.data[pmid] = {}
        self.data[pmid]['sentence_list'] = {}
        sentence_list = nltk.sent_tokenize(text)
        for i in range(len(sentence_list)):
            self.data[pmid]['sentence_list'][i]




class paragraph():

    def __init__(self, text):
        sentence_list = nltk.sent_tokenize(text)


class word():

    def __init__(self, w):
        self.word = w

    def stem(self):
        return self.word


class sentence():
    def __init__(self, sentence, gene_dic, disease_dic, NERsuite_PATH, SETH_PATH):
        self.gene_dic = gene_dic
        self.disease_dic = disease_dic
        self.NERsuite_PATH = NERsuite_PATH
        self.SETH_PATH = SETH_PATH
        self.sentence = sentence


#os.popen('''/usr/local/testlib/NERsuite/bin/nersuite_gtagger -d /home/zili/Desktop/BioNLP/test_lib/models/gtagger 'There is controversy regarding the localization of BRCA1 and BRCA2 proteins to e
#    ...: ither nucleus or cytoplasm and whether the expression is present in premeiotic germ cells or can still be expressed in mitotic spermatogonia.' ''')
#
# nltk.sent_tokenize('BRCA1 and BRCA2 breast cancer susceptibility genes encode proteins, the normal cellular functions of which are complex and multiple, and germ-line mutations in individuals pred
#    ...: ispose both to breast and to ovarian cancer. There is nevertheless substantial evidence linking BRCA1 and BRCA2 to homologous recombination and DNA repair, to transcriptional control and to tissue
#    ...:  proliferation. There is controversy regarding the localization of BRCA1 and BRCA2 proteins to either nucleus or cytoplasm and whether the expression is present in premeiotic germ cells or can sti
#    ...: ll be expressed in mitotic spermatogonia. We report herein an immunohistochemical study of BRCA1 and BRCA2 distribution in a rather unusual tissue (an ovotestis), which addresses this issue.')
# nltk.WordPunctTokenizer().tokenize('There is controversy regarding the localization of BRCA1 and BRCA2 proteins to either nucleus or cytoplasm and whether the expression is present in premeiotic
#     ...: germ cells or can still be expressed in mitotic spermatogonia.')
