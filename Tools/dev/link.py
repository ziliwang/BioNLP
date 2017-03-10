import re
import os
import Base
import pickle


def link_to_transvar(article_ls):
    '''write down a pickle for article ls and transVar'''
    with open('artilce_ls.pkl', 'wb') as f:
        pickle.dump(article_ls, f)
    with open('pass_to_transvar_canno', 'w') as f:
        with open('pass_to_transvar_panno', 'w') as ff:
            for i in article_ls:
                if i.get_abstract_relations():
                    for ii in i.get_abstract_relations():
                        if ii[1][1].type == 'Gene':
                            """ii[0] is mutation, ii[1][1] is gene"""
                            if re.match('c', ii[0].norm):
                                f.write(ii[1][1].norm + ':' + ii[0].norm +
                                        '\n')
                            elif re.match('p', ii[0].norm):
                                ff.write(ii[1][1].norm + ':' + ii[0].norm +
                                         '\n')
                            else:
                                """
                                m level and g level: check the gene in or
                                not
                                """
                                pass
                if i.get_results_relations():
                    for ii in i.get_results_relations():
                        if ii[1][1].type == 'Gene':
                            """ii[0] is mutation, ii[1][1] is gene"""
                            if re.match('c', ii[0].norm):
                                f.write(ii[1][1].norm + ':' + ii[0].norm +
                                        '\n')
                            elif re.match('p', ii[0].norm):
                                ff.write(ii[1][1].norm + ':' + ii[0].norm +
                                         '\n')
                            else:
                                """
                                m level and g level: check the gene in or
                                not
                                """
                                pass


class linker():
    '''a class for call link'''
    def __init__(self, transvar_out_path):
        self.__data = {}
        self.__read(transvar_out_path, 'pass_to_transvar_panno.out')
        self.__read(transvar_out_path, 'pass_to_transvar_canno.out')

    def __read(self, transvar_out_path, part):
        with open(os.path.join(transvar_out_path, part),
                  'r') as f:
            for i in f.readlines():
                items = i.split('\t')
                if items[1] != '.':
                    trans = items[1].split(' ')[0]
                    if items[0] in self.__data.keys():
                        rec_tmp = trans + ':' + items[4]
                        if rec_tmp not in self.__data[items[0]]:
                            self.__data[items[0]].append(trans + ':' +
                                                         items[4])
                    else:
                        self.__data[items[0]] = [trans + ':' + items[4]]

    def is_link(self, string):
        if string in self.__data.keys():
            return True
        return False

    def parser(self, string):
        if string in self.__data.keys():
            return self.__data[string]
        return False
