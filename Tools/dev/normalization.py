import re
from collections import defaultdict


def mutation_string_normalizition(string, lu_nm):
    if not string:
        raise ValueError('void mutation string')
    '''rs number check'''
    if re.search('rs\d+', string):
        tmp = re.search('(rs\d+)', string)
        return tmp.group(1)
    '''mutil describe check, c. XXX (p.xxxx), get protein level'''
    tmp = re.match(r'c\..+?\s+(p\..+)', string)
    if tmp:
        string = tmp.group(1)
    '''if raw string match hgvs rules, return raw string'''
    if re.match(r'^\s*[pc]\s*\.\s*', string):
        return(re.sub(r' ', '', string).replace('X', '*'))
    else:
        '''illegal string trans to hgvs'''
        level, mtype, other_info = lu_nm.split('|', 2)
        if not mtype:
            raise ValueError('illegal string: {0}'.format(lu_nm))
        if level.lower() == 'c':
            if re.search('SUB', mtype):
                ref, loc, ale = other_info.split('|')
                if ref and loc and ale:
                    return '{0}.{1}{2}>{3}'.format(level, loc, ref, ale)
                else:
                    raise ValueError('illegal substitution'
                                     ': {0} {1}'.format(string, lu_nm))
            if re.search('INS|DUP|DEL', mtype):
                loc, ale = other_info.split('|')
                if loc and ale:
                    return '{0}.{1}{3}{2}'.format(
                            level, loc, ale, mtype.lower())
                elif loc and re.search('DEL', mtype):
                    return '{0}.{1}{2}'.format(
                            level, loc, mtype.lower())
                else:
                    raise ValueError('illegal indel(insertion, deletion'
                                     'duplication): {0}'.format(lu_nm))
        else:
            if re.search('SUB', mtype):
                ref, loc, ale = other_info.split('|')
                if ref and loc and ale:
                    return '{0}.{1}{2}{3}'.format(
                                    level, ref, loc, ale.replace('X', '*'))
                else:
                    raise ValueError('illegal substitution'
                                     ': {0}'.format(lu_nm))
            if re.search('INS|DEL|DUP', mtype):
                loc, ale = other_info.split('|')
                if loc and ale:
                    return '{0}.{1}{2}{3}'.format(
                                    level, ale, loc, mtype.lower())
                elif loc and re.search('DEL', mtype):
                    return '{0}.{1}{2}'.format(
                                    level, loc, mtype.lower())
                else:
                    raise ValueError('illegal indel(insertion, deletion, '
                                     'duplication): {0}'.format(lu_nm))
            if re.search('FS', mtype):
                ref, loc, ale, num = other_info.split('|')
                if ref and loc:
                    return '{0}.{1}{2}{3}{4}'.format(
                                level, ref, loc, ale.replace('X', '*'), num)
                else:
                    raise ValueError(
                        'unexpected frameshift format: {0} {1}'.format(
                                            string, lu_nm))


class gene_normor():
    '''gene normlization object. load hgnc file, get a object, which have
    method to covert entry id to symbol'''

    def __init__(self, hgnc):
        with open(hgnc, 'r') as f:
            header = f.readline().strip('\n').split('\t')
            raw = [i.split('\t') for i in f.read().split('\n') if i]
        for i, att in enumerate(header):
            if att == 'symbol':
                sym_index = i
            elif att == 'entrez_id':
                id_index = i
            else:
                pass
        if not id_index or not sym_index:
            raise ValueError('unexpected hgnc file')
        self.__data = {}
        for i in raw:
            en_id, symbol = i[id_index], i[sym_index]
            if not en_id:
                continue  # skip the no entry_id symbol
            if en_id in self.__data.keys():
                raise ValueError('one entrez id maped to more then one'
                                 ' symbol: {0} {1} {2}'.format(
                                                en_id, id_index, sym_index))
            self.__data[i[id_index]] = i[sym_index]

    def norm(self, entry_id):
        if isinstance(entry_id, int):
            entry_id = str(entry_id)
        elif re.match('\d+$', entry_id):
            pass
        else:
            raise ValueError('illegal entrez id: {0}'.format(entry_id))
        return self.__data[entry_id]


class disease_normor():
    '''disease normlization object. load CTD_disease file, get a object, which
    have method to covert disease id to mesh tree node'''

    def __init__(self, CTD_disease, alt=True):
        self.__data = defaultdict(dict)
        with open(CTD_disease, 'r') as f:
            for i in f.readlines():
                if re.match('#', i):
                    continue
                '''
                (name, disease_id, alt_ids, definition, p_ids,
                    nodes, p_nodes, syn, slim) = i.split('\t')
                '''
                (name, disease_id, alt_ids, definition, p_ids,
                    nodes, p_nodes, syn) = i.split('\t')
                if disease_id in self.__data.keys():
                    raise ValueError('one ID mapping to more then one'
                                     'disease: {0}'.format(d_id))
                self.__data[disease_id]['name'] = name
                self.__data[disease_id]['nodes'] = nodes
        if alt:
            with open(CTD_disease, 'r') as f:
                for i in f.readlines():
                    if re.match('#', i):
                        continue
                    '''(name, disease_id, alt_ids, definition, p_ids,
                        nodes, p_nodes, syn, slim) = i.split('\t')'''
                    (name, disease_id, alt_ids, definition, p_ids,
                        nodes, p_nodes, syn) = i.split('\t')
                    for alt_id in alt_ids.split('|'):
                        if alt_id not in self.__data.keys():
                            self.__data[alt_id]['name'] = name
                            self.__data[alt_id]['nodes'] = nodes

    def norm(self, entry_id):
        return self.__data[entry_id]
