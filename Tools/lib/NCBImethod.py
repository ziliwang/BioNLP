#!/usr/bin/env python3
import sys
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
import Base
import re
import os
from pdfminer import encodingdb
import nltk

try:
    import PubMedDB
except:
    sys.path.append('/home/zili/Desktop/BioNLP/Tools/lib')
    import PubMedDB

# previous code for pdf parser
"""
class PDFxmlParser():

    def __init__(self, file):
        import xml.etree.ElementTree as ET
        text = open(file, 'r').read()
        text = re.sub(r'[\x00-\x1f\x7f-\x9f]', '', text)
        text = "".join([i for i in text if ord(i) < 128])
        root = ET.fromstring(text)
        mode_size = self.__get_size(root)
        result_sw = 0  # result switch
        textline_list = []  # result list contains line list
        for page in root.findall('page'):
            for textbox in page.findall('textbox'):
                for textline in textbox.findall('textline'):
                    start_sw = self.__get_result_part(textline)
                    if start_sw:
                        result_sw = start_sw
                        continue
                    else:
                        # check switch
                        if result_sw:
                            # check close
                            if self.__get_result_end(textline, result_sw):
                                result_sw = 0
                                break
                            # check size
                            if self.__size_check(textline, mode_size):
                                text_line_str_list = [i.text for i in
                                                      textline.findall('text')]
                                textline_list.append(text_line_str_list)
        self.text = self.__join_line(textline_list)

    def __get_size(self, root):
        sizels = []
        for page in root.findall('page'):
            for textbox in page.findall('textbox'):
                for textline in textbox.findall('textline'):
                    for text in textline.findall('text'):
                        if 'size' in text.attrib.keys():
                            sizels.append(float(text.attrib['size']))
        most_frequent = self.get_mode(sizels)
        return(most_frequent)

    def __get_result_part(self, textline):
        '''get result return size or return false'''
        sizels = []
        textline_str = ''.join([text.text if text.text else '-'
                                for text in textline.findall('text')])
        if re.match(r'Result|RESULT', textline_str):
            for text in textline.findall('text'):
                if 'size' in text.attrib.keys():
                    sizels.append(float(text.attrib['size']))
        if len(sizels) > 0:
            most_frequent = self.get_mode(sizels)
            return(most_frequent)
        else:
            return(False)

    def __get_result_end(self, textline, result_size):
        for text in textline.findall('text'):
            if 'size' in text.attrib.keys():
                if float(text.attrib['size']) == result_size:
                    return True
        return False

    def __size_check(self, textline, mode_size):
        sizels = []
        for text in textline.findall('text'):
            if 'size' in text.attrib.keys():
                sizels.append(float(text.attrib['size']))
        most_frequent = self.get_mode(sizels)
        if most_frequent == mode_size:
            return True
        return False

    def get_mode(self, numlist):
        from scipy.stats import mode
        most_frequent = mode(numlist)[0][0]  # mode
        # most_frequent = np.argmax(np.bincount(np.array(numlist)))  # mode
        return(most_frequent)

    def __join_line(self, textline_list):
        text = ''
        for i in textline_list:
            '''
            (cid:404) &amp; &lt; (cid:148) (cid:13) (cid:176) (cid:131)
            (cid:147) (cid:150) &quot; (cid:146) &gt; (cid:151)
            '''
            step0 = [j for j in i if j]
            # decode cid
            step1 = [encodingdb.name2unicode(j) if re.search('cid', j)
                     else j for j in step0]
            # decode &amp
            step2 = []
            for j in step1:
                if j == '&amp;':
                    step2.append('&')
                elif j == '&lt;':
                    step2.append('<')
                elif j == '&gt;':
                    step2.append('>')
                elif j == '&quot;':
                    step2.append('"')
                elif j == '&#39;':
                    step2.append("'")
                else:
                    step2.append(j)
            # -
            try:
                if step2[-2] is '-' and step2[-1] is '\n':
                    step3 = step2[0:-2:]
                else:
                    step3 = step2
            except:
                step3 = step2
            # \n
            step4 = [j if j is not '\n' else ' ' for j in step3]
            text += r''.join(step4)
        return text


def Batch_Pdf_xml(dir_path):
    for root, dir, files in os.walk(dir_path):
        for file in files:
            name, suffix = file.split('.')
            name = name.replace(' ', '')
            if suffix == 'pdf':
                if re.search(r' ', file):
                    os.rename(os.path.join(root, file),
                              os.path.join(os.path.join(root, name+'.pdf')))
                os.popen(
                  'pdf2txt.py -M 1.3 -W 0.4 -t xml -o '
                  '{0}.xml {0}.pdf'.format(name)
                  )


"""


def Batch_parse_xml(dir_path):
    '''batch parser pdf into xml, and parser xml into a single string'''
    output = []
    for root, dirs, files in os.walk(dir_path):
        for f in files:
            try:
                name, suffix = f.split('.')
            except:
                continue
            if suffix in ['pdf']:
                filename = os.path.join(root, f)
                res = PDFParser(filename)
                res.parse()
                output.append([name, res.text])
    os.popen('rm {0}/*.xml'.format(dir_path))
    return(output)


def NCBI_input_text_clean(string):
    # replace '|' to ' '
    # replace '\t' to ' '
    # repalce '\n' to ' '
    string = "".join([i for i in string if ord(i) < 128])
    return re.sub(r'\||\t|\n', ' ', string)


def xml2input(ls, output_prefix):
    '''read text and write into input format'''
    input_for_tmvar = output_prefix + '.tmvarI'
    input_for_GNP = output_prefix + '.GNormPlusI'  # dir
    input_for_DN = output_prefix + '.DNormI'
    if not os.path.exists(input_for_GNP):
        os.makedirs(input_for_GNP)
    with open(input_for_tmvar + '.results', 'w') as tmvar:
        with open(input_for_DN + '.results', 'w') as dn:
            for pmid, text in ls:
                if len(text) == 0:
                    continue
                tmp_str = NCBI_input_text_clean(text)
                tmvar.write('{0}|a|{1}\n\n'.format(pmid, tmp_str))
                dn.write('{0}\t{1}\n'.format(pmid, tmp_str))
                with open(input_for_GNP + '/{0}.results.txt'.format(pmid),
                          'w') as gnp:
                    gnp.write('{0}|a|{1}\n\n'.format(pmid, tmp_str))


class DBopt():
    '''operate db. extract info from database'''

    def __init__(self, db):
        engine = create_engine('postgresql://parser:parser@localhost/' + db)
        self.DBSession = sessionmaker(bind=engine)

    def get_text(self, output_prefix):
        input_for_tmvar = output_prefix + '.tmvarI'
        input_for_GNP = output_prefix + '.GNormPlusI'  # dir
        input_for_DN = output_prefix + '.DNormI'
        if not os.path.exists(input_for_GNP):
            os.makedirs(input_for_GNP)
        qurey_session = self.DBSession()
        articles = qurey_session.query(PubMedDB.Citation).all()
        # title and keywords
        with open(input_for_tmvar + '.tk', 'w') as tmvar:
            with open(input_for_DN + '.tk', 'w') as dn:
                for article in articles:
                    keywords = qurey_session.query(PubMedDB.Keyword).filter_by(
                                                    fk_pmid=article.pmid)
                    keyword_ls = [i.keyword for i in keywords if i]
                    tmp_str = article.article_title
                    if len(keyword_ls) > 0:
                        tmp_str += ', '.join(keyword_ls)
                    pmid = article.pmid
                    tmp_str = self.__NCBI_input_text_clean(tmp_str)
                    tmvar.write('{0}|a|{1}\n\n'.format(pmid, tmp_str))
                    dn.write('{0}\t{1}\n'.format(pmid, tmp_str))
                    with open(input_for_GNP + '/{0}.tk.txt'.format(pmid),
                              'w') as gnp:
                        gnp.write('{0}|a|{1}\n\n'.format(pmid, tmp_str))
        with open(input_for_tmvar + '.abstract', 'w') as tmvar:
            with open(input_for_DN + '.abstract', 'w') as dn:
                for article in articles:
                    pmid = article.pmid
                    try:
                        abstract = qurey_session.query(
                                    PubMedDB.Abstract).filter_by(fk_pmid=pmid)
                        tmp_str = abstract[0].abstract_text
                    except:
                        continue
                    tmp_str = self.__NCBI_input_text_clean(tmp_str)
                    tmvar.write('{0}|a|{1}\n\n'.format(pmid, tmp_str))
                    dn.write('{0}\t{1}\n'.format(pmid, tmp_str))
                    with open(input_for_GNP + '/{0}.abstract.txt'.format(pmid),
                              'w') as gnp:
                        gnp.write('{0}|a|{1}\n\n'.format(pmid, tmp_str))

        qurey_session.close()

    def __NCBI_input_text_clean(self, string):
        # replace '|' to ' '
        # replace '\t' to ' '
        # repalce '\n' to ' '
        string = "".join([i for i in string if ord(i) < 128])
        return re.sub(r'\||\t|\n', ' ', string)


def output_parser(outdir, prefix):
    '''output: article obj list'''
    tm_tk = outdir + '/' + prefix + '.tmvarI.tk.PubTator'
    tm_ab = outdir + '/' + prefix + '.tmvarI.abstract.PubTator'
    tm_re = outdir + '/' + prefix + '.tmvarI.results.PubTator'
    dn_tk = outdir + '/' + prefix + '.DNormO.tk'
    dn_ab = outdir + '/' + prefix + '.DNormO.abstract'
    dn_re = outdir + '/' + prefix + '.DNormO.results'
    gn_dir = outdir + '/' + prefix + '.GNormPlusO'
    articles = []
    for pmid, string, mutation_entry_ls in tmvar_parser(tm_tk):
        article_obj = Base.Article(pmid)  # create article_obj
        tk_obj = Base.TK(string)  # create tk_obj
        for e in mutation_entry_ls:
            tk_obj.add_e(e)  # add entry
        article_obj.add_tk(tk_obj)
        articles.append(article_obj)
    # abstract
    for pmid, string, mutation_entry_ls in tmvar_parser(tm_ab):
        article_obj = found_article(pmid, articles)  # found article_obj
        abstract_obj = Base.Abstract(string)  # create abstract_obj
        for e in mutation_entry_ls:
            abstract_obj.add_e(e)  # add entry
        article_obj.add_abstract(abstract_obj)
    # results
    for pmid, string, mutation_entry_ls in tmvar_parser(tm_re):
        article_obj = found_article(pmid, articles)  # found article_obj
        results_obj = Base.Results(string)  # create tk_obj
        for e in mutation_entry_ls:
            results_obj.add_e(e)  # add entry
        article_obj.add_results(results_obj)
    # DNorm
    for pmid, entry in DNorm_parser(dn_tk):
        article_obj = found_article(pmid, articles)
        article_obj.tk.add_e(entry)
    for pmid, entry in DNorm_parser(dn_ab):
        article_obj = found_article(pmid, articles)
        article_obj.abstract.add_e(entry)
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
    return articles


def found_article(pmid, articles):
    for article in articles:
        if pmid == article.pmid:
            return article
    raise ValueError('{0} have no title'.format(pmid))


def tmvar_parser(outputfile):
    '''output: pmid, string, mutation_entry'''
    with open(outputfile, 'r') as f:
        mutation_results = [i.split('\n') for i in f.read().split('\n\n') if i]
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
                mutation_entry_ls.append(
                    Base.BioEntry(start, end, 'mutation', string, norm))
        yield(pmid, text, mutation_entry_ls)


def DNorm_parser(outputfile):
    '''output: pmid, disease_entry'''
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
        yield(pmid, Base.BioEntry(start, end, 'disease', string, id))


def GNormPlus_parser(outputdir):
    '''output: pmid, gene_entry, article_part'''
    for root, dirs, files in os.walk(outputdir):
        for fn in files:
            with open(os.path.join(root, fn), 'r') as f:
                records = [i for i in f.read().split('\n') if i]
                for r in records:
                    is_text = re.match(r'\d+\|a\|(.+)$', r)
                    if is_text:
                        continue
                    else:
                        pmid, start, end, string, type, id = r.split('\t')
                        if not re.match('gene|protein', type, re.IGNORECASE):
                            continue
                        pmid = int(pmid)
                        start = int(start)
                        end = int(end)
                        yield(pmid,
                              Base.BioEntry(start, end, type, string, id),
                              fn)


class PDFParser():
    '''a class for parser pdf based on size and font'''

    def __init__(self, pdffilepath):
        '''default naming rules: pmid + '.pdf' '''
        if not pdffilepath:
            raise ValueError('None type detected')
        if not re.match('pdf', pdffilepath.split('.')[-1], re.IGNORECASE):
            raise ValueError('not pdf file')
        self.pmid = pdffilepath.split('.')[-2]
        self.prefix = '.'.join(pdffilepath.split('.')[:-1])
        self.text = False

    def parse(self):
        '''previous code, for optimize
        # from itertools import product
        # for M, W in product(range(24, 25, 2), range(2, 3)):'''
        M_value = 2.4  # 字母间隔
        W_value = 0.2  # 词间隔
        self.__run_pdf2txt(M_value, W_value)
        self.text = self.__parser_xml()
        '''previous code, for optimize
        # while True:
        #     self.__run_pdf2txt(M_value, W_value)
        #     textline_list = self.__parser_xml()
        #     # if self.__text_check(textline_list):
        #     self.text = self.__join_line(textline_list)
        #     break
        #     # else:
        #     #     M_value += 0.2
        #     #     W_value += 0.1'''

    def __len_check(self, text):
        if not text:
            return False
        if len([i for i in nltk.word_tokenize(text)]) < 50:
            return False
        else:
            return True

    def __run_pdf2txt(self, M_value, W_value):
        import subprocess
        sp = subprocess.Popen(['pdf2txt.py', '-M', str(M_value), '-W',
                               str(W_value), '-t', 'xml', '-o',
                               self.prefix + '.xml', self.prefix + '.pdf'],
                              stdout=subprocess.DEVNULL,
                              stderr=subprocess.DEVNULL)
        out, err = sp.communicate()
        sp.wait()

    def __parser_xml(self):
        '''
        function: extract result part from pdf literature, and return a raw
                  textline list, if extract nothing, extract main body.
        theory: based on the font and size.
        hypothesis: 1. the main body and other part exists distinction in
                       font and size.
                    2. literature contains different part, like abstract and
                       reference
                    3. the first page of pdf contains abstract and last page
                       contains reference which contains less main body.
                    4. the main body have most words.
        '''
        text = ''
        import xml.etree.ElementTree as ET
        with open(self.prefix + '.xml', 'r') as f:
            text = f.read()
        text = re.sub(r'[\x00-\x1f\x7f-\x9f]', '', text)
        text = "".join([i for i in text if ord(i) < 128])
        root = ET.fromstring(text)
        mode_size = self.__get_size(root)
        result_size_font = False  # result switch
        textline_list = []  # result list contains line list
        got_results = False
        fisrt_page_boxnum = len(root.findall('page')[0].findall('textbox'))
        boxnum = 0  # the first 30 line pass
        for page in root.findall('page'):
            if got_results:
                break
            for textbox in page.findall('textbox'):
                if got_results:
                    break
                boxnum += 1
                '''skip fisrt 4/5 text box in first page'''
                if boxnum < 4/5 * fisrt_page_boxnum:
                    continue
                for textline in textbox.findall('textline'):
                    if got_results:
                        break
                    if not result_size_font:
                        start_sw = self.__get_result_part(textline, mode_size)
                        if start_sw:
                            result_size_font = start_sw
                            continue
                    else:
                        # check close
                        if self.__get_result_end(textline, result_size_font):
                            got_results = True
                            break
                        # check size
                        '''with same size is regard as main body'''
                        if self.__size_check(textline, mode_size):
                            text_line_str_list = [i.text for i in
                                                  textline.findall('text')]
                            textline_list.append(text_line_str_list)
        fetched_text = self.__join_line(textline_list)
        if not self.__len_check(fetched_text):
            textline_list = []
            for page in root.findall('page'):
                for textbox in page.findall('textbox'):
                    for textline in textbox.findall('textline'):
                        if self.__size_check(textline, mode_size):
                            text_line_str_list = [i.text for i in
                                                  textline.findall('text')]
                            textline_list.append(text_line_str_list)
            fetched_text = self.__join_line(textline_list)
        return fetched_text

    def __get_size(self, root):
        '''get the size mode'''
        sizels = []
        for page in root.findall('page')[1:-1]:
            for textbox in page.findall('textbox'):
                for textline in textbox.findall('textline'):
                    for text in textline.findall('text'):
                        if 'size' in text.attrib.keys():
                            sizels.append(float(text.attrib['size']))
        most_frequent = self.get_mode(sizels)
        return(most_frequent)

    def __get_result_part(self, textline, mode_size):
        '''get result return size or return false'''
        sizels = []
        fontls = []
        textline_str = ''.join([text.text if text.text else ''
                                for text in textline.findall('text')])
        textline_str = textline_str.replace(' ', '')
        if re.match(r'result', textline_str, re.IGNORECASE) and (
            not re.search('[\.,\-]', textline_str)
        ):
            for text in textline.findall('text'):
                if 'size' in text.attrib.keys():
                    #  print(text.attrib['font'])
                    sizels.append(float(text.attrib['size']))
                    fontls.append(text.attrib['font'])
        if len(sizels) > 0:
            most_frequent = self.get_mode(sizels)
            most_font = self.get_mode(fontls)
            return(most_frequent, most_font)
        return(False)

    def __get_result_end(self, textline, result_size):
        sizels = []
        fontls = []
        string = ''
        for text in textline.findall('text'):
            if 'size' in text.attrib.keys():
                sizels.append(float(text.attrib['size']))
                fontls.append(text.attrib['font'])
                if text.text:
                    string += text.text
        most_size = self.get_mode(sizels)
        most_font = self.get_mode(fontls)
        '''regard as same level with results when with same size and font, and
        also a white list for word'''
        if most_size == result_size[0] and most_font == result_size[1] and (
            re.search('discuss|conclusion|reference|acknowledgement|Material'
                      '|Method|Supporting', string, re.IGNORECASE)
        ):
            return True
        return False

    def __size_check(self, textline, mode_size):
        sizels = []
        for text in textline.findall('text'):
            if 'size' in text.attrib.keys():
                sizels.append(float(text.attrib['size']))
        most_frequent = self.get_mode(sizels)
        if most_frequent == mode_size:
            return True
        return False

    def get_mode(self, ls):
        '''from scipy.stats import mode
        most_frequent = mode(numlist)[0][0]  # mode
        # most_frequent = np.argmax(np.bincount(np.array(numlist)))  # mode'''
        # change to support string
        mode_sth = (0, False)
        for i in set(ls):
            if ls.count(i) > mode_sth[0]:
                mode_sth = (ls.count(i), i)
        return(mode_sth[1])

    def __text_check(self, textline_list):
        cleaned_textline = []
        for i in textline_list:
            cleaned_textline.append(''.join([j for j in i if j]))
        total_line = len(cleaned_textline)
        failed_line = 0
        for i in cleaned_textline:
            tmp = i.split(' ')
            if len(i)/len(tmp) < 4.0:
                failed_line += 1
        if failed_line/total_line < 0.1:
            return True
        else:
            return False

    def __join_line(self, textline_list):
        '''return text, and clean text'''
        text = ''
        for i in textline_list:
            '''
            pdf font change
            (cid:404) &amp; &lt; (cid:148) (cid:13) (cid:176) (cid:131)
            (cid:147) (cid:150) &quot; (cid:146) &gt; (cid:151)
            '''
            step0 = [j for j in i if j]
            # decode cid
            step1 = [encodingdb.name2unicode(j) if re.search('cid', j)
                     else j for j in step0]
            # decode &amp
            step2 = []
            for j in step1:
                if j == '&amp;':
                    step2.append('&')
                elif j == '&lt;':
                    step2.append('<')
                elif j == '&gt;':
                    step2.append('>')
                elif j == '&quot;':
                    step2.append('"')
                elif j == '&#39;':
                    step2.append("'")
                else:
                    step2.append(j)
            # -
            try:
                if step2[-1] == '\n':
                    if step2[-2] is '-':
                        step3 = step2[0:-2:]
                    else:
                        step3 = step2[::]
                else:
                    if step2[-1] is '-':
                        step3 = step2[0:-1:]
                    else:
                        step3 = step2 + ['\n']
            except:
                step3 = step2
            # \n
            step4 = [j if j is not '\n' else ' ' for j in step3]
            text += r''.join(step4)
        return text


def print_sm(article_obj_ls):
    for i in article_obj_ls:
        if i.get_abstract_relations():
            for ii in i.get_abstract_relations():
                yield(i.pmid, 'abstract', ii[0].string, ii[0].id, ii[1][0], ii[1][1].type, ii[1][1].string, ii[1][1].id, ii[1][2])
        if i.get_results_relations():
            for ii in i.get_results_relations():
                yield(i.pmid, 'results', ii[0].string, ii[0].id, ii[1][0], ii[1][1].type, ii[1][1].string, ii[1][1].id, ii[1][2])


def mutation_string_normalizition(string, norm):
    if not string:
        return ''
    if re.match('rs\d+', string):
        return ''
    tmp = re.match(r'c\..+?\s+(p\..+)', string)
    if tmp:
            string = tmp.group(1)
    if re.match(r'^\s*[pc]\s*\.\s*', string):
        return(re.sub(r' ', '', string).replace('X', '*'))
    else:
        level, mtype, other_info = norm.split('|', 2)
        if not mtype:
            return ''
        if level.lower() == 'c':
            if re.search('SUB', mtype):
                ref, loc, ale = other_info.split('|')
                if ref and loc and ale:
                    return '{0}.{1}{2}>{3}'.format(level, loc, ref, ale)
                else:
                    return ''
            if re.search('INS|DEL|DUP', mtype):
                loc, ale = other_info.split('|')
                if loc and ale:
                    return '{0}.{1}{3}{2}'.format(level, loc, ale, mtype.lower())
                else:
                    return ''
        else:
            if re.search('SUB', mtype):
                ref, loc, ale = other_info.split('|')
                if ref and loc and ale:
                    return '{0}.{1}{2}{3}'.format(level, ref, loc, ale.replace('X', '*'))
                else:
                    return ''
            if re.search('INS|DEL|DUP', mtype):
                loc, ale = other_info.split('|')
                if loc and ale:
                    return '{0}.{1}{2}{3}'.format(level, ale, loc, mtype.lower())
                else:
                    return ''
            if re.search('FS', mtype):
                ref, loc, ale, num = other_info.split('|')
                if ref and loc:
                    return '{0}.{1}{2}{3}{4}'.format(level, ref, loc, ale.replace('X', '*'), num)
                else:
                    return ''


def write_in_file(a):
    with open('tmp', 'w') as f:
        for i in print_sm(a):
            p = [j if isinstance(j, str) else str(j) for j in i]
            p.insert(4, mutation_string_normalizition(p[2], p[3]))
            f.write('\t'.join(p) + '\n')
