#!/usr/bin/env python3
from pdfminer import encodingdb
import re
import os


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


def Batch_parse_xml(dir_path):
    '''
    batch parser pdf into xml, and parser xml into a single string.
    return: [[pmid, string], ...]
    warn: pdf format need as pmid.pdf
    '''
    output = []
    for root, dirs, files in os.walk(dir_path):
        for f in files:
            try:
                name, suffix = f.split('.')
            except BaseException as e:
                print('{0} passed: {1}'.format(f, e))
                continue
            if suffix in ['pdf']:
                filename = os.path.join(root, f)
                res = PDFParser(filename)
                res.parse()
                output.append([name, res.text])
    os.popen('rm {0}/*.xml'.format(dir_path))
    return(output)
