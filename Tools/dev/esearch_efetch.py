import urllib.request
import re
import xml.etree.ElementTree as ET


class Esearch_Efetch():

    def __init__(self, db='pubmed', terms='deafness+or+%22hearing+loss%22',
                 retmax=1):
        base_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?'
        db_url = 'db=' + db
        terms_url = 'term=' + terms
        type_url = 'datetype=edat'
        retmax_url = 'retmax=' + str(retmax)
        add_url = 'usehistory=y'
        self.q_url = base_url + '&'.join([db_url, terms_url, type_url,
                                         retmax_url, add_url])
        self.__webenv = False

    def esearch(self):
        try:
            _xml = urllib.request.urlopen(self.q_url).read().decode('utf-8')
            root = ET.fromstring(_xml)
            self.__webenv = root.find('WebEnv').text
            self.__count = root.find('Count').text
        except BaseException as e:
            print('Esearch err: {0}'.format(e))

    def efecth(self, output_prefix, start=0, retmax=1000):
        base_url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?'
        db_url = 'db=pubmed&query_key=1'
        retmax_url = 'retmax=' + str(retmax)
        self.esearch()
        for i in range(start, int(self.__count), retmax):
            _write = False
            _try = 0
            while not _write and _try < 3:
                start_url = 'retstart=' + str(i)
                add_url = 'retmode=xml&WebEnv=' + self.__webenv
                q_url = base_url + '&'.join([db_url, start_url, retmax_url,
                                             add_url])
                try:
                    urllib.request.urlretrieve(q_url, output_prefix + '_' +
                                               str(i) + '.xml')
                    _write = True
                except BaseException as e:
                    _try += 1
            if not _write:
                raise ValueError('too many times')
