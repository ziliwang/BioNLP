## pubmed下载数据
###### 1. pubmed检索相关关键字，记录下URL中term=后的字符串(query_terms)。
![QQ截图20170309170759](/assets/QQ截图20170309170759.png)
###### 2. 打开ipython, 调用esearch_efetch包。
```python
>>> from esearch_efetch import Esearch_Efetch
>>> d = Esearch_Efetch(terms='query_terms')  # query_terms为第一步中的字符串
>>> d.efecth('output_prefix')  # output_prefix为输出文件前缀
```
如果中途网络中断了，可以通过查了中断处，可以通过以下进行恢复。
```python
>>> d.efecth('output_prefix', start=1000)  # start参数为中断位置
```
## 导入到postgresql数据库
```shell
$ python PubMedParser.py -i path_to_xml_dir -d psql_db
```
## 生成NER输入文件
```python
>>> from preprocess import db
>>> dbobj = db('database_name')  # database_name: postgresql database name
>>> dbobj.get_text('project_name')  # 生成NER输入文件, project_name为前缀
```
## NER
### tmVar
将`tmvarI.*`结尾的文件复制到tmVar的input目录下，执行
```
$ bash tmVar.sh
```
output目录输出结果文件。
### EMU
对`EMU.*`文件逐个执行。
```
$ perl path/to/EMUv1.0.19.pl -i 文件  
```
`EMU_1.19_HUGO_*.EMU.*`为输出文件。
### GNormPlus
将`GNormPlusI.*`结尾的文件复制到GNormPlus的input目录下，执行
```
$ bash GNormPlus.sh
```
output目录输出结果文件。
### DNorm
对`DNormI.*`文件，逐个执行
```
$ ./RunDNorm.sh config/banner_NCBIDisease_TEST.xml data/CTD_disease.tsv output/simmatrix_NCBIDisease_e4.bin 输入 输出
```
## NER输出文件解析
将所有输出文件放在一个目录下，打开ipython，执行：
```python
>>> from ner import output_parser
>>> article_obj_list = output_parser('./', 'project_name')
```
## transvar排错
### 产生transVar输入文件
```python
>>> from link import link_to_transvar
>>> link_to_transvar(article_obj_list)
```
这将生成`pass_to_transvar_canno`和`pass_to_transvar_panno`
### transVar 注释
```
$ sort -u pass_to_transvar_canno > u.pass_to_transvar_canno
$ sort -u pass_to_transvar_panno > u.pass_to_transvar_panno
$ transvar canno -l u.pass_to_transvar_canno --noheader --refseq --ccds --ensembl > pass_to_transvar_canno.out
$ transvar panno -l u.pass_to_transvar_panno --noheader --refseq --ccds --ensembl > pass_to_transvar_panno.out
```
`pass_to_transvar_panno.out`和`pass_to_transvar_canno.out`为输出文件。
## mutation-gene-disease关系提取
ipython下执行：
```python
>>> from mutation_gene import rank_with_dt
>>> rank_with_dt(article_obj_list, 'transvar_output_dir')
>>> from mutation_to_disease import topic_anno
>>> topic_anno(article_obj_list, 'gene-mutation.csv')
```
