# API Reference
version 1.1
## package **esearch_efetch**
### class **Esearch_Efetch**
#### example
```python
>>> from esearch_efetch import Esearch_Efetch
>>> d = Esearch_Efetch(db='pubmed', terms='deafness+or+%22hearing+loss%22', remax=1)
>>> d.esearch()
>>> d.efecth('deafness', start=0, retmax=1000)
```
#### method
**\_\_init\_\_**(*db=*, *terms=*, *retmax=*)
parameters|returns
---|---
**db=**: string, default `'pubmed'`, NCBI database; **terms=**: string, default `'deafness+or+%22hearing+loss%22'`, query string; **retmax=**: int, default `1`, return max limit|instance
**esearch**()
parameters|returns
---|---
*null* | *null*
**efecth**(*output_prefix*, *start=*, *retmax=*)
parameters|returns
---|---
**output_prefix**: string, required, prefix of output files; **start=**: int, default 0, start record;**retmax=**: int, default 1000, max record numbers per batch.|*null*
#### attribute
**q_url**
|returns|
|---|
|string, query url string|
-------
## package **Base**
数据结构
### class **BioEntry**
#### example
```python
>>> from Base import BioEntry
>>> BioEntry('1', '11', 'gene', 'cdh23', '2312', 'CDH23')
<Base.BioEntry at 0x4bdf2b0>
>>> a = BioEntry('1', '11', 'gene', 'cdh23', '2312', 'CDH23')
>>> a.start
'1'
>>> a.end
'11'
>>> a.type
'gene'
>>> a.string
'cdh23'
>>> a.id
'2312'
>>> a.norm
'CDH23'
```
#### method
**\_\_init\_\_**(*start*, *end*, *type*, *string*, *id*, *norm*)
parameters|returns
---|---
**start**: string, required, start position of entity string in text. **end**: string, required, end position of entity string in text. **type**: string, required, entity type. **string**: string, required, entity string in text. **id**: string, required, entity id. **norm**: string, required, normalization entity. | instance, have several attribute.
#### attribute
**start**
|returns|
|---|
|string, start position of entity|
**end**
|returns|
|---|
|string, end position of entity|
**type**
|returns|
|---|
|string, entity type|
**string**
|returns|
|---|
|string, raw string of entity|
**id**
|returns|
|---|
|string, entity id|
**norm**
|returns|
|---|
|string, entity normalization|
### class **TK**
#### example
```python
>>> from Base import TK
>>> TK('the tittle and the keywords')
<Base.TK at 0x67878f0>
>>> tk = TK('a mutation c.216A > T')
>>> tk.add_e(a)
```
#### method
**\_\_init\_\_**(*text*)
parameters|returns
---|---
**text**: string, required, text of title and keywords|instance
**add_e**(*bioentry*)
parameters|returns
---|---
**bioentry**: BioEntry instance, required, add BioEntry into title and keywords.| *null*
#### attribute
**text**
|returns|
|---|
|string, text|
**bioe**
|returns|
|---|
|list, BioEntry instance list|
### class **Abstract**
#### example
```python
>>> from Base import Abstract
>>> Abstract('this is abstract, there is a mutation c.213A>T in P53.')
<Base.Abstract at 0xfb1c10>
>>> abstr = Abstract('this is abstract, there is a mutation c.213A>T in P53.')
>>> abstr.add_e(a)
>>> abstr.sort_entry()
```
#### method
**\_\_init\_\_**(*text*)
parameters|returns
---|---
**text**: string, required, text of abstract|instance
**add_e**(*bioentry*)
parameters|returns
---|---
**bioentry**: BioEntry instance, required, add BioEntry into abstract.| *null*
**sort_entry**()
sorted the entities.
parameters|returns
---|---
*null*| *null*
#### attribute
**bioe**
|returns|
|---|
|list, BioEntry instances list|
**sorted**
|returns|
|---|
|bool|
**text**
|returns|
|---|
|string, text|
### class **Results**
#### example
```python
>>> from Base import Results
>>> Results('This is results, there is a mutation c.213A>T in P53.')
<Base.Results at 0xfb3650>
>>> res = Results('This is results, there is a mutation c.213A>T in P53.')
>>> res.add_e(a)
>>> res.sort_entry()
```
#### method
**\_\_init\_\_**(*text*)
parameters|returns
---|---
**text**: string, required, text of results|instance
**add_e**(*bioentry*)
parameters|returns
---|---
**bioentry**: BioEntry instance, required, add BioEntry into results.| *null*
**sort_entry**()
sorted the entities.
parameters|returns
---|---
*null*| *null*
#### attribute
**bioe**
|returns|
|---|
|list, BioEntry instances list|
**sorted**
|returns|
|---|
|bool|
**text**
|returns|
|---|
|string, text|
### class **Article**
#### example
```python
>>> from Base import Article
>>> Article('11121112')
<Base.Article at 0xfb5c50>
>>> article = Article('11121112')
>>> article.add_tk(tk)
>>> article.add_abstract(abstr)
>>> article.add_results(res)
>>> event_abstract = article.get_abstract_relations()
>>> event_result =  article.get_results_relations()
```
#### method
**\_\_init\_\_**(*pmid*)
parameters|returns
---|---
**pmid**: string, required, pmid|instance
**add_tk**(*tk*)
parameters|returns
---|---
**tk**: instance, required, title and keywords instance. | *null*
**add_abstract**(*abstract*)
parameters|returns
---|---
**abstract**: instance, required, abstract instance. | *null*
**add_results**(*results*)
parameters|returns
---|---
**results**: instance, required, results instance. | *null*
**get_abstract_relations**()
get the mutation-gene relations in abstract.
parameters|returns
---|---
*null* | list, by `(mutation_entity, [relation_level, other_entity, evidence])`
**get_results_relations**()
get the mutation-gene relations in results.
parameters|returns
---|---
*null* | list, by `(mutation_entity, [relation_level, other_entity, evidence])`
### class **w_tree**
#### example
```python
>>> from Base import w_tree
>>> mytree = w_tree('root')
>>> child = w_tree('child')
>>> mytree.add_child(child, 32)
>>> mytree.children_set()
(<Base.w_tree at 0xfb5c51>)
>>> mytree.get('child')
<Base.w_tree at 0xfb5c51>
>>> mytree.add_child_weight('child', '3.2')
```
#### method
**\_\_init\_\_**(*value*)
parameters|returns
---|---
**value**: string, required, node value|instance, a node
**add_child**(*node*, *weight*)
parameters|returns
---|---
**node**: instance, required, node instance. **weight**: float, required, the weight of edge.|*null*
**children_set**()
parameters|returns
---|---
*null* |list, the value of children nodes.
**get_child**(*value*)
parameters|returns
---|---
**value**: string, required, node value.| instance/False, node instance.
**add_child_weight**(*value*, *weight*)
parameters|returns
---|---
**value**: string, required, node value. **weight**: float, required, addtional weight.|*null*
#### attribute
**value**
|returns|
|---|
|string, node value|
**children**
|returns|
|---|
|list, a list of node instance|
### class **mesh_tree**
#### example
```python
>>> from Base import mesh_tree
>>> my_mesh_tree = mesh_tree()
>>> my_mesh_tree.add('C06.12.34.45', 1)
>>> my_mesh_tree.add('C06.12.34.47', 2)
>>> my_mesh_tree.add('C06.12.58.77', 1)
>>> my_mesh_tree.show()
graph TB
ROOT--1.0-->C06
C06--1.0-->C06.12
C06.12--0.75-->C06.12.34
C06.12.34--0.25-->C06.12.34.45
C06.12.34--0.5-->C06.12.34.47
C06.12--0.25-->C06.12.58
C06.12.58--0.25-->C06.12.58.77
>>> my_mesh_tree.get_weight_path()
<generator object mesh_tree.get_weight_path at 0x00FE4F30>
>>> list(my_mesh_tree.get_weight_path())
[('C06.12.34.47', 3.25)]
```
#### method
**\_\_init\_\_**()
parameters|returns
---|---
*null* | *null*
**add**(*node*, *weight*)
parameters|returns
---|---
**node**: string, required, node in mesh tree. **weight**: float, required, the edge weight.|*null*
**show**
parameters|returns
---|---
*null* | string
**get_weight_path**
parameters|returns
---|---
*null* | list， weightest path
#### attribute
**root**
|returns|
|---|
|instance, root node instance|
----
## package **pdf_parser**
### class **PDFParser**
解析PDF文件类，提取PDF文件中result部分信息，若无result信息，提取全文。
#### example
```python
>>> from pdf_parser import PDFParser
>>> pdf_obj = PDFParser('PMID.pdf')
>>> pdf_obj.parse()
>>> pdf_obj.text
'Textual information in the world can be broadly categorized into two main types: facts and opinions. Facts are objective expressions about entities, events and their properties. Opinions are usually subjective expressions that describe people’s sentiments, appraisals or feelings toward entities, events and their properties. The concept of opinion is very broad.'
```
#### method
**\_\_init\_\_**(*pdf_file*)
parameters|returns
---|---
**pdf_file**: string, required, pdf file path|a pdf class which have **parse()** method
**parse**()
parsing the pdf.
parameters|returns
---|---
*null* |*null*
#### attribute
**text**
parameters|returns
---|---
required called **parse**() | the text from pdf
### function **Batch_parse_pdf**
#### example
```python
>>> from pdf_parser import Batch_parse_pdf
>>> pdf_text_ls = Batch_parse_pdf('/path/to/pdf/dir/')
>>> from preprocess import pdf_ls_to_input
>>> pdf_ls_to_input(pdf_text_ls, '/output/prefix')
```
**Batch_parse_pdf**(*dir_path*)
parameters|returns
---|---
**dir_path**: string, required, the path of folder which stored pdfs | a list, `[<pmid(string), text(string)>...]`
-----
## package **preprocess**
预处理模块，导出NER输入文件
### class **db**
从postgresql数据库重导出tmVar，EMU， DNorm，GNormPlus输入文件。
#### example
```python
>>> from preprocess import db
>>> dbobj = db('database_name')
>>> dbobj.get_text('project_name')
```
#### method
**\_\_init\_\_**(*database_name*)
parameters|returns
---|---
**database_name**: string, required, the name of database in postgresql.|A db object
**get_text**(*output_file_prefix*)
parameters|returns
---|---
**output_file_prefix**: string, required, a prefix for output files, suggest seting by project name.|tmVar input: *output_file_prefix.tmvarI*; EMU input: *output_file_prefix.EMU*; DNorm input: *output_file_prefix.DNormI*; GNormPlus input: *output_file_prefix.GNormPlusI*
### function **pdf_ls_to_input**
#### example
```python
>>> from pdf_parser import Batch_parse_pdf
>>> pdf_text_ls = Batch_parse_pdf('/path/to/pdf/dir/')
>>> from preprocess import pdf_ls_to_input
>>> pdf_ls_to_input(pdf_text_ls, 'project_name')
```
**pdf_ls_to_input**(*pdf_text_list*, *output_prefix*)
parameters|returns
---|---
**pdf_text_list**: list, required, the return of pdf_parser.Batch_parse_pdf; **output_prefix**: string, required, a prefix for output files, suggest seting by project name.|tmVar input: *output_file_prefix.tmvarI*; EMU input: *output_file_prefix.EMU*; DNorm input: *output_file_prefix.DNormI*; GNormPlus input: *output_file_prefix.GNormPlusI*
-----
## package **ner**
解析NER
### function **output_parser**
#### example
```python
>>> from ner import output_parser
>>> article_obj_list = output_parser('NER_out_dir', 'project_name')
>>> article_obj_list
[<Base.Article at 0x7f23267254bd9b0>,
<Base.Article at 0x7f23267254bd9b1>,
<Base.Article at 0x7f23267254bd9b2>,
<Base.Article at 0x7f23267254bd9b3>,
<Base.Article at 0x7f23267254bd9b4>,
<Base.Article at 0x7f23267254bd9b5>]
```
**output_parser**(*NER_output_dir*,*prefix*)
parameters|returns
---|---
**NER_output_dir**: string, required, the folder which NER output in. **prefix**: string, required, also named the prefix of the files. | list of Base.Article instances
----
## package **normalization**
归一化
### function **mutation_string_normalizition**
#### example
```python
>>> from normalization import mutation_string_normalizition
>>> mutation_string_normalizition('c.1752 G > A', 'c|SUB|G|1752|A')
'c.1752G>A'
```
**mutation_string_normalizition**(*raw_text*, *tmVar_normal_mutation*)
parameters|returns
---|---
**raw_text**: string, required, raw string from the original text; **tmVar_normal_mutation**: string, required, tmVar normalization string;|string, HGVS format
### class **gene_normor**
#### example
```python
>>> from normalization import gene_normor
>>> gn = gene_normor('path/to/hgnc.txt')
>>> gn.norm('673')
'BRAF'
```
#### method
**\_\_init\_\_**(*hgnc_file*)
parameters|returns
---|---
**hgnc_file**: string, required, hgnc file | instances, which can parse entrez ID to gene symbol
**norm**(*entrez_id*)
parameters|returns
---|---
**entrez_id**: string, required, entrez_id| string, gene symbol
### class **disease_normor**
#### example
```python
>>> from normalization import disease_normor
>>> dn = disease_normor('CTD_disease.csv')
>>> dn.norm('MESH:C535857')
{'name': 'Hecht syndrome',
'node': 'C05.550.150|C05.651.102|C16.131.621.077|C23.888.592.608.750.700'}
```
#### method
**\_\_init\_\_**(*CTD_disease_file*)
parameters|returns
---|---
**CTD_disease_file**: string, required, CTD disease file, csv format. | instances, which can parse disease ID to name and its node in MESH tree
**norm**(*disease_ID*)
parameters|returns
---|---
**disease_ID**: string, required, disease ID|dict, `{'name':name, 'node':nodes}`
----
## package **link**
transvar排除假阳性
### function **link_to_transvar**
#### example
```python
>>> from link import link_to_transvar
>>> link_to_transvar(article_obj_list)
```
**link_to_transvar**(*article_obj_list*)
parameters|returns
---|---
**article_obj_list**: list, required, a list of article instances|*null*, but ouput transvar input file
### class **linker**()
#### example
```python
>>> from link import linker
>>> l = liner('transvar_output_dir')
>>> l.is_link('CDH23:c.123A>T')
True
>>> l.is_link('CDH23:p.A213T')
False
>>> l.parse('CDH23:c.123A>T')
'XM_011540046:chr10:g.71566840A>T/c.123A>T/p.Q41H|XM_011540047:chr10:g.71638501A>T/c.123A>T/p.E41D'
>>> l.parse('CDH23:p.A213T')
False
```
#### method
**__init__**(*transvar_output_dir*)
parameters|returns
---|---
**transvar_output_dir**: string, required, the transvar output folder path.|instance, which can parse `gene:hgvs` into standard hgvs.
**is_link**(*gene_hgvs_string*)
parameters|returns
---|---
**gene_hgvs_string**: string, required, format:`gene:hgvs`, i.e. `CDH23:c.234A>T`|Bool
**parser**(*gene_hgvs_string*)
parameters|returns
---|---
**gene_hgvs_string**: string, required. |string, HGVS format
----
## package **mutation_gene**
基因突变关系挖掘
### function **rank_with_dt**
#### example
```python
>>> from mutation_gene import rank_with_dt
>>> rank_with_dt(article_obj_list, 'transvar_output_dir')
```
**rank_with_dt**(*article_obj_list*, *transvar_output_dir*)
parameters|returns
---|---
**article_obj_list**: list, required, article instances. **transvar_output_dir**: string, required, transvar output folder|*null*, but output the gene - mutation - evidence document.
----
## package **mutation_to_disease**
突变疾病关系挖掘
### function **topic_anno**
#### example
```python
>>> from mutation_to_disease import topic_anno
>>> topic_anno(article_obj_list, 'gene-mutation.tsv')
```
**topic_anno**(*article_obj_list*, *gene_mutation_file*)
parameters|returns
---|---
**article_obj_list**: list, required, article instances list; **gene_mutation_file**: string, required, gene-mutation file|*null*, but output disease-gene-mutation tsv file.
## package **dependency_tree**
依存文法解析
### function **dependency_tree_shortest_path**
#### example
```python
>>> from dependency_tree import dependency_tree_shortest_path
>>> qdic = dependency_tree_shortest_path(qdic)
```
**dependency_tree_shortest_path**(*qdic*)
parameters|returns
---|---
**qdic**: dict, required, preform a data structure. | dict, qdic
