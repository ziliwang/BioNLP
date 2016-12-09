## 储存位置
~/BioNLP/Tools/lib/MutationNorm.py
## 任务
通过基因名加cDNA层面上的碱基改变，通过transVar进行查询，并将transVar的输出转化为VCF格式，再将结果与原始记录对比。
## 流程
#### 1. xlsx中提取碱基改变信息
调用MutationNorm模块extract_mutation_string获取目标文件夹所有gene.xlsx文件，并提取其中第五列中的cDNA层面上的碱基改变，输出（文件前缀， 行， 剪辑改变），将此存储至json文件中。下一步调用。
#### 2. transVar批量查询
将上步结构作为参数输入至transvar_bath_query函数，批量进行transVar，输出为（文件前缀，行号， transVar结果）
#### 3. 解析transVar结果
调用parser_transver_output，将上步结果作为输入，输出为filename, rownum, [mutation_type, position, ref_ale, alt_ale, chr],
           raw_transcript_name, mutation_description， i.e.:
```python
('SLC26A5', '3656', ('snv', '103054757', 'C', 'G', 'chr7'),
          'NM_206885 (protein_coding)', 'chr7:g.103054757C>G/c.293-1198G>C/.')
```
### 4. samtools校正
将上步结果作为samtools_check的输入，通过samtools补充部分位点信息，输出为上步补充后的的结果，数据结构未变
### 5. 转化为字典
将上步结果的数据结构转换为字典
### 6. 对比原始数据，并输出对比结果
将上步结果输入 processdic 函数， 得到对比后的txt文件
## 结构说明
### PASS
原始记录与transVar结果一致
### FAILED
原始记录与transVar结果不一致，主要原因包括：1、突变位点不一致；2、原始记录unalign；3、原始记录部分信息缺失;4、突变位点在反链上；
### Invalid
原始记录中使用的转录本未在transVar中找到。
### NoData
  1. 原始记录中未描述cDNA的突变形式，或是其他原因未能匹配到相应信息（6条）
  2. transVar中没有相应记录，原因推测是原始描述有误。

## 注意
原始文件名更换：

更换前  |  更换后
---------|----------
ADGRV1|GPR98
PD2D7| PDZD7
POU3F4| OTC
SIPR2|S1PR2
PAX3| SOX10

有俩文件格式与预期格式不一致：SOX10, PAX3, POU3F4
