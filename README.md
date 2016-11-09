# BioNLP-GDV
## literature resource
- PubMeb ( the literature frame)
- PMC (the free text)

### 1. PubMeb XML to PostgreSQL
PubMedPortable Toolkit: PubMedDB.py and PubMedParser.py
```bash
sudo su - postgres
psql -f PATH/TO/auth.sql
exit
python PubMedParser.py -c -d DATABASE -i XML
```
### 2. Padding PMCID by NCBI API
### 3. Download free full text from PMC by NCBI API
### 4. Parse the xml of free full text

# test set
*  Cytoplasmic Mislocalization of POU3F4 Due to Novel Mutations Leads to Deafness in Humans and Mice
*  Novel domain-specific POU3F4 mutations are associated with X-linked deafness: examples from different populations

-----------------
# Plan
## 11/3
target | status | note
-------|--------|-------
PubMebParser convert python2 to python3 | pass | hand over in engineering
convert PMID to PMCID | finished  |
Download full text from PMCID | finished | Some text un-available in NCBI OA API
## 11/4
target | status | note
-------|--------|------
fengci | | |
NER | | SNP:http://rockt.github.io/SETH/; phenotypes(disease):http://phenominer.mml.cam.ac.uk/ |
relation recognition | | |
