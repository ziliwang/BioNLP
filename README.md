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





# Plan
## 11/3
*  PubMebParser convert python2 to python3
*  convert PMID to PMCID
*  Download full text from PMCID
