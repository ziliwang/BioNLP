# NLP tools
part|tool
---|---
Pos|HarvardNLPcore-2016-12
word tokenizer|NLTK
sentence tokenizer|myself-defined
dependency tree|HarvardNLPcore-2016-12

# text preprocess
1. Irregular/Imbalanced variant forms including unex- pected space in between variant forms.
2. Special symbols: tmVar and EMU misidentified minus (-), plus (+) and asterisk ( * ) symbols.
3. PDF parse error: 2365A.T which should be 2365A>T

reference:
[Lee, K., Lee, S., Park, S., Kim, S., Kim, S., Choi, K., … Kang, J. (2016). BRONCO: Biomedical entity Relation ONcology COrpus for extracting gene-variant-disease-drug relations. Database, 2016, 1–13.](https://doi.org/10.1093/database/baw043)

# NER
entry|tools
---|---
gene|GNormPlus
disease|DNorm
variants|EMU and tmVar

# Normalization
entry|database
---|---
gene|hgnc(2017.1)
disease|CTD_disease_2013
variants|hgvs

# IR
## mutation-gene
### evidence rank
```
for article in articles:
    r3_candidate = [genes entries in titile and keywords]
    r4_candidate = most_common([genes in abstract])
    r5_candidate = most_common([genes in result])
    for mutation in article:
        r0_candidate = [genes in same sentence]
        r1_candidate = [genes in previous sentence and next sentence]
        r6_candidate = [previous gene]
        for level in [r0_candidate, r1_candidate, ... r6_candidate]:
            if transvar_verify(mutation, level.genes):
                if level == 'r0_candidate' and level.varified_gene_number > 1:
                    return mutation, shortest_in_dt_tree(level.verified_genes)
                else:
                    return mutation, level.verified_genes
```
## mutation-disease
### topic inference
```
for article in articles:
    for diseases in [tittle or keywords, abstract]:
        next if diseases == 0
        most_common_parent_diseases = integrat diseases base mesh tree
        for disease in most_common_parent_diseases:
            for mutation, gene in mutation_gene_pairs:
                yied disease, mutation, gene
        break
```

# Pipeline
```@mermaid
graph TB;
    A[PubMed];
    B[local database];
    C[GNorm input];
    D[DNorm input];
    E[tmVar input];
    F[EMU input];
    AA(GNormPlus);
    AB(DNorm);
    AC(tmVar);
    AD(EMU);
    BA[Article_obs];
    BB(normalization);
    BC(IR: mutation-gene evidence rank);
    BD(IR:topic disease);
    CA[pmid, M, G, evidence]
    CB[pmid, D]
    CC[pmid, D, M, G, evidence]
    A-->B;
    B-->C;
    B-->D;
    B-->E;
    B-->F;
    C-->AA;
    D-->AB;
    E-->AC;
    F-->AD;
    AA-->BA;
    AB-->BA;
    AC-->BA;
    AD-->BA;
    BA-->BB;
    BB-->BC;
    BB-->BD;
    BC-->CA;
    BD-->CB;
    CA-->CC;
    CB-->CC;
```
