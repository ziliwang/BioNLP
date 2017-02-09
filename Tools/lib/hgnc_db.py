#!/usr/bin/env python3
from sqlalchemy import Column, Integer, String
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import create_engine, engine, MetaData
from sqlalchemy.orm import sessionmaker
import sqlalchemy
Base = declarative_base()


class HGNC(Base):
    """hgnc complete data"""
    __tablename__ = 'hgnc'

    hgnc_id = Column(Integer, primary_key=True)  # col 1
    symbol = Column(String, unique=True, index=True)         # col 2
    name = Column(String, unique=True)           # col 3
    locus_group = Column(String)                 # col 4
    locus_type = Column(String)                  # col 5
    # status = Column(String)                      # col 6
    # location = Column(String)                    # col 7
    location_sortable = Column(String)           # col 8
    alias_symbol = Column(String)                # col 9
    alias_name = Column(String)                  # col 10
    prev_symbol = Column(String)                 # col 11
    prev_name = Column(String)                   # col 12
    gene_family = Column(String)                 # col 13
    gene_family_id = Column(String)              # col 14 | split
    # date_approved_reserved = Column(String)      # col 15, time
    # date_symbol_changed = Column(String)         # col 16, time
    # date_name_changed = Column(String)           # col 17, time
    # date_modified = Column(String)               # col 18, time
    entrez_id = Column(Integer, index=True, unique=True)      # col 19, ncbi gene id
    ensembl_gene_idvega_id = Column(String, index=True)
    # vega_id = Column(String)
    ucsc_id = Column(String, index=True, unique=True)
    # ena = Column(String)
    refseq_accession = Column(String, index=True)
    ccds_id = Column(String)
    uniprot_ids = Column(String, index=True)  # col 26
    pubmed_id = Column(String)
    # mgd_id = Column(String)
    # rgd_id = Column(String)
    # lsdb = Column(String)
    # cosmic = Column(String)
    omim_id = Column(String, unique=True, index=True)

    def __repr__(self):
        print('hgnc item remove tail column')

if __name__ == '__main__':
    db = 'kownledge_db'
    con = 'postgresql://s:123456@localhost/'+db
    my_engine = create_engine(con)
    try:
        HGNC.__table__.drop(my_engine)
    except:
        pass
    Base.metadata.create_all(my_engine)
    fl = 'hgnc_complete_set.txt'
    Session = sessionmaker(my_engine)
    session = Session()
    with open(fl, 'r') as f:
        for i in f.readlines():
            items = i.split('\t')
            if items[5] != 'Approved':
                continue
            record = HGNC(hgnc_id=items[0].split(':')[1],
                          symbol=items[1],
                          name=items[2],
                          locus_group=items[3],
                          locus_type=items[4],
                          # status=items[5],
                          # location=items[6],
                          location_sortable=items[7],
                          alias_symbol=items[8],
                          alias_name=items[9],
                          prev_symbol=items[10],
                          prev_name=items[11],
                          gene_family=items[12],
                          gene_family_id=items[13],
                          #  date_approved_reserved=items[14],
                          #  date_symbol_changed=items[15],
                          #  date_name_changed=items[16],
                          #  date_modified=items[17],
                          entrez_id=items[18] if items[18] else None,
                          ensembl_gene_idvega_id=items[19] if items[19] else None,
                          #  vega_id=items[20],
                          ucsc_id=items[21] if items[21] else None,
                          #  ena=items[22],
                          refseq_accession=items[23] if items[23] else None,
                          ccds_id=items[24] if items[24] else None,
                          uniprot_ids=items[25] if items[25] else None,
                          pubmed_id=items[26] if items[26] else None,
                          #  mgd_id=items[27],
                          #  rgd_id=items[28],
                          #  lsdb=items[29],
                          #  cosmic=items[30],
                          omim_id=items[31] if items[31] else None)
            session.add(record)
            session.commit()
