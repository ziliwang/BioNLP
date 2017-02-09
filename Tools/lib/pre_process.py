#!/usr/bin/env python3
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
import re
from hgnc_db import HGNC
import pickle


def query(id):
    db = 'kownledge_db'
    con = 'postgresql://s:123456@localhost/'+db
    engine = create_engine(con)
    Session = sessionmaker(engine)
    session = Session()
    result = session.query(HGNC.symbol).filter_by(entrez_id=id).one()
    session.commit()
    print(result[0])
    return(result[0])


def gene_norm(file):
    normgene = []
    with open(file, 'r') as f:
        for i in f.readlines():
            items = i.split('\t')
            if items[6] == 'Gene':
                if re.match('\d+$', items[8]):
                    try:
                        items.insert(9, query(items[8]))
                    except:
                        items.insert(9, None)
                else:
                    items.insert(9, None)
            else:
                items.insert(9, 'wait')
            normgene.append(items)
    with open('gene_norm', 'wb') as f:
        pickle.dump(normgene, f)
