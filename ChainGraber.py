# rabid-juniper

import sqlite3
import pypdb
import pandas as pd

conn = sqlite3.connect('test_proteins2.db')
c = conn.cursor()
first_k = c.execute("SELECT s_m_title,s_bioluminate_Antigen_Type, s_bioluminate_Antigen_Chain, s_bioluminate_Antigen_Seq FROM 'PrimeStructureDB_DataTable'")
k = c.fetchall()
k = [x for x in k if x[1] is not None]
    
def get_pdb():
    k02 = [k0[2] for k0 in k]
    k00 = [k0[0] for k0 in k]
    pdbs = [pypdb.get_pdb_file(k_i, filetype='pdb', compression=False).split('\n') for k_i in k00[:4]] # проверила на первых 4 шт.
    return pdbs

get_pdb1 = get_pdb()

def pdb_info(get_pdb1):
    keys = []
    values = []
    for i_ in range(len(get_pdb1)):
        pdbs2 = get_pdb1[i_]
        ind = [pdbs2.index(line) for line in pdbs2 if line[18:19] == k02[i_]]
        for i in range(ind[0]-2, ind[0]+1):
            j = pdbs2[i]
            keys += [j[10:j.find(':')]]
            values += [j[j.find(':')+1:j.find(';')]]
    values = [values[i:i+3] for i in range(0,len(values),3)]
    return values, keys[:3]

pdb_inf = pdb_info(get_pdb1)

frame = pd.DataFrame([i for i in pdb_inf[0]], columns=pdb_inf[1])
print(frame)
