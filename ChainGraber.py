# Этот код должен выводить датафрейм с антигенами, полученными из pdb, в порядке убывания их встречаемости в базе

import sqlite3
import pypdb
import pandas as pd
    
def get_pdb(): # функция извлекает pdb файлы для каждого pdb_id, содержащегося в k00
    conn = sqlite3.connect('test_proteins2.db')
    c = conn.cursor()
    first_k = c.execute("SELECT s_m_title,s_bioluminate_Antigen_Type, s_bioluminate_Antigen_Chain, s_bioluminate_Antigen_Seq FROM 'PrimeStructureDB_DataTable'")
    k = c.fetchall()
    k = [x for x in k if x[1] is not None]
    k02 = [k0[2] for k0 in k] # символ цепи, к которой относится антиген
    k00 = [k0[0] for k0 in k] # pdb_id из test_proteins2.db
    pdbs = [pypdb.get_pdb_file(k_i, filetype='pdb', compression=False).split('\n') for k_i in k00[:4]] # проверила на первых 4 шт.
    return pdbs
get_pdb1 = get_pdb()

def pdb_info(pdb): # функция возвращает списки для построения датафрейма
    keys = []
    values = []
    for i_ in range(len(pdb)):
        pdbs2 = pdb[i_]
        ind = [pdbs2.index(line) for line in pdbs2 if line[18:19] == k02[i_]]
        for i in range(ind[0]-2, ind[0]+1):
            j = pdbs2[i]
            keys += [j[10:j.find(':')]]
            values += [j[j.find(':')+1:j.find(';')]]
    values = [values[i:i+3] for i in range(0,len(values),3)]
    return values, keys[:3]
pdb_inf = pdb_info(get_pdb1)

# содержит инфрмацию об id молекулы, pdb_id, типе цепочки, и название антигена
def antigen_frame(pdb_inf):
    frame = pd.DataFrame([i for i in pdb_inf[0]], columns=pdb_inf[1])
    frame['PDB ID'] = k00[:4] # 4 взялось из 15 строчки
    return frame
frame = antigen_frame(pdb_inf)

def sorted_frame(frame):
    molecule_counts = frame[' MOLECULE'].value_counts()
    for_sorting = [molecule_counts[i] for i in frame[' MOLECULE']]
    frame['MOLECULE COUNTS'] = for_sorting
    after_sorting = frame.sort_values(['MOLECULE COUNTS'], ascending=[False])
    return after_sorting
prevailing_antigens = sorted_frame(frame)

