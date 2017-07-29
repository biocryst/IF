import urllib
import os

f = open('cc_dataset_mmol_all').read().splitlines()

pdb_ids = [id_str[0:4] for id_str in f]
filenames = [id_str.split()[0] for id_str in f]


link_pref = "http://coiledcoils.chm.bris.ac.uk/ccplus/files/data/{0}/{1}/structures/{2}"
file_pref = "MMOL/"
pdb_pref = "http://www.rcsb.org/pdb/files/{0}"

if not os.path.exists(file_pref):
        os.makedirs(file_pref)

for pdb_id,filename in zip(pdb_ids,filenames):

    link_pdb = link_pref.format(pdb_id[1:3],pdb_id,filename)
    file_pdb = file_pref + filename

    if os.path.exists(file_pdb):
        continue

    if file_pdb.endswith('.pdb'):
        link_pdb = pdb_pref.format(filename)
        urllib.urlretrieve(link_pdb, file_pdb)
        continue

    try:
        urllib.urlretrieve(link_pdb, file_pdb)
        print pdb_id
    except:        
        print pdb_id, ' error'

    
