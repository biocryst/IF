import os
from Bio import PDB
from Bio.SeqUtils import seq1
import warnings
import numpy as np
import cPickle
from collections import defaultdict
from difflib import SequenceMatcher
from sklearn.neighbors.kde import KernelDensity


def get_CA_distances(model, selection):
    (pdb_id, chains, start_res, end_res, sequence1, sequence2,_,_) = selection

    chain1 = model[chains[0]]
    chain2 = model[chains[1]]
    CA_distances = []
    for res1_id,res2_id,aminoacid1,aminoacid2 in zip(range(start_res[0],end_res[0]+1), range(start_res[1],end_res[1]+1),sequence1,sequence2):
        try:
            res1 = chain1[res1_id]
        except:
            res1 = chain1[('H_MSE', res1_id, ' ')]
        try:
            res2 = chain2[res2_id]
        except:
            res2 = chain2[('H_MSE', res2_id, ' ')]
        if aminoacid1 != 'X':
            assert aminoacid1 == seq1(res1.resname)
        if aminoacid2 != 'X':
            assert aminoacid2 == seq1(res2.resname)
        dist = np.linalg.norm(res1['CA'].coord-res2['CA'].coord)
        CA_distances.append(dist)

    return CA_distances

def similar(a, b):
    return SequenceMatcher(None, a, b).ratio()


f = open('cc_dataset_mmol_all').read().splitlines()

records = [line.split() for line in f]
ccoils = []

for c1,c2 in zip(records[0::2], records[1::2]):
    assert c1[0] == c2[0]
    filename = c1[0]
    chains = (c1[1],c2[1])
    start_res = (int(c1[2]),int(c2[2]))
    end_res = (int(c1[3]),int(c2[3]))
    if len(c1[4]) != len(c2[4]):
        continue
    sequence1 = c1[4]
    sequence2 = c2[4]
    ccoils.append((filename, chains, start_res, end_res, sequence1,sequence2))

print len(ccoils)

ccoils.sort(key = lambda s: len(s[4]),reverse=True)

pdb_path = "MMOL"

warnings.filterwarnings("ignore")

ccoils_filt = []

for item1 in ccoils:
    filename1, chains1, start_res1, end_res1, sequence11, sequence12 = item1
    unique1 = True
    unique2 = True
    for item2 in ccoils_filt:
        filename2, chains2, start_res2, end_res2, sequence21, sequence22,_,_ = item2
        if sequence11 in sequence21 or sequence11 in sequence22:
            unique1 = False
        if sequence12 in sequence21 or sequence12 in sequence22:
            unique2 = False
        if unique1 and (similar(sequence11,sequence21)>0.8 or similar(sequence11,sequence22)>0.8):
            unique1 = False
        if unique2 and (similar(sequence12,sequence21)>0.8 or similar(sequence12,sequence22)>0.8):
            unique2 = False
        if not unique1 and not unique2:
            break

    if unique1 or unique2:
        ccoils_filt.append((filename1, chains1, start_res1, end_res1, sequence11, sequence12,unique1,unique2))

print len(ccoils_filt)
print np.sum([len(s[4]) for s in ccoils_filt if s[6]])
print np.sum([len(s[5]) for s in ccoils_filt if s[7]])

ca_distances_all = []
for item in ccoils_filt:
    (filename, chains, start_res, end_res, sequence1, sequence2, unique1, unique2) = item
    structure = PDB.PDBParser().get_structure(filename, os.path.join(pdb_path,filename))
    model = structure[0]
    try:
        CA_distances = get_CA_distances(model, item)
        ca_distances_all.append(CA_distances)
    except:
        print 'problem with', item
        continue


aa_distr = defaultdict(list)

for cc_item, cc_dists in zip(ccoils_filt,ca_distances_all):
    filename, chains, start_res, end_res, sequence1, sequence2,unique1,unique2 = cc_item
    if unique1:
        for aa,dist in zip(sequence1,cc_dists):
            if aa =='X':
                aa = 'M'
            aa_distr[aa].append(dist)
    if unique2 and sequence2 != sequence1:
        for aa,dist in zip(sequence2,cc_dists):
            if aa =='X':
                aa = 'M'
            aa_distr[aa].append(dist)

aa_distr_kde = defaultdict(list)
for aa,dists in aa_distr.iteritems():
    data = np.array(dists).reshape(-1, 1)
    kde = KernelDensity(kernel='gaussian', bandwidth=1).fit(data)
    aa_distr_kde[aa]=kde


cPickle.dump(aa_distr_kde, open('cc_aa_distance_distr.pkl', "wb"))
