import os
from Bio import PDB
from Bio.SVDSuperimposer import SVDSuperimposer
from Bio.SeqUtils import seq3
import numpy as np
import multiprocessing
from copy import deepcopy
from pyRMSD import RMSDCalculator

def new_cc(sequence, coords):
    seq_len = min(len(sequence)*5,coords.shape[0]/2)/5
    new_model = PDB.Model.Model(0)
    segid = '    '
    atomname = ['N', 'CA', 'C', 'O', 'CB']
    bfactor = 0
    occupancy = 1
    altloc = ' '
    fullname = [' N  ', ' CA ',' C  ', ' O  ', ' CB ']
    element = ['N', 'C', 'C', 'O', 'C']
    serial_number = 1
    for chain_id in ['A','B']:
        new_chain = PDB.Chain.Chain(chain_id)
        for res_i in range(1,seq_len+1):
            res_id = (' ', res_i, ' ')
            resname = seq3(sequence[res_i-1]).upper()
            new_res = PDB.Residue.Residue(res_id, resname, segid)
            for atom_i in range(5):
                new_atom = PDB.Atom.Atom(atomname[atom_i], coords[serial_number-1], bfactor, occupancy, altloc, fullname[atom_i], serial_number, element[atom_i])
                new_res.add(new_atom)
                serial_number +=1
            new_chain.add(new_res)
        new_model.add(new_chain)
    return new_model


def merge_dimer(coords_list, res_overlap):

    ref_coords = coords_list[0]
    aligned_coords = [deepcopy(coords_list[0])]
    n_atoms_per_res = 5
    n_atoms_mono = int(ref_coords.shape[0]/2)
    msds = []
    for coords, cc_overlap in zip(coords_list[1:],res_overlap):
        n_atoms_overlap = cc_overlap*n_atoms_per_res
        h1_ref = ref_coords[n_atoms_mono-n_atoms_overlap:n_atoms_mono]
        h2_ref = ref_coords[n_atoms_mono:n_atoms_mono+n_atoms_overlap]
        ref_atoms = np.append(h1_ref, h2_ref, axis=0)

        h1 = coords[:n_atoms_overlap]
        h2 = coords[-n_atoms_overlap:]

        sup_atoms = np.append(h1, h2, axis=0)

        sup=SVDSuperimposer()
        sup.set(ref_atoms,sup_atoms)
        sup.run()
        msds.append(sup.get_rms()**2)
        rot,tran = sup.get_rotran()
        coord_new = np.dot(coords, rot) + tran
        aligned_coords.append(coord_new)
        ref_coords = coord_new

    rmsd = np.sqrt(np.sum(msds))

    h1_all = aligned_coords[0][:n_atoms_mono]
    h2_all = aligned_coords[0][n_atoms_mono:]

    for coords,cc_overlap in zip(aligned_coords[1:],res_overlap):
        h1 = coords[:n_atoms_mono]
        h2 = coords[n_atoms_mono:]
        n_atoms_overlap = cc_overlap*n_atoms_per_res
        for ind_overlap in range(cc_overlap):
            weight = (ind_overlap+1)/float(cc_overlap+1)
            for ind_atom in range(n_atoms_per_res):
                ind_shift = ind_overlap*n_atoms_per_res+ind_atom
                coord1_prev = h1_all[-n_atoms_overlap+ind_shift]
                coord1_next = h1[ind_shift]
                h1_all[-n_atoms_overlap+ind_shift] = (1-weight)*coord1_prev+weight*coord1_next
                coord2_prev = h2_all[ind_shift]
                coord2_next = h2[-n_atoms_overlap+ind_shift]
                h2_all[ind_shift] = (1-weight)*coord2_next+weight*coord2_prev
        h1_rest = h1[n_atoms_overlap:]
        h2_rest = h2[:-n_atoms_overlap]

        h1_all = np.append(h1_all,h1_rest,axis=0)
        h2_all = np.append(h2_rest,h2_all,axis=0)

    res_dimer = np.append(h1_all,h2_all,axis=0)
    return res_dimer, rmsd

def best_model(args):
    variants,rms_dict = args
    n_classes = len(variants)
    solutions = [[0]*len(variants[0])]
    track = [[0]*len(variants[0])]

    for i in range(1, n_classes):
        cur_solutions = []
        cur_track = []
        var_from, var_to, rms_matrix = rms_dict[i-1]
        for str_ind in variants[i]:
            scores = []
            for sol_ind,prev_str_ind in enumerate(variants[i-1]):
                s1 = var_from.index(prev_str_ind)
                s2 = var_to.index(str_ind)
                rmsd = rms_matrix[s1,s2]
                scores.append(rmsd+solutions[i-1][sol_ind])

            min_score = min(scores)
            cur_solutions.append(min_score)
            cur_track.append(np.argmin(scores))

        solutions.append(cur_solutions)
        track.append(cur_track)

    sorted_solutions = sorted(solutions[-1])
    min_score =sorted_solutions[0]
    best_trace = []
    last_best = np.argmin(solutions[-1])

    for level_struct, back_pointers,sol_level in reversed(zip(variants, track,solutions)):
        best_trace.append(level_struct[last_best])
        last_best = back_pointers[last_best]
    best_trace.reverse()
    return best_trace, np.sqrt(min_score)

def get_all_coords(cc_coords_all,cc_inds,dimer_overlap, return_rms = False):
    cc_coords_sel = []
    n_classes = len(cc_inds)

    for i in range(n_classes):
        cc_coords_sel.append(cc_coords_all[cc_inds[i]])

    coords, rms = merge_dimer(cc_coords_sel,dimer_overlap)

    if return_rms:
        return coords, rms
    return coords

def write_model(coords,target_seq,filename):
    cc_model = new_cc(target_seq,coords)

    io=PDB.PDBIO()
    io.set_structure(cc_model)
    io.save(filename+".pdb")
    return

def est_logprob_cached(seq1,seq2,aas_list,aa_struct_prob):
    sum_logprob = np.zeros(aa_struct_prob.shape[2])
    for i_window,aa in enumerate(seq1):
        i_aa = aas_list.index(aa)
        sum_logprob += aa_struct_prob[i_aa,i_window,:]
    for i_window,aa in enumerate(seq2):
        i_aa = aas_list.index(aa)
        sum_logprob += aa_struct_prob[i_aa,i_window,:]
    return sum_logprob

def pyrmsd_table(args):
    cc_coords_level1,cc_coords_level2,dimer_overlap = args

    n_models = cc_coords_level1.shape[0]
    n_atoms_total = cc_coords_level1.shape[1]
    n_atoms_mono = int(n_atoms_total/2)
    n_atoms_per_res = 5
    n_atoms_overlap = dimer_overlap*n_atoms_per_res

    range1 = range(n_atoms_mono-n_atoms_overlap,n_atoms_mono)
    range2 = range(n_atoms_mono,n_atoms_mono+n_atoms_overlap)
    ref_atoms = cc_coords_level1[:,range1+range2]

    range1 = range(n_atoms_overlap)
    range2 = range(n_atoms_total-n_atoms_overlap,n_atoms_total)
    sup_atoms = cc_coords_level2[:,range1+range2]

    rms_matrix = np.zeros((n_models,n_models))

    for j in range(n_models):
        CC_cur = ref_atoms[j]
        aln_models = np.insert(sup_atoms, 0, CC_cur, 0)
        calculator = RMSDCalculator.RMSDCalculator("QCP_OMP_CALCULATOR", aln_models)
        dist = calculator.oneVsTheOthers(0)
        rms_matrix[j,:]=dist

    return np.square(rms_matrix)

def belief_propagation(variants,rms_matrix,n_select_bp):
    n_variables = len(variants)
    n_models = [len(v) for v in variants]
    forward_acc = [np.ones(n_i, dtype=np.float128) / n_i for n_i in n_models]
    backward_acc = [np.ones(n_i, dtype=np.float128) / n_i for n_i in n_models]
    forward_acc_log = [np.log(np.ones(n_i, dtype=np.float128) / n_i) for n_i in n_models]
    backward_acc_log = [np.log(np.ones(n_i, dtype=np.float128) / n_i) for n_i in n_models]
    forward_trans = []
    for level_rmsd in rms_matrix:
        level = np.copy(level_rmsd)+0.01
        sum_vec = np.sum(level, axis=1)
        level = level / sum_vec[:, None]
        forward_trans.append(level)

    backward_trans = []
    for level_rmsd in rms_matrix:
        level = np.copy(level_rmsd.T)+0.01
        sum_vec = np.sum(level, axis=1)
        level = level / sum_vec[:, None]
        backward_trans.append(level)

    for i in range(1, n_variables):
        for j in range(0, n_models[i]):
            forward_acc[i][j] = np.sum(np.multiply(forward_acc[i - 1], forward_trans[i - 1][:, j]))
            products = forward_acc_log[i - 1] + np.log(forward_trans[i - 1][:, j])
            max_log = np.max(products)
            forward_acc_log[i][j] = np.log(np.sum(np.exp(products-max_log)))+max_log
    for i in reversed(range(0, n_variables - 1)):
        for j in range(0, n_models[i]):
            backward_acc[i][j] = np.sum(np.multiply(backward_acc[i + 1], backward_trans[i][:, j]))
            products = backward_acc_log[i + 1] + np.log(backward_trans[i][:, j])
            max_log = np.max(products)
            backward_acc_log[i][j] = np.log(np.sum(np.exp(products - max_log))) + max_log

    variants_bp = []
    ind_rearrs = []
    for i in range(0, n_variables):
        solutions = []
        for j in range(0, n_models[i]):
            total_lg_ij = forward_acc_log[i][j] + backward_acc_log[i][j]
            solutions.append(total_lg_ij)

        ind_rearr = np.argsort(solutions)[:n_select_bp]
        ind_rearrs.append(ind_rearr)
        var_level = np.array(variants[i])[ind_rearr].tolist()
        variants_bp.append(var_level)

    rms_dict = []
    for i,(row_idx, col_idx,var_from,var_to) in enumerate(zip(ind_rearrs[:-1],ind_rearrs[1:],variants_bp[:-1],variants_bp[1:])):
        rms_matrix_red = rms_matrix[i,row_idx[:, None], col_idx]
        rms_dict.append((var_from,var_to,rms_matrix_red))

    return variants_bp, rms_dict

def run():

    target_seq = 'ATTFARLCQQVDMTQKHLEEEIARLSKEIDQLEKMQNNSKLLRNKAVQLESELENFSKQFLHAAAAA'
    dimer_overlap = [4, 5, 4, 6, 4]
    start_inds = [0, 11, 21, 32, 41, 52]
    pdb_id = '1t3j'

    # target_seq = 'ALKKHHENEISHHAKEIERLQKEIERHKQSIKKLKQSEDDDAA'
    # dimer_overlap = [5, 7, 5]
    # start_inds = [0, 10, 18, 28]
    # pdb_id = '1hf9'

    # target_seq = 'CGGDNIEQKIDDIDHEIADLQAKRTRLVQQHPRAAA'
    # dimer_overlap = [4, 5]
    # start_inds = [0, 11, 21]
    # pdb_id = '1r48'

    out_path = pdb_id

    if not os.path.exists(out_path):
        os.makedirs(out_path)

    aas_sorted = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    n_window = 15

    data_dir = os.path.dirname(os.path.realpath(__file__))
    loaded = np.load(os.path.join(data_dir,'cc_data_anti.npz'))
    cc_coords_all = loaded['cc_coords_all']
    aa_struct_prob = loaded['aa_struct_prob']

    variants = []

    target_seq_reversed = target_seq[::-1]
    n_select = 350

    cc_coords_levels = []
    for start_ind in start_inds:
        seq_window = target_seq[start_ind:start_ind+n_window]
        seq_window_rev = target_seq_reversed[start_ind:start_ind+n_window]

        dists = est_logprob_cached(seq_window,seq_window_rev,aas_sorted,aa_struct_prob)
        sorted_inds = np.argsort(dists)[::-1].tolist()
        var_level = sorted_inds[:n_select]
        cc_coords_levels.append(cc_coords_all[var_level, :, :])
        print seq_window
        variants.append(var_level)

    print 'logprob done'

    n_proc = min(max(1, multiprocessing.cpu_count()), len(variants)-1)
    pool = multiprocessing.Pool(n_proc)

    from_to = zip(cc_coords_levels[:-1], cc_coords_levels[1:], dimer_overlap)

    rms_matrix = pool.map(pyrmsd_table, from_to, 1)
    pool.close()
    pool.join()

    rms_matrix = np.array(rms_matrix)
    print 'rms matrix done'

    variants_bp,rms_dict = belief_propagation(variants,rms_matrix,50)

    print 'belief propagation done'

    cc_inds_logprob = [v[0] for v in variants_bp]
    most_probable_coords,rmsd = get_all_coords(cc_coords_all, cc_inds_logprob, dimer_overlap, True)
    write_model(most_probable_coords, target_seq, os.path.join(out_path, pdb_id+'_bp_{0:0.3f}'.format(rmsd)))

if __name__ == "__main__":
    run()
