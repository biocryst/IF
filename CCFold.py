import os, sys, argparse, itertools, multiprocessing
import numpy as np
from copy import deepcopy
from collections import defaultdict
from Bio import PDB, SeqIO
from Bio.SVDSuperimposer import SVDSuperimposer
from Bio.SeqUtils import seq1,seq3
from pyRMSD import RMSDCalculator


def new_cc(sequences, coords):
    n_cc_helices = len(sequences)
    seq_len = int(min(len(sequences[0])*5,coords.shape[0]/n_cc_helices)/5)
    new_model = PDB.Model.Model(0)
    segid = '    '
    atomname = ['N', 'CA', 'C', 'O', 'CB']
    bfactor = 30
    occupancy = 1
    altloc = ' '
    fullname = [' N  ', ' CA ',' C  ', ' O  ', ' CB ']
    element = ['N', 'C', 'C', 'O', 'C']
    serial_number = 1
    chain_id_base = ord('A')
    for chain_i, sequence in enumerate(sequences):
        chain_id = chr(chain_id_base+chain_i)
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


def merge_cc(coords_list, res_overlap,n_cc_helices):

    ref_coords = coords_list[0]
    aligned_coords = [deepcopy(coords_list[0])]
    n_atoms_per_res = 5
    n_atoms_mono = int(ref_coords.shape[0]/n_cc_helices)
    msds = []
    for coords, cc_overlap in zip(coords_list[1:],res_overlap):

        n_atoms_overlap = cc_overlap*n_atoms_per_res

        for i in range(n_cc_helices):
            hi_ref = ref_coords[(i+1)*n_atoms_mono-n_atoms_overlap:(i+1)*n_atoms_mono]
            if i==0:
                ref_atoms = hi_ref
            else:
                ref_atoms = np.append(ref_atoms, hi_ref, axis=0)

        for i in range(n_cc_helices):
            hi = coords[i*n_atoms_mono:i*n_atoms_mono+n_atoms_overlap]
            if i==0:
                sup_atoms = hi
            else:
                sup_atoms = np.append(sup_atoms, hi, axis=0)

        sup=SVDSuperimposer()
        sup.set(ref_atoms,sup_atoms)
        sup.run()
        msds.append(sup.get_rms()**2)
        rot,tran = sup.get_rotran()
        coord_new = np.dot(coords, rot) + tran
        aligned_coords.append(coord_new)
        ref_coords = coord_new

    rmsd = np.sqrt(np.sum(msds))

    hi_all = []
    for i in range(n_cc_helices):
        hi_all.append(aligned_coords[0][i*n_atoms_mono:(i+1)*n_atoms_mono])

    for coords,cc_overlap in zip(aligned_coords[1:],res_overlap):
        hi = []
        for i in range(n_cc_helices):
            hi.append(coords[i * n_atoms_mono:(i + 1) * n_atoms_mono])

        n_atoms_overlap = cc_overlap*n_atoms_per_res
        for ind_overlap in range(cc_overlap):
            weight = (ind_overlap+1)/float(cc_overlap+1)
            for ind_atom in range(n_atoms_per_res):
                ind_shift = ind_overlap*n_atoms_per_res+ind_atom

                for i in range(n_cc_helices):
                    coordi_prev = hi_all[i][-n_atoms_overlap+ind_shift]
                    coordi_next = hi[i][ind_shift]
                    hi_all[i][-n_atoms_overlap+ind_shift] = (1-weight)*coordi_prev+weight*coordi_next

        for i in range(n_cc_helices):
            hi_rest = hi[i][n_atoms_overlap:]
            hi_all[i] = np.append(hi_all[i],hi_rest,axis=0)

    res_dimer =hi_all[0]
    for i in range(1,n_cc_helices):
        res_dimer = np.append(res_dimer,hi_all[i],axis=0)

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

def get_all_coords(cc_coords_all,cc_inds,dimer_overlap, n_cc_helices, return_rms = False):
    cc_coords_sel = []
    n_classes = len(cc_inds)

    for i in range(n_classes):
        cc_coords_sel.append(cc_coords_all[cc_inds[i]])

    coords, rms = merge_cc(cc_coords_sel,dimer_overlap, n_cc_helices)

    if return_rms:
        return coords, rms
    return coords

def write_model(coords,target_seq,filename):
    cc_model = new_cc(target_seq,coords)

    io=PDB.PDBIO()
    io.set_structure(cc_model)
    io.save(filename+".pdb")
    return

def est_logprob_cached(seq_windows, aas_list,aa_struct_prob):
    sum_logprob = np.zeros(aa_struct_prob.shape[2])
    for seq_window in seq_windows:
        for i_window,aa in enumerate(seq_window):
            i_aa = aas_list.index(aa)
            sum_logprob += aa_struct_prob[i_aa,i_window,:]
    return sum_logprob

def pyrmsd_table(args):
    cc_coords_level1,cc_coords_level2,dimer_overlap, n_cc_helices = args

    n_models = cc_coords_level1.shape[0]
    n_atoms_total = cc_coords_level1.shape[1]
    n_atoms_mono = int(n_atoms_total/n_cc_helices)
    n_atoms_per_res = 5
    n_atoms_overlap = dimer_overlap*n_atoms_per_res

    range_ref = []
    for i in range(n_cc_helices):
        range_i = range((i + 1) * n_atoms_mono - n_atoms_overlap, (i + 1) * n_atoms_mono)
        range_ref += range_i

    ref_atoms = cc_coords_level1[:,range_ref]

    range_sup = []
    for i in range(n_cc_helices):
        range_i = range(i * n_atoms_mono, i * n_atoms_mono + n_atoms_overlap)
        range_sup += range_i

    sup_atoms = cc_coords_level2[:,range_sup]

    rms_matrix = np.zeros((n_models,n_models))

    for j in range(n_models):
        CC_cur = ref_atoms[j]
        aln_models = np.insert(sup_atoms, 0, CC_cur, 0)
        calculator = RMSDCalculator.RMSDCalculator("QCP_SERIAL_CALCULATOR", aln_models)
        dist = calculator.oneVsFollowing(0)
        rms_matrix[j,:]=dist

    return np.square(rms_matrix)


def sel_straight(coords_arr, n_cc_helices):
    n_atoms_mono = int(coords_arr[0].shape[0] / n_cc_helices)
    chain_rmss = []
    for coords in coords_arr:

        hi_all = []
        for i in range(n_cc_helices):
            hi_all.append(coords[i * n_atoms_mono:(i + 1) * n_atoms_mono])

        rmss = []
        for i in range(n_cc_helices-1):
            sup=SVDSuperimposer()
            sup.set(hi_all[i],hi_all[i+1])
            sup.run()
            rms = sup.get_rms()
            rmss.append(rms)
        chain_rmss.append(np.mean(rmss))

    return np.argmin(chain_rmss), np.min(chain_rmss)


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


def straighten(cc_coords_all, dimer_overlap, rms_dict, variants,n_cc_helices):
    variants_sym = []
    change_sym = range(len(variants))
    for change_sym_ind in change_sym:
        for model_ind in range(len(variants[change_sym_ind])):
            model_var = variants[change_sym_ind][model_ind]
            new_var = deepcopy(variants[:change_sym_ind] + [[model_var]] + variants[change_sym_ind + 1:])
            variants_sym.append(new_var)

    n_proc = min(len(variants_sym), multiprocessing.cpu_count())
    pool = multiprocessing.Pool(n_proc)
    params = itertools.product(variants_sym, [rms_dict])
    results = pool.map(best_model, params, 1)
    pool.close()
    pool.join()

    median_score = np.median([r[1] for r in results])
    results_filt = [r for r in results if r[1] < median_score]
    cand_coords_arr = []
    for result in results_filt:
        cc_inds, score = result
        cand_coords = get_all_coords(cc_coords_all, cc_inds, dimer_overlap,n_cc_helices)
        cand_coords_arr.append(cand_coords)

    straight_coords_ind, straight_rms = sel_straight(cand_coords_arr, n_cc_helices)
    cc_inds, score = results_filt[straight_coords_ind]
    coords = cand_coords_arr[straight_coords_ind]
    return coords, score, straight_rms

class split_sequence():
    def __init__(self, target_seq, aas_sorted,n_select,aa_struct_prob,cc_coords_all):
        self.target_seq = target_seq
        self.aa_struct_prob = aa_struct_prob
        self.cc_coords_all = cc_coords_all
        self.n_window = 15
        self.n_select = n_select
        self.aas_sorted = aas_sorted
        self.n_cc_helices = len(target_seq)
        self.n_res = len(target_seq[0])

        self.n_levels = self.n_res - self.n_window + 1
        self.cc_coords_levels = np.zeros((self.n_levels, self.n_select, self.cc_coords_all.shape[-2], self.cc_coords_all.shape[-1]))

        self.window_rms_table = np.zeros((self.n_levels, self.n_levels))
        self.solutions_rms_dic = defaultdict(list)
        self.solutions_overlap_dic = defaultdict(list)

    def overlaps2start(self,path):
        ind_windows = [0]
        cur_start = 0
        for overlap in path:
            cur_start = cur_start + 15 - overlap
            ind_windows.append(cur_start)
        return ind_windows

    def dyn_split_sequence(self,cur_start, cur_end, minmax_rms=0):
        if cur_start > cur_end:
            return [float("inf")], [0]
        if cur_start == cur_end:
            return [], []
        if self.solutions_rms_dic[cur_start]:
            return self.solutions_rms_dic[cur_start], self.solutions_overlap_dic[cur_start]

        rms_paths = []
        overlap_paths = []
        for overlap in range(4, 9):
            next_start = cur_start + 15 - overlap
            rms_path, overlap_path = self.dyn_split_sequence(next_start, cur_end, minmax_rms)
            rms_path = [self.get_window_rms(cur_start, next_start)] + rms_path
            overlap_path = [overlap] + overlap_path
            rms_paths.append(rms_path)
            overlap_paths.append(overlap_path)

        max_rmss = np.zeros((len(rms_paths)))
        sum_rmss = np.zeros((len(rms_paths)))
        for i, rms_path in enumerate(rms_paths):
            max_rmss[i] = np.max(rms_path)
            sum_rmss[i] = np.sum(rms_path)
        if minmax_rms > 0:
            cand_indices = np.where(max_rmss <= minmax_rms)
        else:
            cand_indices = np.where(np.isclose(max_rmss, max_rmss.min()))

        if len(cand_indices) == 1:
            cand_indices = cand_indices[0]
        if len(cand_indices) == 0:
            cand_indices = [0]

        ind_min = np.argmin(sum_rmss[cand_indices])
        solution_rms = rms_paths[cand_indices[ind_min]]
        solution_overlap = overlap_paths[cand_indices[ind_min]]
        self.solutions_rms_dic[cur_start] = solution_rms
        self.solutions_overlap_dic[cur_start] = solution_overlap

        return solution_rms, solution_overlap

    def get_window_rms(self,start1, start2):
        if start2 >= self.window_rms_table.shape[1]:
            return float('inf')
        if self.window_rms_table[start1, start2] > 0:
            return self.window_rms_table[start1, start2]
        res_overlap = start1 + 15 - start2

        n_models = self.cc_coords_levels.shape[1]
        n_atoms_total = self.cc_coords_levels.shape[2]
        n_atoms_mono = int(n_atoms_total / self.n_cc_helices)
        n_atoms_per_res = 5
        n_atoms_overlap = res_overlap * n_atoms_per_res

        range_ref = []
        for i in range(self.n_cc_helices):
            range_i = range((i+1)*n_atoms_mono - n_atoms_overlap, (i+1)*n_atoms_mono)
            range_ref += range_i

        ref_atoms = self.cc_coords_levels[start1][:, range_ref]


        range_sup = []
        for i in range(self.n_cc_helices):
            range_i = range(i*n_atoms_mono, i*n_atoms_mono+n_atoms_overlap)
            range_sup += range_i

        sup_atoms = self.cc_coords_levels[start2][:, range_sup]

        dists = []
        for j in range(n_models):
            CC_cur = ref_atoms[j]
            aln_models = np.insert(sup_atoms, 0, CC_cur, 0)
            calculator = RMSDCalculator.RMSDCalculator("QCP_OMP_CALCULATOR", aln_models)
            dist = calculator.oneVsFollowing(0)
            dists.append(dist)

        res_dist = np.mean(dists)
        self.window_rms_table[start1, start2] = res_dist
        return res_dist

    def run(self):
        for start_ind in range(0, self.n_levels, 1):
            seq_window = [t[start_ind:start_ind + self.n_window] for t in self.target_seq]
            dists = est_logprob_cached(seq_window, self.aas_sorted, self.aa_struct_prob)
            var_level = np.argsort(dists)[-self.n_select:][::-1]
            self.cc_coords_levels[start_ind, :, :, :] = self.cc_coords_all[var_level, :, :]

        opt_rms, opt_overlaps = self.dyn_split_sequence(0, self.n_levels - 1)

        self.solutions_rms_dic.clear()
        self.solutions_overlap_dic.clear()
        opt_rms, opt_overlaps = self.dyn_split_sequence(0, self.n_levels - 1, np.max(opt_rms))
        return opt_overlaps, self.overlaps2start(opt_overlaps),opt_rms

def run():

    n_window = 15
    aas_sorted = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

    parser = argparse.ArgumentParser(description='Coiled coils folding')
    parser.add_argument('-n', action='store', default=2, choices=['2','3'],
                    dest='n_cc_helices',metavar='2|3',help='Number of helices in a coiled coil. 2 (default) or 3 are supported at the moment.')

    parser.add_argument('-s', action='store_true', default=False,
                    dest='symmetric',help='Straighten the model. Prefer symmetric output at some cost to smoothness. '
                                          'Number of fragments tried for each segment is defined by the -spp parameter.')

    parser.add_argument('-t', action='store', default='one', choices=['one','termini','all'],
                    dest='vary',metavar='one|termini|all',help='Produce single model (default), '
                                                               'vary termini or vary fragments in all windows. '
                                                               'Number of fragments to choose from is defined by the -spp parameter. ')
    parser.add_argument('-sp', action='store', default=100,
                    dest='segmentation_pool',metavar='100',help='Number of fragments to consider for initial segmentation.')

    parser.add_argument('-bp', action='store', default=350,
                    dest='bp_pool',metavar='350',help='Number of fragments to pick for mutual alignment.')

    parser.add_argument('-spp', action='store', default=50,
                    dest='spp_pool',metavar='50',help='Number of fragments to use in shortest path outputs.')

    parser.add_argument('sequence', help='Fasta file with the target sequence (hmodimer), or a pair of aligned sequences (heterodimer).')
    parser.add_argument('output_dir', nargs='?',default='',help='Directory to output the model files.')
    args = parser.parse_args()

    seq_name = os.path.splitext(os.path.basename(args.sequence))[0]
    if args.output_dir == '':
        args.output_dir = seq_name

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    n_cc_helices = int(args.n_cc_helices)
    data_dir = os.path.dirname(os.path.realpath(__file__))
    loaded = np.load(os.path.join(data_dir,'cc_data'+str(n_cc_helices)+'.npz'))
    cc_coords_all = loaded['cc_coords_all']
    aa_struct_prob = loaded['aa_struct_prob']

    seq_file = list(SeqIO.parse(args.sequence, 'fasta'))
    target_seq = []
    homooligomer = True

    if len(seq_file) == 1:
        for i in range(n_cc_helices):
            target_seq.append(seq_file[0].seq)
    else:
        assert n_cc_helices == len(seq_file)
        for i in range(n_cc_helices):
            target_seq.append(seq_file[i].seq)
            assert len(seq_file[i].seq) == len(target_seq[0])
        homooligomer = False

    print 'Segmenting the sequence...'
    ss = split_sequence(target_seq, aas_sorted, args.segmentation_pool,aa_struct_prob, cc_coords_all)
    dimer_overlap, start_inds, opt_rms = ss.run()
    print 'Sequence windows:'
    variants = []
    cc_coords_levels = []
    for start_ind in start_inds:
        seq_window = [t[start_ind:start_ind+n_window] for t in target_seq]

        dists = est_logprob_cached(seq_window, aas_sorted,aa_struct_prob)
        sorted_inds = np.argsort(dists)[::-1].tolist()
        var_level = sorted_inds[:int(args.bp_pool)]
        cc_coords_levels.append(cc_coords_all[var_level, :, :])
        variants.append(var_level)
        if homooligomer:
            print seq_window[0]
        else:
            print [s for s in seq_window]

    n_proc = min(multiprocessing.cpu_count(), len(variants)-1)
    pool = multiprocessing.Pool(n_proc)

    from_to = zip(cc_coords_levels[:-1], cc_coords_levels[1:], dimer_overlap,itertools.repeat(n_cc_helices))

    print 'Calculating RMSD matrix...'

    rms_matrix = pool.map(pyrmsd_table, from_to, 1)
    pool.close()
    pool.join()

    rms_matrix = np.array(rms_matrix)

    print 'Belief propagation...'
    variants_bp,rms_dict = belief_propagation(variants,rms_matrix,int(args.spp_pool))

    print 'Fusing the model...'
    if args.vary == 'one':
        if args.symmetric:
            coords, score, straight_rms = straighten(cc_coords_all, dimer_overlap, rms_dict, variants_bp,n_cc_helices)
            fname = seq_name+"_{0:0.3f}_{1:0.3f}".format(score, straight_rms)
            write_model(coords, target_seq,os.path.join(args.output_dir, fname))
        else:
            cc_inds_logprob = [v[0] for v in variants_bp]
            most_probable_coords,rmsd = get_all_coords(cc_coords_all, cc_inds_logprob, dimer_overlap, n_cc_helices,True)
            write_model(most_probable_coords, target_seq,os.path.join(args.output_dir, seq_name+'_bp_{0:0.3f}'.format(rmsd)))
        sys.exit()

    if args.vary == 'termini':
        change_arr = [0,len(variants_bp)-1]
        prefixes = ['nterm','cterm']
    else:
        change_arr = range(len(variants_bp))
        prefixes = ['W'+str(c) for c in change_arr]

    for change_ind1,prefix in zip(change_arr,prefixes):
        print 'Varying sequence window', prefix
        variants_arr1 = []
        for model_ind in range(len(variants_bp[change_ind1])):
            model_var = variants_bp[change_ind1][model_ind]
            new_var = deepcopy(variants_bp[:change_ind1]+[[model_var]]+variants_bp[change_ind1+1:])
            variants_arr1.append(new_var)

        if args.symmetric:
            print 'Looking for straight models, this will take a while...'
            for ind_var1,var1 in enumerate(variants_arr1):
                coords, score, straight_rms = straighten(cc_coords_all, dimer_overlap, rms_dict, var1, n_cc_helices)
                fname = prefix+"_{0:0.3f}_{1:0.3f}_{2:04d}".format(score, straight_rms,ind_var1)
                write_model(coords, target_seq, os.path.join(args.output_dir, fname))

        else:
            n_proc = min(len(variants_arr1), n_proc)
            pool = multiprocessing.Pool(n_proc)
            params = itertools.product(variants_arr1, [rms_dict])
            results = pool.map(best_model, params, 1)
            pool.close()
            pool.join()

            for ind_var1,result in enumerate(results):
                cc_inds, score = result
                coords = get_all_coords(cc_coords_all, cc_inds, dimer_overlap,n_cc_helices)
                tmp_var, straight_rms = sel_straight([coords], n_cc_helices)
                fname = prefix + "_{0:0.3f}_{1:0.3f}_{2:04d}".format(score, straight_rms, ind_var1)
                write_model(coords,target_seq,os.path.join(args.output_dir,fname))


if __name__ == "__main__":
    run()
