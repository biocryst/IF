import numpy as np
import cPickle

def est_logprob(seq,kde_dict,ca_dists):
    sum_logprob = 0
    for aa,ca_dist in zip(seq,ca_dists):
        sum_logprob += kde_dict[aa].score_samples([[ca_dist]])[0]
    return sum_logprob

def est_logprob_cached(seq,aas_list,aa_struct_prob):
    sum_logprob = np.zeros(aa_struct_prob.shape[2])
    for i_window,aa in enumerate(seq):
        i_aa = aas_list.index(aa)
        sum_logprob += aa_struct_prob[i_aa,i_window,:]
    return sum_logprob

def run():

    n_window = 15
    aa_distr_kde = cPickle.load(open('cc_aa_distance_distr.pkl', "rb"))
    cc_coords_all, cc_ca_dists_all = cPickle.load(open('cc_coords_unique_0.2.pkl', "rb"))

    n_structures = len(cc_coords_all)
    n_aas = len(aa_distr_kde.keys())

    aa_struct_prob = np.zeros((n_aas,n_window,n_structures))

    aas_sorted = sorted(aa_distr_kde.keys())

    for i_aa in range(n_aas):
        aa = aas_sorted[i_aa]
        print aa
        for i_window in range(n_window):
            for i_structure in range(n_structures):
                ca_dist = cc_ca_dists_all[i_structure][i_window]
                aa_struct_prob[i_aa,i_window,i_structure]=aa_distr_kde[aa].score_samples([[ca_dist]])[0]

    cPickle.dump((aa_distr_kde,aa_struct_prob), open('cc_aa_distance_distr_cached.pkl', "wb"))


if __name__ == "__main__":
    run()
