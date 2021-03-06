import numpy as np
import itertools
from collections import namedtuple

import WR.Operators.Chebyshev as Chebyshev
import WR.Operators.operator_from_matrix as ofm

import time

import os.path


__author__ = ["Benjamin, Bykowski", "Jean-Luc Bouchot"]
__copyright__ = "Copyright 2019, Chair C for Mathematics (Analysis), RWTH Aachen and Seminar for Applied Mathematics, ETH Zurich and School of Mathematics and Statistics, Beijing Institute of Technology"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Falk Pulsmeyer", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "jlbouchot@gmail.com"
__status__ = "Development"
__lastmodified__ = "2019/06/05"

from time import sleep

CSPDEResult = namedtuple('CSPDEResult', ['J_s', 'N', 's', 'm', 'd', 'Z', 'y', 'A', 'w', 'result', 't_samples', 't_matrix', 't_recovery'])

# def CSPDE_ML(spde_model, wr_model, unscaledNbIter, epsilon, L_first = 1, L=1, dat_constant = 5, ansatz_space = 0, cspde_result = None, sampling_fname = None, datamtx_fname = None, t = 1, tprime = 1, p0 = 1./3., p = 2./3., energy_constant = 10): # the filenames are early only if we already computed quite a few solutions and don't want to have to recompute the whole matrix. In theory, we wouldn't need this for all practical purposes.
def CSPDE_ML(spde_model, wr_model, dict_config, cspde_result = None, sampling_fname = None, datamtx_fname = None): 
    """
    Parameters
    ----------
    spde_model : SPDEModel 
        an SPDEModel containing all the FEniCS models required for the computation
    wr_model : WRModel 
        an algorithm used for the recovery and containing the sampling matrix 
    unscaledNbIter : int / float
        an (absolute) number of iterations for the iterative recovrery algo 
    epsilon : float 
        a tolerance on the residual for the recovery algorithms 
    L : int 
        number of levels being studied -- Default = 1
    dat_constant : int 
        a very poorly chosen name to describe the proportionality constant in the number of samples -- Default = 10 
    ansatz_space : int, >= 0 
        Describes how the Ansatz space of polynomial is chosen. 0 is the one describe in the papers, while anything > 0 defines the total degree Ansatz space -- Default = 0 
    cspde_result = None --- This is reminiscent from an old version and should be deleted /!\ /!\ /!\
    sampling_fname : string 
        a string containing the filename to which the samples will be saved -- Default = None
    datamtx_fname : string 
        a file containing the precomputed sensing matrix -- Default = None
    """
    # First set up some basic things needed for the rest of the computations
    lvl_by_lvl_result = [] # This will keep the results

    # Load all the important things from the dictionary 
    unscaledNbIter = dict_config["iter"] 
    epsilon = dict_config["tolres"]
    L_first = dict_config["J"] 
    L = dict_config["L"] 
    dat_constant = dict_config["s_L"] 
    p = dict_config["p"]
    p0 = dict_config["p0"] 
    t = dict_config["t"]
    tprime = dict_config["tprime"]
    energy_constant = dict_config["s_J"]
    ansatz_space = dict_config["ansatz"] 

    s_L = np.ceil((dat_constant*(L-L_first))**(p/(1-p))) # This is basically the multiplicative constant in front of the sparsity at the finest level


    # Approximate the Jth level with a single level CSPG 
    energy_constant = np.max((energy_constant, dat_constant**(p/p0*(1-p0)/(1-p)) * (L-L_first)**(p/p0*(1-p0)/(1-p)) * 2**((L-L_first-1)*p/p0*(1-p0)/(1-p)*(t+tprime)) * 2**(-L*(t+tprime))+1)) # This ensures that the Jth level has more samples than the J+1
    s_J = np.ceil(energy_constant**(p0/(1-p0))*2**(L*p0*(t+tprime)/(1-p0)))
    print("Computing level {0} (this is a Single Level approximation) from a total of {1}. Current sparsity = {2}".format(L_first,L,s_J))
    ## 1. Create index set and draw random samples
    print("Generating J_s ...")
    
    # Compute "active index set" J_s
    if ansatz_space == 0: 
        J_s = J(s_J, wr_model.operator.theta, wr_model.weights)
    else:
        J_s = J_tot_degree(wr_model.weights, ansatz_space)

    # Get total number of coefficients in tensorized chebyshev polynomial base
    N = len(J_s)

    # Calculate number of samples
    m = wr_model.get_m_from_s_N(s_J, N)
        
    # Get sample dimension
    d = len(J_s[0])

    # Check whether this even an interesting case
    print("   It is N={0}, m={1} and d={2} ... ".format(N, m, d))
    wr_model.check(N, m)

    y_new, y_old, Z, t_samples = get_samples(spde_model, wr_model, m, d, L_first, L_first, L, s_J, sampling_fname, datamtx_fname)
    A, t_matrix = get_mtx(wr_model, J_s, Z, d, L_first, L_first, L, s_J, sampling_fname, datamtx_fname)

    print("   Computing weights ...")
    w = calculate_weights(wr_model.operator.theta, np.array(wr_model.weights), J_s)

    print("   Weighted minimization ...")
    t_start = time.time()
    # print("Shape of y_new is {} versus the shape of y_old {}".format(y_new.shape,y_old.shape))
    result = wr_model.method(A, y_new-y_old, w, s_J, epsilon, unscaledNbIter) # note that if we decide to not have a general framework, but only a single recovery algo, we can deal with a much better scaling: i.e. 13s for omp, 3s for HTP, and so on...
    t_stop = time.time()
    t_recovery = t_stop-t_start
    # result = wr_model.method(A, y_new-y_old, w, sl, np.sqrt(m) *epsilon, unscaledNbIter) # note that if we decide to not have a general framework, but only a single recovery algo, we can deal with a much better scaling: i.e. 13s for omp, 3s for HTP, and so on...
    lvl_by_lvl_result.append(CSPDEResult(J_s, N, s_J, m, d, Z, y_new-y_old, 0, w, result, t_samples, t_matrix, t_recovery))
    print("\n\tRecovery time: {0} \t Building the Matrix: {1} \t Computing the samples {2}\n".format(t_recovery, t_matrix, t_samples))



    # Deal with the approximation of the details
    for oneLvl in range(L_first+1,L+1):
        # sl = 10+np.max([2**(L-oneLvl),1])

        # sl = np.floor(dat_constant*2**(L-oneLvl))
        sl = np.ceil(s_L*2**((L-oneLvl)*(t+tprime)*p/(1-p)))
        print("Computing level {0} from a total of {1}. Current sparsity = {2}".format(oneLvl,L,sl))
        ## 1. Create index set and draw random samples
        print("Generating J_s ...")
        
        # Compute "active index set" J_s
        if ansatz_space == 0: 
            J_s = J(sl, wr_model.operator.theta, wr_model.weights)
        else:
            J_s = J_tot_degree(wr_model.weights, ansatz_space)

        # Get total number of coefficients in tensorized chebyshev polynomial base
        N = len(J_s)

        # Calculate number of samples
        m = wr_model.get_m_from_s_N(sl, N)
        
        # Get sample dimension
        d = len(J_s[0])

        # if not cspde_result is None:
        #     assert d == cspde_result.d, "New sample space dimension is different from old sample space dimension."

        # Check whether this even an interesting case
        print("   It is N={0}, m={1} and d={2} ... ".format(N, m, d))
        wr_model.check(N, m)


        y_new, y_old, Z, t_samples = get_samples(spde_model, wr_model, m, d, oneLvl, J, L, sl, sampling_fname, datamtx_fname)
        A, t_matrix = get_mtx(wr_model, J_s, Z, d, oneLvl, J, L, sl, sampling_fname, datamtx_fname)

        
        print("   Computing weights ...")
        w = calculate_weights(wr_model.operator.theta, np.array(wr_model.weights), J_s)    

        print("   Weighted minimization ...")
        t_start = time.time()
        result = wr_model.method(A, y_new-y_old, w, sl, epsilon, unscaledNbIter) # note that if we decide to not have a general framework, but only a single recovery algo, we can deal with a much better scaling: i.e. 13s for omp, 3s for HTP, and so on...
        t_stop = time.time()
        t_recovery = t_stop-t_start
        # result = wr_model.method(A, y_new-y_old, w, sl, np.sqrt(m) *epsilon, unscaledNbIter) # note that if we decide to not have a general framework, but only a single recovery algo, we can deal with a much better scaling: i.e. 13s for omp, 3s for HTP, and so on...
        lvl_by_lvl_result.append(CSPDEResult(J_s, N, sl, m, d, Z, y_new-y_old, 0, w, result, t_samples, t_matrix, t_recovery))
        print("\n\tRecovery time: {0} \t Building the Matrix: {1} \t Computing the samples {2}\n".format(t_recovery, t_matrix, t_samples))
    
    
    return lvl_by_lvl_result


def get_mtx(wr_model, J_s, Z, d, oneLvl, J, L, sl, sampling_fname, datamtx_fname): 
    if (sampling_fname is None) or (datamtx_fname is None):  
        # Create sampling matrix and weights
        print("   Creating sample operator ...")
        t_start = time.time()
        A = wr_model.operator.create(J_s, Z)
        t_stop = time.time()
        t_matrix = t_stop-t_start
    else:
        # Create sampling matrix and weights
        print("   Creating sample operator ...")
        t_start = time.time()
        mtx_file = datamtx_fname + '_d' + str(d) + '_l' + str(oneLvl) + '_s_' + str(sl) + '.npy'
        if os.path.isfile(mtx_file):
            # A = ofm(Chebyshev, np.load(mtx_file))
            A = wr_model.operator.load(mtx_file)
        else: 
            A = wr_model.operator.create(J_s, Z)
            A.save(mtx_file)
            # np.save(mtx_file, A.A)
        t_stop = time.time()
        t_matrix = t_stop-t_start

    return A, t_matrix



def get_samples(spde_model, wr_model, m, d, oneLvl, J, L, sl, sampling_fname, datamtx_fname):
    if (sampling_fname is None) or (datamtx_fname is None):  
        Z = wr_model.operator.apply_precondition_measure(np.random.uniform(-1, 1, (m, d)))
        print("\nComputing {0} SPDE sample approximations ...".format(m))
        # Get samples
        t_start = time.time()
        if onelvl != J:
            y_old = spde_model.samples(Z)
            # print("Here are the samples computed: {}".format(y_old))
            # print("Here are the size of the samples: {}".format(y_old.len))
            spde_model.refine_mesh()
        else:
            y_old = np.zeros([m,441])
        y_new = spde_model.samples(Z)
        t_stop = time.time()
        t_samples = t_stop-t_start


    else:
        sampling_file = sampling_fname + '_d' + str(d) + '_l' + str(oneLvl) + '_s_' + str(sl) + '.npy'
        # sampling_file = sampling_fname + '_d' + str(d) + '_l' + str(oneLvl) + '.npy' ## Really HAVE to do this better one day!
        if os.path.isfile(sampling_file):
            Z = np.load(sampling_file)
        else: 
            Z = wr_model.operator.apply_precondition_measure(np.random.uniform(-1, 1, (m, d)))
            np.save(sampling_file, Z)
        print("\nComputing {0} SPDE sample approximations ...".format(m))
        # Get samples
        t_start = time.time()
        if oneLvl != J:
            y_old = spde_model.samples(Z)
            spde_model.refine_mesh()
        else:
            y_old = np.zeros([m,1681])
        y_new = spde_model.samples(Z)
        t_stop = time.time()
        t_samples = t_stop-t_start


    return y_new, y_old, Z, t_samples






def J_tot_degree(v, max_degree = 2, threshold = np.inf):
    print("Generating an Ansatz space of multiindices what have total degree <= {}".format(max_degree))
    # Remember v contains the weights associated to the operators in the expansion. 
    # We assume that above a certain weight, it can simply be discarded, the associated coefficient can be discarded. 
    aux = np.array([list(x) for x in itertools.product(range(max_degree+1), repeat=len([v_i for v_i in v if v_i < threshold]))]) # This creates a set of multi-indices with max norm max_degree
    J = [one_multi_index for one_multi_index in aux if one_multi_index.sum() <= max_degree]
    return J


def J(s, theta, v):
    print("s = {0}, theta = {1}, v = {2}".format(s,theta,v))
    # Function for generating all admissible indices over given index set S
    def iterate(M, a, B, S):
        def iterate_(B, S, p):
            L = []

            if len(S):
                while True:
                    # Take the highest index in S
                    r     = S[-1]

                    # Substract weight from B
                    B    -= a[r]

                    # Increase nu_r by one
                    p[r] += 1

                    
                    if len(S) > 1:
                        # If there is more than one index left, recurse with remaining indices and 'remaining weight' B
                        e = iterate_(B, S[0:-1], p.copy())
                    else:
                        # If B - sum over a_j is larger than 0, add this multiindex
                        if B >= 0:
                            e = [p.copy()]
                        else:
                            e = []

                    # If the list of new multiindices is empty, there is nothing left to be done
                    # thanks to the monotonicity of a
                    if not e:
                        break

                    # Add found multiindices
                    L += e

            return L

        return iterate_(B, S, np.zeros(M, dtype='int'))


    # Set A and a as in Theorem 5.2
    A = np.log2(s/2.)
    a = 2 * np.log2(v)
    T = 2 * np.log2(theta)

    # Determine maximal M s.t. for j = 0 ... M-1 is a_j <= A - T
    # M is also the maximal support size
    # M = np.argmin(a <= A - T)
    M = np.argmin(a <= A )
    assert 0 != M, "Weight array too short. (Last element: {0}. Threshold: {1})".format(a[-1], A-T)

    # If A is non-negative the zero vector is always admissible
    assert A >= 0, "Negative A, i.e. sparsity less than 2."
    L = [np.zeros(M, dtype='int')]

    # Iterate through support sets of cardinality k = 1 ... M
    for k in range(1, M + 1):
        new_indices = []
    
        for S in itertools.combinations(range(M), k):
            new_indices += iterate(M, a, A - k*T, list(S))

        if [] == new_indices:
            break

        L += new_indices
    
    return L


def calculate_weights(theta, v, J_s):
    return np.array([theta**np.count_nonzero(nu) * np.product(v[np.where(nu > 0)]**nu[np.where(nu > 0)]) for nu in J_s])
    # return np.array([theta**np.count_nonzero(nu) * np.product(v[nu > 0]**nu[nu > 0]) for nu in J_s]) old version for earlier python distributions
