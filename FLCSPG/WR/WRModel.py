import numpy as np
import WR.Algorithms as Algorithms

__author__ = ["Benjamin, Bykowski", "Jean-Luc Bouchot"]
__copyright__ = "Copyright 2015-2020, Chair C for Mathematics (Analysis), RWTH Aachen; Seminar for Applied Mathematics, ETH Zurich; School of Mathematics and Statistics, Beijing Institute of Technology"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPLv3"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "jlbmathit@gmail.com"
__status__ = "Development"
__lastmodified__ = "2020/04/02"

class WRModel:
    def __init__(self, method, operator, weights, m_from_s_N, check):
        self.method         = get_recovery_algo_from_string(method)
        self.operator       = operator
        self.weights        = weights
        self.get_m_from_s_N = m_from_s_N
        self.check          = check

    def estimate_ML_samples(self, cspde_result, Z_cross):
        nbLevel = len(cspde_result)
        nb_nodes = cspde_result[0].result.x.shape[1]
        # cspde_result[l-1].d contains the number of dimensions d required at level l. 
        # The Z_cross contains max_{for 0 \leq l \leq L-1} {cspde_result[l-1].d} dimensions!
        y_recon = np.zeros([len(Z_cross),nb_nodes])

        for oneLevel in range(0,nbLevel):
            # F_recon = self.operator.create(np.array(cspde_result[oneLevel].J_s)[cspde_result[oneLevel].result.x != 0], Z_cross[:,:cspde_result[oneLevel].d])

            # Keep track of the interesting rows
            row_norms = np.linalg.norm(cspde_result[oneLevel].result.x, axis=1)

            if np.sum(row_norms  > 1e-10) > 0 : # Avoids running into problems in case the resulting sparse vector is the 0 vector!
#                F_recon = self.operator.create(np.array(cspde_result[oneLevel].J_s)[cspde_result[oneLevel].result.x != 0], Z_cross[:,:cspde_result[oneLevel].d], normalization=np.sqrt(cspde_result[oneLevel].m))
                F_recon = self.operator.create(np.array(cspde_result[oneLevel].J_s)[row_norms > 1e10], Z_cross[:,:cspde_result[oneLevel].d], normalization=np.sqrt(cspde_result[oneLevel].m))
                y_recon = y_recon+F_recon.apply(cspde_result[oneLevel].result.x[row_norms > 1e10, :])
            #print cspde_result[oneLevel].result.x[cspde_result[oneLevel].result.x != 0]

        return y_recon


def check_cs(N, m):
    if N < m: 
        print("!!! WARNING !!! It holds N < m ({0} < {1}). You might want to use other L2 minimization procedures instead of weighted l1".format(N, m))
    # assert N >= m, "It is N < m ({0} < {1}). Aborting ...".format(N, m)

def cs_theoretic_m(s, N):
    # TODO: Constant factor is known to be 2?
    return int(np.ceil(2 * np.log(s)**3 * s * np.log(N)))

def cs_pragmatic_m(s, N):
    # Calculate number of samples using the "pragmatic" bound. That is:
    #  Factor np.log(s)**3 probably ignoreable in practice
    #  Constant factor "C" is about 2
    return int(np.ceil(2 * s * np.log(N)))
	

def cs_theoretic_m_new(s, N):
    # TODO: Constant factor is known to be 2?
    return int(np.ceil(2 * np.log(s)**2 * s * np.log(N)))
	
def get_recovery_algo_from_string(algo_name):
    switcher = {
        "whtp": Algorithms.whtp,
        "wiht": Algorithms.wiht,
        "womp": Algorithms.womp,
#        "bp": Algorithms.exact_wbp_cvx,
#        "bpdn": Algorithms.qc_wbp_cvx,
    }
    return switcher.get(algo_name, Algorithms.whtp)
