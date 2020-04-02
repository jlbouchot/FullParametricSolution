import numpy as np

from WR import eps
from .Result import *
from .weighted_quasi_abslargest import *

__author__ = ["Benjamin, Bykowski", "Jean-Luc Bouchot"]
__copyright__ = "Copyright 2015-2020, Chair C for Mathematics (Analysis), RWTH Aachen; Seminar for Applied Mathematics, ETH Zurich; School of Mathematics and Statistics, Beijing Institute of Technology"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPLv3"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "jlbmathit@gmail.com"
__status__ = "Development"
__lastmodified__ = "2020/04/02"

def wiht(Operator, y, w, s, eta, maxiter):
    x         = np.zeros((Operator.n, y.shape[1]))
    last_norm = 0
    k         = 0
    print(y)

    # print("Size of x is {} to be compared to the size of Ax {}".format(x.shape,Operator.apply(x).shape))

    while np.linalg.norm(Operator.apply(x) - y) > eta*y.shape[1]: # XXX Check if this is the appropriate norm
        print("Iteration number {} -- Current norm on the residual {}".format(k,np.linalg.norm(Operator.apply(x) - y)))
        residuum = y - Operator.apply(x)
        cur_norm = np.linalg.norm(residuum) # See above

        # print("A times x is {}".format(Operator.apply(x).tolist()))

        # print("A star residual is {}".format(Operator.apply_adj(residuum).tolist()))

        x, dummy, row_norms  = weighted_quasi_abslargest(x + Operator.apply_adj(residuum), 3 * s, w)
        last_norm = cur_norm
        k         = k + 1
        # print(x)
        # print(np.sort(dummy))
        print(row_norms)

        if k > maxiter:
            print('WIHT did not converge after {0} iterations.'.format(k))
            break

    return Result(x, k, 'Weighted Iterative Hard Thresholding')



