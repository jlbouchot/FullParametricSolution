import numpy as np

from WR import eps
from .Result import Result
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

def whtp(Operator, y, w, s, eta, maxiter):
    x     = np.zeros_like(Operator.apply_adj(y))
    S_old = np.array([])
    k     = 0

    while True:
        if k > maxiter:
            print('WHTP did not converge after {0} iterations.'.format(k))
            break

        dummy, S,row_norms = weighted_quasi_abslargest(x + Operator.apply_adj(y - Operator.apply(x)), 3 * s, w)
        print('Norm of the error {0}, while the objective norm is {1} and the ratio is {2}'.format(np.linalg.norm(Operator.apply(x) - y),np.linalg.norm(y),np.linalg.norm(Operator.apply(x) - y)/np.linalg.norm(y)))
        if set(S) == set(S_old) or np.linalg.norm(Operator.apply(x) - y)/np.linalg.norm(y) <= eta:
            print('WHTP Converged after {0} iterations. Norm of residual {1}'.format(k,np.linalg.norm(Operator.apply(x) - y)))
            break

        SubMatrix = Operator.genSubMatrix(S)

        x    = np.zeros_like(x)
        x[S] = np.linalg.lstsq(SubMatrix, y)[0]

        S_old = S
        k     = k + 1

    return Result(x, k, 'Weighted HTP')
