import numpy as np

__author__ = ["Benjamin, Bykowski", "Jean-Luc Bouchot"]
__copyright__ = "Copyright 2015-2020, Chair C for Mathematics (Analysis), RWTH Aachen; Seminar for Applied Mathematics, ETH Zurich; School of Mathematics and Statistics, Beijing Institute of Technology"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPLv3"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "jlbmathit@gmail.com"
__status__ = "Development"
__lastmodified__ = "2020/04/02"

def weighted_quasi_abslargest(x, s, w):
    # Get indices of elements sorted in descending order
    row_norms = np.linalg.norm(x, axis=1) # have to check what type of norm this is 
    # print("Row norms is {} while the weights are {}".format(row_norms,w))
    # This would equally work for the case of a quantity of interest!
    sortIndex = np.argsort((row_norms) * (w**(-1)))[::-1]

    k      = 0
    w_test = 0
    w_cur  = 0

    while True:
        w_test = w_cur + w[sortIndex[k]]**2

        if s < w_test:
            break

        w_cur = w_test
        k     = k + 1


    i    = sortIndex[0:k]

    r    = np.zeros_like(x)
    r[i,:] = x[i,:]

    return r, i, row_norms

# def get_group_norms(x): 
#     # Simply compute the l2 norm for all the rows of the data matrix
    
