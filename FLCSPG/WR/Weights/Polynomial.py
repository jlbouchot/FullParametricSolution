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

class Polynomial:
    def __init__(self, c, alpha, stepsize=None):
        if stepsize is None:
            stepsize = 1

        if c == 1:
            self.name = 'Polynomial $j^{{{0}}}$'.format(alpha)
        else:
            self.name = 'Polynomial ${0} \cdot j^{{{1}}}$'.format(c, alpha)
            
        self.w    = c * np.repeat(np.arange(1,100000), stepsize)**alpha
