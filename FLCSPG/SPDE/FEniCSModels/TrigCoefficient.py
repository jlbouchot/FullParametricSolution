from dolfin import *

import numpy as np
from scipy.special import zeta

__author__ = ["Benjamin, Bykowski", "Jean-Luc Bouchot"]
__copyright__ = "Copyright 2015-2020, Chair C for Mathematics (Analysis), RWTH Aachen; Seminar for Applied Mathematics, ETH Zurich; School of Mathematics and Statistics, Beijing Institute of Technology"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPLv3"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "jlbmathit@gmail.com"
__status__ = "Development"
__lastmodified__ = "2020/04/02"

class TrigCoefficient:
    def __init__(self, d, alpha, c = None):
        self.num_params = 2*d

        if c is None:
            self.c = zeta(alpha, 1) * 2 + 0.1
        else:
            self.c = c

        self.alpha = alpha

    def __call__(self, x, z):
        r = 0
        for k in range(int(len(z)/2)):
            r +=     z[2*k]   * cos(np.pi * (k+1) * x[0]) / (k+1)**self.alpha

            if 2*k+1 < len(z):
                r += z[2*k+1] * sin(np.pi * (k+1) * x[0]) / (k+1)**self.alpha

        return self.c + r
