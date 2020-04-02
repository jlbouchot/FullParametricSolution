import numpy as np
from progressbar import Bar, ETA, Percentage, ProgressBar

__author__ = ["Benjamin, Bykowski", "Jean-Luc Bouchot"]
__copyright__ = "Copyright 2015-2020, Chair C for Mathematics (Analysis), RWTH Aachen; Seminar for Applied Mathematics, ETH Zurich; School of Mathematics and Statistics, Beijing Institute of Technology"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPLv3"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "jlbmathit@gmail.com"
__status__ = "Development"
__lastmodified__ = "2020/04/02"

class SPDEModel:
    def samples(self, Z, ratio=None):

        if ratio is not None:
            self.refine_mesh(ratio)

        m, _ = Z.shape

        # Show a progressbar
        widgets = [Percentage(), ' ', Bar(), ' ', ETA()]
        pbar    = ProgressBar(widgets=widgets)

        # Compute all functional evalutations with given precision
        y0 = self.sample(Z[0])
        nb_basis_fct = y0.size # Dunno how to do this quite yet. 
        y = np.zeros([m,nb_basis_fct])
        y[0,:] = y0
        for k in pbar(range(m-1)):
            y[k+1,:] = self.sample(Z[k+1]) # Make sure this is the correct format

        return y
