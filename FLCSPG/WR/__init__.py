__author__ = ["Benjamin, Bykowski", "Jean-Luc Bouchot"]
__copyright__ = "Copyright 2015-2020, Chair C for Mathematics (Analysis), RWTH Aachen; Seminar for Applied Mathematics, ETH Zurich; School of Mathematics and Statistics, Beijing Institute of Technology"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPLv3"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "jlbmathit@gmail.com"
__status__ = "Development"
__lastmodified__ = "2020/04/05"

__all__ = [
    'Algorithms',
    'Operators',
    'Weights',
    'WRModel'
]

import numpy as np
eps = np.finfo(float).eps

from .Algorithms import *
from .Operators import *
from .Weights import *
from .WRModel import *
