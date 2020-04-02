# FullParametricSolution
Compressed sensing approximation of the full solution to parametric PDEs
========================================================================

This repository is the prototype python implementation of the works related to compressed sensing approximation for parametric PDEs. 
Specifically, this repository deals with the recovery of the parametric solution and not only of a quantity of interest. 

## Remarks

-----------------------------------
This toolbox is provided as is so far. There is absolutely no guarantees that it will work on any of your problems!
Most of the parametric PDE computing still require the appropriate formulation to be passed to this tool. 
-----------------------------------

### Anaconda environment
You will need some of the following packages to run: 
* FEniCS: http://fenicsproject.org/ ; used as the black box solver for the PDEs
* ProgressBar: https://pypi.python.org/pypi/progressbar ; used to keep track of the advancement of the PDE solves
* Scientific Linux (including in particular scipy, numpy, etc... ) - Ubuntu works just fine too!
* cvxpy: http://www.cvxpy.org/en/latest/install/index.html ; used for the convex minimizations needed for the compressed sensing part

Therefore, I strongly recommend to work with a Docker image or a conda environment, associated with [this file](env-fenics.yml). 
This file is the one I am currently using for presenting numerical results.

Note that this package has been developped with Ubuntu 18.04 and Python 3.6. 
Some functions may not work with other versions of Python (some issues were encountered with Python 2, mainly when using printing functions)

Package Installation
====================

We recommend using anaconda for most of the installation. 
Be aware that these packages are evolving very fast and you might need to adjust here and there some versions of the various packages (most recently (2019/06/04): FEniCS and the latest version of cvxpy are not compatible). 
Try the following sequence of commands: 

First, create a separate environment for the project (note that you could use the [env-fenics.yml file](env-fenics.yml))
```
conda create --name env-fenics 
conda activate env-fenics
```

Make sure you have `FEniCS` installed
```
conda install -c conda-forge fenics
```

Add `cvxpy` with the right version
```
conda config --add channels oxfordcontrol
conda install -c conda-forge lapack
conda install -c cvxgrp cvxpy=1.0.11
```

(Note that the specified version is not the latest one!)
```
conda install nose
nosetests cvxpy
```

`numba` can be useful for certain speed-ups (not implemented yet)
```
conda install numba
```

You also need to install ```progressbar``` separately using -- otherwise, you might be standing in front of your computer, not knowing what is happening
```
python setup.py install 
```
in the folder where you downloaded and extracted the archive. 






### Installation
Nothing difficult here: simply clone this repository: 
```
git clone https://github.com/jlbouchot/FullParametricSolution.git 
cd FullParametricSolution
```
and you are ready to start playing with things. 

### Testing
Of course the tests or computations should be run on your own data / problems. 
However, you may run a small script to make sure everything works fine with all the packages: 
```
bash myFirstTest.sh
```


### Literature
This repository is created in order to numerically verify the theoretical results introduced in [Weighted block compressed sensing for parametrized function approximation](Papers/FullSolution.pdf). 
Some related and relevant literature can be found in [Multi-level compressed sensing Petrov-Galerkin](https://github.com/jlbouchot/CSPDEs). 

Should you have more questions, you may direct them to me: jlbmathit@gmail.com .
