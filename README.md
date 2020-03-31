# FullParametricSolution
Compressed sensing approximation of the full solution to parametric PDEs
========================================================================

This repository is the prototype python implementation of the works related to compressed sensing approximation for parametric PDEs. 
Specifically, this repository deals with the recovery of the parametric solution and not only of a quantity of interest. 

## Remarks

### Anaconda environment
You will need some of the following packages to run: 

Therefore, I strongly recommend to work with a Docker image or a conda environment, associated with the following file. 
This file is the one I am currently using for presenting numerical results.

Note that this package has been developped with Ubuntu 18.04 and Python 3.6. 
Some functions may not work with other versions of Python (some issues were encountered with Python 2, mainly when using printing functions)

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

