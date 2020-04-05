#!/bin/bash
#title           :testsForFullSolution.sh
#description     :Bunch of tests for the full solution scenario.
#author		 :Jean-Luc Bouchot
#date            :2020/04/05
#version         :0.1    
#usage		 :bash testsForFullSolution.sh
#notes           :Install FEniCS, CVXPY, progressbar before using. See README.md
#===============================================================================

# Script to test if the full solution works fine 
#
# It deals with a 2 (spatial) dimensional problem in 20 (parametric) dimensions. 

 
# List of parameters
d=5 # number of cosines
start_h0=20 # will be used as the discretization step for the first level. 
Lmax=1 # Target number of discretization steps -> Let this vary! 
# No need in this example !! Lmin=3 # First level of discretization (this one will use a single level approximation)
algo=whtp # Which algorithm should be used
vj=1.05 # Value of the constant coefficients
nbSamples=pragmatic # What should be the number of samples 
nbtests=100 # A few tests at the end to make sure it somewhat worked
powerTrig=2.5 # Power of the trigonometric decay
abar=10 # Constant mean field
flucImportance=1 # Importance of the fluctuations
sL=40 # Constant appearing in front of the number of samples of the target discretization
dotensor=FALSE # Use a tensor-based computation instead of building the whole sensing matrix
wCosine=0.25 # How important the 'j' component in the decay is
p0=0.25 # Compressibility in the original space
p=0.3 # Compressibility in the smoothness scale
#sJ=40 # Constant used for the first level of approximation
sJ0=10 # Constant appearing in front of the number of samples of the target discretization
J=$Lmax
nbIter=50

expBasename="FullSolutionResults"

for ((sJ=$sJ0; sJ<=100; sJ=sJ+10))
do
	echo "RUNNING THE EXPERIMENT WITH sJ = $sJ"
	folder=$expBasename$sL
        python test_FullSolution_2D.py -d $d -o WeightedCosine2D -L $Lmax -s $J -x $start_h0 -y $start_h0 -t $nbSamples -r $algo -g $vj -n $nbtests -p $powerTrig -a $abar -c $sL -b $dotensor -i $flucImportance -w $wCosine --smooth_0 $p0 --smooth_t $p --const_sJ $sJ -f $folder -N $nbIter #--onb ONB.npy 
done


#for ((L=$J; L<=$Lmax; L++))
#do
#	echo "RUNNING THE EXPERIMENT WITH L = $L"
#	folder=$expBasename$L
#	python test_wCosine_2D_avg_v_ML.py -d $d -o WeightedCosine2D -L $L -s $J -x $start_h0 -y $start_h0 -t $nbSamples -r $algo -g $vj -n $nbtests -p $powerTrig -a $abar -c $sL -b $dotensor -i $flucImportance -w $wCosine --smooth_0 $p0 --smooth_t $p --const_sJ $sJ -f $folder
#done
#
#python test_FullSolution_2D.py -d $d -o WeightedCosine2D -L $Lmax -s $J -x $start_h0 -y $start_h0 -t $nbSamples -r $algo -g $vj -n $nbtests -p $powerTrig -a $abar -c $sL -b $dotensor -i $flucImportance -w #$wCosine --smooth_0 $p0 --smooth_t $p --const_sJ $sJ -f $folder -N $nbIter #--onb ONB.npy 
