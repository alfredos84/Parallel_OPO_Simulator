#!/bin/bash

# This file contains a set of instructions to run the main file cuOPO.cu.
# Use this file to perform simulations in bulk. This means, for example, 
# systematically varying the input power, the cavity reflectivity, etc.

TIMEFORMAT='It took %0R seconds.' 
time { 

clear      # Clear screen
rm *.dat   # This removes all .dat files in the current folder. Comment this line for safe.
rm cuOPO   # This removes a previuos executable file (if it exist)


#####################################################################
# -----Compilation-----
# There are two compilation modes: with and without thermal calculations. To set the thermal calculations
# define the preprocessor variable -DTHERMAL in the compilation line.

# Uncomment the compilation line in line 22 for the inclusion of thermal calculations and comment the 
# compilation line in line 27.
nvcc cuOPO.cu -diag-suppress 177 --gpu-architecture=sm_60 -lcufftw -lcufft -o cuOPO
FOLDERSIM="GDD"

# There are three flags specific for CUDA compiler:
# --gpu-architecture=sm_60: please check your GPU card architecture (Ampere, Fermi, Tesla, Pascal, Turing, Kepler, etc) 
#                           to set the correct number sm_XX. This code was tested using a Nvidia GeForce GTX 1650 card (Turing
#                           architecture). 
# 				    Please visit https://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/index.html#gpu-feature-list
#                           to set the proper flag according to your GPU card.
# -lcufftw -lcufft        : flags needed to perform cuFFT (the CUDA Fourier transform)
#####################################################################

# The variables defined below (ARGX) will be passed as arguments to the main file 
# cuOPO.cu on each execution via the argv[X] instruction.

ARG1=(16)     	# Pump power                (ARG1)

# x=$(awk 'BEGIN{for(i=20; i<=40; i=i+1)print i}')
# n=$(awk 'BEGIN{for(i=1; i<=31; i=i+3)print i}')
# t1=$(awk 'BEGIN{for(i=54; i<=59; i=i+0.1)print i}')
# Each for-loop span over one or more values defined in the previous arguments. 
# 		for N in $x
# 		do


for (( p=0; p<${#ARG1[@]}; p++ ))
do
	N=${ARG1[$p]}
	printf "\nPower			= ${N} W\n" 

	printf "\nMaking directory...\n"
	FOLDER="sPPLT_N_${N}"
	FILE="sPPLT_N_${N}.txt"
	
	printf "Bash execution and writing output file...\n\n"
	./cuOPO $N | tee -a $FILE
	
	printf "Bash finished!!\n\n" 
	mkdir $FOLDER
	mv *.dat $FOLDER"/"
	mv FILE.txt $FOLDER"/"
done


if [ -d "$FOLDERSIM" ]; then
	echo "Moving simulations in ${FOLDERSIM}..."
	mv sPPLT* $FOLDERSIM"/" 
	mv *.txt $FOLDERSIM"/" 
else

	mkdir $FOLDERSIM
	echo "Creating and moving simulations in ${FOLDERSIM}..."
	mv sPPLT* $FOLDERSIM"/" 
	mv *.txt $FOLDERSIM"/" 
fi

}
