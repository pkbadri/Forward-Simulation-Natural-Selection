# Forward-Simulation-Natural-Selection

The program given below is based on the algorithm described in Padhukasahasram et al. 2008.


https://github.com/pkbadri/Forward-Simulation-Natural-Selection/blob/main/newsel.c

https://github.com/pkbadri/Forward-Simulation-Natural-Selection/blob/main/mtrand.h

https://github.com/pkbadri/Forward-Simulation-Natural-Selection/blob/main/mtrand.cpp



Download the 3 program files in the links above. newsel.c can simulate a Wright-Fisher population of constant size undergoing mutation, recombination and natural selection at multiple sites. First create a file named newselinput.txt in the same directory in which you want to run this program. Fix the number and positions of selected sites and selection coefficients in newselinput.txt. This file shud have the following:


1st line: Number of sites under selection.
2nd line: Positions of selected sites in ascending order separated by spaces. Positions vary from 0 to len - 1, where len is the length of the DNA sequence in base pairs.
3rd line: Respective selection coefficients for the above positions for the genotype heterozygous for alleles (10 or 01). These should be separated by spaces.
4th line: Respective selection coefficients for the above positions for the genotype homozygous for one allele (11). These should be separated by spaces.
5th line: Respective selection coefficients for the above positions for the genotype homozygous for other allele (00). These should be separated by spaces.

For example:
5
0 100 200 300 499
0.5 0.5 0.5 0.5 0.5
0 0 0 0 0
0 0 0 0 0

The fitness for each genotype is equal to 1 + selection coefficient. Fitness effects are either added or multiplied across the sites. Compile using g++ -O3 -o sel newsel.c mtrand.cpp -Wno-deprecated. There are 10 different command line parameters:

-samples (Number of samples to output from final population. Choose < pop/2)
-pop (Total number of chromsomes in the diploid population. Should be a even number.)
-len (Length of the DNA sequence in basepairs. Choose larger than del*u*pop)
-r (Per-generation per-sequence rate of recombination)
-u (Per-generation per-sequence rate of mutation)
-s (Probability of self-fertilization)
-gen (Total number of generations)
-del (Generations after which fixed mutations are removed)
-reps (Number of replicates)
-fitness (1 to add fitness effects. 2 to multiply effects)

To run type ./sel along with all the command line parameters in the same order. For example: ./sel -samples 100 -pop 1000 -len 1000000 -r 0.50 -u 0.50 -s 0.0 -gen 10000 -del 500 -reps 10 -fitness 2. Outputs positions of biallelic polymorphisms and haplotypes in 0-1 format for samples collected from the final population.
