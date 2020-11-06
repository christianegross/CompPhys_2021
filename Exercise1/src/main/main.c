//Exercise 1 for Computational Physics: Simulation of the 1D Ising-Model
//Nico dichter, Christiane Groß
//06. November 2020

//computation would be faster with char instead of int in the arrays, however use int at first to get human readable output

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>//random number generator
#include "math.h"//exp-Function
#include <sys/time.h>//measuring of wallclock execution time

void generate_random_state(int *lattice, int length, gsl_rng *generator){
	//uses generator to generate a random state of 1 and -1 in lattice with length length
	double randomstate;
	for(int i=0;i<length; i+=1){
		//generate one random number for each element in lattice, decide if lattice point positive or negative
		randomstate=gsl_rng_uniform(generator);
		if (randomstate<0.5){
			lattice[i]=-1;
		}
		else{
			lattice[i]=1;
		}
	}
}

/**
 * @fn double hamiltonian(int *lattice, int length,double h, double J);
 * @brief Calculates the hamiltonian of the given configuration
 *
 *
 * @param lattice	Array of Spins
 * @param length	Length of lattice
 * @param h			Magnetic field strength
 * @param J			Coppling strength
 * 
 * @return hamiltonian
 */
double hamiltonian(int *lattice, int length,double h, double J){
	double interaction=0;
	for(int i=0;i<length-1; i++){
		interaction-=J*lattice[i]*lattice[i+1]+h*lattice[i];
	}
	interaction-=J*lattice[length-1]*lattice[0]+h*lattice[length-1];
	return interaction;
}

/**
 * @fn double mean_spin(int *lattice, int length);
 * @brief Calculates the mean over all spins in the given configuration
 *
 * @param lattice	Array of Spins
 * @param length	Length of lattice
 * 
 * @return mean over all spins
 *
 */
double mean_spin(int *lattice, int length){
	double sum_of_spins=0;
	for(int i=0;i<length-1; i++){
		sum_of_spins+=lattice[i];
	}
	return sum_of_spins/length;
}

int main(int argc, char **argv){
	/**
	 * Declarations:
	 *
	 * @var length		length of each configuartion
	 * @var h			Magnetic field strength
	 * @var conf		#configurations
	 * @var Temp		Temperatur of the heat bath
	 * @var J			Coppling strength
	 * @var lattice		Array of Spins=single configuration
	 *
	 */
	const int length=atoi(argv[1]);
	const double h=atof(argv[2]);
	const int conf=pow (2, length-1)*2;
	const double Temp=1;
	const double J=1;
	int lattice[length];
	
	double part_fct=0;
	double boltz_weight=0;
	double magnetization=0;
	//set and allocate random number generator
	int seed=2;//use fixed seed: result should be exactly reproduced using the same seed
	gsl_rng *generator;
	
	generator=gsl_rng_alloc(gsl_rng_mt19937);//use mersenne-twister
	gsl_rng_set(generator, seed);
	//try if allocations from command line, randeóm state generation work as intended
	printf("l=%d\th=%f\t#configurations=%d\n", length, h,conf);
	
	/**
	 * @note	Generate a new configuration "conf"-times.
	 * 			Calculate the boltzman weight and the weighted mean-spin 
	 * 			of each configuration.
	 * 			At the end the expected value of magnetization gets calculated.
	 *
	 */
	for(int i=0;i<conf; i+=1){
		generate_random_state(lattice,length,generator);
		boltz_weight=exp (-hamiltonian (lattice, length, h,  J)/Temp);
		part_fct+=boltz_weight;
		magnetization+=mean_spin(lattice,length)*boltz_weight;
	}
	magnetization/=part_fct;
	printf("<m>=%f\n",magnetization);
	gsl_rng_free(generator);
	return 0;
}
