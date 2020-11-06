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
		randomstate=gsl_rng_uniform_pos(generator);
		if (randomstate<0.5){
			lattice[i]=-1;
		}
		else{
			lattice[i]=1;
		}
	}
}

int main(int argc, char **argv){
	const int length=atoi(argv[1]);
	const double h=atof(argv[2]);
	int lattice[length];
	//set and allocate random number generator
	int seed=2;//use fixed seed: result should be exactly reproduced using the same seed
	gsl_rng *generator;
	generator=gsl_rng_alloc(gsl_rng_mt19937);//use mersenne-twister
	gsl_rng_set(generator, seed);
	//try if allocations from command line, randeóm state generation work as intended
	printf("l=%d\th=%f\n", length, h);
	generate_random_state(lattice, length, generator);
	for(int i=0;i<length; i+=1){
		printf(" %d ", lattice[i]);
	}
	printf("\n");
	gsl_rng_free(generator);
	return 0;
}
