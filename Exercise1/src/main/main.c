//Exercise 1 for Computational Physics: Simulation of the 1D Ising-Model
//Nico dichter, Christiane Groß
//06. November 2020

//computation would be faster with char instead of int in the arrays, however use int at first to get human readable output

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_block.h>
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
	for(int i=0;i<length; i++){
		sum_of_spins+=lattice[i];
	}
	return sum_of_spins/length;
}

int amount_conf(int length){
	return pow (2, length-1+5);
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
	 * @var length_max	Maximum length to be simulated
	 */
	int length=2;
	double h=-1;
	int conf=20;
	double Temp=1;
	const double J=1;
	int lattice[length];
	int length_max=20;
	double part_fct=0;
	double magnetization=0;
	double magnet_var=0;
	
	//set and allocate random number generator
	int seed=2;//use fixed seed: result should be exactly reproduced using the same seed
	gsl_rng *generator;
	generator=gsl_rng_alloc(gsl_rng_mt19937);//use mersenne-twister
	gsl_rng_set(generator, seed);
	
	/**
	 * @note	Allocate memory for all the boltzmann weights and mean 
	 * 			of the spins for a given parameter set. 
	 * 			And open the needed file stream.
	 *
	 */
	gsl_block *b_weights=gsl_block_alloc(amount_conf (length_max));
	gsl_block *means_spin=gsl_block_alloc(amount_conf (length_max));
	char filename [100];
	snprintf (filename, 100, "data/%d.dat",length);
	FILE * savedata=fopen (filename, "w");
	
	/**
	 * @note	Itterate through each paramter set.
	 */
	for(length=2;length<=length_max;length*=2){
		printf ("Starting calculation for N=%d....\n",length);
		snprintf (filename, 100, "data/%d.dat",length);
		if(freopen (filename, "w", savedata)==NULL){
			fprintf(stderr, "Error in freopen.\n");
			exit (EXIT_FAILURE);
		};
		conf=amount_conf (length);
		fprintf(savedata,"#h\tT\t<m>\t<m>_err\n");
		for(h=-1;h<1.1;h+=0.5){
			for(Temp=0.2;Temp<8.1;Temp+=0.2){
				/**
				 * @note	Generate a new configuration "conf"-times.
				 * 			Calculate the boltzman weight and the weighted mean-spin 
				 * 			of each configuration.
				 * 			At the end the expected value of magnetization gets calculated.
				 *
				 */
				part_fct=0;
				magnetization=0;
				magnet_var=0;
				for(int i=0;i<conf; i+=1){
					generate_random_state(lattice,length,generator);
					b_weights->data[i]=exp (-hamiltonian (lattice, length, h,  J)/Temp);
					means_spin->data[i]=mean_spin(lattice,length);
					part_fct+=b_weights->data[i];
					magnetization+=(means_spin->data[i])*(b_weights->data[i]);
				}
				magnetization/=part_fct;
				
				/**
				 * @note	Calculate the variance of the magnetization.
				 * 			And save the data.
				 */
				for(int i=0;i<conf; i+=1){
					magnet_var+=(magnetization-means_spin->data[i])*(magnetization-means_spin->data[i])*(b_weights->data[i]);
				}
				magnet_var/=part_fct;
				fprintf(savedata,"%f\t%f\t%f\t%f\n",h,Temp,magnetization,sqrt (magnet_var));
			}
		}
	}
	
	/**
	 * @note	Cleanup.
	 */
	gsl_rng_free(generator);
	gsl_block_free (b_weights);
	gsl_block_free (means_spin);
	fclose (savedata);
	return 0;
}
