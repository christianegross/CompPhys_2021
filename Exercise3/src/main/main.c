//Exercise 1 for Computational Physics: Simulation of the 1D Ising-Model
//Nico Dichter, Christiane Groß
//06. November 2020-11. Novemebr 2020
//13. November 2020-
//Further development: 2D, Metropolis-Hastings-algorithm instead of simple sampling

//computation would be faster with char instead of int in the arrays, however use int at first to get human readable output

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_block.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>//random number generator
#include "math.h"//exp-Function
#include <sys/time.h>//measuring of wallclock execution time

/**
 * @fn void generate_random_state(gsl_matrix_int* lattice, gsl_rng *generator);
 * @brief uses generator to generate a random state of 1 and -1 in lattice
 *
 *
 * @param lattice	Matrix of Spins
 * @param generator	used for random states
 */
void generate_random_state(gsl_matrix_int* lattice, gsl_rng *generator){
	double randomstate;
	for(int i=0;i<lattice->size1; i+=1){
		for(int j=0;j<lattice->size2; j+=1){
			//generate one random number for each element in lattice, decide if lattice point positive or negative
			randomstate=gsl_rng_uniform(generator);
			if (randomstate<0.5){
				gsl_matrix_int_set (lattice, i, j, -1);
			}
			else{
				gsl_matrix_int_set (lattice, i, j, 1);
			}
		}
	}
}



/**
 * @fn double b_hamiltonian(gsl_matrix_int* lattice, double h_T, double J_T);
 * @brief Calculates the beta*hamiltonian of the given configuration, implements periodic boundary conditions with modulus
 * @note at the moment, beta=1 for alle measurements, so not considered
 *
 *
 * @param lattice	Matrix of Spins
 * @param h_T		Magnetic field strength
 * @param J_T		Coppling strength
 * 
 * @return beta*hamiltonian
 */
double b_hamiltonian(gsl_matrix_int* lattice, double h_T, double J_T){
	double interaction=0;
	for(int i=0;i<lattice->size1; i++){//all rows
		for(int j=0;j<lattice->size2; j++){//each column in row
			interaction-=gsl_matrix_int_get (lattice, i, j)*(h_T//external magnetic field
												+J_T*(gsl_matrix_int_get (lattice, i, (j+1)%lattice->size2)//right neighbour with periodic boundary conditions
												+gsl_matrix_int_get (lattice, (i+1)%lattice->size1, j)));//lower neighbour with periodic boundary conditions
		}
	}
	return interaction;
}


/**
 * @fn double mean_spin(gsl_matrix_int* lattice);
 * @brief Calculates the mean over all spins in the given configuration
 *
 * @param lattice	Matrix of Spins
 * 
 * @return mean over all spins
 *
 */
double mean_spin(gsl_matrix_int* lattice){
	double sum_of_spins=0;
	for(int i=0;i<lattice->size1; i+=1){
		for(int j=0;j<lattice->size2; j+=1){
			sum_of_spins+=gsl_matrix_int_get (lattice, i, j);
		}
	}
	return sum_of_spins/lattice->size1/lattice->size2;
}



int main(int argc, char **argv){
	/**
	 * Declarations:
	 *
	 * @var N_x_max,N_y_max		Maximum dimensions of the lattice
	 * @var N					Current size of the lattice (side length)
	 * @var lambda				#sites
	 * @var h_T					Externel magnetic field/Temperatur
	 * @var J_T					Coppling strength/Temperatur
	 * @var therm_sweeps		#thermalization sweeps
	 * @var amount_meas			#measurements
	 */
	int N_x_max=20;
	int N_y_max=20;
	int N=4;
	int lambda=N*N;
	double h_T=-1;
	double J_T=0.25;
	double h_T_max=1.;
	double J_T_max=2.;
	int therm_sweeps=1000;
	int amount_meas=2000;
	
	double magnetization=0;
	double squared_mean=0;
	double abs_magnetization=0;
	double avr_energy_ps=0;
	double avr_squared_energy_ps=0;
	

	
	//set and allocate random number generator
	int seed=2;//use fixed seed: result should be exactly reproduced using the same seed
	gsl_rng *generator;
	generator=gsl_rng_alloc(gsl_rng_mt19937);//use mersenne-twister
	gsl_rng_set(generator, seed);
	
	/**
	 * @note	Allocate blocks for several measurements
	 * 			and allocate memory for the largest lattice and set its initial state randomly
	 */
	gsl_block *means_spin=gsl_block_alloc(amount_meas);
	gsl_block *energy_ps=gsl_block_alloc(amount_meas);
	gsl_matrix_int* lattice_mem=gsl_matrix_int_alloc (N_x_max, N_y_max);
	generate_random_state(lattice_mem,generator);
	gsl_matrix_int_view lattice;
	
	/**
	 * @note	Open file streams to save data into
	 */
	FILE * savedata=fopen ("data/N_J_h.dat", "w");
	fprintf(savedata,"#N\tJ\th\t<m>\t<m>_err\t<e>\t<e>_err\n");
	FILE * savedata_abs_m=fopen ("data/abs_m.dat", "w");
	fprintf(savedata_abs_m,"#N\tJ\t<|m|>\t<|m|>_err\n");
	
	
	/**
	 * @note	Cleanup
	 */
	gsl_block_free (means_spin);
	gsl_block_free (energy_ps);
	gsl_rng_free (generator);
	gsl_matrix_int_free (lattice_mem);
	fclose (savedata);
	fclose (savedata_abs_m);
	return 0;
}
