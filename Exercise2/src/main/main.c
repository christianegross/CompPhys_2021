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

void generate_random_state(gsl_matrix_int* lattice, gsl_rng *generator){
	//uses generator to generate a random state of 1 and -1 in lattice with dimension lengthx*lengthy
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
 * @fn double hamiltonian(int *lattice, int length,double h, double J);
 * @brief Calculates the beta*hamiltonian of the given configuration, implements periodic boundary conditions with modulus
 *
 *
 * @param lattice	Array of Spins
 * @param lengthx	Length of lattice in x-direction
 * @param lengthy	Length of lattice in y-direction
 * @param h			Magnetic field strength
 * @param J			Coppling strength
 * 
 * @return hamiltonian
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
 * @fn double mean_spin(int *lattice, int lengthx, int lengthy);
 * @brief Calculates the mean over all spins in the given configuration
 *
 * @param lattice	Array of Spins
 * @param lengthx	Length of lattice in x-direction
 * @param lengthy	Length of lattice in y-direction
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

/**
 * @fn double energy_change(int *lattice, int lengthx, int lengthy, int posx, int posy, double J, double h);
 * @brief calculates the change in energy if the sppin at position (posx, posy) in the lattice is flipped
 * 
 * @param lattice	Array of Spins
 * @param lengthx	Length of lattice in x-direction
 * @param lengthy	Length of lattice in y-direction
 * @param posx		x-position of given spin(row)
 * @param posy		y-position of given spin(column)
 * @param J			coupling constant
 * @param h			strength of external magnetic field
 * 
 * @return change in energy if spin is flipped
 */
double energy_change(gsl_matrix_int* lattice, int i, int j, double h_T, double J_T){
	double change=0;
	change+=2*gsl_matrix_int_get (lattice, i, j)*(h_T//external field
										+J_T*(gsl_matrix_int_get (lattice, i, (j+1)%lattice->size2)//right neighbour with periodic boundary conditions
										+gsl_matrix_int_get (lattice, (i+1)%lattice->size1, j)//lower neighbour with periodic boundary conditions
										+gsl_matrix_int_get (lattice, i, (j-1+lattice->size2)%lattice->size2)//left neighbour with periodic boundary conditions
										+gsl_matrix_int_get (lattice, (i-1+lattice->size1)%lattice->size1, j)));//upper neighbour with periodic boundary conditions
	return change;
}
	

int amount_conf(int length){
	return pow (2, 21);//length-1+8); 
	//amount of configurations can be set to be different per length, because phase space has different sizes per length
	//measurements show not needed
}

void complete_sweep(gsl_matrix_int* lattice, gsl_rng *generator,double h_T, double J_T){
	double delta_S=0;
	double randomstate=0;
	int amount_accepts=0;
	for(int i=0;i<lattice->size1; i+=1){
		for(int j=0;j<lattice->size2; j+=1){
			delta_S=energy_change(lattice,i,j,h_T,J_T);
			randomstate=gsl_rng_uniform(generator);
			if(randomstate<=exp (-delta_S)){
				gsl_matrix_int_set (lattice, i, j, -1*gsl_matrix_int_get (lattice, i, j));
				amount_accepts++;
			}
		}
	}
	printf ("%f\n",((double)amount_accepts)/lattice->size1/lattice->size2);
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
	int N_x=4;
	int N_y=4;
	double h=-1;
	double T=1;
	double J=1;
	double h_T=h/T;
	double J_T=J/T;
	int therm_sweeps=200;
	int amount_meas=10;
	

	
	//set and allocate random number generator
	int seed=42;//use fixed seed: result should be exactly reproduced using the same seed
	gsl_rng *generator;
	generator=gsl_rng_alloc(gsl_rng_mt19937);//use mersenne-twister
	gsl_rng_set(generator, seed);
	
	gsl_block *means_spin=gsl_block_alloc(amount_meas);
	gsl_matrix_int* lattice=gsl_matrix_int_alloc (N_x, N_y);
	generate_random_state(lattice,generator);
	
	/**
	 * @note	Thermalization
	 */
	for(int k=0;k<therm_sweeps-1;k++){
		complete_sweep (lattice, generator, h_T, J_T);
		//printf ("%d\n",k);
	}
	/**
	 * @note	Measurements
	 */
	for(int k=0;k<amount_meas;k++){
		complete_sweep (lattice, generator, h_T, J_T);
		means_spin->data[k]=mean_spin(lattice);
		printf ("%f\n",means_spin->data[k]);
	}
	
	
	return 0;
}
