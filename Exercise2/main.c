//Exercise 1 for Computational Physics: Simulation of the 1D Ising-Model
//Nico Dichter, Christiane Groß
//06. November 2020-11. Novemebr 2020
//13. November 2020-
//Further development: 2D, Metropolis-Hastings-algorithm instead of simple sampling

//computation would be faster with char instead of int in the arrays, however use int at first to get human readable output

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_block.h>
#include <gsl/gsl_rng.h>//random number generator
#include "math.h"//exp-Function
#include <sys/time.h>//measuring of wallclock execution time

void generate_random_state(int *lattice, int lengthx, int lengthy, gsl_rng *generator){
	//uses generator to generate a random state of 1 and -1 in lattice with dimension lengthx*lengthy
	double randomstate;
	for(int i=0;i<lengthx; i+=1){
		for(int k=0;k<lengthy; k+=1){
			//generate one random number for each element in lattice, decide if lattice point positive or negative
			randomstate=gsl_rng_uniform(generator);
			if (randomstate<0.5){
				lattice[lengthx*i+k]=-1;
			}
			else{
				lattice[lengthx*i+k]=1;
			}
		}
	}
}

/**
 * @fn double hamiltonian(int *lattice, int length,double h, double J);
 * @brief Calculates the hamiltonian of the given configuration, implements periodic boundary conditions with modulus
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
double hamiltonian(int *lattice, int lengthx, int lengthy, double h, double J){
	double interaction=0;
	for(int i=0;i<lengthx; i++){//all rows
		for(int k=0;k<lengthy; k++){//each column in row
			interaction-=lattice[(i*lengthx)+k]*(h//external magnetic field
												+lattice[(i*lengthx)+(k+1)%lengthy]//right neighbour with periodic boundary conditions
												+lattice[(((i+1)%lengthx)*lengthx)+k]);//lower neighbour with periodic boundary conditions
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
double mean_spin(int *lattice, int lengthx, int lengthy){
	double sum_of_spins=0;
	for(int i=0;i<lengthx; i++){
		for(int k=0;k<lengthy; k++){
			sum_of_spins+=lattice[i*lengthx+k];
		}
	}
	return sum_of_spins/lengthx/lengthy;
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
double energy_change(int *lattice, int lengthx, int lengthy, int posx, int posy, double J, double h){
	double change=0;
	change+=2*lattice[lengthx*posx+posy]*(h//external field
										+lattice[(posx*lengthx)+(posy+1)%lengthy]//right neighbour with periodic boundary conditions
										+lattice[(((posx+1)%lengthx)*lengthx)+posy]//lower neighbour with periodic boundary conditions
										+lattice[(posx*lengthx)+(posy-1+lengthy)%lengthy]//left neighbour with periodic boundary conditions
										+lattice[(((posx-1+lengthx)%lengthx)*lengthx)+posy]);//upper neighbour with periodic boundary conditions
	return change;
}
	

int amount_conf(int length){
	return pow (2, 21);//length-1+8); 
	//amount of configurations can be set to be different per length, because phase space has different sizes per length
	//measurements show not needed
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
	int length_max=20;
	int lattice[length_max];
	const int amount_N=10;
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
	gsl_block_int *lengths=gsl_block_int_alloc(amount_N);
	FILE * savedata=fopen ("data/N_h.dat", "w");
	fprintf(savedata,"#N\th\t<m>\t<m>_err\n");
	
	/**
	 * @note	Set the wanted values for N
	 * 			Itterate through each paramter set: N, h
	 */
	for(int i=0;i<amount_N;i++){
		lengths->data[i]=2*i+2;
	}
	for(int j=0;j<amount_N;j++){
		length=lengths->data[j];
		printf ("Starting calculation for N=%d....\n",length);
		conf=amount_conf (length);
		
		for(h=-1;h<1.1;h+=0.2){
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
			fprintf(savedata,"%d\t%f\t%f\t%f\n",length,h,magnetization,sqrt (magnet_var));
		}
	}
	
	/**
	 * @note	Cleanup.
	 */
	gsl_rng_free(generator);
	gsl_block_free (b_weights);
	gsl_block_free (means_spin);
	gsl_block_int_free (lengths);
	fclose (savedata);
	return 0;
}
