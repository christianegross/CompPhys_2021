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
#include <gsl/gsl_randist.h>
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

void leapfrog(double p_0,double phi_0,double* p_f,double* phi_f, int N_md,
			   double J_hat_T, double h_T, int N){
	double epsilon=1/((double)N_md);
	*p_f=p_0;
	*phi_f=phi_0;
	*phi_f=*phi_f+*p_f*epsilon/2;
	for(int i=1;i<N_md;i++){
		*p_f=*p_f-epsilon*(*phi_f/(J_hat_T)-N*tanh (h_T+*phi_f));
		*phi_f=*phi_f+*p_f*epsilon;
	}
	*p_f=*p_f-epsilon*(*phi_f/(J_hat_T)-N*tanh (h_T+*phi_f));
	*phi_f=*phi_f+*p_f*epsilon/2;
}

double art_hamiltonian(double p,double phi, double J_hat_T, double h_T, int N){
	return p*p/2.+phi*phi/(2.*J_hat_T)-log (2.*cosh (h_T+phi))*N;
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
	 * @var amount_conf			#measurements
	 */
	
	int therm_steps=500;
	int amount_conf=5000;

	int N=5;
	double J_T=0.2;
	double J_hat_T=J_T/N;
	double h_T=0.5;
	
	
	//set and allocate random number generator
	int seed=2;//use fixed seed: result should be exactly reproduced using the same seed
	gsl_rng *generator;
	generator=gsl_rng_alloc(gsl_rng_mt19937);//use mersenne-twister
	gsl_rng_set(generator, seed);
	
	/**
	 * @note	Allocate blocks for several measurements
	 * 			and allocate memory for the largest lattice and set its initial state randomly
	 */
	gsl_block * set_of_phis=gsl_block_alloc (amount_conf);
	gsl_block * magnetizations=gsl_block_alloc (amount_conf);
	
	/**
	 * @note	Open file streams to save data into
	 */
	FILE * converge_data=fopen ("data/converge.dat", "w");
	fprintf(converge_data,"#N_md\tH_rel_delta\n");
	FILE * raw_data=fopen ("data/raw.dat", "w");
	fprintf(raw_data,"#N\tJ\t<m>\n");
	
	double p_0=1;
	double phi_0=2;
	double p_f=0;
	double phi_f=0;
	double delta=0;
	double prob=1;
	int N_md=4;
	double mean_magnetization=0;
	double accept_rate=0;
	int amount_ar=0;
	double H_0=art_hamiltonian (p_0, phi_0,  J_hat_T,  h_T,  N);
	for(int i=2;i<101;i++){
		leapfrog (p_0, phi_0, &p_f, &phi_f, i,  J_hat_T, h_T,  N);
		delta=(art_hamiltonian (p_f, phi_f,  J_hat_T,  h_T,  N)-H_0)/H_0;
		fprintf (converge_data,"%d\t%.10f\n",i , fabs(delta));
		
	}
	for(N=5;N<21;N+=5){
		for(J_T=0.2;J_T<2.0001;J_T+=0.01){
			J_hat_T=J_T/((double)N);
			mean_magnetization=0;
			for(int i=0;i<therm_steps;i++){
				p_0=gsl_ran_ugaussian (generator);
				H_0=art_hamiltonian (p_0, phi_0,  J_hat_T,  h_T,  N);
				leapfrog (p_0, phi_0, &p_f, &phi_f, N_md,  J_hat_T, h_T,  N);
				prob=exp(H_0-art_hamiltonian (p_f, phi_f,  J_hat_T,  h_T,  N));
				if(gsl_rng_uniform(generator)<prob){
					phi_0=phi_f;
				}
			}
			
			for(int i=0;i<amount_conf;i++){
				p_0=gsl_ran_ugaussian (generator);
				H_0=art_hamiltonian (p_0, phi_0,  J_hat_T,  h_T,  N);
				leapfrog (p_0, phi_0, &p_f, &phi_f, N_md,  J_hat_T, h_T,  N);
				prob=exp(H_0-art_hamiltonian (p_f, phi_f,  J_hat_T,  h_T,  N));
				if(gsl_rng_uniform(generator)<prob){
					phi_0=phi_f;
					accept_rate+=1.;
				}
				amount_ar++;
				set_of_phis->data[i]=phi_0;
				magnetizations->data[i]=tanh (h_T+set_of_phis->data[i]);
				mean_magnetization+=magnetizations->data[i];
			}
			mean_magnetization/=amount_conf;
			fprintf (raw_data,"%d\t%f\t%f\n",N,J_T,mean_magnetization);
		}
	}
	
	accept_rate/=amount_ar;
	printf ("Acceptance Rate:%f\n",accept_rate);
	
	
	
	/**
	 * @note	Cleanup
	 */
	fclose (converge_data);
	fclose (raw_data);
	gsl_rng_free (generator);
	gsl_block_free (set_of_phis);
	return 0;
}
