//Exercise 3 for Computational Physics: Simulation of the 1D Ising-Model
//Nico Dichter, Christiane Groß
//16. November 2020-

//computation would be faster with char instead of int in the arrays, however use int at first to get human readable output

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_block.h>
#include <gsl/gsl_rng.h>//random number generator
#include <gsl/gsl_randist.h>
#include "math.h"//exp-Function


/**
 * @fn double mean_spin(gsl_matrix_int* lattice);
 * @brief Calculates the mean over all spins in the given configuration
 *
 * @param lattice	Matrix of Spins
 * 
 * @return mean over all spins
 *
 */
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
	 * @var therm_steps			#thermalization steps
	 * @var amount_conf			#configurations
	 * @var N					Size of the lattice
	 * @var J_T=J/T				Coppling strength
	 * @var h_T=h/T				Strength of external field
	 * @var p					Conjugated momentum
	 * @var phi					Continuous paramater of H
	 * @var delta				Relative error when going fom phi,p->phi',p'
	 * @var prob				Probability to accept a proposed phi'
	 * @var N_md				#steps in leap frog
	 * @var mean_magnetization	mean magnetization
	 * @var accept_rate			acceptance rate
	 * @var amount_ar			amount of accept/reject steps
	 */
	
	int therm_steps=500;
	int amount_conf=5000;

	int N=5;
	double J_T=0.2;
	double J_hat_T=J_T/N;
	double h_T=0.5;
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
	
	//set and allocate random number generator
	int seed=2;//use fixed seed: result should be exactly reproduced using the same seed
	gsl_rng *generator;
	generator=gsl_rng_alloc(gsl_rng_mt19937);//use mersenne-twister
	gsl_rng_set(generator, seed);
	
	/**
	 * @note	Allocate blocks for several measurements
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
	
	
	/**
	 * @note	Test rather the convergence conditions for leap frog are met
	 */
	double H_0=art_hamiltonian (p_0, phi_0,  J_hat_T,  h_T,  N);
	for(int i=2;i<101;i++){
		leapfrog (p_0, phi_0, &p_f, &phi_f, i,  J_hat_T, h_T,  N);
		delta=(art_hamiltonian (p_f, phi_f,  J_hat_T,  h_T,  N)-H_0)/H_0;
		fprintf (converge_data,"%d\t%.10f\n",i , fabs(delta));
		
	}
	/**
	 * @note	Itterate through different sets of N/J
	 */
	for(N=5;N<21;N+=5){
		for(J_T=0.2;J_T<2.0001;J_T+=0.01){
			J_hat_T=J_T/((double)N);
			mean_magnetization=0;
			
			/**
			 * @note	Thermalization:
			 * 			1. sample p
			 * 			2. leapfrog
			 * 			3. accept/reject
			 */
			for(int i=0;i<therm_steps;i++){
				p_0=gsl_ran_ugaussian (generator);
				H_0=art_hamiltonian (p_0, phi_0,  J_hat_T,  h_T,  N);
				leapfrog (p_0, phi_0, &p_f, &phi_f, N_md,  J_hat_T, h_T,  N);
				prob=exp(H_0-art_hamiltonian (p_f, phi_f,  J_hat_T,  h_T,  N));
				if(gsl_rng_uniform(generator)<prob){
					phi_0=phi_f;
				}
			}
			/**
			 * @note	Generate ensemble of phis:
			 * 			1. sample p
			 * 			2. leapfrog
			 * 			3. accept/reject
			 * 			And take measurements for <m> and <epsilon>
			 */
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
	
	/**
	 * @note	Print out the overall acceptance rate
	 */
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
