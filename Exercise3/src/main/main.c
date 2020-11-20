//Exercise 3 for Computational Physics: Simulation of the long range Ising-Model
//Nico Dichter, Christiane Groß
//16. November 2020-

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_block.h>
#include <gsl/gsl_rng.h>//random number generator
#include <gsl/gsl_randist.h>//pull random number from gaussian
#include "math.h"//exp-Function
#include <gsl/gsl_vector.h>
double art_hamiltonian(double p,double phi, double J_hat_T, double h_T, int N);

/**
 * @fn void leapfrog(double p_0,double phi_0,double* p_f,double* phi_f, int N_md,
			   double J_hat_T, double h_T, int N);
 * @brief Integrates using the force equations of the hamiltonian, by N_md steps
 * 			(using leap frog algorithm)
 *
 * @param p_0		Intial momentum
 * @param phi_0		Intial "position"
 * @param p_f		Final momentum
 * @param phi_f		Final "position"
 * @param N_md		#steps in leap frog
 * @param J_hat_T	J_T/N
 * @param h_T=h/T	Strength of external field
 * @param N			Size of the lattice
 * 
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

inline double art_hamiltonian(double p,double phi, double J_hat_T, double h_T, int N){
	return p*p/2.+phi*phi/(2.*J_hat_T)-log (2.*cosh (h_T+phi))*N;
}


/**
 * @fn 		inline void binning(gsl_block *measurements, gsl_block *binneddata, int lengthofbin)
 * @brief 	takes data points in measurements and writes bins over lengthofbin in binneddata
 * 
 * @note 	inline makes execution time faster (proposes different implementation to the compiler)
 * 
 * @param measurements 	determind through algorithm, correlated
 * @param binneddata 	holds results from binning
 * @param lengthofbin	determins how many elements are in one bin, so also has a role in how many elements there are in binneddata
 */
inline void binning(gsl_block *measurements, gsl_vector *binneddata, int lengthofbin){
	double bincontent=0;
	//assertion if lengths of bins, measurements fit together
	int numberofbins=measurements->size/lengthofbin;
	//printf("numberofbins=%d\tnumberofelemtns=%lu\tlengthofbin=%d\n", numberofbins, binneddata->size, lengthofbin);
	if (numberofbins!=binneddata->size){fprintf(stderr, "Gave wrong lengths! %d number of bins to be calculated, but storage allocated for %lu!", numberofbins, binneddata->size);}
	//calculate content of single bin as arithmetic mean oover lengthofbin datapoints
	for (int bin=0; bin<numberofbins; bin+=1){
		bincontent=0;
		for (int datapoint=0; datapoint<lengthofbin; datapoint+=1){
			bincontent+=measurements->data[bin*lengthofbin+datapoint];
		}
		gsl_vector_set(binneddata, bin, bincontent/lengthofbin);
	}
}

/**
 * @fn makebootstrapreplica(gsl_block * measurements, gsl_rng * generator)
 * @brief makes one bootstrapreplica out of the data in measurements
 * 
 * @param measurements contains all the measurements which can be used to calculate the replica
 * @param generator random number generator which determines which elements of measurements are used to calaculate the replica
 * 
 * @return replica arithmetic mean of chosen measurements, one bootstrapreplica
 */
double makebootstrapreplica(gsl_vector * measurements, gsl_rng * generator){
	double replica=0;
	int randomnumber;
	for (int datapoint=0; datapoint<measurements->size; datapoint+=1){
		randomnumber=gsl_rng_uniform_int(generator, measurements->size);
		replica+=gsl_vector_get(measurements, randomnumber);
	}
	return replica/measurements->size;
}
/**
 * @fn bootstrap(gsl_block *measurements, gsl_rng *generator, int R, double *mean, double *variance)
 * @brief uses the bootstrapmethod to calculate mean and variance of the data in measurements
 * 
 * @param measurements contains all the datappoints which are used
 * @param generator random number generator which determines which elements of measurements are used to calaculate the replica
 * @param R number of replicas which are used to calculate mean and variance
 * @param calculated mean over all replicas
 * @param variance calculated variance over all replica
 */ 
void bootstrap(gsl_vector *measurements, gsl_rng *generator, int R, double *mean, double *variance){
	*mean=0;
	*variance=0;
	gsl_block *replicalist=gsl_block_alloc(R);
	for (int i=0; i<R; i+=1){
		replicalist->data[i]=makebootstrapreplica(measurements, generator);
		*mean+=replicalist->data[i];
	}
	//calculate arithmetic mean over replicas
	*mean/=R;
	//calculate standard deviation of replicas
	for (int i=0; i<R; i+=1){
		*variance+=(*mean-replicalist->data[i])*(*mean-replicalist->data[i]);
	}
	*variance/=R-1;
	gsl_block_free(replicalist);
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
	 * @var var_magnetization	Variance of the magnetization
	 * @var mean_energy_p_site	mean energy per site
	 * @var var_energy_p_site	Variance of the energy per site
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
	double var_magnetization=0;
	double mean_energy_p_site=0;
	double var_energy_p_site=0;
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
	gsl_block * energies=gsl_block_alloc (amount_conf);
	gsl_vector * binnedmagnetization_mem=gsl_vector_alloc(amount_conf);
	gsl_vector_view binnedmagnetization;
	gsl_vector * binnedenergy_mem=gsl_vector_alloc(amount_conf);
	gsl_vector_view binnedenergy;
	/**
	 * @note	Open file streams to save data into
	 */
	FILE * converge_data=fopen ("data/converge.dat", "w");
	fprintf(converge_data,"#N_md\tH_rel_delta\n");
	FILE * raw_data=fopen ("data/raw.dat", "w");
	fprintf(raw_data,"#N\tJ\t<m>\t<m>_err\t<epsilon>\t<epsilon>_err\tlofb\n");
	
	
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
		for(J_T=0.2;J_T<2.0001;J_T+=0.05){
			J_hat_T=J_T/((double)N);
			mean_magnetization=0;
			var_magnetization=0;
			mean_energy_p_site=0;
			var_energy_p_site=0;
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
				energies->data[i]=-(set_of_phis->data[i])*(set_of_phis->data[i])/2/N/J_hat_T-h_T*tanh (h_T+set_of_phis->data[i])-1;
			}
			/**
			 * @note	analyse measured data for magnetization and energy:
			 * 			binning to reduce correltations
			 * 			bootstrapping to get a good error estimate
			 * 	
			 * @note	write binned data in vector instead of block, so size can be changed with subvector-view. 
			 * 			block was tried, but led to problems with memory in reallocing.
			 */ 
			for (int lengthofbin=2; lengthofbin<1023; lengthofbin*=2){
				binnedmagnetization=gsl_vector_subvector(binnedmagnetization_mem, 0, magnetizations->size/lengthofbin);
				binning(magnetizations, &binnedmagnetization.vector, lengthofbin);
				bootstrap(&binnedmagnetization.vector, generator, 4*amount_conf, &mean_magnetization, &var_magnetization);
				binnedenergy=gsl_vector_subvector(binnedenergy_mem, 0, energies->size/lengthofbin);
				binning(energies, &binnedenergy.vector, lengthofbin);
				bootstrap(&binnedenergy.vector, generator, 4*amount_conf, &mean_energy_p_site, &var_energy_p_site);
				fprintf (raw_data,"%d\t%f\t%f\t%f\t%f\t%f\t%d\n", N, J_T, mean_magnetization,sqrt (var_magnetization), mean_energy_p_site ,sqrt (var_energy_p_site), lengthofbin);
			}
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
	gsl_block_free (magnetizations);
	gsl_vector_free(binnedmagnetization_mem);
	gsl_block_free (energies);
	gsl_vector_free(binnedenergy_mem);
	return 0;
}
