//Exercise 4 for Computational Physics: Simulation of the long range Ising-Model
//Nico Dichter, Christiane Groß
//23. November 2020-

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
inline void binning(gsl_vector *measurements, gsl_vector *binneddata, int lengthofbin){
	double bincontent=0;
	//assertion if lengths of bins, measurements fit together
	int numberofbins=measurements->size/lengthofbin;
	//printf("numberofbins=%d\tnumberofelemtns=%lu\tlengthofbin=%d\n", numberofbins, binneddata->size, lengthofbin);
	if (numberofbins!=binneddata->size){fprintf(stderr, "Gave wrong lengths! %d number of bins to be calculated, but storage allocated for %lu!", numberofbins, binneddata->size);}
	//calculate content of single bin as arithmetic mean oover lengthofbin datapoints
	for (int bin=0; bin<numberofbins; bin+=1){
		bincontent=0;
		for (int datapoint=0; datapoint<lengthofbin; datapoint+=1){
			bincontent+=gsl_vector_get(measurements,bin*lengthofbin+datapoint);
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
		
/**
 * Maybe make functions thermalize, measure, to make sure main is not too long?
 * 
 */

void autocorrelation(gsl_vector *measurements, gsl_vector *results, double mean){
	/**
	 * make block, calculate gamma with tau=0, write (0,1) to results
	 * for range of tau: calculate gamma, divide by gamma(0), write (tau, gamma/gamma) to reaults
	 */
	double czero=0;
	double ctime=0;
	for (int i=0; i<measurements->size; i+=1){
		 czero+=(gsl_vector_get(measurements,i)-mean)*(gsl_vector_get(measurements,i)-mean);
	}
	czero/=measurements->size;
	gsl_vector_set (results, 0,1);
	for (int time=1; time<results->size; time+=1){
		ctime=0;
		for (int i=0; i<measurements->size-time; i+=1){
		 ctime+=(gsl_vector_get(measurements,i)-mean)*(gsl_vector_get(measurements,i+time)-mean);
		}
		ctime/=(measurements->size-time);
		gsl_vector_set (results, time, ctime/czero);
	}
}

int main(int argc, char **argv){
	/**
	 * Declarations:
	 *
	 * @var therm_steps			#thermalization steps
	 * @var amount_conf			#configurations
	 * @var N					Size of the lattice
	 * @var N_bs				#bootstrap replicas
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
	
	int therm_steps=2000;
	int amount_conf=12800;

	int N=5;
	int N_bs=4*amount_conf;
	double J_T=0.1;
	double J_hat_T=J_T/N;
	double h_T=0.5;
	double p_0=1;
	double phi_0=0;
	double p_f=0;
	double phi_f=0;
	double H_0;
	double prob=1;
	int N_md=4;
	double mean_magnetization=0;
	double var_magnetization=0;
	double var_magnetization_binned=0;
	double simple_mean_magnetization=0;
	double simple_mean_magnetization_binned=0;
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
	gsl_vector * set_of_phis=gsl_vector_alloc (amount_conf);
	gsl_vector * magnetizations=gsl_vector_alloc (amount_conf);
	gsl_vector * binnedmagnetization_mem=gsl_vector_alloc(amount_conf);
	gsl_vector_view binnedmagnetization;
	gsl_vector * magnetizationcorrelation =gsl_vector_alloc((int)amount_conf/4);
	gsl_vector * binnedcorrelation =gsl_vector_alloc((int)amount_conf/128);
	/**
	 * @note	Open file streams to save data into
	 */
	FILE * raw_data=fopen ("data/raw4.dat", "w");
	fprintf(raw_data,"#N\tJ\t<m>\t<m>_err\t<m>_err_bin\t<epsilon>\t<epsilon>_err\tlofb\n");
	FILE * magnetization_nmd4=fopen ("data/magnetizationnmd4.dat", "w");
	FILE * magnetizationcorrelation_nmd4=fopen("data/magnetizationcorrelationnmd4.dat", "w");
	

	mean_magnetization=0;
	var_magnetization=0;
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
		gsl_vector_set (set_of_phis, i, phi_0);
		gsl_vector_set(magnetizations, i, tanh (h_T+phi_0));
		simple_mean_magnetization+=gsl_vector_get (magnetizations, i);
	}
	simple_mean_magnetization/=amount_conf;
	/**
	 * @note	Print out the raw data
	 * 			Calculate and print out its correlation
	 */
	printf ("simple_mean_magnetization:%f\n",simple_mean_magnetization);
	autocorrelation(magnetizations, magnetizationcorrelation, simple_mean_magnetization);
	gsl_vector_fprintf(magnetization_nmd4, magnetizations, "%f");
	fprintf (magnetizationcorrelation_nmd4, "#1\n");
	gsl_vector_fprintf(magnetizationcorrelation_nmd4, magnetizationcorrelation, "%f");
	
	/**
	 * @note	analyse measured data for magnetization and energy:
	 * 			binning to reduce correltations
	 * 			bootstrapping to get a good error estimate
	 * 	
	 * @note	write binned data in vector instead of block, so size can be changed with subvector-view. 
	 * 			block was tried, but led to problems with memory in reallocing.
	 */ 
	for (int lengthofbin=2; lengthofbin<65; lengthofbin*=2){
		fprintf (magnetizationcorrelation_nmd4, "\n#%d\n",lengthofbin);
		binnedmagnetization=gsl_vector_subvector(binnedmagnetization_mem, 0, magnetizations->size/lengthofbin);
		binning(magnetizations, &binnedmagnetization.vector, lengthofbin);
		
		/**
		 * @note	Analyse blocked data
		 */
		var_magnetization_binned=0;
		for(int i=0;i<(&binnedmagnetization.vector)->size;i++){
			simple_mean_magnetization_binned+=gsl_vector_get (&binnedmagnetization.vector, i);
		}
		simple_mean_magnetization_binned/=(&binnedmagnetization.vector)->size;
		for(int i=0;i<(&binnedmagnetization.vector)->size;i++){
			var_magnetization_binned+=(simple_mean_magnetization_binned-gsl_vector_get (&binnedmagnetization.vector, i))*(simple_mean_magnetization_binned-gsl_vector_get (&binnedmagnetization.vector, i));
		}
		var_magnetization_binned/=((&binnedmagnetization.vector)->size-1)*(magnetizations->size/lengthofbin);

		autocorrelation (&binnedmagnetization.vector, binnedcorrelation, simple_mean_magnetization_binned);
		gsl_vector_fprintf(magnetizationcorrelation_nmd4, binnedcorrelation, "%f");
		bootstrap(&binnedmagnetization.vector, generator, N_bs, &mean_magnetization, &var_magnetization);
		
		fprintf (raw_data,"%d\t%f\t%f\t%f\t%f\t%f\t%d\n", N, J_T, mean_magnetization,sqrt (var_magnetization), simple_mean_magnetization_binned, sqrt (var_magnetization_binned), lengthofbin);
	}

	
	/**
	 * @note	Print out the overall acceptance rate
	 */
	accept_rate/=amount_ar;
	printf ("Acceptance Rate:%f\n",accept_rate);
	
	
	
	
	/**
	 * @note	Cleanup
	 */
	fclose (raw_data);
	fclose (magnetization_nmd4);
	fclose (magnetizationcorrelation_nmd4);
	gsl_rng_free (generator);
	gsl_vector_free (set_of_phis);
	gsl_vector_free (magnetizations);
	gsl_vector_free(binnedmagnetization_mem);
	gsl_vector_free(magnetizationcorrelation);
	return 0;
}
