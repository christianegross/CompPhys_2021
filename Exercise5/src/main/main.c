//Exercise 5 for Computational Physics: Simulation of the Gaussian Model with Multigrid Monte Carlo
//Nico Dichter, Christiane Groß
//02. December 2020-

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_block.h>
#include <gsl/gsl_rng.h>//random number generator
#include <gsl/gsl_randist.h>//pull random number from gaussian
#include "math.h"//exp-Function
#include <gsl/gsl_vector.h>



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
void binning(gsl_vector *measurements, gsl_vector *binneddata, int lengthofbin){
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
		

void autocorrelation(gsl_vector *measurements, gsl_vector *results, double mean){
	/**
	 * @note implements normalized autocorrelation by first calculation c(0) with the definition from the lectures, then calculating c(tau) and then storing c(tau)/c(0) in results
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

/**
 * @brief implements fine-to-coarse restricion given in eq. 6 on the sheet
 * 
 * @param ufine field wih smaller distances between sites, used to populate coarser field
 * @param ucoarse field with larger distances
 */
void finetocoarserestriction(gsl_vector *ufine, gsl_vector *ucoarse){
	for (int i=0; i<ucoarse->size-1; i+=1){
		gsl_vector_set(ucoarse, i, gsl_vector_get(ufine,2*i));
	}
	gsl_vector_set(ucoarse,ucoarse->size-1, 0);
}

/**
 * @brief adds coarse-to-fine interpolaton to current u
 * 
 * @param ufine field wih smaller distances between sites
 * @param ucoarse field with larger distances, used to update finer field
 * */
 void addcoarsetofineinterpolation(gsl_vector *ufine, gsl_vector *ucoarse){
	 for (int i=0; i<ucoarse->size-1; i+=1){
		 gsl_vector_set(ufine,2*i, gsl_vector_get(ufine,2*i)+gsl_vector_get(ucoarse,i));
		 gsl_vector_set(ufine,2*i+1, gsl_vector_get(ufine,2*i+1)+0.5*(gsl_vector_get(ucoarse,i)+gsl_vector_get(ucoarse,i+1)));
	 }
	 gsl_vector_set(ufine,ufine->size-1, 0);
}


/**
 * @brief implements Hamiltonian as given on the sheet
 * 
 * @param u field to be analysed
 * @param phi external field
 * @param a lattice spacing
 * */	
double hamiltonian(gsl_vector *u, gsl_vector *phi, double a){
	double H=0;
	for (int i=1; i<u->size-1; i+=1){
		H+=(gsl_vector_get(u, i)-gsl_vector_get(u, i-1))*(gsl_vector_get(u, i)-gsl_vector_get(u, i-1))+a*a*gsl_vector_get(phi, i)*gsl_vector_get(u, i);
	}
	H+=(gsl_vector_get(u, u->size-1)-gsl_vector_get(u, u->size-2))*(gsl_vector_get(u, u->size-1)-gsl_vector_get(u, u->size-2));
	return H/a;
}

/**
 * @brief fills phi of coarser field with the values of the finer field, relation determined in the report
 * 
 * @param phifine external field at the finer level
 * @param phicoarse external field at the coarser level
 * @param ufine field at the finer level
 * @param a lattice spacing
 * */
void interpolatephi(gsl_vector* phifine, gsl_vector *ufine, gsl_vector* phicoarse, double a){
	gsl_vector_set(phicoarse, 0, 0);
	//insert phicoarse(phifine, ufine)
	for (int i=1; i<phicoarse->size-1;i+=1){
		//gsl_vector_set(phicoarse, i,1);
		gsl_vector_set(phicoarse, i, 
		1/a/a*(gsl_vector_get(ufine, 2*i)-gsl_vector_get(ufine, 2*i-1)-gsl_vector_get(ufine, 2*i+2)+gsl_vector_get(ufine, 2*i+1))
		+0.5*gsl_vector_get(phifine, 2*i)+0.25*(gsl_vector_get(phifine, 2*i+1)+gsl_vector_get(phifine, 2*i-1)));
	}
	gsl_vector_set(phicoarse, phicoarse->size-1, 0);
}

/**
 * @brief does one sweep of the lattice as explained on the sheet: select random site, propose change, metropolis a/r step as often as there are free lattice points
 * 
 * @param u simulated lattice
 * @param phi external field
 * @param delta maximum step at proposed change
 * @param a lattice spacing
 * @param generator random numbers for selecting sit, proposing change, a/r
 * */
void onesweep(gsl_vector *u, gsl_vector *phi, double delta, double a, gsl_rng *generator){
	int randomsite, randomstep;
	double deltah;
	for (int l=1; l<u->size-1; l+=1){
		randomsite=gsl_rng_uniform_int(generator, u->size-2)+1; //ensures only values between 1 and N-1 are selected
		randomstep=gsl_ran_flat(generator, -delta, delta);
		deltah=a*gsl_vector_get(phi, randomsite)*randomstep+2/a*randomstep*(randomstep+2*gsl_vector_get(u, randomsite)-gsl_vector_get(u, randomsite+1)-gsl_vector_get(u, randomsite-1));
		if (gsl_rng_uniform(generator)<exp(-deltah)){
			gsl_vector_set(u, randomsite, gsl_vector_get(u, randomsite)+randomstep);
		}
	}
}
/**
 * @brief returns magnetisation of given lattice according to the definition on the sheet
 * */
double magnetisation(gsl_vector *u){
	double magnet=0;
	for (int i=1; i<u->size-1; i+=1){
		magnet+=gsl_vector_get(u, i);
	}
	return magnet/(u->size-1);
}

/**
 * @brief returns number of pre and postsweeps to be done per level
 * */
inline int nu(int level){
	return pow(2, level-1);
}

/**
 * @brief implementation of the multigrid algorithm as given on the sheet: presweeps, recursion, postsweeps
 * 
 * @param u simulated lattice
 * @param phi external field
 * @param generator for random numbers for sweep
 * @param a lattice spacing
 * @param delta maximum step in sweep
 * @param level current level of recursion
 * @param levelmax maximum number of recursions
 * @param gamma how often should multigrid algorithm be done in recursion step
 * */
void multigrid(gsl_vector* u, gsl_vector *phi, gsl_rng *generator, double a, double delta, int level, int levelmax, int gamma){
	//presweeps
	for (int i=0; i<nu(level); i+=1){
		onesweep(u, phi, delta, a, generator);
	}
	//recursive steps: assign coarser field, do multigrid, evaluate finer fields
	if (level<levelmax){
		gsl_vector *ucoarse=gsl_vector_alloc((u->size-1)/2+1);
		gsl_vector *phicoarse=gsl_vector_alloc((phi->size-1)/2+1);
		finetocoarserestriction(u, ucoarse);
		interpolatephi(phi, u, phicoarse, a);
		for (int i=0; i<gamma; i+=1){
			multigrid(ucoarse, phicoarse, generator, 2*a, delta, level+1, levelmax, gamma);
		}
		addcoarsetofineinterpolation(u, ucoarse);
		gsl_vector_free(ucoarse);
		gsl_vector_free(phicoarse);
	}
	//postsweeps
	for (int i=0; i<nu(level); i+=1){
		onesweep(u, phi, delta, a, generator);
	}
}
	

int main(int argc, char **argv){
	/**
	 * Declarations:
	 *
	 * @var delta
	 * @var N_meas	#measurements
	 * 
	 */
	double delta=2.;
	int N=64;
	double magnet, magnet_mean=0, magnet_var=0,  magnet_sqr_mean=0, hamilton, hamilton_mean=0, hamilton_var=0; 
	double a=1.0;
	int startlevel=1;
	int maxlevel=1;//TODO test maxlevel=3, yields inf atm
	int gamma=1;
	int N_meas=1000;
	int N_therm=100;
	int lengthofbin=10;
	int N_bs=4*N_meas;
	
	//set and allocate random number generator
	int seed=2;//use fixed seed: result should be exactly reproduced using the same seed
	gsl_rng *generator;
	generator=gsl_rng_alloc(gsl_rng_mt19937);//use mersenne-twister
	gsl_rng_set(generator, seed);
	
	/**
	 * @note	Allocate memory for u and phi and initialize u
	 */
	gsl_vector *u=gsl_vector_calloc(N+1);
	gsl_vector *phi=gsl_vector_calloc(N+1);
	for (int i=1; i<u->size-1; i+=1){
		gsl_vector_set(u, i, gsl_rng_uniform(generator)+5);
	}
	gsl_vector_set(u, 0, 0);
	gsl_vector_set(u, u->size-1, 0);
	
	/**
	 * @note	Open output streams and allocate memory for measurements
	 */
	FILE * meas_stream=fopen ("data/meas.dat", "w");
	FILE * vcycle_stream=fopen ("data/vcycle.dat", "w");
	FILE * wcycle_stream=fopen ("data/wcycle.dat", "w");
	gsl_vector *magnet_vec=gsl_vector_alloc(N_meas);
	gsl_vector *magnet_vec_sqr=gsl_vector_alloc(N_meas);
	gsl_vector *hamilton_vec=gsl_vector_alloc(N_meas);
	gsl_vector *correlation_magn_sqr=gsl_vector_alloc(N_meas/2);
	gsl_vector * binnedmagnetization_mem=gsl_vector_alloc(N_meas);
	gsl_vector_view binnedmagnetization;
	gsl_vector * binnedhamilton_mem=gsl_vector_alloc(N_meas);
	gsl_vector_view binnedhamilton;
	
	
	
	
	
	
	/**
	 * @note	Standard Metropolis-Hastings:
	 * 			Thermalization and measurement
	 */
	for(int i=0;i<N_therm;i++){
		onesweep(u, phi, delta, a, generator);
	}
	for(int i=0;i<N_meas;i++){
		onesweep(u, phi, delta, a, generator);
		magnet=magnetisation(u);
		hamilton=hamiltonian(u, phi, a);
		gsl_vector_set (magnet_vec, i, magnet);
		gsl_vector_set (hamilton_vec, i, hamilton);
	}
	
	/**
	 * @note	Binning and bootstrap
	 */
	binnedmagnetization=gsl_vector_subvector(binnedmagnetization_mem, 0, magnet_vec->size/lengthofbin);
	binnedhamilton=gsl_vector_subvector(binnedhamilton_mem, 0, hamilton_vec->size/lengthofbin);
	binning(magnet_vec, &binnedmagnetization.vector, lengthofbin);
	binning(hamilton_vec, &binnedhamilton.vector, lengthofbin);
	bootstrap(&binnedmagnetization.vector, generator, N_bs, &magnet_mean, &magnet_var);
	bootstrap(&binnedhamilton.vector, generator, N_bs, &hamilton_mean, &hamilton_var);
	fprintf (meas_stream, "#simplehaste\n%e\t%e\t%e\t%e\n",magnet_mean,magnet_var,hamilton_mean,hamilton_var);
	
	/**
	 * @note	Reset u
	 */
	for (int i=1; i<u->size-1; i+=1){
		gsl_vector_set(u, i, gsl_rng_uniform(generator)+5);
	}
	gsl_vector_set(u, 0, 0);
	gsl_vector_set(u, u->size-1, 0);
	magnet_mean=0;
	hamilton_mean=0;
	magnet_sqr_mean=0;
	magnet_var=0;
	hamilton_var=0;
	
	
	
	
	
	
	
	/**
	 * @note	v-cycle
	 * 			Thermalization and measurement
	 */
	gamma=1;
	for(int i=0;i<N_therm;i++){
		multigrid(u, phi, generator, a, delta, startlevel, maxlevel, gamma);
	}
	for (int i=0;i<N_meas;i+=1){
		multigrid(u, phi, generator, a, delta, startlevel, maxlevel, gamma);
		magnet=magnetisation(u);
		hamilton=hamiltonian(u, phi, a);
		gsl_vector_set (magnet_vec, i, magnet);
		gsl_vector_set (magnet_vec_sqr, i, magnet*magnet);
		magnet_sqr_mean+=gsl_vector_get (magnet_vec_sqr, i);
		gsl_vector_set (hamilton_vec, i, hamilton);
	}
	magnet_sqr_mean/=N_meas;
	
	/**
	 * @note	Binning, bootstrap and autocorrelation
	 */
	binnedmagnetization=gsl_vector_subvector(binnedmagnetization_mem, 0, magnet_vec->size/lengthofbin);
	binnedhamilton=gsl_vector_subvector(binnedhamilton_mem, 0, hamilton_vec->size/lengthofbin);
	binning(magnet_vec, &binnedmagnetization.vector, lengthofbin);
	binning(hamilton_vec, &binnedhamilton.vector, lengthofbin);
	bootstrap(&binnedmagnetization.vector, generator, N_bs, &magnet_mean, &magnet_var);
	bootstrap(&binnedhamilton.vector, generator, N_bs, &hamilton_mean, &hamilton_var);
	autocorrelation (magnet_vec_sqr, correlation_magn_sqr, magnet_sqr_mean);
	
	fprintf (meas_stream, "#v-cycle\n%e\t%e\t%e\t%e\n",magnet_mean,magnet_var,hamilton_mean,hamilton_var);
	gsl_vector_fprintf (vcycle_stream, correlation_magn_sqr, "%f");
	
	/**
	 * @note	Reset u
	 */
	for (int i=1; i<u->size-1; i+=1){
		gsl_vector_set(u, i, gsl_rng_uniform(generator)+5);
	}
	gsl_vector_set(u, 0, 0);
	gsl_vector_set(u, u->size-1, 0);
	magnet_mean=0;
	hamilton_mean=0;
	magnet_sqr_mean=0;
	magnet_var=0;
	hamilton_var=0;
	
	
	
	
	
	
	
	
	/**
	 * @note	w-cycle
	 * 			Thermalization and measurement
	 */
	gamma=2;
	for(int i=0;i<N_therm;i++){
		multigrid(u, phi, generator, a, delta, startlevel, maxlevel, gamma);
	}
	for (int i=0;i<N_meas;i+=1){
		multigrid(u, phi, generator, a, delta, startlevel, maxlevel, gamma);
		magnet=magnetisation(u);
		hamilton=hamiltonian(u, phi, a);
		gsl_vector_set (magnet_vec, i, magnet);
		gsl_vector_set (magnet_vec_sqr, i, magnet*magnet);
		magnet_sqr_mean+=gsl_vector_get (magnet_vec_sqr, i);
		gsl_vector_set (hamilton_vec, i, hamilton);
	}
	magnet_sqr_mean/=N_meas;
	
	/**
	 * @note	Binning, bootstrap and autocorrelation
	 */
	binnedmagnetization=gsl_vector_subvector(binnedmagnetization_mem, 0, magnet_vec->size/lengthofbin);
	binnedhamilton=gsl_vector_subvector(binnedhamilton_mem, 0, hamilton_vec->size/lengthofbin);
	binning(magnet_vec, &binnedmagnetization.vector, lengthofbin);
	binning(hamilton_vec, &binnedhamilton.vector, lengthofbin);
	bootstrap(&binnedmagnetization.vector, generator, N_bs, &magnet_mean, &magnet_var);
	bootstrap(&binnedhamilton.vector, generator, N_bs, &hamilton_mean, &hamilton_var);
	autocorrelation (magnet_vec_sqr, correlation_magn_sqr, magnet_sqr_mean);
	
	fprintf (meas_stream, "#w-cycle\n%e\t%e\t%e\t%e\n",magnet_mean,magnet_var,hamilton_mean,hamilton_var);
	gsl_vector_fprintf (wcycle_stream, correlation_magn_sqr, "%f");
	
	/**
	 * @note	Cleanup
	 */
	fclose (meas_stream);
	fclose (vcycle_stream);
	fclose (wcycle_stream);
	gsl_rng_free(generator);
	gsl_vector_free(hamilton_vec);
	gsl_vector_free(magnet_vec);
	gsl_vector_free(u);
	gsl_vector_free(phi);
	return 0;
}


