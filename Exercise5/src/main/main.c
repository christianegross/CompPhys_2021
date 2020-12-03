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
 * @note implements fine-to-coarse restricion given ineq. 6 on the sheet
 */
void finetocoarserestriction(gsl_vector *ufine, gsl_vector *ucoarse){
	for (int i=0; i<ucoarse->size-1; i+=1){
		gsl_vector_set(ucoarse, i, gsl_vector_get(ufine,2*i));
	}
	gsl_vector_set(ucoarse,ucoarse->size-1, 0);
}

/**
 * @note adds coarse-to-fine interpolaton to current u
 * */
 void addcoarsetofineinterpolation(gsl_vector *ufine, gsl_vector *ucoarse){
	 for (int i=0; i<ucoarse->size-1; i+=1){
		 gsl_vector_set(ufine,2*i, gsl_vector_get(ufine,2*i)+gsl_vector_get(ucoarse,i));
		 gsl_vector_set(ufine,2*i+1, gsl_vector_get(ufine,2*i+1)+0.5*(gsl_vector_get(ucoarse,i)+gsl_vector_get(ucoarse,i+1)));
	 }
	 gsl_vector_set(ufine,ufine->size-1, 0);
}
	
double hamiltonian(gsl_vector *u, gsl_vector *phi, double a){
	double H=0;
	for (int i=1; i<u->size-1; i+=1){
		H+=(gsl_vector_get(u, i)-gsl_vector_get(u, i-1))*(gsl_vector_get(u, i)-gsl_vector_get(u, i-1))+a*a*gsl_vector_get(phi, i)*gsl_vector_get(u, i);
	}
	H+=(gsl_vector_get(u, u->size-1)-gsl_vector_get(u, u->size-2))*(gsl_vector_get(u, u->size-1)-gsl_vector_get(u, u->size-2));
	return H/a;
}

void interpolatephi(gsl_vector* phifine, gsl_vector *ufine, gsl_vector* phicoarse, double a){
	gsl_vector_set(phicoarse, 0, 0);
	//insert phicoarse(phifine, ufine)
	for (int i=1; i<phicoarse->size-1;i+=1){
		//gsl_vector_set(phicoarse, i,1);
		gsl_vector_set(phicoarse, i, 
		1/a/a*(gsl_vector_get(ufine, 2*i)-gsl_vector_get(ufine, 2*i-1)-gsl_vector_get(ufine, 2*i+2)+gsl_vector_get(ufine, 2*i+1))
		+gsl_vector_get(phifine, 2*i)+0.5*(gsl_vector_get(phifine, 2*i+1)+gsl_vector_get(phifine, 2*i-1)));
	}
	gsl_vector_set(phicoarse, phicoarse->size-1, 0);
}

void onesweep(gsl_vector *u, gsl_vector *phi, double delta, double a, gsl_rng *generator){
	int randomsite, randomstep;
	double deltah;
	for (int l=1; l<u->size-1; l+=1){
		randomsite=gsl_rng_uniform_int(generator, u->size-2)+1; //ensures only values between 1 and N-1 are selected
		randomstep=gsl_ran_flat(generator, -delta, delta);
		deltah=-2*a*gsl_vector_get(phi, randomsite)*randomstep+2/a*randomstep*(-randomstep+2*gsl_vector_get(u, randomsite)+gsl_vector_get(u, randomsite+1)+gsl_vector_get(u, randomsite-1));
		if (gsl_rng_uniform(generator)>exp(deltah)){
			gsl_vector_set(u, randomsite, gsl_vector_get(u, randomsite)+randomstep);
		}
	}
}

double magnetisation(gsl_vector *u){
	double magnet=0;
	for (int i=1; i<u->size-1; i+=1){
		magnet+=gsl_vector_get(u, i);
	}
	return magnet/(u->size-1);
}

inline int nu(int level){
	return pow(2, level-1);
}

void multigrid(gsl_vector* u, gsl_vector *phi, gsl_rng *generator, double a, double delta, int level, int levelmax, int gamma){
	//printf("level=%d\n", level);
	for (int i=0; i<nu(level); i+=1){
		onesweep(u, phi, delta, a, generator);
	}
	//printf("done presweep\n");
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
	for (int i=0; i<nu(level); i+=1){
		onesweep(u, phi, delta, a, generator);
	}
	
	//printf("done postsweep\n");
}
	

int main(int argc, char **argv){	
	//set and allocate random number generator
	int seed=2;//use fixed seed: result should be exactly reproduced using the same seed
	gsl_rng *generator;
	generator=gsl_rng_alloc(gsl_rng_mt19937);//use mersenne-twister
	gsl_rng_set(generator, seed);
	
	gsl_vector *u=gsl_vector_calloc(65);
	gsl_vector *uold=gsl_vector_calloc(65);
	gsl_vector *phi=gsl_vector_calloc(65);
	for (int i=1; i<u->size-1; i+=1){
		gsl_vector_set(u, i, gsl_rng_uniform(generator));
	}
	double magnet, hamilton, a=1.0/64, deltah;
	for (int i=0;i<10000;i+=1){
		onesweep(u, phi, 2, a, generator);
		//~ gsl_vector_memcpy(uold, u);
		//~ multigrid(u, phi, generator, a, 2, 1, 1, 1);
		//~ deltah=hamiltonian(uold, phi, a)-hamiltonian(u, phi, a);
		//~ if (gsl_rng_uniform(generator)>exp(-deltah)){
			//~ gsl_vector_memcpy(u, uold);
			//~ //printf("not accepted\t");
		//~ }
		magnet=magnetisation(u);
		hamilton=hamiltonian(u, phi, a);
		printf("%d\t%e\t%e\t%e\n", i, magnet, hamilton, deltah);
	}
	for (int i=0; i<u->size; i+=1){
	//	printf("%e\t", gsl_vector_get(u, i));
	}
	printf("\n");
	gsl_rng_free(generator);
	gsl_vector_free(u);
	gsl_vector_free(phi);
	return 0;
}


