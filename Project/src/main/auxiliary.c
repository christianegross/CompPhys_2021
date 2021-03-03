//Project in Computational Physics: Yang-Mills/QQbar
//Nico Dichter, Christiane Groß
//02. February 2021-
//auxiliary functions so mainfile is not as big



#include <stdio.h>
#include <stdlib.h>
#include "math.h"						//exp-Function etc
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>		//needed for LU-decomp
#include <gsl/gsl_rng.h>				//random number generator
#include <gsl/gsl_randist.h>			//pull random number from gaussian


#define GSL_COMPLEX_NAN (gsl_complex_rect(GSL_NAN,0.0))

/**
 * @brief returns trace of complex matrix
 * */
gsl_complex trace(gsl_matrix_complex *matrix){
	if (matrix->size1!=matrix->size2){fprintf(stderr, "argument has to be square matrix!\n");return GSL_COMPLEX_NAN;}
	gsl_complex number=GSL_COMPLEX_ZERO;
	for (int i=0; i<matrix->size1; i+=1){
		number=gsl_complex_add(number, gsl_matrix_complex_get(matrix, i, i));
	}
	return number;
}


/**
 * @brief returns det of complex matrix
 * */
gsl_complex det(gsl_matrix_complex *matrix){
	if (matrix->size1!=matrix->size2||matrix->size1!=3){fprintf(stderr, "argument has to be 3x3 matrix!\n");return GSL_COMPLEX_NAN;}
	gsl_complex number=GSL_COMPLEX_ZERO;
	number=gsl_complex_add(number,gsl_complex_mul(gsl_complex_mul(gsl_matrix_complex_get(matrix,0,0),gsl_matrix_complex_get(matrix,1,1)),gsl_matrix_complex_get(matrix,2,2)));
	number=gsl_complex_add(number,gsl_complex_mul(gsl_complex_mul(gsl_matrix_complex_get(matrix,0,1),gsl_matrix_complex_get(matrix,1,2)),gsl_matrix_complex_get(matrix,2,0)));
	number=gsl_complex_add(number,gsl_complex_mul(gsl_complex_mul(gsl_matrix_complex_get(matrix,0,2),gsl_matrix_complex_get(matrix,1,0)),gsl_matrix_complex_get(matrix,2,1)));

	number=gsl_complex_sub(number,gsl_complex_mul(gsl_complex_mul(gsl_matrix_complex_get(matrix,0,2),gsl_matrix_complex_get(matrix,1,1)),gsl_matrix_complex_get(matrix,2,0)));
	number=gsl_complex_sub(number,gsl_complex_mul(gsl_complex_mul(gsl_matrix_complex_get(matrix,0,1),gsl_matrix_complex_get(matrix,1,0)),gsl_matrix_complex_get(matrix,2,2)));
	number=gsl_complex_sub(number,gsl_complex_mul(gsl_complex_mul(gsl_matrix_complex_get(matrix,0,0),gsl_matrix_complex_get(matrix,1,2)),gsl_matrix_complex_get(matrix,2,1)));
	return number;
}

/** 
 * @brief returns the determinant of a complex 2x2 matrix {{a,b},{c,d}} by returning ad-bc */
extern inline gsl_complex det2(gsl_matrix_complex *matrix){
	if (matrix->size1!=matrix->size2||matrix->size1!=2){fprintf(stderr, "argument has to be 2x2 matrix!\n");return GSL_COMPLEX_NAN;}
	return gsl_complex_sub(gsl_complex_mul(gsl_matrix_complex_get(matrix, 0,0), gsl_matrix_complex_get(matrix, 1,1)), gsl_complex_mul(gsl_matrix_complex_get(matrix, 1,0), gsl_matrix_complex_get(matrix, 0,1)));
	}

/**
 * @brief transforms square matrix into its conjugate transposed. 
 * @note loops over i, conjugates diagonal elements, saves ij, sets ij=(ji)*, sets ji=(ij, old)*
 * */
void conjugatetranspose(gsl_matrix_complex *matrix){
	if (matrix->size1!=matrix->size2){fprintf(stderr, "argument has to be square matrix!\n");return;}
	gsl_complex number=GSL_COMPLEX_ZERO;
	for (int i=0;i<matrix->size1;i+=1){
		gsl_matrix_complex_set(matrix, i,i, gsl_complex_conjugate(gsl_matrix_complex_get(matrix, i,i)));
		for(int j=i+1;j<matrix->size1;j+=1){
			number=gsl_matrix_complex_get(matrix, i, j);
			gsl_matrix_complex_set(matrix, i,j, gsl_complex_conjugate(gsl_matrix_complex_get(matrix, j,i)));
			gsl_matrix_complex_set(matrix, j,i, gsl_complex_conjugate(number));
		}
	}
}

/** sets all matrix elements to zero **/
void settozero(gsl_matrix_complex *matrix){
	if (matrix->size1!=matrix->size2){fprintf(stderr, "argument has to be square matrix!\n");return;}
	for (int i=0;i<matrix->size1;i+=1){
		gsl_matrix_complex_set(matrix, i,i, GSL_COMPLEX_ZERO);
		for(int j=i+1;j<matrix->size1;j+=1){
			gsl_matrix_complex_set(matrix, i,j, GSL_COMPLEX_ZERO);
			gsl_matrix_complex_set(matrix, j,i, GSL_COMPLEX_ZERO);
	
		}
	}
}

/** sets matrix to unity matrix **/
void settounity(gsl_matrix_complex *matrix){
	if (matrix->size1!=matrix->size2){fprintf(stderr, "argument has to be square matrix!\n");return;}
	for (int i=0;i<matrix->size1;i+=1){
		gsl_matrix_complex_set(matrix, i,i, GSL_COMPLEX_ONE);
		for(int j=i+1;j<matrix->size1;j+=1){
			gsl_matrix_complex_set(matrix, i,j, GSL_COMPLEX_ZERO);
			gsl_matrix_complex_set(matrix, j,i, GSL_COMPLEX_ZERO);
		}
	}
}

/**
 * @brief
 * @note
 */
void crossproduct(gsl_vector_complex* a,gsl_vector_complex* b,gsl_vector_complex* result){
	gsl_vector_complex_set (result, 0, gsl_complex_conjugate(gsl_complex_sub(gsl_complex_mul (gsl_vector_complex_get(a,1),  (gsl_vector_complex_get(b,2))),gsl_complex_mul (gsl_vector_complex_get(a,2),  (gsl_vector_complex_get(b,1))))));
	gsl_vector_complex_set (result, 1, gsl_complex_conjugate(gsl_complex_sub(gsl_complex_mul (gsl_vector_complex_get(a,2),  (gsl_vector_complex_get(b,0))),gsl_complex_mul (gsl_vector_complex_get(a,0),  (gsl_vector_complex_get(b,2))))));
	gsl_vector_complex_set (result, 2, gsl_complex_conjugate(gsl_complex_sub(gsl_complex_mul (gsl_vector_complex_get(a,0),  (gsl_vector_complex_get(b,1))),gsl_complex_mul (gsl_vector_complex_get(a,1),  (gsl_vector_complex_get(b,0))))));
}


/**
 * @fn 		inline void binning(gsl_block *measurements, gsl_block *binneddata, int lengthofbin)
 * @brief 	takes data points in measurements and writes bins over lengthofbin in binneddata
 */
extern inline void binning(gsl_vector *measurements, gsl_vector *binneddata, int lengthofbin){
	double bincontent=0;
	//assertion if lengths of bins, measurements fit together
	int numberofbins=measurements->size/lengthofbin;
	//printf("numberofbins=%d\tnumberofelemtns=%lu\tlengthofbin=%d\n", numberofbins, binneddata->size, lengthofbin);
	if (numberofbins!=binneddata->size){fprintf(stderr, "Gave wrong lengths! %d number of bins to be calculated, but storage allocated for %lu!", numberofbins, binneddata->size); return;}
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

