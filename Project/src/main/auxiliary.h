//header-Datei zu den auxiliary functions
#pragma once

#define GSL_COMPLEX_NAN (gsl_complex_rect(GSL_NAN,0.0))

/**
 * @brief returns trace of complex matrix
 * */
gsl_complex trace(gsl_matrix_complex *matrix);

/**
 * @brief returns det of complex matrix
 * */
gsl_complex det(gsl_matrix_complex *matrix);

/** 
 * @brief returns the determinant of a complex 2x2 matrix {{a,b},{c,d}} by returning ad-bc */
gsl_complex det2(gsl_matrix_complex *matrix);

/**
 * @brief transforms square matrix into its conjugate transposed. 
 * @note loops over i, conjugates diagonal elements, saves ij, sets ij=(ji)*, sets ji=(ij, old)*
 * */
void conjugatetranspose(gsl_matrix_complex *matrix);

/** sets all matrix elements to zero **/
void settozero(gsl_matrix_complex *matrix);

/** sets matrix to unity matrix **/
void settounity(gsl_matrix_complex *matrix);

/**
 * @brief
 * @note
 */
void crossproduct(gsl_vector_complex* a,gsl_vector_complex* b,gsl_vector_complex* result);

/**
 * @brief 	takes data points in measurements and writes bins over lengthofbin in binneddata
 */
void binning(gsl_vector *measurements, gsl_vector *binneddata, int lengthofbin);

/**
 * @brief makes one bootstrapreplica out of the data in measurements
 */
double makebootstrapreplica(gsl_vector * measurements, gsl_rng * generator);

/**
 * @brief uses the bootstrapmethod to calculate mean and variance of the data in measurements
 */ 
void bootstrap(gsl_vector *measurements, gsl_rng *generator, int R, double *mean, double *variance);


void autocorrelation(gsl_vector *measurements, gsl_vector *results, double mean);
