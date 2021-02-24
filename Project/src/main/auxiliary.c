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
	gsl_vector_complex_set (result, 0, gsl_complex_sub(gsl_complex_mul (gsl_vector_complex_get(a,1), gsl_vector_complex_get(b,2)),gsl_complex_mul (gsl_vector_complex_get(a,2), gsl_vector_complex_get(b,1))));
	gsl_vector_complex_set (result, 1, gsl_complex_sub(gsl_complex_mul (gsl_vector_complex_get(a,2), gsl_vector_complex_get(b,0)),gsl_complex_mul (gsl_vector_complex_get(a,0), gsl_vector_complex_get(b,2))));
	gsl_vector_complex_set (result, 2, gsl_complex_sub(gsl_complex_mul (gsl_vector_complex_get(a,0), gsl_vector_complex_get(b,1)),gsl_complex_mul (gsl_vector_complex_get(a,1), gsl_vector_complex_get(b,0))));
}
