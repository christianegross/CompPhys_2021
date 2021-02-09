//Project in Computational Physics: Yang-Mills/QQbar
//Nico Dichter, Christiane Groß
//02. February 2021-

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
 * @brief transforms matrix into its conjugate transposed. 
 * @note loops over i, conjugates diagonal elements, saves ij, sets ij=(ji)*, sets ji=(ij, old)*
 * */
void conjugatetranspose(gsl_matrix_complex *matrix){
	if (matrix->size1!=matrix->size2){fprintf(stderr, "argument has to be square matrix!\n");return GSL_NAN;}
	gsl_complex number=GSL_COMPLEX_ZERO;
	for (int i=0;i<matrix->size1;i+=1){
		gsl_matrix_complex_set(matrix, i,i, gsl_complex_conjugate(gsl_matrix_complex_get(matrix, i,i)));
		for(int j=i+1;j<matrix->size1;j+=1){
			number=gsl_matrix_complex_get(matrix, i, j);
			gsl_matrix_complex_set(matrix, i,j, gsl_complex_conjugate(gsl_matrix_complex_get(matrix, j,i)));
			gsl_matrix_complex_set(matrix, j,i, gsl_complex_conjugate(number))
			;
		}
	}
}

void generatesu2(gsl_matrix_complex * matrix, double epsilon, gsl_rng * generator){
	//error handling
	if(matrix->size1!=2||matrix->size2!=2){fprintf(stderr, "argument has to be 2x2 matrix!\n");return GSL_NAN;}
	//generate random numbers
	double e,f,g, det;
	gsl_vector_complex_view columnzero;
	e=2.0*gsl_rng_uniform(generator)-1.0;
	f=2.0*gsl_rng_uniform(generator)-1.0;
	g=2.0*gsl_rng_uniform(generator)-1.0;
	//set matrix up to be hermitian
	gsl_matrix_complex_set(matrix, 0, 0, gsl_complex_rect(1, epsilon*e));
	gsl_matrix_complex_set(matrix, 1, 0, gsl_complex_rect(epsilon*g, epsilon*f));
	//normalize vector zero to make sure determinant is one
	//rescale by 1/norm
	columnzero=gsl_matrix_complex_column(matrix, 0);
	det=gsl_blas_dznrm2(&columnzero.vector);		
	gsl_blas_zdscal(1.0/det, &columnzero.vector);
	//set up rest of matrix as 11=(00)*, 01=-(10)*
	gsl_matrix_complex_set(matrix, 1,1,gsl_complex_conjugate(gsl_matrix_complex_get(matrix, 0,0)));
	gsl_matrix_complex_set(matrix, 0,1,gsl_complex_mul_real(gsl_complex_conjugate(gsl_matrix_complex_get(matrix, 1,0)), -1.0));	
}



int main(int argc, char **argv){
	//~ double a,b,c,d,det;
	double epsilon=1;
	//~ double e,f,g,h;
	gsl_complex complexproduct=GSL_COMPLEX_ZERO;
	gsl_matrix_complex *one=gsl_matrix_complex_alloc(2,2);
	gsl_matrix_complex *two=gsl_matrix_complex_alloc(2,2);
	//~ gsl_matrix_complex *three=gsl_matrix_complex_alloc(2,2);
	//~ gsl_matrix_complex *four=gsl_matrix_complex_alloc(2,2);
	gsl_matrix_complex *multiplied=gsl_matrix_complex_alloc(2,2);	
	gsl_matrix_complex *five=gsl_matrix_complex_calloc(3,2);	
	//~ gsl_vector_complex_view columnzero, columnone;
	int seed=2;//use fixed seed: result should be exactly reproduced using the same seed
	gsl_rng *generator;
	generator=gsl_rng_alloc(gsl_rng_mt19937);//use mersenne-twister
	gsl_rng_set(generator, seed);
	
	//~ for(int i=0;i<1;i+=1){	
		 //~ /**
		  //~ * choose four random doubles between -1 and 1, normalise through determinant of future matrix, so determinant is 1
		  //~ * */
		//~ a=2.0*gsl_rng_uniform(generator)-1.0;
		//~ b=2.0*gsl_rng_uniform(generator)-1.0;
		//~ c=2.0*gsl_rng_uniform(generator)-1.0;
		//~ d=2.0*gsl_rng_uniform(generator)-1.0;
		//~ det=a*a+b*b+c*c+d*d;
		//~ a/=sqrt(det);
		//~ b/=sqrt(det);
		//~ c/=sqrt(det);
		//~ d/=sqrt(det);
		
		//~ /** set complex hermitian matrix*/
		//~ gsl_matrix_complex_set(one, 0, 0, gsl_complex_rect(a, b));
		//~ gsl_matrix_complex_set(one, 1, 1, gsl_complex_rect(a,-b));
		//~ gsl_matrix_complex_set(one, 0, 1, gsl_complex_rect(c, d));
		//~ gsl_matrix_complex_set(one, 1, 0, gsl_complex_rect(-c,d));
		
		//~ gsl_matrix_complex_set(two, 0, 0, gsl_complex_rect(a,-b));
		//~ gsl_matrix_complex_set(two, 1, 1, gsl_complex_rect(a, b));
		//~ gsl_matrix_complex_set(two, 0, 1, gsl_complex_rect(-c,-d));
		//~ gsl_matrix_complex_set(two, 1, 0, gsl_complex_rect(c,-d));
		//~ /** see if uudagger =1*/
		//~ gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(1,0), one, two, gsl_complex_rect(0,0), multiplied);
		
		//~ printf("run no. %d\n matrix one \n", i);
		//~ gsl_matrix_complex_fprintf(stdout, one, "%f");
		//~ printf("\n matrix two\n");
		//~ gsl_matrix_complex_fprintf(stdout, two, "%f");
		//~ printf("\n multiplied\n");
		//~ gsl_matrix_complex_fprintf(stdout, multiplied, "%f");
		//~ printf("\n old and new determinant\n");
		//~ printf("%f\t%f\n\n", det, a*a+b*b+c*c+d*d);
		
		//~ /**
		 //~ * Make hermitian Matrix 3 with entries between -1 and 1 {{e, f+ig},{f-ig, h}}
		 //~ * Fill Matrix four =1+iepsilon*three={{1+iepsilon e, -epsilon g+iepsilon f}, {epsilong+iepsilon f, 1+iepsilon h}}
		 //~ * */  
		
		//~ e=2.0*gsl_rng_uniform(generator)-1.0;
		//~ f=2.0*gsl_rng_uniform(generator)-1.0;
		//~ g=2.0*gsl_rng_uniform(generator)-1.0;
		//~ h=2.0*gsl_rng_uniform(generator)-1.0;
		
		//~ gsl_matrix_complex_set(three, 0, 0, gsl_complex_rect(e, 0));
		//~ gsl_matrix_complex_set(three, 1, 1, gsl_complex_rect(-e, 0));
		//~ gsl_matrix_complex_set(three, 0, 1, gsl_complex_rect(f, g));
		//~ gsl_matrix_complex_set(three, 1, 0, gsl_complex_rect(f,-g));
		
		//~ gsl_matrix_complex_set(four, 0, 0, gsl_complex_rect(1, epsilon*e));
		//~ gsl_matrix_complex_set(four, 1, 1, gsl_complex_rect(1, -epsilon*e));
		//~ gsl_matrix_complex_set(four, 0, 1, gsl_complex_rect(-epsilon*g, epsilon*f));
		//~ gsl_matrix_complex_set(four, 1, 0, gsl_complex_rect(epsilon*g, epsilon*f));
		
		
		//~ printf("\n four before normalization\n");
		//~ gsl_matrix_complex_fprintf(stdout, four, "%f");
		
		//~ columnzero=gsl_matrix_complex_column(four, 0);	//set up columns as vectors
		//~ columnone=gsl_matrix_complex_column(four, 1);
		
		
		//~ gsl_blas_zdotu(&columnzero.vector, &columnone.vector, &complexproduct); //check orthogonality
		//~ printf("\northogonal?\n%f+%f*i\n\n", GSL_REAL(complexproduct), GSL_IMAG(complexproduct));
		
		//~ det=gsl_blas_dznrm2(&columnzero.vector);		//normalize vector zero
		//~ gsl_blas_zdscal(1.0/det, &columnzero.vector);
		
		//~ gsl_blas_zdotu(&columnzero.vector, &columnone.vector, &complexproduct);		//do Gram-Schmidt orthonormalization as in Script from T. Raesch
		//~ gsl_blas_zaxpy(gsl_complex_mul_real(complexproduct, -1.0), &columnzero.vector, &columnone.vector); // should be one=-1*complexprodukt*zero+one
		
		//~ det=gsl_blas_dznrm2(&columnone.vector);		//normalize vector one
		//~ gsl_blas_zdscal(1.0/det, &columnone.vector);
		
		//~ //set 11=(00)*, 01=-(10)*
		//~ gsl_matrix_complex_set(four, 1,1,gsl_complex_conjugate(gsl_matrix_complex_get(four, 0,0)));
		//~ gsl_matrix_complex_set(four, 0,1,gsl_complex_mul_real(gsl_complex_conjugate(gsl_matrix_complex_get(four, 1,0)), -1.0));
		
		//~ printf("\n vectors normalized?\n");
		//~ printf("\n%f\n", gsl_blas_dznrm2(&columnzero.vector));		//check normalization and Matrix
		//~ printf("%f\n", gsl_blas_dznrm2(&columnone.vector));
		
		//~ printf("\n four after normalization\n");
		//~ gsl_matrix_complex_fprintf(stdout, four, "%f");
		
		//~ complexproduct=gsl_complex_sub(gsl_complex_mul(gsl_matrix_complex_get(four, 0,0), gsl_matrix_complex_get(four, 1,1)), gsl_complex_mul(gsl_matrix_complex_get(four, 1,0), gsl_matrix_complex_get(four, 0,1)));
		//~ printf("\ndeterminant\n%f+%f*i\t%f\n\n", GSL_REAL(complexproduct), GSL_IMAG(complexproduct), gsl_complex_abs(complexproduct));
		
		//~ gsl_blas_zdotu(&columnzero.vector, &columnone.vector, &complexproduct); //check determinant and orthogonality
		//~ printf("\northogonal\n%f+%f*i\t%f\n\n", GSL_REAL(complexproduct), GSL_IMAG(complexproduct), gsl_complex_abs(complexproduct));
		
		
		//~ gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, gsl_complex_rect(1,0), four, four, gsl_complex_rect(0,0), multiplied);
		//~ printf("\n multiplied\n");
		//~ gsl_matrix_complex_fprintf(stdout, multiplied, "%f");
		//~ complexproduct=trace(four);
		
		//~ printf("\ntrace\n%f+%f*i\t%f\n\n", GSL_REAL(complexproduct), GSL_IMAG(complexproduct), gsl_complex_abs(complexproduct));
	//~ }
	
	for(int i=0; i<100; i+=1){
		epsilon=0.01*i;
		generatesu2(one, epsilon, generator);
		gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, gsl_complex_rect(1,0), one, one, gsl_complex_rect(0,0), multiplied);
		printf("\n multiplied in run %d\n", i);
		//~ gsl_matrix_complex_fprintf(stdout, multiplied, "%f");
		gsl_matrix_complex_fprintf(stdout,one, "%f");
		complexproduct=trace(one);
		printf("\ntrace\n%f+%f*i\t%f\n\n", GSL_REAL(complexproduct), GSL_IMAG(complexproduct), gsl_complex_abs(complexproduct));
	}
	
		complexproduct=trace(five);
		conjugatetranspose(five);
		generatesu2(five, epsilon, generator);
	
	gsl_rng_free(generator);
	gsl_matrix_complex_free(one);
	gsl_matrix_complex_free(two);
	//~ gsl_matrix_complex_free(three);
	//~ gsl_matrix_complex_free(four);
	gsl_matrix_complex_free(multiplied);
	gsl_matrix_complex_free(five);
	return 0;
}


