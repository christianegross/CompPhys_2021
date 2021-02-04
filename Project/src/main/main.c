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



int main(int argc, char **argv){
	double a,b,c,d,det;
	double epsilon=1;
	double e,f,g,h;
	gsl_complex complexproduct=GSL_COMPLEX_ZERO;
	gsl_matrix_complex *one=gsl_matrix_complex_alloc(2,2);
	gsl_matrix_complex *two=gsl_matrix_complex_alloc(2,2);
	gsl_matrix_complex *three=gsl_matrix_complex_alloc(2,2);
	gsl_matrix_complex *four=gsl_matrix_complex_alloc(2,2);
	gsl_matrix_complex *multiplied=gsl_matrix_complex_alloc(2,2);	
	gsl_vector_complex_view columnzero, columnone;
	int seed=2;//use fixed seed: result should be exactly reproduced using the same seed
	gsl_rng *generator;
	generator=gsl_rng_alloc(gsl_rng_mt19937);//use mersenne-twister
	gsl_rng_set(generator, seed);
	
	for(int i=0;i<1;i+=1){	
		 /**
		  * choose four random doubles between -1 and 1, normalise through determinant of future matrix, so determinant is 1
		  * */
		a=2.0*gsl_rng_uniform(generator)-1.0;
		b=2.0*gsl_rng_uniform(generator)-1.0;
		c=2.0*gsl_rng_uniform(generator)-1.0;
		d=2.0*gsl_rng_uniform(generator)-1.0;
		det=a*a+b*b+c*c+d*d;
		a/=sqrt(det);
		b/=sqrt(det);
		c/=sqrt(det);
		d/=sqrt(det);
		
		/** set complex hermitian matrix*/
		gsl_matrix_complex_set(one, 0, 0, gsl_complex_rect(a, b));
		gsl_matrix_complex_set(one, 1, 1, gsl_complex_rect(a,-b));
		gsl_matrix_complex_set(one, 0, 1, gsl_complex_rect(c, d));
		gsl_matrix_complex_set(one, 1, 0, gsl_complex_rect(-c,d));
		
		gsl_matrix_complex_set(two, 0, 0, gsl_complex_rect(a,-b));
		gsl_matrix_complex_set(two, 1, 1, gsl_complex_rect(a, b));
		gsl_matrix_complex_set(two, 0, 1, gsl_complex_rect(-c,-d));
		gsl_matrix_complex_set(two, 1, 0, gsl_complex_rect(c,-d));
		/** see if uudagger =1*/
		gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(1,0), one, two, gsl_complex_rect(0,0), multiplied);
		
		printf("run no. %d\n matrix one \n", i);
		gsl_matrix_complex_fprintf(stdout, one, "%f");
		printf("\n matrix two\n");
		gsl_matrix_complex_fprintf(stdout, two, "%f");
		printf("\n multiplied\n");
		gsl_matrix_complex_fprintf(stdout, multiplied, "%f");
		printf("\n old and new determinant\n");
		printf("%f\t%f\n\n", det, a*a+b*b+c*c+d*d);
		
		/**
		 * Make hermitian Matrix 3 with entries between -1 and 1 {{e, f+ig},{f-ig, h}}
		 * Fill Matrix four =1+iepsilon*three={{1+iepsilon e, -epsilon g+iepsilon f}, {epsilong+iepsilon f, 1+iepsilon h}}
		 * */  
		
		e=2.0*gsl_rng_uniform(generator)-1.0;
		f=2.0*gsl_rng_uniform(generator)-1.0;
		g=2.0*gsl_rng_uniform(generator)-1.0;
		h=2.0*gsl_rng_uniform(generator)-1.0;
		
		gsl_matrix_complex_set(three, 0, 0, gsl_complex_rect(e, 0));
		gsl_matrix_complex_set(three, 1, 1, gsl_complex_rect(h, 0));
		gsl_matrix_complex_set(three, 0, 1, gsl_complex_rect(f, g));
		gsl_matrix_complex_set(three, 1, 0, gsl_complex_rect(f,-g));
		
		gsl_matrix_complex_set(four, 0, 0, gsl_complex_rect(1, epsilon*e));
		gsl_matrix_complex_set(four, 1, 1, gsl_complex_rect(1, epsilon*h));
		gsl_matrix_complex_set(four, 0, 1, gsl_complex_rect(-epsilon*g, epsilon*f));
		gsl_matrix_complex_set(four, 1, 0, gsl_complex_rect(epsilon*g, epsilon*f));
		
		
		gsl_matrix_complex_fprintf(stdout, four, "%f");
		
		columnzero=gsl_matrix_complex_column(four, 0);	//set up columns as vectors
		columnone=gsl_matrix_complex_column(four, 1);
		
		det=gsl_blas_dznrm2(&columnzero.vector);		//normalize vector zero
		gsl_blas_zdscal(1.0/det, &columnzero.vector);
		
		gsl_blas_zdotu(&columnzero.vector, &columnone.vector, &complexproduct);		//do Gram-Schmidt orthonormalization as in Script from T. Raesch
		gsl_blas_zaxpy(gsl_complex_mul_real(complexproduct, -1.0), &columnzero.vector, &columnone.vector); // should be one=-1*complexprodukt*zero+one
		
		det=gsl_blas_dznrm2(&columnone.vector);		//normalize vector one
		gsl_blas_zdscal(1.0/det, &columnone.vector);
		
		printf("\n%f\n", gsl_blas_dznrm2(&columnzero.vector));		//check normalization and Matrix
		printf("%f\n", gsl_blas_dznrm2(&columnone.vector));
		gsl_matrix_complex_fprintf(stdout, four, "%f");
		
		complexproduct=gsl_complex_add(gsl_complex_mul(gsl_matrix_complex_get(four, 0,0), gsl_matrix_complex_get(four, 1,1)), gsl_complex_mul_real(gsl_complex_mul(gsl_matrix_complex_get(four, 1,0), gsl_matrix_complex_get(four, 1,0)), -1.0));
		printf("\n%f+%f*i\n", GSL_REAL(complexproduct), GSL_IMAG(complexproduct));
		
		gsl_blas_zdotu(&columnzero.vector, &columnone.vector, &complexproduct); //check determinant and orthogonality
		printf("\n%f+%f*i\n\n", GSL_REAL(complexproduct), GSL_IMAG(complexproduct));
	}
	
	gsl_rng_free(generator);
	gsl_matrix_complex_free(one);
	gsl_matrix_complex_free(two);
	gsl_matrix_complex_free(multiplied);
	return 0;
}


