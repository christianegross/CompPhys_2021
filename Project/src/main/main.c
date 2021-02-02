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
	gsl_matrix_complex *one=gsl_matrix_complex_alloc(2,2);
	gsl_matrix_complex *two=gsl_matrix_complex_alloc(2,2);
	gsl_matrix_complex *multiplied=gsl_matrix_complex_alloc(2,2);	
	int seed=2;//use fixed seed: result should be exactly reproduced using the same seed
	gsl_rng *generator;
	generator=gsl_rng_alloc(gsl_rng_mt19937);//use mersenne-twister
	gsl_rng_set(generator, seed);
	
	for(int i=1;i<50;i+=1){	
		a=2.0*gsl_rng_uniform(generator)-1.0;
		b=2.0*gsl_rng_uniform(generator)-1.0;
		c=2.0*gsl_rng_uniform(generator)-1.0;
		d=2.0*gsl_rng_uniform(generator)-1.0;
		det=a*a+b*b+c*c+d*d;
		a/=sqrt(det);
		b/=sqrt(det);
		c/=sqrt(det);
		d/=sqrt(det);
		
		gsl_matrix_complex_set(one, 0, 0, gsl_complex_rect(a, b));
		gsl_matrix_complex_set(one, 1, 1, gsl_complex_rect(a,-b));
		gsl_matrix_complex_set(one, 0, 1, gsl_complex_rect(c, d));
		gsl_matrix_complex_set(one, 1, 0, gsl_complex_rect(-c,d));
		
		gsl_matrix_complex_set(two, 0, 0, gsl_complex_rect(a,-b));
		gsl_matrix_complex_set(two, 1, 1, gsl_complex_rect(a, b));
		gsl_matrix_complex_set(two, 0, 1, gsl_complex_rect(-c,-d));
		gsl_matrix_complex_set(two, 1, 0, gsl_complex_rect(c,-d));
		
		gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(1,0), one, two, gsl_complex_rect(0,0), multiplied);
		
		printf("run no. %d\n matrix one \n", i);
		gsl_matrix_complex_fprintf(stdout, one, "%f");
		printf("\n matrix two\n");
		gsl_matrix_complex_fprintf(stdout, two, "%f");
		printf("\n multiplied\n");
		gsl_matrix_complex_fprintf(stdout, multiplied, "%f");
		printf("\n old and new determinant\n");
		printf("%f\t%f\n\n", det, a*a+b*b+c*c+d*d); 
	}
	
	gsl_rng_free(generator);
	gsl_matrix_complex_free(one);
	gsl_matrix_complex_free(two);
	gsl_matrix_complex_free(multiplied);
	return 0;
}


