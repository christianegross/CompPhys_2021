//Exercise 6 for Computational Physics: Form factor of a boson bound state
//Nico Dichter, Christiane Groß
//18. December 2020-

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_block.h>
#include <gsl/gsl_rng.h>//random number generator
#include <gsl/gsl_randist.h>//pull random number from gaussian
#include "math.h"//exp-Function
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_integration.h>

/**
 * @brief determines Gauss-Legendre integration points and weights for given grid size, also sets additional point pN=q with weight 1
 * */
void getgridpoints(gsl_vector *momenta, gsl_vector *weights, double q, double pmax, int sizeofgrid){
	gsl_integration_glfixed_table *table=gsl_integration_glfixed_table_alloc(sizeofgrid);
	double momentum, weight;
	int errorcode;
	for (int k=0; k<sizeofgrid; k+=1){
		errorcode=gsl_integration_glfixed_point(0, pmax, k, &momentum, &weight, table);
		gsl_vector_set(momenta, k, momentum);
		gsl_vector_set(weights, k, weight);
	}
	gsl_vector_set(momenta, momenta->size-1, q);
	gsl_vector_set(weights, weights->size-1, 1);
	gsl_integration_glfixed_table_free(table);	
}

/**
 * @brief uses Gauss-Legendre integration to calculate the potential, form of the potential given in lecture 7, 
 * fills matrix with the potential
 * */
inline void fillpotentialmatrix(gsl_matrix *pot, gsl_vector *momenta, int l, int sizeofangulargrid, double mu){
	gsl_integration_glfixed_table *table=gsl_integration_glfixed_table_alloc(sizeofangulargrid);
	double result, x, weight, pi, pj;
	int errorcode;
	for (int i=0; i<pot->size1; i+=1){
		pi=gsl_vector_get(momenta, i);
		for (int j=0; j<pot->size2; j+=1){
			pj=gsl_vector_get(momenta, j);
			result=0;
			for (int k=0; k<sizeofangulargrid; k+=1){
				errorcode=gsl_integration_glfixed_point(-1, 1, k, &x, &weight, table);
				result+=weight*gsl_sf_legendre_Pl(l, x)/(pi*pi+pj*pj-2*pi*pj*x+mu*mu);
			}
			gsl_matrix_set(pot, i, j, 2*M_PI*result);
		}
	}
	gsl_integration_glfixed_table_free(table);	
}

int main(int argc, char **argv){
	/**
	 * @note set up parameters
	 * */
	int sizeofgrid=20;
	int sizeofangulargrid=20;
	int l=5;
	double E=1;
	double mu=938.92;
	double q=sqrt(2*mu*E);
	double pmax=100;
	
	/**
	 * @note set up gsl vectors and matrices, streams for toring results
	 * */
	gsl_vector *momenta=gsl_vector_alloc(sizeofgrid+1);
	gsl_vector *weights=gsl_vector_alloc(sizeofgrid+1);
	
	gsl_matrix *pot=gsl_matrix_alloc(sizeofgrid+1, sizeofgrid+1);
	
	FILE *test=fopen("data/test.dat", "w");
	
	/**
	 * @note test functions
	 * */
	
	getgridpoints(momenta, weights, q, pmax, sizeofgrid);
	
	gsl_vector_fprintf(test, momenta, "%e");
	fprintf(test, "\n\n");
	gsl_vector_fprintf(test, weights, "%e");
	fprintf(test, "\n\n");
	
	fillpotentialmatrix(pot, momenta, l, sizeofangulargrid, mu);
	
	gsl_matrix_fprintf(test, pot, "%e");
	fprintf(test, "\n\n");
	
	/**
	 * @note cleanup
	 * */
	
	gsl_vector_free(momenta);
	gsl_vector_free(weights);
	gsl_matrix_free(pot);
	
	return 0;
}


