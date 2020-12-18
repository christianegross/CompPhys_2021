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
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_blas.h>

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
inline void fillpotentialmatrix(gsl_matrix_complex *pot, gsl_vector *momenta, int l, int sizeofangulargrid, double mu){
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
			gsl_matrix_complex_set(pot, i, j, gsl_complex_rect(2*M_PI*result, 0));
		}
	}
	gsl_integration_glfixed_table_free(table);	
}

/**
 * @brief determines elements of matrix A as given on sheet
 * */
inline void fillmatrixa(gsl_matrix_complex *a, gsl_matrix_complex *pot, gsl_vector *momenta, gsl_vector *weights, double mu, double pmax){
	double result, resultimaginary, pk, pm;
	double q=gsl_vector_get(momenta, momenta->size-1);
	for (int i=0; i<a->size1; i+=1){
		//k!=N
		for (int k=0; k<a->size2-1; k+=1){
			result=0;
			pk=gsl_vector_get(momenta, k);
			if(i==k){result+=1;}
			result-=2*mu*GSL_REAL(gsl_matrix_complex_get(pot, i, k))*pk*pk*gsl_vector_get(weights, k)/(q*q-pk*pk);
			gsl_matrix_complex_set(a, i, k, gsl_complex_rect(result, 0));
		}
		//k==N
		result=0;
		if(i==a->size2-1){result+=1;}
		for (int m=0; m<a->size2-1; m+=1){
			pm=gsl_vector_get(momenta, m);
			result+=2*mu*GSL_REAL(gsl_matrix_complex_get(pot, i, m))*pm*pm*gsl_vector_get(weights, m)/(q*q-pm*pm);
		}
		result-=mu*q*GSL_REAL(gsl_matrix_complex_get(pot, i, pot->size2-1))*log((pmax+q)/(pmax-q));
		resultimaginary=M_PI*mu*q*GSL_REAL(gsl_matrix_complex_get(pot, i, pot->size2-1));
		gsl_matrix_complex_set(a, i, a->size2-1, gsl_complex_rect(result, resultimaginary));
	}
} 

/**
 * @brief determines elements of matrix A as given on sheet, neglects the contributions from pmax
 * */
inline void fillmatrixawopm(gsl_matrix_complex *a, gsl_matrix_complex *pot, gsl_vector *momenta, gsl_vector *weights, double mu){
	double result, resultimaginary, pk, pm;
	double q=gsl_vector_get(momenta, momenta->size-1);
	for (int i=0; i<a->size1; i+=1){
		//k!=N
		for (int k=0; k<a->size2-1; k+=1){
			result=0;
			pk=gsl_vector_get(momenta, k);
			if(i==k){result+=1;}
			result-=2*mu*GSL_REAL(gsl_matrix_complex_get(pot, i, k))*pk*pk*gsl_vector_get(weights, k)/(q*q-pk*pk);
			gsl_matrix_complex_set(a, i, k, gsl_complex_rect(result, 0));
		}
		//k==N
		result=0;
		if(i==a->size2-1){result+=1;}
		for (int m=0; m<a->size2-1; m+=1){
			pm=gsl_vector_get(momenta, m);
			result+=2*mu*GSL_REAL(gsl_matrix_complex_get(pot, i, m))*pm*pm*gsl_vector_get(weights, m)/(q*q-pm*pm);
		}
		resultimaginary=M_PI*mu*q*GSL_REAL(gsl_matrix_complex_get(pot, i, pot->size2-1));
		gsl_matrix_complex_set(a, i, a->size2-1, gsl_complex_rect(result, resultimaginary));
	}
} 

/**
 * @brief calculates the modulus of s
 * */
inline double absofs(gsl_complex tnn, double mu, double q){
	gsl_complex result=gsl_complex_add_real(gsl_complex_mul_imag(gsl_complex_mul_real(tnn, 2*M_PI*mu*q), -1), 1);
	return gsl_complex_abs(result);
}

int main(int argc, char **argv){
	/**
	 * @note set up parameters
	 * */
	int sizeofgrid=40;
	int sizeofangulargrid=40;
	int l=0;
	double E=1;
	double mu=938.92;
	double q=sqrt(2*mu*E);
	double pmax=100;
	gsl_complex tnn;
	
	/**
	 * @note set up gsl vectors and matrices, streams for toring results
	 * */
	gsl_vector *momenta=gsl_vector_alloc(sizeofgrid+1);
	gsl_vector *weights=gsl_vector_alloc(sizeofgrid+1);
	
	gsl_matrix_complex *pot=gsl_matrix_complex_alloc(sizeofgrid+1, sizeofgrid+1);
	gsl_matrix_complex *a=gsl_matrix_complex_alloc(sizeofgrid+1, sizeofgrid+1);
	gsl_matrix_complex *awopm=gsl_matrix_complex_alloc(sizeofgrid+1, sizeofgrid+1); 		//a without pmax
	gsl_matrix_complex *t=gsl_matrix_complex_alloc(sizeofgrid+1, sizeofgrid+1);
	
	gsl_vector_view p, w;
	gsl_matrix_complex_view V, A, Awopm, T;
	
	FILE *test=fopen("data/test.dat", "w");
	FILE *result3=fopen("data/result3.dat", "w");
	
	/**
	 * @note test functions
	 * */
	
	getgridpoints(momenta, weights, q, pmax, sizeofgrid);
	
	//~ gsl_vector_fprintf(test, momenta, "%e");
	//~ fprintf(test, "\n\n");
	//~ gsl_vector_fprintf(test, weights, "%e");
	//~ fprintf(test, "\n\n");
	
	fillpotentialmatrix(pot, momenta, l, sizeofangulargrid, mu);
	
	//~ gsl_matrix_complex_fprintf(test, pot, "%e");
	//~ fprintf(test, "\n\n");
	
	fillmatrixa(a, pot, momenta, weights, mu, pmax);
	
	//~ gsl_matrix_complex_fprintf(test, a, "%e");
	//~ fprintf(test, "\n\n");
	
	/**
	 * @note determine t by using t=A^-1*V, copy V in empty t so V does not get lost
	 * */
	gsl_matrix_complex_memcpy(t, pot);
	
	//~ gsl_matrix_complex_fprintf(test, t, "%e");
	//~ fprintf(test, "\n\n");
	
	gsl_blas_ztrsm(CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, gsl_complex_rect(1.0, 0), a, t);
	
	gsl_matrix_complex_fprintf(test, t, "%e");
	fprintf(test, "\n\n");
	
	/**
	 * @note determine t_NN, check if |S|==1
	 * */
	tnn=gsl_matrix_complex_get(t, t->size1-1, t->size2-1);
	printf("%e\n", GSL_REAL(tnn));
	//~ s=gsl_complex_sub_real(gsl_complex_mul_imag(gsl_complex_mul_real(gsl_matrix_complex_get(t, t->size1-1, t->size2-1), 2*M_PI*mu*q), -1), 1);
	//~ printf("%e\n", gsl_complex_abs(s));
	printf("%e\n", absofs(tnn, mu, q));
	
	/**
	 * @note measurements for exercise 3
	 * */
	 fprintf(result3, "size\tangularsize\tpmax\tR(tnn)\tI(tnn)\tAbs(tnn)\tR(tnnwopm)\tI(tnnwopm)\tAbs(tnnwopm)\n");
	for (int size=4; size<sizeofgrid; size+=4){
		p=gsl_vector_subvector(momenta, 0, size+1);
		w=gsl_vector_subvector(weights, 0, size+1);
		V=gsl_matrix_complex_submatrix(pot, 0, 0, size+1, size+1);
		A=gsl_matrix_complex_submatrix(a, 0, 0, size+1, size+1);
		Awopm=gsl_matrix_complex_submatrix(awopm, 0, 0, size+1, size+1);
		T=gsl_matrix_complex_submatrix(t, 0, 0, size+1, size+1);
		for (int angularsize=4; angularsize<sizeofangulargrid; angularsize+=4){
			for (int maxp=1; maxp<20; maxp+=1){
				pmax=100.0*maxp;
				getgridpoints(&p.vector, &w.vector, q, pmax, size);
				fillpotentialmatrix(&V.matrix, &p.vector, l, angularsize, mu);
				fillmatrixa(&A.matrix, &V.matrix, &p.vector, &w.vector, mu, pmax);
				fillmatrixawopm(&Awopm.matrix, &V.matrix, &p.vector, &w.vector, mu);
				fprintf(result3, "%d\t%d\t%e\t", size, angularsize, pmax);
				gsl_matrix_complex_memcpy(&T.matrix, &V.matrix);
				gsl_blas_ztrsm(CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, gsl_complex_rect(1.0, 0), &A.matrix, &T.matrix);
				tnn=gsl_matrix_complex_get(&T.matrix, size, size);
				fprintf(result3, "%e\t%e\t%e\t", GSL_REAL(tnn), GSL_IMAG(tnn), gsl_complex_abs(tnn));
				gsl_matrix_complex_memcpy(&T.matrix, &V.matrix);
				gsl_blas_ztrsm(CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, gsl_complex_rect(1.0, 0), &Awopm.matrix, &T.matrix);
				tnn=gsl_matrix_complex_get(&T.matrix, size, size);
				fprintf(result3, "%e\t%e\t%e\n", GSL_REAL(tnn), GSL_IMAG(tnn), gsl_complex_abs(tnn));
			}
		}
	}
	/**
	 * @note cleanup
	 * */
	
	gsl_vector_free(momenta);
	gsl_vector_free(weights);
	gsl_matrix_complex_free(pot);
	gsl_matrix_complex_free(a);
	fclose(test);
	fclose(result3);
	
	return 0;
}


