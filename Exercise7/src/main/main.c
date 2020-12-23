//Exercise 7 for Computational Physics: Two-body scattering
//Nico Dichter, Christiane Groß
//18. December 2020-

#include <stdio.h>
#include <stdlib.h>
#include "math.h"						//exp-Function etc
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_legendre.h>		//Legendre-polynomial
#include <gsl/gsl_integration.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>		//needed for LU-decomp

/**
 * @brief determines Gauss-Legendre integration points and weights for given grid size, also sets additional point pN=q with weight 1
 * */
void getgridpoints(gsl_vector *momenta, gsl_vector *weights, double q, double pmax, int sizeofgrid){
	gsl_integration_glfixed_table *table=gsl_integration_glfixed_table_alloc(sizeofgrid);
	double momentum, weight;
	for (int k=0; k<sizeofgrid; k+=1){
		gsl_integration_glfixed_point(0, pmax, k, &momentum, &weight, table);
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
inline void fillpotentialmatrix(gsl_matrix_complex *pot, gsl_vector *momenta, int l, int sizeofangulargrid, double mb, double Apref, double C_0, double lambda){
	gsl_integration_glfixed_table *table=gsl_integration_glfixed_table_alloc(sizeofangulargrid);
	double result, x, weight, pi, pj, qvalsquare;
	for (int i=0; i<pot->size1; i+=1){
		pi=gsl_vector_get(momenta, i);
		for (int j=0; j<pot->size2; j+=1){
			pj=gsl_vector_get(momenta, j);
			result=0;
			for (int k=0; k<sizeofangulargrid; k+=1){
				gsl_integration_glfixed_point(-1, 1, k, &x, &weight, table);
				qvalsquare=pi*pi+pj*pj-2*pi*pj*x;
				result+=weight*gsl_sf_legendre_Pl(l, x)/(qvalsquare+mb*mb)*exp(-1.0*(qvalsquare+mb*mb)/lambda/lambda);
			}
			result*=Apref;
			if (l==0){
				result+=C_0*exp(-1.0*(pi*pi+pj*pj)/lambda/lambda);
			}
			gsl_matrix_complex_set(pot, i, j, gsl_complex_rect(result, 0));
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
			result-=2.0*mu*GSL_REAL(gsl_matrix_complex_get(pot, i, k))*pk*pk*gsl_vector_get(weights, k)/(q*q-pk*pk);
			gsl_matrix_complex_set(a, i, k, gsl_complex_rect(result, 0));
		}
		//k==N
		result=0;
		if(i==a->size2-1){result+=1;}
		for (int m=0; m<a->size2-1; m+=1){
			pm=gsl_vector_get(momenta, m);
			result+=2.0*mu*GSL_REAL(gsl_matrix_complex_get(pot, i, pot->size2-1))*q*q*gsl_vector_get(weights, m)/(q*q-pm*pm);
		}
		result-=1.0*mu*q*GSL_REAL(gsl_matrix_complex_get(pot, i, pot->size2-1))*log(1.0*(pmax+q)/(pmax-q));
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
			result-=2.0*mu*GSL_REAL(gsl_matrix_complex_get(pot, i, k))*pk*pk*gsl_vector_get(weights, k)/(q*q-pk*pk);
			gsl_matrix_complex_set(a, i, k, gsl_complex_rect(result, 0));
		}
		//k==N
		result=0;
		if(i==a->size2-1){result+=1;}
		for (int m=0; m<a->size2-1; m+=1){
			pm=gsl_vector_get(momenta, m);
			result+=2.0*mu*GSL_REAL(gsl_matrix_complex_get(pot, i, pot->size2-1))*q*q*gsl_vector_get(weights, m)/(q*q-pm*pm);
		}
		resultimaginary=M_PI*mu*q*GSL_REAL(gsl_matrix_complex_get(pot, i, pot->size2-1));
		gsl_matrix_complex_set(a, i, a->size2-1, gsl_complex_rect(result, resultimaginary));
	}
} 

/**
 * @brief calculates s
 * */
inline gsl_complex calculate_s(gsl_complex tnn, double mu, double q){
	gsl_complex result=gsl_complex_add_real(gsl_complex_mul_imag(gsl_complex_mul_real(tnn, 2*M_PI*mu*q), -1), 1);
	return result;
}

int main(int argc, char **argv){
	/**
	 * @note set up parameters
	 * */
	double hbarc=197.3;
	double lambda=800.0/hbarc;
	double Apref=-0.1544435;
	double C_0=2.470795e-2;
	int sizeofgrid=60;
	int sizeofangulargrid=60;
	int l=0;
	int lmax=7;
	double E=1.0/hbarc;
	double mu=938.92/hbarc;
	double mb=138.0/hbarc;
	double q=sqrt(2.0*mu*E);
	double pmax=100;
	gsl_complex tnn, s;
	int size, angularsize;
	
	/**
	 * @note set up gsl vectors and matrices, streams for storing results
	 * */
	gsl_vector *momenta=gsl_vector_alloc(sizeofgrid+1);
	gsl_vector *weights=gsl_vector_alloc(sizeofgrid+1);
	gsl_vector_complex *t_l=gsl_vector_complex_alloc (lmax);
	gsl_vector_complex *tknvector=gsl_vector_complex_alloc(sizeofgrid+1);
	gsl_vector_complex *potvector=gsl_vector_complex_alloc(sizeofgrid+1);
	
	gsl_matrix_complex *pot=gsl_matrix_complex_alloc(sizeofgrid+1, sizeofgrid+1);
	gsl_matrix_complex *a=gsl_matrix_complex_alloc(sizeofgrid+1, sizeofgrid+1);
	gsl_matrix_complex *awopm=gsl_matrix_complex_alloc(sizeofgrid+1, sizeofgrid+1); 		//a without pmax
	
	gsl_vector_view p, w;
	gsl_vector_complex_view tkn, vin;
	gsl_matrix_complex_view V, A, Awopm;

	
	FILE *test=fopen("data/test.dat", "w");
	FILE *result3=fopen("data/result3.dat", "w");
	FILE *result4=fopen("data/result4.dat", "w");
	FILE *result5=fopen("data/result5.dat", "w");

	/**
	 * @note measurements for exercise 3
	 * @note get grids, fill matrices, do LU decomposition, solve system Aik*tkN=ViN, do the same for Awopm
	 * */
	fprintf(result3, "size\tangularsize\tpmax\tR(tnn)\tI(tnn)\tAbs(tnn)\tR(tnnwopm)\tI(tnnwopm)\tAbs(tnnwopm)\n");
	for (size=4; size<sizeofgrid; size+=4){
		/**
		 * @note set up subvectors, submatrices, permutation
		 * */
		p=gsl_vector_subvector(momenta, 0, size+1);
		w=gsl_vector_subvector(weights, 0, size+1);
		tkn=gsl_vector_complex_subvector(tknvector, 0, size+1);
		vin=gsl_vector_complex_subvector(potvector, 0, size+1);
		V=gsl_matrix_complex_submatrix(pot, 0, 0, size+1, size+1);
		A=gsl_matrix_complex_submatrix(a, 0, 0, size+1, size+1);
		Awopm=gsl_matrix_complex_submatrix(awopm, 0, 0, size+1, size+1);
		gsl_permutation *permutation= gsl_permutation_alloc(size+1);
		int signum=1;
		
		for (angularsize=4; angularsize<=sizeofangulargrid; angularsize+=4){
			for (int maxp=1; maxp<=200; maxp+=1){
				pmax=75.0*maxp;
				getgridpoints(&p.vector, &w.vector, q, pmax, size);
				fillpotentialmatrix(&V.matrix, &p.vector, l, angularsize, mb, Apref, C_0, lambda);
				fillmatrixa(&A.matrix, &V.matrix, &p.vector, &w.vector, mu, pmax);
				fillmatrixawopm(&Awopm.matrix, &V.matrix, &p.vector, &w.vector, mu);
				vin=gsl_matrix_complex_column(&V.matrix, size);
				fprintf(result3, "%3d\t%3d\t%e\t", size, angularsize, pmax);
				
				gsl_linalg_complex_LU_decomp(&A.matrix, permutation, &signum);														//make LU-decomposition of A
				gsl_linalg_complex_LU_solve(&A.matrix, permutation, &vin.vector, &tkn.vector);										//solve Aik*tkN=Vin
				tnn=gsl_vector_complex_get(&tkn.vector, size);
				fprintf(result3, "%e\t%e\t%e\t", GSL_REAL(tnn), GSL_IMAG(tnn), gsl_complex_abs(tnn));
				
				gsl_linalg_complex_LU_decomp(&Awopm.matrix, permutation, &signum);													//make LU-decomposition of A
				gsl_linalg_complex_LU_solve(&Awopm.matrix, permutation, &vin.vector, &tkn.vector);									//solve Aik*tkN=Vin
				tnn=gsl_vector_complex_get(&tkn.vector, size);
				fprintf(result3, "%e\t%e\t%e\n", GSL_REAL(tnn), GSL_IMAG(tnn), gsl_complex_abs(tnn));
			}
		}
		gsl_permutation_free(permutation);
	}
	
	
	/**
	 * @note measurements for exercise 4
	 * */
	size=50;
	angularsize=sizeofangulargrid;
	p=gsl_vector_subvector(momenta, 0, size+1);
	w=gsl_vector_subvector(weights, 0, size+1);
	tkn=gsl_vector_complex_subvector(tknvector, 0, size+1);
	vin=gsl_vector_complex_subvector(potvector, 0, size+1);
	V=gsl_matrix_complex_submatrix(pot, 0, 0, size+1, size+1);
	A=gsl_matrix_complex_submatrix(a, 0, 0, size+1, size+1);
	Awopm=gsl_matrix_complex_submatrix(awopm, 0, 0, size+1, size+1);
	gsl_permutation *permutation= gsl_permutation_calloc(size+1);
	int signum=1;
	fprintf(result4, "energy\tq\tAbs(s)\tArg(s)\tAbs(swopm)\tArg(swopm)\n");
	for (int energy=0; energy<=200; energy +=1){p=gsl_vector_subvector(momenta, 0, size+1);
		q=sqrt(2.0*mu*energy/hbarc);
		pmax=100000;
		
		getgridpoints(&p.vector, &w.vector, q, pmax, size);
		fillpotentialmatrix(&V.matrix, &p.vector, l, angularsize, mb, Apref, C_0, lambda);
		fillmatrixa(&A.matrix, &V.matrix, &p.vector, &w.vector, mu, pmax);
		fillmatrixawopm(&Awopm.matrix, &V.matrix, &p.vector, &w.vector, mu);
		vin=gsl_matrix_complex_column(&V.matrix, size);
		fprintf(result4, "%3d\t%e\t", energy, q);
		
		gsl_linalg_complex_LU_decomp(&A.matrix, permutation, &signum);														//make LU-decomposition of A
		gsl_linalg_complex_LU_solve(&A.matrix, permutation, &vin.vector, &tkn.vector);										//solve Aik*tkN=Vin
		tnn=gsl_vector_complex_get(&tkn.vector, size);
		s=calculate_s(tnn, mu, q);
		fprintf(result4, "%e\t%e\t",gsl_complex_abs(s)-1.0, gsl_complex_arg(s));
		
		gsl_linalg_complex_LU_decomp(&Awopm.matrix, permutation, &signum);													//make LU-decomposition of A
		gsl_linalg_complex_LU_solve(&Awopm.matrix, permutation, &vin.vector, &tkn.vector);									//solve Aik*tkN=Vin
		tnn=gsl_vector_complex_get(&tkn.vector, size);
		s=calculate_s(tnn, mu, q);
		fprintf(result4, "%e\t%e\n", gsl_complex_abs(s)-1.0, gsl_complex_arg(s));
	}
	gsl_permutation_free(permutation);
	
	/**
	 * @note measurements for exercise 5
	 * */
	E=10.0/hbarc;
	q=sqrt(2.0*mu*E);
	size=50;
	angularsize=sizeofangulargrid;
	pmax=100000;
	fprintf(result5, "x\tl\tdc\n");
	
	/**
	 * @note set up subvectors, submatrices, permutation
	 * */
	p=gsl_vector_subvector(momenta, 0, size+1);
	w=gsl_vector_subvector(weights, 0, size+1);
	tkn=gsl_vector_complex_subvector(tknvector, 0, size+1);
	vin=gsl_vector_complex_subvector(potvector, 0, size+1);
	V=gsl_matrix_complex_submatrix(pot, 0, 0, size+1, size+1);
	A=gsl_matrix_complex_submatrix(a, 0, 0, size+1, size+1);
	getgridpoints(&p.vector, &w.vector, q, pmax, size);
	
	/**
	 * @note	Save all required t_l
	 */
	for (l=0;l<lmax;l++){
		signum=1;
		permutation= gsl_permutation_alloc(size+1);
		fillpotentialmatrix(&V.matrix, &p.vector, l, angularsize, mb, Apref, C_0, lambda);
		fillmatrixa(&A.matrix, &V.matrix, &p.vector, &w.vector, mu, pmax);
		vin=gsl_matrix_complex_column(&V.matrix, size);
	
		gsl_linalg_complex_LU_decomp(&A.matrix, permutation, &signum);														//make LU-decomposition of A
		gsl_linalg_complex_LU_solve(&A.matrix, permutation, &vin.vector, &tkn.vector);										//solve Aik*tkN=Vin
		gsl_vector_complex_set (t_l, l, gsl_vector_complex_get(&tkn.vector, size));
		gsl_permutation_free(permutation);
	}
	
	/**
	 * @note Calculate the differential cross sction depending on x
	 */
	gsl_complex t=gsl_complex_rect (0, 0);
	for(double x=-0.99;x<=0.99;x+=0.01){
		t=gsl_complex_rect (0, 0);
		for(l=0;l<lmax;l++){
			t=gsl_complex_add (t, gsl_complex_mul_real(gsl_vector_complex_get(t_l,l),gsl_sf_legendre_Pl(l, x)*(2*l+1)/(4*M_PI)));
			fprintf(result5, "%e\t%d\t%e\n", x, l, gsl_pow_4 (2*M_PI)*mu*mu*gsl_complex_abs2 (t));
		}
	}
	
	
	/**
	 * @note cleanup
	 * */
	
	gsl_vector_free(momenta);
	gsl_vector_free(weights);
	gsl_vector_complex_free (t_l);
	gsl_vector_complex_free(tknvector);
	gsl_vector_complex_free(potvector);
	gsl_matrix_complex_free(pot);
	gsl_matrix_complex_free(a);
	gsl_matrix_complex_free(awopm);
	fclose(test);
	fclose(result3);
	fclose(result4);
	fclose(result5);
	
	return 0;
}


