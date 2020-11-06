// Die folgenden Befehle müssen im Verzeichnis, 
// in dem die "Makefile" liegt ausgeführt werden:
// make
// make calc
// Nico Dichter, Florian Giefer
// Die Plots können mit "make plot" erstellt werden.
/** 
 *  @file global_structs.h 
 *  @brief Provides structs,to be use in all general functions
 *  @author Dichter nd@nocoffeetech.de 
 * 
 *   
 * 
 *  @date 14.06.2020	First Implementation of: data_sets_t,parameters_t
 *  @date 08.07.2020	Rewrote some macros and added more parameters
 * 
 * 
 *  @bug No known bugs  
 * 
 *  @version 0.2
 */

#ifndef GLOBAL_STRUCTS_H
#define GLOBAL_STRUCTS_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <complex.h>

#define GLOBAL_STRUCTS_H
#define alpha (paraset->a_d)
#define t_max (paraset->b_d)
#define t_l (paraset->c_d)
#define amp_a (paraset->d_d)
#define amp_b (paraset->e_d)
#define beta (paraset->f_d)
#define n_step (paraset->a_i)
#define m_pow (paraset->b_i)
#define fp_len (paraset->c_i)
#define a_h (paraset->h1)
#define t_int (paraset->z1)
#define alp_0 (paraset->var1)
#define alp_1 (paraset->var2)

#define h_e (dataset->a_c)
#define h_s (dataset->b_c)
#define filepath (dataset->s1)

typedef struct data_sets{
	double *a;
	double *b;
	double *c;
	double *d;
	double complex *a_c;
	double complex *b_c;
	char *s1;
	gsl_matrix* m1;
	gsl_matrix* m2;
	gsl_matrix* m3;
	gsl_vector* v1;
	gsl_vector* v2;
	gsl_vector* v3;
}data_sets_t;

typedef struct parameters{
	double a_d;
	double b_d;
	double c_d;
	double d_d;
	double e_d;
	double f_d;
	int    a_i;
	int    b_i;
	int    c_i;
	int    d_i;
	int    e_i;
	double h1;
	double h2;
	double h3;
	double z1;
	double var1;
	double var2;
	int check1;
}parameters_t;

#endif // GLOBAL_STRUCTS_H