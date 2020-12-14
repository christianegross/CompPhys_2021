//Exercise 6 for Computational Physics: Form factor of a boson bound state
//Nico Dichter, Christiane Groß
//12. December 2020-

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_block.h>
#include <gsl/gsl_rng.h>//random number generator
#include <gsl/gsl_randist.h>//pull random number from gaussian
#include "math.h"//exp-Function
#include <gsl/gsl_vector.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_integration.h>

/**
 * @brief reads content of the given file into the vectors, ignores other constants given in the file
 * 
 * @param wavefunctionfile file containing the data
 * @param p, w, wf stores contents of columns: momentum, weight for integration, value of the ewavefunction at the given p
 * */
void readinwavefunction(FILE * wavefunctionfile, gsl_vector *p, gsl_vector *w, gsl_vector* wf){
	/**
	 * @note cannot use built-in function because there are lines with other data, and several columns
	 * @note use read variables to store values from file and assign to columns, empty for values from file that are not needed
	 * */
	double pread, wread, wfread, empty;
	rewind(wavefunctionfile);
	empty=fscanf(wavefunctionfile, "# Lam=   %le  A=  %le  mb=   %le C0=  %le \n# p [fm-1]       w             wf [fm3/2]\n", &empty, &empty, &empty, &empty);
	for (int line=0; line<p->size; line+=1){
		if(3== fscanf(wavefunctionfile, "%le%le%le\n", &pread, &wread, &wfread)){
		gsl_vector_set(p, line, pread);
		gsl_vector_set(w, line, wread);
		gsl_vector_set(wf, line, wfread);
		}
	}
}

/**
 * @brief calculates form factor according to equation given on sheet
 * @param q value of q for which F is calculated
 * @param p, p_weights, wf momnetum, weight for momentum integration, value of wavefunction for given p
 * @param interpolater, accelerator used to calculate wavefunction for values not in p-vector
 * @param nx grid points used for angular integration
 * */
double formfactor(double q, gsl_vector *p, gsl_vector* p_weights, gsl_vector *wf, gsl_interp* interpolater, gsl_interp_accel* accelerator, int nx ){
	double x=0.5;			//angle
	double p_p=1;			//momentum
	double x_weight=1;		//weight for angle integration
	double p_p_z=0;			//z-component of momentum
	double F=0;				//Form factor
	double F_add=0;			//Form factor contribution for one given p_p for all angles
	double p_2=0;			//modulus of p-0.5q
	int l_z=5, l=5;			//angular momentum
	
	gsl_integration_glfixed_table* table= gsl_integration_glfixed_table_alloc (nx);		//used for integrating over angles
	
	for(int j=1;j<p->size-1;j++){	//for all p_p: first sum over all angles for one p_p, then sum over these contributions for all p_p
		F_add=0;
		p_p=gsl_vector_get (p, j);
		for(int i=0;i<nx;i++){
			gsl_integration_glfixed_point (-1, 1, i, &x, &x_weight, table);		//determine weight for x-integration
			p_p_z=p_p*x;
			p_2=sqrt (p_p*p_p-p_p_z*q+q*q/4.);
			//modulus of p-0.5q has to be between p_min and p_max									
			if (p_2>1.03e-03&&p_2<34.9){
				//Y_lm(p)*Y_lm(p-0.5q)*psi(p)*psi(p-0.5q), assume psi is real everywhere
				F_add+=x_weight*gsl_sf_legendre_sphPlm (l, l_z, x)*gsl_sf_legendre_sphPlm (l, l_z, (p_p_z-0.5*q)/p_2)
					*gsl_interp_eval (interpolater, p->data, wf->data, p_p, accelerator)*gsl_interp_eval (interpolater, p->data, wf->data, p_2, accelerator);
				//printf ("F_add=%e,p_p=%e,x=%e\n",F_add,p_p,x);
			}
			else{printf("p not accepted\n");}
		}
		F+=F_add*p_p*p_p*gsl_vector_get (p_weights, j);
	}
	
	gsl_integration_glfixed_table_free (table);
	
	return F*2*M_PI;
}

int main(int argc, char **argv){
	int datasetsize=60;
	double result=0;
	
	/**
	 * @note	General allocations and opening of streams
	 */
	gsl_vector *p=gsl_vector_alloc(datasetsize);
	gsl_vector *w=gsl_vector_alloc(datasetsize);
	gsl_vector *wf=gsl_vector_alloc(datasetsize);
	gsl_vector *norm=gsl_vector_alloc(datasetsize);
	FILE *wavefunctionfile=fopen("data/wavefunctions/wf-obe-lam=1200.00.dat", "r");
	gsl_interp* interpolate_test=gsl_interp_alloc(gsl_interp_cspline ,datasetsize);
	gsl_interp_accel* acc_test=gsl_interp_accel_alloc ();
	
	
	/**
	 * @note	Reading of dataset and init of interpolater
	 */
	readinwavefunction(wavefunctionfile, p, w, wf);
	gsl_interp_init (interpolate_test, p->data, wf->data, datasetsize);
	
	/**
	 * @note Test of normalization
	 */
	gsl_vector_memcpy (norm, wf);
	gsl_vector_mul (norm, wf);
	gsl_vector_mul (norm, p);
	gsl_vector_mul (norm, p);
	gsl_vector_mul (norm, w);
	for(int i=0;i<norm->size;i++){
		result+=gsl_vector_get (norm, i);
	}
	
	
	printf ("norm=%e, F=%e\n",result,formfactor (0, p, w, wf, interpolate_test, acc_test, 100));
	
	/**
	 * @note	Cleanup after each dataset
	 */
	gsl_interp_accel_reset (acc_test);
	
	
	
	
	/**
	 * @note	Complete cleanup
	 */
	fclose(wavefunctionfile);
	gsl_vector_free(p);
	gsl_vector_free(w);
	gsl_vector_free(wf);
	gsl_vector_free(norm);
	gsl_interp_free (interpolate_test);
	gsl_interp_accel_free (acc_test);

	return 0;
}


