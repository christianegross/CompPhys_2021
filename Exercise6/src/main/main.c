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

void readinwavefunction(FILE * wavefunctionfile, gsl_vector *p, gsl_vector *w, gsl_vector* wf){
	/**
	 * @note cannot use built-in function because there are lines with other data, and several columns
	 * */
	 double pread, wread, wfread, empty;
	 rewind(wavefunctionfile);
	 empty=fscanf(wavefunctionfile, "# Lam=   %le  A=  %le  mb=   %le C0=  %le \n# p [fm-1]       w             wf [fm3/2]\n", &empty, &empty, &empty, &empty);
	 for (int line=0; line<p->size; line+=1){
		 if(3== fscanf(wavefunctionfile, "%le%le%le\n", &pread, &wread, &wfread)){
		 gsl_vector_set(p, line, pread);
		 gsl_vector_set(w, line, wread);
		 gsl_vector_set(wf, line, wfread);
		 //printf("sucessfully read\n");
	 }
	 }
}

double formfactor(double q, gsl_vector *p,gsl_vector *wf, gsl_interp* interpolater, gsl_interp_accel* accelerator ){
	double x=0.5;
	double p_p=3;
	double p_p_z=0;
	double F=0;
	double p_2=0;
	int l_z=2, l=2;
	
	p_p_z=p_p*x;
	p_2=sqrt (p_p*p_p-p_p_z*q+q*q/4.);
	F=p_p*p_p*gsl_sf_legendre_sphPlm (l, l_z, x)*gsl_sf_legendre_sphPlm (l, l_z, (p_p_z-0.5*q)/p_2)
				*gsl_interp_eval (interpolater, p->data, wf->data, p_p, accelerator)*gsl_interp_eval (interpolater, p->data, wf->data, p_2, accelerator);
	return F;
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
	//gsl_vector_mul (norm, w);
	for(int i=0;i<norm->size;i++){
		result+=gsl_vector_get (norm, i);
	}
	
	
	printf ("norm=%e, F=%e\n",result,formfactor (0, p, wf, interpolate_test, acc_test));
	
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


