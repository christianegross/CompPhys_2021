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

void readinwavefunction(FILE * wavefunctionfile, gsl_vector *p, gsl_vector *w, gsl_vector* wf){
	/**
	 * @note cannot use built-in function because there are lines with other data, and several columns
	 * */
	 double pread, wread, wfread, empty;
	 rewind(wavefunctionfile);
	 fscanf(wavefunctionfile, "# Lam=   %le  A=  %le  mb=   %le C0=  %le \n# p [fm-1]       w             wf [fm3/2]\n", &empty, &empty, &empty, &empty);
	 for (int line=0; line<p->size; line+=1){
		 if(3== fscanf(wavefunctionfile, "%le%le%le\n", &pread, &wread, &wfread)){
		 gsl_vector_set(p, line, pread);
		 gsl_vector_set(w, line, wread);
		 gsl_vector_set(wf, line, wfread);
		 //printf("sucessfully read\n");
	 }
	 }
}

int main(int argc, char **argv){
	int datasetsize=60;
	FILE *wavefunctionfile=fopen("data/wavefunctions/wf-obe-lam=1200.00.dat", "r");
	FILE *testfile=fopen("data/test.dat", "w");
	gsl_vector *p=gsl_vector_alloc(datasetsize);
	gsl_vector *w=gsl_vector_alloc(datasetsize);
	gsl_vector *wf=gsl_vector_alloc(datasetsize);
	gsl_vector *norm=gsl_vector_alloc(datasetsize);
	double result=0;
	readinwavefunction(wavefunctionfile, p, w, wf);
	gsl_vector_memcpy (norm, wf);
	gsl_vector_mul (norm, wf);
	gsl_vector_mul (norm, w);
	for(int i=0;i<norm->size;i++){
		result+=gsl_vector_get (norm, i);
	}
	
	gsl_vector_fprintf(testfile, p, "%e");
	
	gsl_interp* interpolate_test=gsl_interp_alloc(gsl_interp_cspline ,datasetsize);
	gsl_interp_accel* acc_test=gsl_interp_accel_alloc ();
	gsl_interp_init (interpolate_test, p->data, wf->data, datasetsize);
	printf ("%e,%e\n",result,gsl_interp_eval (interpolate_test, p->data, wf->data, 30, acc_test));
	fclose(wavefunctionfile);
	fclose(testfile);
	gsl_vector_free(p);
	gsl_vector_free(w);
	gsl_vector_free(wf);
	gsl_vector_free(norm);
	gsl_interp_free (interpolate_test);
	
	return 0;
}


