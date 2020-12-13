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
	FILE *wavefunctionfile=fopen("data/wavefunctions/wf-obe-lam=0300.00.dat", "r");
	FILE *testfile=fopen("data/test.dat", "w");
	gsl_vector *p=gsl_vector_alloc(60);
	gsl_vector *w=gsl_vector_alloc(60);
	gsl_vector *wf=gsl_vector_alloc(60);
	readinwavefunction(wavefunctionfile, p, w, wf);
	gsl_vector_fprintf(testfile, p, "%e");
	fclose(wavefunctionfile);
	fclose(testfile);
	gsl_vector_free(p);
	gsl_vector_free(w);
	gsl_vector_free(wf);
	
	
	return 0;
}


