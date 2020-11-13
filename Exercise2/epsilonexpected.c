//Christiane Groß, Nico Dichter
//Computational Physics WS2021
//Exercise sheet 2 - 2D Ising model
//additional programm to determine the literature result for the average energy per site using the math.h and gsl libraries

#include <stdio.h>
#include <stdlib.h>
#include "math.h"//exp-Function
#include <gsl/gsl_sf.h>//Kcomp
#include <gsl/gsl_math.h>//M_PI
#include <gsl/gsl_monte_vegas.h>//Evaluation of Kcomp

double epsilon_expected(double J){
	double epsilon=gsl_sf_ellint_Kcomp(4*sinh(2*J)*sinh(2*J)/cosh(2*J)/cosh(2*J)/cosh(2*J)/cosh(2*J),GSL_VEGAS_MODE_IMPORTANCE);
	epsilon*=2/M_PI*(2*tanh(2*J)*tanh(2*J)-1);
	epsilon+=1;
	epsilon*=-J*sinh(2*J)/cosh(2*J);
	return epsilon;
}

int main(int argc, char **argv){
	double J;
	FILE *epsilonvalues=fopen("epsilonexpected.txt", "w");
	for (int i=25; i<=200; i+=1){
		J=i/100.0;
		fprintf(epsilonvalues, "%f\t%f\n", J, epsilon_expected(J));
	}
	fclose(epsilonvalues);
}
