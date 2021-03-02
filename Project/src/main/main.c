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
#include "auxiliary.h"					//own functions, e.g. trace, det,..., analysis by binning/bootstrapping...

#define GSL_COMPLEX_NAN (gsl_complex_rect(GSL_NAN,0.0))

/** 
 * @brief returns a SU(2) matrix in matrix, generates three free paramters needed, normalizes first column {z1, z2} to 1, sets second column as {-z2*, z1*} */
void generatesu2(gsl_matrix_complex * matrix, double epsilon, gsl_rng * generator){
	//error handling
	if(matrix->size1!=2||matrix->size2!=2){fprintf(stderr, "argument has to be 2x2 matrix!\n");return;}
	//generate random numbers
	double e,f,g, det;
	gsl_vector_complex_view columnzero;
	e=2.0*gsl_rng_uniform(generator)-1.0;
	f=2.0*gsl_rng_uniform(generator)-1.0;
	g=2.0*gsl_rng_uniform(generator)-1.0;
	//set matrix 0 column of 1+i*epsilon*H
	gsl_matrix_complex_set(matrix, 0, 0, gsl_complex_rect(1, epsilon*e));
	gsl_matrix_complex_set(matrix, 1, 0, gsl_complex_rect(epsilon*g, epsilon*f));
	//normalize vector zero to make sure determinant is one
	//rescale by 1/norm
	columnzero=gsl_matrix_complex_column(matrix, 0);
	det=gsl_blas_dznrm2(&columnzero.vector);		
	gsl_blas_zdscal(1.0/det, &columnzero.vector);
	//set up rest of matrix as 11=(00)*, 01=-(10)*
	gsl_matrix_complex_set(matrix, 1,1,gsl_complex_conjugate(gsl_matrix_complex_get(matrix, 0,0)));
	gsl_matrix_complex_set(matrix, 0,1,gsl_complex_mul_real(gsl_complex_conjugate(gsl_matrix_complex_get(matrix, 1,0)), -1.0));	
}

void generatesu3(gsl_matrix_complex * matrix, double epsilon, gsl_rng * generator){
	//error handling
	if(matrix->size1!=3||matrix->size2!=3){fprintf(stderr, "argument has to be 3x3 matrix!\n");return;}
	//generate random numbers
	double e,f,g,h,j,k,l,m,r, norm;
	gsl_vector_complex_view columnzero,columnone,columntwo;
	gsl_complex complexscalarproduct=GSL_COMPLEX_ZERO, complexscalarproduct2=GSL_COMPLEX_ZERO;
	e=2.0*gsl_rng_uniform(generator)-1.0;
	f=2.0*gsl_rng_uniform(generator)-1.0;
	g=2.0*gsl_rng_uniform(generator)-1.0;
	h=2.0*gsl_rng_uniform(generator)-1.0;
	j=2.0*gsl_rng_uniform(generator)-1.0;
	k=2.0*gsl_rng_uniform(generator)-1.0;
	l=2.0*gsl_rng_uniform(generator)-1.0;
	m=2.0*gsl_rng_uniform(generator)-1.0;
	r=2.0*gsl_rng_uniform(generator)-1.0;
	//set matrix 0 column of 1+i*epsilon*H
	gsl_matrix_complex_set(matrix, 0, 0, gsl_complex_rect(1, epsilon*e));
	gsl_matrix_complex_set(matrix, 1, 0, gsl_complex_rect(-epsilon*g, epsilon*f));
	gsl_matrix_complex_set(matrix, 2, 0, gsl_complex_rect(-epsilon*j, epsilon*h));
	//set matrix 1 column of 1+i*epsilon*H
	gsl_matrix_complex_set(matrix, 0, 1, gsl_complex_rect(epsilon*g, epsilon*f));
	gsl_matrix_complex_set(matrix, 1, 1, gsl_complex_rect(1, epsilon*k));
	gsl_matrix_complex_set(matrix, 2, 1, gsl_complex_rect(-epsilon*m, epsilon*l));
	//set matrix 2 column of 1+i*epsilon*H
	gsl_matrix_complex_set(matrix, 0, 2, gsl_complex_rect(epsilon*j, epsilon*h));
	gsl_matrix_complex_set(matrix, 1, 2, gsl_complex_rect(epsilon*m, epsilon*l));
	gsl_matrix_complex_set(matrix, 2, 2, gsl_complex_rect(1, epsilon*r));
	//normalize vector zero to make sure determinant is one
	//rescale by 1/norm
	columnzero=gsl_matrix_complex_column(matrix, 0);
	norm=gsl_blas_dznrm2(&columnzero.vector);
	gsl_blas_zdscal(1.0/norm, &columnzero.vector);

	columnone=gsl_matrix_complex_column(matrix, 1);
	gsl_blas_zdotc (&columnzero.vector, &columnone.vector, &complexscalarproduct);
	gsl_blas_zaxpy (gsl_complex_mul_real(complexscalarproduct, -1.0), &columnzero.vector, &columnone.vector);
	norm=gsl_blas_dznrm2(&columnone.vector);
	gsl_blas_zdscal(1.0/norm, &columnone.vector);
	columntwo=gsl_matrix_complex_column(matrix, 2);
	//crossproduct(&columnzero.vector,&columnone.vector,&columntwo.vector);
	gsl_blas_zdotc (&columnzero.vector, &columntwo.vector, &complexscalarproduct);
	gsl_blas_zdotc (&columnone.vector, &columntwo.vector, &complexscalarproduct2);
	gsl_blas_zaxpy (gsl_complex_mul_real(complexscalarproduct, -1.0), &columnzero.vector, &columntwo.vector);
	gsl_blas_zaxpy (gsl_complex_mul_real(complexscalarproduct2, -1.0), &columnone.vector, &columntwo.vector);
	norm=gsl_blas_dznrm2(&columntwo.vector);
	gsl_blas_zdscal(1.0/norm, &columntwo.vector);
	//~ printf("%f\t%f\n", GSL_REAL(complexscalarproduct), GSL_REAL(complexscalarproduct2));
}

/** @brief returns the difference in action before and after change of one link
 * @note prototype: returns a random number between -1 and 1 to test MH-setup*/
inline double deltaS1(gsl_rng * generator){
return 2.0*gsl_rng_uniform(generator)-1.0;
}

/** @brief returns the difference in action before and after change of one link
 * @note takes sum over unchanged matrices as argument, same computation results can be used for several hits, precalculation in loop 
 * @note calculates (newmatrix-oldmatrix)*sum over contributions from different directions 
 * **/
inline double deltaS(gsl_matrix_complex *deltacontribution, gsl_matrix_complex* newmatrix, gsl_matrix_complex *oldmatrix, gsl_matrix_complex *helpone, gsl_matrix_complex *helptwo,int dim){
	gsl_matrix_complex_memcpy(helpone, newmatrix);
	gsl_matrix_complex_sub(helpone, oldmatrix);
	gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(1,0), helpone, deltacontribution, gsl_complex_rect(0,0), helptwo);
	//possibly wrong numerical prefactor
	return GSL_REAL(trace(helptwo))/dim;
}

/**@fn calculateGamma(int * neighbour, int counter, int mu, int size, gsl_matrix_complex *plaquettecontribution, gsl_matrix_complex *deltacontribution, gsl_matrix_complex *helpone, gsl_matrix_complex *helptwo, gsl_matrix_complex **matrixarray)
 * @brief calculates the contribution to all plaquettes U_mu(x) is involved in
 * @note U-mu(x) is involved in six different plaquettes. calculating all of them every time for the change in the action would take too long
 * **/ 
void calculateGamma(int * neighbour, int counter, int mu, int size, gsl_matrix_complex *plaquettecontribution, gsl_matrix_complex *deltacontribution, gsl_matrix_complex *helpone, gsl_matrix_complex *helptwo, gsl_matrix_complex **matrixarray){
	/**
	*forward U_nu(x+amu)*U^dagger_mu(x+anu)*U^dagger_nu(x)
	*backward U^dagger_nu(x+amu-anu)*U^dagger_mu(x-anu)*U_nu(x-anu) -> gives P_munu(x-anu) when multiplied by U_mu(x)
	* */
	for (int nu =0;nu<4;nu+=1){ if(mu!=nu){
		//forward plaquette
		gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, GSL_COMPLEX_ONE, matrixarray[counter+neighbour[2*mu]-mu+nu], matrixarray[counter+neighbour[2*nu]], GSL_COMPLEX_ZERO, helpone);
		gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, GSL_COMPLEX_ONE, helpone, matrixarray[counter-mu+nu] , GSL_COMPLEX_ZERO, helptwo);
		gsl_matrix_complex_add(deltacontribution, helptwo);
		if(mu>nu){gsl_matrix_complex_add(plaquettecontribution, helptwo);}
		//backward plaquette
		gsl_blas_zgemm(CblasConjTrans, CblasConjTrans, GSL_COMPLEX_ONE, matrixarray[counter+neighbour[2*mu]+neighbour[2*nu+1]-mu+nu], matrixarray[counter+neighbour[2*nu+1]], GSL_COMPLEX_ZERO, helpone);
		gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, helpone, matrixarray[counter+neighbour[2*nu+1]-mu+nu] , GSL_COMPLEX_ZERO, helptwo);
		gsl_matrix_complex_add(deltacontribution, helptwo);
	}}
}

/**
 * @fn calculateplaquette(gsl_matrix_complex ** matrixarray, int counter, int* neighbour, int mu, int nu, gsl_matrix_complex *helpone, gsl_matrix_complex *helptwo, gsl_matrix_complex *helpthree)
 * @brief calculates the plaquette at the position counter in the direction mu and nu
 * **/
double calculateplaquette(gsl_matrix_complex ** matrixarray, int counter, int* neighbour, int mu, int nu, gsl_matrix_complex *helpone, gsl_matrix_complex *helptwo, gsl_matrix_complex *helpthree,int dim){
	gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, matrixarray[counter+mu], matrixarray[counter+neighbour[2*mu]+nu], GSL_COMPLEX_ZERO, helpone);
	gsl_blas_zgemm(CblasConjTrans, CblasConjTrans, GSL_COMPLEX_ONE, matrixarray[counter+neighbour[2*nu]+mu], matrixarray[counter+nu], GSL_COMPLEX_ZERO, helptwo);
	gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, helpone, helptwo, GSL_COMPLEX_ZERO, helpthree);
	return GSL_REAL(trace(helpthree))/((double)dim);
}

double calculatewilsonloop(gsl_matrix_complex ** matrixarray, gsl_matrix_complex ** helparray, int x, int y, int z, int t, int r1, int r2, int r3, int tdistance, int size, int dim){
/** see two problems:
 * 1. What order to use for x,y,z contributions? average over all possible permutations (xyz, xzy, yxz, yzx, zxy, zyx) oder always do order xyz?
 * 2. How to implement boundary conditions? How to see, where to insert boundaries?
 * ->probably made easier by using x, y, z, t seperately, maybe not put counter as an argument, but individual values for x,y,z,t?
 * **/

 /** first calculation: go in x-direction for r1 steps, y for r2 steps, z for r3 steps, t for tdistance steps, and back with xdagger, ydagger, zdagger, tdagger
  first multiplication: U_x(pos)U_x(pos+ax)**/
	gsl_matrix_complex* m1,*m2, *intermediate;
	int len1init=0,len2init=0,len3init=0,lentinit=0;
	if(r1>=1){
		gsl_matrix_complex_memcpy (helparray[0], matrixarray[x*size*size*size*4+y*size*size*4+z*size*4+t*4+3]);
		len1init=1;
	}else if(r2>=1){
		gsl_matrix_complex_memcpy (helparray[0], matrixarray[x*size*size*size*4+y*size*size*4+z*size*4+t*4+2]);
		len2init=1;
	}else if(r3>=1){
		gsl_matrix_complex_memcpy (helparray[0], matrixarray[x*size*size*size*4+y*size*size*4+z*size*4+t*4+1]);
		len3init=1;
	}else{
		gsl_matrix_complex_memcpy (helparray[0], matrixarray[x*size*size*size*4+y*size*size*4+z*size*4+t*4+0]);
		lentinit=1;
	}
	m1=helparray[0];
	m2=helparray[1];
	//forwards:
	for (int len1=len1init;len1<r1;len1+=1){
		gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, m1, matrixarray[((x+len1)%size)*size*size*size*4+y*size*size*4+z*size*4+t*4+3], GSL_COMPLEX_ZERO, m2);
		intermediate=m1;
		m1=m2;
		m2=intermediate;
	}
	for (int len2=len2init;len2<r2;len2+=1){
		gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, m1, matrixarray[((x+r1)%size)*size*size*size*4+((y+len2)%size)*size*size*4+z*size*4+t*4+2], GSL_COMPLEX_ZERO, m2);
		intermediate=m1;
		m1=m2;
		m2=intermediate;
	}
	for (int len3=len3init;len3<r3;len3+=1){
		gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, m1, matrixarray[((x+r1)%size)*size*size*size*4+((y+r2)%size)*size*size*4+((z+len3)%size)*size*4+t*4+1], GSL_COMPLEX_ZERO, m2);
		intermediate=m1;
		m1=m2;
		m2=intermediate;
	}
	for (int lent=lentinit;lent<tdistance;lent+=1){
		gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, m1, matrixarray[((x+r1)%size)*size*size*size*4+((y+r2)%size)*size*size*4+((z+r3)%size)*size*4+((t+lent)%size)*4+0], GSL_COMPLEX_ZERO, m2);
		intermediate=m1;
		m1=m2;
		m2=intermediate;
	}
	//backwards:
	for (int len1=1;len1<=r1;len1+=1){
		gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, GSL_COMPLEX_ONE, m1, matrixarray[((x+r1-len1)%size)*size*size*size*4+((y+r2)%size)*size*size*4+((z+r3)%size)*size*4+((t+tdistance)%size)*4+3], GSL_COMPLEX_ZERO, m2);
		intermediate=m1;
		m1=m2;
		m2=intermediate;
	}
	for (int len2=1;len2<=r2;len2+=1){
		gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, GSL_COMPLEX_ONE, m1, matrixarray[x*size*size*size*4+((y+r2-len2)%size)*size*size*4+((z+r3)%size)*size*4+((t+tdistance)%size)*4+2], GSL_COMPLEX_ZERO, m2);
		intermediate=m1;
		m1=m2;
		m2=intermediate;
	}
	for (int len3=1;len3<=r3;len3+=1){
		gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, GSL_COMPLEX_ONE, m1, matrixarray[x*size*size*size*4+y*size*size*4+((z+r3-len3)%size)*size*4+((t+tdistance)%size)*4+1], GSL_COMPLEX_ZERO, m2);
		intermediate=m1;
		m1=m2;
		m2=intermediate;
	}
	for (int lent=1;lent<=tdistance;lent+=1){
		gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, GSL_COMPLEX_ONE, m1, matrixarray[x*size*size*size*4+y*size*size*4+z*size*4+((t+tdistance-lent)%size)*4+0], GSL_COMPLEX_ZERO, m2);
		intermediate=m1;
		m1=m2;
		m2=intermediate;
	}
/** possible additional runs: x-r1-z-r2-y-r3, ...**/
 return GSL_REAL(trace(m1))/dim;
 }
	

int main(int argc, char **argv){
	//set up constants, matrices, generator
	int dim=3; //switches between SU2 and SU3
	double epsilon=0.2;
	int hotstart=1;
	int size=8;
	double beta=2.3;
	int numberofthermalizations=100;

	int numberofmeasurements=2048; //=pow(2, 13)
	

	gsl_matrix_complex *newmatrix=gsl_matrix_complex_alloc(dim,dim);
	gsl_matrix_complex *multiplier=gsl_matrix_complex_alloc(dim,dim);
	gsl_matrix_complex *plaquettecontribution=gsl_matrix_complex_alloc(dim,dim);
	gsl_matrix_complex *deltacontribution=gsl_matrix_complex_alloc(dim,dim);
	//~ gsl_matrix_complex *helpone=gsl_matrix_complex_alloc(2,2);
	//~ gsl_matrix_complex *helpthree=gsl_matrix_complex_alloc(2,2);
	//~ gsl_matrix_complex *helptwo=gsl_matrix_complex_alloc(2,2);
	int seed=2;//use fixed seed: result should be exactly reproduced using the same seed
	gsl_rng *generator;
	generator=gsl_rng_alloc(gsl_rng_mt19937);//use mersenne-twister
	gsl_rng_set(generator, seed);
	
	
	/**set up array to store link matrices, put SU(2) matrix in every link
	every matrix can be initialised to unity or a random SU(2) matrix**/
	gsl_matrix_complex* matrixarray[size*size*size*size*4];
	for (int i=0;i<size*size*size*size*4;i+=1){
		matrixarray[i]=gsl_matrix_complex_alloc(dim,dim);
		//~ generatesu2(matrixarray[i], epsilon, generator);  //set up as random or unity to test different configurations
		if(hotstart){
			if(dim==2){
				generatesu2(matrixarray[i], epsilon, generator);
			}else if (dim==3){
		 		generatesu3(matrixarray[i], epsilon, generator);
			}else {fprintf (stderr, "Invalid value for dim.");return -1;}
		}else{
			settounity(matrixarray[i]);
		}
	}
	/** set up array with helpmatrices to store intermediate results
	 * size: max(3, 2*(r1max+r2max+r3max+tmax)) **/
	 gsl_matrix_complex* helparray[30];
	 for (int i=0;i<30;i+=1){
		 helparray[i]=gsl_matrix_complex_alloc(dim,dim);
	 }
	
	/** set up vectors and streams for analysis of results **/
	double mean_plaquette, var_plaquette;
	gsl_vector * plaquette=gsl_vector_alloc(numberofmeasurements);
	gsl_vector * binned_plaquette_mem=gsl_vector_alloc(numberofmeasurements);
	gsl_vector_view binned_plaquette;
	//~ gsl_vector * plaquette_correlation=gsl_vector_alloc(numberofmeasurements/4);
	gsl_vector * plaquette_correlation_binned=gsl_vector_alloc(numberofmeasurements/32);
	FILE * plaquette_data=fopen("data/plaquettedata.dat", "w");
	FILE * plaquette_autocorrelation=fopen("data/plaquetteautocorrelation.dat", "w");
	FILE * plaquette_analysis=fopen("data/plaquette.dat", "w");
	
	//test of MH: go through lattice, perform 10 accept/reject steps at every link, measure acceptance rate
	/**thermalizations**/
	/** counter defined as position with direction**/
	/**neighbours: need two neighbours for every direction, one forward and one backward
	 * 0: t-forward	1:t-backward
	 * 2: z-forward	3:z-backward
	 * 4: y-forward	5:y-backward
	 * 6: x-forward	7:x-backward
	 * forward: 2*mu, backward: 2*mu+1
	 * relative to position
	 * mu			0	1	2	3
	 * direction	t	z	y	x
	 * **/
	int counter, acceptance;
	int neighbour[8]; //for implementing (periodic) boundary conditions
	double plaquetteexpectation, plaquetteafter, wilsonexpectation[4], wilsonexpectationmedium;
	for(int runs=0;runs<numberofthermalizations;runs+=1){
		acceptance=0;
		plaquetteexpectation=0;
		wilsonexpectation[0]=0;
		wilsonexpectation[1]=0;
		wilsonexpectation[2]=0;
		wilsonexpectation[3]=0;
		wilsonexpectationmedium=0;
		plaquetteafter=0;
		for (int x=0;x<size;x+=1){
			neighbour[6]=(x==size-1)?-(size-1)*pow(size, 3)*4:pow(size,3)*4;
			neighbour[7]=(x==0)?(size-1)*pow(size, 3)*4:-pow(size,3)*4;
			for (int y=0;y<size;y+=1){
				neighbour[4]=(y==size-1)?-(size-1)*pow(size, 2)*4:pow(size,2)*4;
				neighbour[5]=(y==0)?(size-1)*pow(size, 2)*4:-pow(size,2)*4;
				for (int z=0;z<size;z+=1){
					neighbour[2]=(z==size-1)?-(size-1)*size*4:size*4;
					neighbour[3]=(z==0)?(size-1)*size*4:-size*4;
					for (int t=0;t<size;t+=1){
						neighbour[0]=(t==size-1)?-(size-1)*4:4;
						neighbour[1]=(t==0)?(size-1)*4:-4;
						for (int mu=0;mu<4;mu+=1){
							
							counter=x*size*size*size*4+y*size*size*4+z*size*4+t*4+mu;
							settozero(plaquettecontribution); settozero(deltacontribution);
							calculateGamma(neighbour, counter, mu, size, plaquettecontribution, deltacontribution, helparray[0], helparray[1], matrixarray);
							for (int attempts=0;attempts<10;attempts+=1){
								if(dim==2){
									generatesu2(multiplier, epsilon, generator);
								}else if (dim==3){
								 	generatesu3(multiplier, epsilon, generator);
								}else {fprintf (stderr, "Invalid value for dim.");return -1;}

								gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(1,0), multiplier, matrixarray[counter], gsl_complex_rect(0,0), newmatrix);
								if (exp(beta*deltaS(deltacontribution, newmatrix, matrixarray[counter], helparray[0], helparray[1],dim))>gsl_rng_uniform(generator)){
									acceptance+=1;
									gsl_matrix_complex_memcpy(matrixarray[counter], newmatrix); 
								}
								/**what to use for plaquette: sum over (mu>nu) or 1/2*sum over (mu)? Also include backwards plaquettes?**/
								gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(1,0), matrixarray[counter], plaquettecontribution, gsl_complex_rect(0,0), helparray[0]);
								plaquetteexpectation+=GSL_REAL(trace(helparray[0]));
								
							}
						}

						wilsonexpectation[0]+=calculatewilsonloop(matrixarray, helparray, x, y, z, t, 4, 0,0,1,size,dim);
						wilsonexpectation[1]+=calculatewilsonloop(matrixarray, helparray, x, y, z, t, 4, 0,0,2,size,dim);
						wilsonexpectation[2]+=calculatewilsonloop(matrixarray, helparray, x, y, z, t, 4, 0,0,3,size,dim);
						wilsonexpectation[3]+=calculatewilsonloop(matrixarray, helparray, x, y, z, t, 4, 0,0,4,size,dim);
					}
				}
			}
		}	
		/**measure plaquette after one sweep is complete**/
		/** counter defined as position without direction**/
		/**for (int x=0;x<size;x+=1){
			neighbour[6]=(x==size-1)?-(size-1)*pow(size, 3)*4:pow(size,3)*4;
			neighbour[7]=(x==0)?(size-1)*pow(size, 3)*4:-pow(size,3)*4;
			for (int y=0;y<size;y+=1){
				neighbour[4]=(y==size-1)?-(size-1)*pow(size, 2)*4:pow(size,2)*4;
				neighbour[5]=(y==0)?(size-1)*pow(size, 2)*4:-pow(size,2)*4;
				for (int z=0;z<size;z+=1){
					neighbour[2]=(z==size-1)?-(size-1)*size*4:size*4;
					neighbour[3]=(z==0)?(size-1)*size*4:-size*4;
					for (int t=0;t<size;t+=1){
						neighbour[0]=(t==size-1)?-(size-1)*4:4;
						neighbour[1]=(t==0)?(size-1)*4:-4;
						counter=x*size*size*size*4+y*size*size*4+z*size*4+t*4;
						for (int mu=0;mu<4;mu+=1){
							for (int nu=mu+1;nu<4;nu+=1){
								plaquetteafter+=calculateplaquette(matrixarray, counter, neighbour, mu, nu, helparray[0], helparray[1], helparray[2],dim);
							}
						}
					}
				}
			}
		}**/
		/**where to measure plaquette? measure directly after one link is switched, and get contributions from links that are changed in the next step, or loop over entire lattice after every sweep and take longer? 
		Or maybe not longer, since matrix links are looked at ten times per sweep? Maybe look during sweep, but only after ten attempts have ben made?**/
		/**factors for plaquette and acceptance rate: both have to be 1.0 when filled with unity matrices and epsilon=0**/
		fprintf(stdout, "%d\tacc=%f\tplaq=%f\tw1=%f\tw2=%f\tw3=%f\tw4=%f\n", runs, (double)acceptance/((double)10*size*size*size*size*4),plaquetteexpectation/((double)10*size*size*size*size*4*3*dim/2),
				/**plaquetteafter/((double)size*size*size*size*4*3*0.5),**/ wilsonexpectation[0]/((double)size*size*size*size),wilsonexpectation[1]/((double)size*size*size*size),wilsonexpectation[2]/((double)size*size*size*size),wilsonexpectation[3]/((double)size*size*size*size));
	}
	
	
	//~ /**measurements **/
	//~ for(int runs=0;runs<numberofmeasurements;runs+=1){
		//~ acceptance=0;
		//~ plaquetteexpectation=0;
/**copy from thermalization-> still to be done, but set up functions for gamma, plaquette, wilson-loop, etc first**/
		//~ gsl_vector_set(plaquette, runs, plaquetteexpectation/((double)10*size*size*size*size*4*3*0.5));
	//~ }
	
	/** analysis by binning and bootstrapping**/
	//~ fprintf(plaquette_analysis, "bin\t<plaq>\tvar(plaq)\n", binsize, mean_plaquette, var_plaquette);
	//~ for (int binsize=1;binsize<33;binsize*=2){
		//~ binned_plaquette=gsl_vector_subvector(binned_plaquette_mem, 0, plaquette->size/binsize);
		//~ binning(plaquette, &binned_plaquette.vector, binsize);
		//~ bootstrap(&binned_plaquette.vector, generator, 2, &mean_plaquette, &var_plaquette);
		//~ autocorrelation(&binned_plaquette.vector, plaquette_correlation_binned, mean_plaquette);
		
		//~ fprintf(plaquette_data, "\nbinsize %d\n", binsize);
		//~ fprintf(plaquette_autocorrelation, "\nbinsize %d\n", binsize);
		//~ gsl_vector_fprintf(plaquette_data, &binned_plaquette.vector, "%f");
		//~ gsl_vector_fprintf(plaquette_autocorrelation, plaquette_correlation_binned, "%f");
		
		//~ fprintf(plaquette_analysis, "%.2d\t%f\t%f\n", binsize, mean_plaquette, var_plaquette);
	//~ }
	
	//Test SU3:
	gsl_matrix_complex *su3=gsl_matrix_complex_alloc(3,3);
	gsl_vector_complex_view columnzero,columnone,columntwo;
	double norm;
	gsl_complex dot01, dot02, dot12;
	for(int i=0;i<50;i++){
		generatesu3 (su3, i/50., generator);
		//gsl_matrix_complex_fprintf (stdout, su3, "%e");
		columnzero=gsl_matrix_complex_column(su3, 0);
		columnone=gsl_matrix_complex_column(su3, 1);
		columntwo=gsl_matrix_complex_column(su3, 2);
		gsl_blas_zdotc (&columnzero.vector, &columnone.vector, &dot01);
		gsl_blas_zdotc (&columnzero.vector, &columntwo.vector, &dot02);
		gsl_blas_zdotc (&columnone.vector, &columntwo.vector, &dot12);
		norm=gsl_blas_dznrm2(&columntwo.vector);
		printf ("eps: %e\tnorm: %e\tdot01: %e\tdot02: %e\tdot12: %e\n",i/50.,norm,gsl_complex_abs (dot01),gsl_complex_abs (dot02),gsl_complex_abs (dot12));
		gsl_matrix_complex_set_zero (su3);
	}





	//cleanup
	gsl_rng_free(generator);
	gsl_matrix_complex_free(multiplier);
	gsl_matrix_complex_free(newmatrix);
	gsl_matrix_complex_free(deltacontribution);
	gsl_matrix_complex_free(plaquettecontribution);
	//~ gsl_matrix_complex_free(helpone);
	//~ gsl_matrix_complex_free(helptwo);
	//~ gsl_matrix_complex_free(helpthree);
	for (int i=0;i<size*size*size*size*4;i+=1){
		gsl_matrix_complex_free(matrixarray[i]);
	}
	for (int i=0;i<30;i+=1){
		gsl_matrix_complex_free(helparray[i]);
	 }
	gsl_vector_free(plaquette);
	gsl_vector_free(binned_plaquette_mem);
	//~ gsl_vector_free(plaquette_correlation);
	gsl_vector_free(plaquette_correlation_binned);
	fclose(plaquette_data);
	fclose(plaquette_autocorrelation);
	fclose(plaquette_analysis);
	return 0;
}


