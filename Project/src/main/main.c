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


#define GSL_COMPLEX_NAN (gsl_complex_rect(GSL_NAN,0.0))

/**
 * @brief returns trace of complex matrix
 * */
gsl_complex trace(gsl_matrix_complex *matrix){
	if (matrix->size1!=matrix->size2){fprintf(stderr, "argument has to be square matrix!\n");return GSL_COMPLEX_NAN;}
	gsl_complex number=GSL_COMPLEX_ZERO;
	for (int i=0; i<matrix->size1; i+=1){
		number=gsl_complex_add(number, gsl_matrix_complex_get(matrix, i, i));
	}
	return number;
}


/**
 * @brief returns det of complex matrix
 * */
gsl_complex det(gsl_matrix_complex *matrix){
	if (matrix->size1!=matrix->size2||matrix->size1!=3){fprintf(stderr, "argument has to be 3x3 matrix!\n");return GSL_COMPLEX_NAN;}
	gsl_complex number=GSL_COMPLEX_ZERO;
	number=gsl_complex_add(number,gsl_complex_mul(gsl_complex_mul(gsl_matrix_complex_get(matrix,0,0),gsl_matrix_complex_get(matrix,1,1)),gsl_matrix_complex_get(matrix,2,2)));
	number=gsl_complex_add(number,gsl_complex_mul(gsl_complex_mul(gsl_matrix_complex_get(matrix,0,1),gsl_matrix_complex_get(matrix,1,2)),gsl_matrix_complex_get(matrix,2,0)));
	number=gsl_complex_add(number,gsl_complex_mul(gsl_complex_mul(gsl_matrix_complex_get(matrix,0,2),gsl_matrix_complex_get(matrix,1,0)),gsl_matrix_complex_get(matrix,2,1)));

	number=gsl_complex_sub(number,gsl_complex_mul(gsl_complex_mul(gsl_matrix_complex_get(matrix,0,2),gsl_matrix_complex_get(matrix,1,1)),gsl_matrix_complex_get(matrix,2,0)));
	number=gsl_complex_sub(number,gsl_complex_mul(gsl_complex_mul(gsl_matrix_complex_get(matrix,0,1),gsl_matrix_complex_get(matrix,1,0)),gsl_matrix_complex_get(matrix,2,2)));
	number=gsl_complex_sub(number,gsl_complex_mul(gsl_complex_mul(gsl_matrix_complex_get(matrix,0,0),gsl_matrix_complex_get(matrix,1,2)),gsl_matrix_complex_get(matrix,2,1)));
	return number;
}

/** 
 * @brief returns the determinant of a complex 2x2 matrix {{a,b},{c,d}} by returning ad-bc */
inline gsl_complex det2(gsl_matrix_complex *matrix){
	if (matrix->size1!=matrix->size2||matrix->size1!=2){fprintf(stderr, "argument has to be 2x2 matrix!\n");return GSL_COMPLEX_NAN;}
	return gsl_complex_sub(gsl_complex_mul(gsl_matrix_complex_get(matrix, 0,0), gsl_matrix_complex_get(matrix, 1,1)), gsl_complex_mul(gsl_matrix_complex_get(matrix, 1,0), gsl_matrix_complex_get(matrix, 0,1)));
	}

/**
 * @brief transforms square matrix into its conjugate transposed. 
 * @note loops over i, conjugates diagonal elements, saves ij, sets ij=(ji)*, sets ji=(ij, old)*
 * */
void conjugatetranspose(gsl_matrix_complex *matrix){
	if (matrix->size1!=matrix->size2){fprintf(stderr, "argument has to be square matrix!\n");return;}
	gsl_complex number=GSL_COMPLEX_ZERO;
	for (int i=0;i<matrix->size1;i+=1){
		gsl_matrix_complex_set(matrix, i,i, gsl_complex_conjugate(gsl_matrix_complex_get(matrix, i,i)));
		for(int j=i+1;j<matrix->size1;j+=1){
			number=gsl_matrix_complex_get(matrix, i, j);
			gsl_matrix_complex_set(matrix, i,j, gsl_complex_conjugate(gsl_matrix_complex_get(matrix, j,i)));
			gsl_matrix_complex_set(matrix, j,i, gsl_complex_conjugate(number));
		}
	}
}

/** sets all matrix elements to zero **/
void settozero(gsl_matrix_complex *matrix){
	if (matrix->size1!=matrix->size2){fprintf(stderr, "argument has to be square matrix!\n");return;}
	for (int i=0;i<matrix->size1;i+=1){
		gsl_matrix_complex_set(matrix, i,i, GSL_COMPLEX_ZERO);
		for(int j=i+1;j<matrix->size1;j+=1){
			gsl_matrix_complex_set(matrix, i,j, GSL_COMPLEX_ZERO);
			gsl_matrix_complex_set(matrix, j,i, GSL_COMPLEX_ZERO);
	
		}
	}
}

/** sets matrix to unity matrix **/
void settounity(gsl_matrix_complex *matrix){
	if (matrix->size1!=matrix->size2){fprintf(stderr, "argument has to be square matrix!\n");return;}
	for (int i=0;i<matrix->size1;i+=1){
		gsl_matrix_complex_set(matrix, i,i, GSL_COMPLEX_ONE);
		for(int j=i+1;j<matrix->size1;j+=1){
			gsl_matrix_complex_set(matrix, i,j, GSL_COMPLEX_ZERO);
			gsl_matrix_complex_set(matrix, j,i, GSL_COMPLEX_ZERO);
		}
	}
}

/**
 * @brief
 * @note
 */
void crossproduct(gsl_vector_complex* a,gsl_vector_complex* b,gsl_vector_complex* result){
	gsl_vector_complex_set (result, 0, gsl_complex_sub(gsl_complex_mul (gsl_vector_complex_get(a,1), gsl_vector_complex_get(b,2)),gsl_complex_mul (gsl_vector_complex_get(a,2), gsl_vector_complex_get(b,1))));
	gsl_vector_complex_set (result, 1, gsl_complex_sub(gsl_complex_mul (gsl_vector_complex_get(a,2), gsl_vector_complex_get(b,0)),gsl_complex_mul (gsl_vector_complex_get(a,0), gsl_vector_complex_get(b,2))));
	gsl_vector_complex_set (result, 2, gsl_complex_sub(gsl_complex_mul (gsl_vector_complex_get(a,0), gsl_vector_complex_get(b,1)),gsl_complex_mul (gsl_vector_complex_get(a,1), gsl_vector_complex_get(b,0))));
}

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
	double e,f,g,h,j,k,l,m, norm;
	gsl_vector_complex_view columnzero,columnone,columntwo;
	gsl_complex complexscalarproduct=GSL_COMPLEX_ZERO;
	e=2.0*gsl_rng_uniform(generator)-1.0;
	f=2.0*gsl_rng_uniform(generator)-1.0;
	g=2.0*gsl_rng_uniform(generator)-1.0;
	h=2.0*gsl_rng_uniform(generator)-1.0;
	j=2.0*gsl_rng_uniform(generator)-1.0;
	k=2.0*gsl_rng_uniform(generator)-1.0;
	l=2.0*gsl_rng_uniform(generator)-1.0;
	m=2.0*gsl_rng_uniform(generator)-1.0;
	//set matrix 0 column of 1+i*epsilon*H
	gsl_matrix_complex_set(matrix, 0, 0, gsl_complex_rect(1, epsilon*e));
	gsl_matrix_complex_set(matrix, 1, 0, gsl_complex_rect(-epsilon*g, epsilon*f));
	gsl_matrix_complex_set(matrix, 2, 0, gsl_complex_rect(-epsilon*j, epsilon*h));
	//set matrix 1 column of 1+i*epsilon*H
	gsl_matrix_complex_set(matrix, 0, 1, gsl_complex_rect(epsilon*g, epsilon*f));
	gsl_matrix_complex_set(matrix, 1, 1, gsl_complex_rect(1, epsilon*k));
	gsl_matrix_complex_set(matrix, 2, 1, gsl_complex_rect(-epsilon*m, epsilon*l));
	//normalize vector zero to make sure determinant is one
	//rescale by 1/norm
	columnzero=gsl_matrix_complex_column(matrix, 0);
	norm=gsl_blas_dznrm2(&columnzero.vector);
	gsl_blas_zdscal(1.0/norm, &columnzero.vector);

	columnone=gsl_matrix_complex_column(matrix, 1);
	gsl_blas_zdotu (&columnone.vector, &columnzero.vector, &complexscalarproduct);
	gsl_blas_zaxpy (gsl_complex_mul_real(complexscalarproduct, -1.0), &columnzero.vector, &columnone.vector);
	norm=gsl_blas_dznrm2(&columnone.vector);
	gsl_blas_zdscal(1.0/norm, &columnone.vector);
	//TODO: add cross product
	columntwo=gsl_matrix_complex_column(matrix, 2);
	crossproduct(&columnzero.vector,&columnone.vector,&columntwo.vector);
}

/** @brief returns the difference in action before and after change of one link
 * @note prototype: returns a random number between -1 and 1 to test MH-setup*/
inline double deltaS1(gsl_rng * generator){
return 2.0*gsl_rng_uniform(generator)-1.0;
}

/** @brief returns the difference in action before and after change of one link
 * @note takes sum over unchanged matrices as argument, same computation results can be used for several hits, precalculation in loop 
 * @note caqlculates (newmatrix-oldmatrix)*sum over contributions from different directions **/
inline double deltaS(gsl_matrix_complex *deltacontribution, gsl_matrix_complex* newmatrix, gsl_matrix_complex *oldmatrix, gsl_matrix_complex *helpone, gsl_matrix_complex *helptwo){
	gsl_matrix_complex_memcpy(helpone, newmatrix);
	gsl_matrix_complex_sub(helpone, oldmatrix);
	gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(1,0), helpone, deltacontribution, gsl_complex_rect(0,0), helptwo);
	return 0.5*GSL_REAL(trace(helptwo));
}

int main(int argc, char **argv){
	//~ double a,b,c,d,det;
	//set up constants, matrices, generator
	double epsilon=0.2;
	int size=8;
	//~ double e,f,g,h;
	gsl_complex complexproduct=GSL_COMPLEX_ZERO;
	gsl_matrix_complex *one=gsl_matrix_complex_alloc(2,2);
	gsl_matrix_complex *two=gsl_matrix_complex_alloc(2,2);
	gsl_matrix_complex *multiplier=gsl_matrix_complex_alloc(2,2);
	gsl_matrix_complex *result=gsl_matrix_complex_alloc(2,2);
	gsl_matrix_complex *newmatrix=gsl_matrix_complex_alloc(2,2);
	gsl_matrix_complex *plaquettecontribution=gsl_matrix_complex_alloc(2,2);
	gsl_matrix_complex *deltacontribution=gsl_matrix_complex_alloc(2,2);
	gsl_matrix_complex *helpone=gsl_matrix_complex_alloc(2,2);
	gsl_matrix_complex *helpthree=gsl_matrix_complex_alloc(2,2);
	gsl_matrix_complex *helptwo=gsl_matrix_complex_alloc(2,2);
	gsl_matrix_complex *three=gsl_matrix_complex_alloc(3,3);
	//~ gsl_matrix_complex *four=gsl_matrix_complex_alloc(2,2);
	gsl_matrix_complex *multiplied=gsl_matrix_complex_alloc(2,2);	
	gsl_matrix_complex *multiplied3=gsl_matrix_complex_alloc(3,3);
	gsl_matrix_complex *five=gsl_matrix_complex_calloc(3,2);	
	//~ gsl_vector_complex_view columnzero, columnone;
	int seed=2;//use fixed seed: result should be exactly reproduced using the same seed
	gsl_rng *generator;
	generator=gsl_rng_alloc(gsl_rng_mt19937);//use mersenne-twister
	gsl_rng_set(generator, seed);
	
	
	//set up array to store link matrices, put SU(2) matrix in every link
	gsl_matrix_complex* matrixarray[size*size*size*size*4];
	for (int i=0;i<size*size*size*size*4;i+=1){
		matrixarray[i]=gsl_matrix_complex_alloc(2,2);
		//~ generatesu2(matrixarray[i], epsilon, generator);  //set up as random or unity to test different configurations
		settounity(matrixarray[i]);
		complexproduct=det2(matrixarray[i]);
		//~ fprintf(stdout, "det of %d:\t%f+%fi\n", i, GSL_REAL(complexproduct), GSL_IMAG(complexproduct));
	}
	
	//test of MH: go through lattice, perform 10 accept/reject steps at every link, measure acceptance rate
	int counter;
	int neighbour[4]; //for implementing periodic boundary conditions
	double plaquetteexpectation;
	for(int runs=0;runs<100;runs+=1){
		int acceptance=0;
		plaquetteexpectation=0;
		for (int x=0;x<size;x+=1){
			neighbour[0]=(x==size-1)?-(size-1)*pow(size, 3)*4:pow(size,3)*4;
			for (int y=0;y<size;y+=1){
				neighbour[1]=(y==size-1)?-(size-1)*pow(size, 2)*4:pow(size,2)*4;
				for (int z=0;z<size;z+=1){
					neighbour[2]=(z==size-1)?-(size-1)*size*4:size*4;
					for (int t=0;t<size;t+=1){
						neighbour[3]=(x==size-1)?-(size-1)*4:4;
						for (int mu=0;mu<4;mu+=1){
							counter=x*size*size*size*4+y*size*size*4+z*size*4+t*4+mu;
							settozero(plaquettecontribution); settozero(deltacontribution);
							//calculate contribution of different plaquettes to action
							//U_nu(x+amu)*U^dagger_mu(x+anu)*U^dagger_nu(x)
							for (int nu =0;nu<4;nu+=1){
								if(mu!=nu){
								gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, GSL_COMPLEX_ONE, matrixarray[counter-mu+neighbour[mu]+nu], matrixarray[counter+neighbour[nu]], GSL_COMPLEX_ZERO, helpone);
								gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, GSL_COMPLEX_ONE, helpone, matrixarray[counter-mu+nu] , GSL_COMPLEX_ZERO, helptwo);
								gsl_matrix_complex_add(deltacontribution, helptwo);
								//~ fprintf(stdout, "%.1f\t", GSL_REAL(det2(helptwo)));
								if(mu>nu){gsl_matrix_complex_add(plaquettecontribution, helptwo);}
							}}
							for (int attempts=0;attempts<10;attempts+=1){
								generatesu2(multiplier, epsilon, generator);
								gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(1,0), multiplier, matrixarray[counter], gsl_complex_rect(0,0), newmatrix);
								//~ fprintf(stdout, "%.1f\t", GSL_REAL(det2(result)));
								//beta in exp(-betaDeltaS)
								if (exp(/*-1.0* */2.3*deltaS(plaquettecontribution, newmatrix, matrixarray[counter], helpone, helptwo))>gsl_rng_uniform(generator)){
									acceptance+=1;
									gsl_matrix_complex_memcpy(matrixarray[counter], newmatrix); 
								}
								gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(1,0), matrixarray[counter], plaquettecontribution, gsl_complex_rect(0,0), helpone);
								plaquetteexpectation+=0.5*GSL_REAL(trace(helpone));
							}
						}
					}
				}
			}
		}
		//~ fprintf(stdout, "\nacceptance rate: %f\n", (double)acceptance/((double)10*size*size*size*size*4));
		//~ fprintf(stdout, "\nplaquette expectation: %f\n", plaquetteexpectation/((double)10*size*size*size*size*4*3*0.5));
		fprintf(stdout, "%f\t%f\n", (double)acceptance/((double)10*size*size*size*size*4),plaquetteexpectation/((double)10*size*size*size*size*4*3*0.5));
	}
	
	
	//test behaviour of SU(2), SU(3) matrices
	for(int i=0; i<50; i+=10){
		epsilon=0.01*i;
		generatesu2(one, epsilon, generator);
		gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, gsl_complex_rect(1,0), one, one, gsl_complex_rect(0,0), multiplied);
		//printf("\n multiplied in run %d\n", i);
		//~ gsl_matrix_complex_fprintf(stdout, multiplied, "%f");
		//gsl_matrix_complex_fprintf(stdout,one, "%f");
		complexproduct=trace(one);
		//printf("\ntrace\n%f+%f*i\t%f\n\n", GSL_REAL(complexproduct), GSL_IMAG(complexproduct), gsl_complex_abs(complexproduct));
		complexproduct=gsl_complex_sub(gsl_complex_mul(gsl_matrix_complex_get(one, 0,0), gsl_matrix_complex_get(one, 1,1)), gsl_complex_mul(gsl_matrix_complex_get(one, 1,0), gsl_matrix_complex_get(one, 0,1)));
		//~ printf("\ndet\n%f+%f*i\t%f\n\n", GSL_REAL(complexproduct), GSL_IMAG(complexproduct), gsl_complex_abs(complexproduct));
		complexproduct=det2(one);
		//~ printf("\ndet\n%f+%f*i\t%f\n\n", GSL_REAL(complexproduct), GSL_IMAG(complexproduct), gsl_complex_abs(complexproduct));


		generatesu3(three, epsilon, generator);
		gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, gsl_complex_rect(1,0), three, three, gsl_complex_rect(0,0), multiplied3);
		//~ printf("\n multiplied3 in run %d\n", i);
		//~ gsl_matrix_complex_fprintf(stdout, multiplied3, "%f");
		//~ gsl_matrix_complex_fprintf(stdout,three, "%f");
		complexproduct=trace(three);
		//~ printf("\ntrace\n%f+%f*i\t%f\n\n", GSL_REAL(complexproduct), GSL_IMAG(complexproduct), gsl_complex_abs(complexproduct));
		complexproduct=det(three);
		//~ printf("\ndet\n%f+%f*i\t%f\n\n", GSL_REAL(complexproduct), GSL_IMAG(complexproduct), gsl_complex_abs(complexproduct));
	}
	
	//test error handling
		//~ complexproduct=trace(five);
		//~ conjugatetranspose(five);
		//~ generatesu2(five, epsilon, generator);
	
	//cleanup
	gsl_rng_free(generator);
	gsl_matrix_complex_free(one);
	gsl_matrix_complex_free(two);
	gsl_matrix_complex_free(multiplier);
	gsl_matrix_complex_free(result);
	gsl_matrix_complex_free(newmatrix);
	gsl_matrix_complex_free(deltacontribution);
	gsl_matrix_complex_free(plaquettecontribution);
	gsl_matrix_complex_free(helpone);
	gsl_matrix_complex_free(helptwo);
	gsl_matrix_complex_free(helpthree);
	gsl_matrix_complex_free(three);
	//~ gsl_matrix_complex_free(four);
	gsl_matrix_complex_free(multiplied);
	gsl_matrix_complex_free(multiplied3);
	gsl_matrix_complex_free(five);
	for (int i=0;i<size*size*size*size*4;i+=1){
		gsl_matrix_complex_free(matrixarray[i]);
	}
	return 0;
}


