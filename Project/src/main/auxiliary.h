//header-Datei zu den auxiliary functions
#pragma once

#define GSL_COMPLEX_NAN (gsl_complex_rect(GSL_NAN,0.0))

/**
 * @brief returns trace of complex matrix
 * */
gsl_complex trace(gsl_matrix_complex *matrix);

/**
 * @brief returns det of complex matrix
 * */
gsl_complex det(gsl_matrix_complex *matrix);

/** 
 * @brief returns the determinant of a complex 2x2 matrix {{a,b},{c,d}} by returning ad-bc */
extern inline gsl_complex det2(gsl_matrix_complex *matrix);

/**
 * @brief transforms square matrix into its conjugate transposed. 
 * @note loops over i, conjugates diagonal elements, saves ij, sets ij=(ji)*, sets ji=(ij, old)*
 * */
void conjugatetranspose(gsl_matrix_complex *matrix);

/** sets all matrix elements to zero **/
void settozero(gsl_matrix_complex *matrix);

/** sets matrix to unity matrix **/
void settounity(gsl_matrix_complex *matrix);

/**
 * @brief
 * @note
 */
void crossproduct(gsl_vector_complex* a,gsl_vector_complex* b,gsl_vector_complex* result);

