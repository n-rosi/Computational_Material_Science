#include <stdint.h>

#include "gsl/gsl_eigen.h"
#include "gsl/gsl_vector_complex_double.h"
#include "gsl/gsl_vector_double.h"
#include "gsl/gsl_complex.h"
#include "gsl/gsl_complex_math.h"

#ifndef LINALG_H
#define LINALG_H

typedef uint16_t sizeVectorType;
typedef uint16_t sizeMatrixType[2];

void copyVector(VectorType *destination, const VectorType *source, size_t size);
void DiscoverVectorSize(VectorType *v, sizeVectorType *sizeV);
void DiscoverMatrixSize(MatrixType *M, sizeMatrixType *sizeM);
double* MatrixVectorProduct(VectorType *v, MatrixType *M);
void FillGSLMatrix(gsl_matrix *M_gsl, MatrixType *M);
void FillGSLMatrixComplex (gsl_matrix_complex *M_gsl, MatrixType *M);
void RealSymmetricDefinite_GeneralizedEigenvalueProblem(MatrixType *S, MatrixType *F, VectorType *new_c, double *E); 
void copyVector(VectorType *destination, const VectorType *source, size_t size);
void copyVector_gls(gsl_vector *fillingVector, double *newVec);

/* Functions to solve generalized eigenvalue problem with complex eigenvalue and eigenvectors */
// void complexv_ToRealv(gsl_vector_complex *complexVec, double *realVec);
// void ComplexHermitianDefinite_GeneralizedEigenvalueProblem(VectorType *c, MatrixType *S, MatrixType *F, VectorType *new_c, double *E);

#endif // LINALG_H