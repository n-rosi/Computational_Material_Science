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
void FillGSLMatrix(gsl_matrix *M_gsl, MatrixType *M);
void RealSymmetricDefinite_GeneralizedEigenvalueProblem(MatrixType *S, MatrixType *F, VectorType *c, double *E); 
void copyVector(VectorType *destination, const VectorType *source, size_t size);
void fillC(gsl_vector *fillingVector, double *newVec);

#endif // LINALG_H