#include <stdint.h>

#ifndef LINALG_H
#define LINALG_H

typedef uint16_t sizeVectorType;
typedef uint16_t sizeMatrixType[2];

void copyVector(VectorType *destination, const VectorType *source, size_t size);
void DiscoverVectorSize(VectorType *v, sizeVectorType *sizeV);
void DiscoverMatrixSize(MatrixType *M, sizeMatrixType *sizeM);
double* MatrixVectorProduct(VectorType *v, MatrixType *M);
void FillGSLMatrix(gsl_matrix *M_gsl, MatrixType *M);
void GeneralizedEigenvalueProblem(VectorType *c, MatrixType *S, MatrixType *F, VectorType *new_c, double *Eg);
void copyVector(VectorType *destination, const VectorType *source, size_t size);

#endif // LINALG_H