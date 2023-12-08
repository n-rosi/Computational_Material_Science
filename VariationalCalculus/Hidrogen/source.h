#ifndef SOURCE_H
#define SOURCE_H

#define PI 3.14159265358979323846

#define dim1 4
#define dim2 4

typedef double VectorType[dim1];
typedef double MatrixType[dim1][dim2]; 

void S_pq(const VectorType *alpha, MatrixType *S);
void T_pq(const VectorType *alpha, MatrixType *T);
void A_pq(const VectorType *alpha, MatrixType *A);
void h_pq(const VectorType *alpha, MatrixType *h, MatrixType *T, MatrixType *A);

#endif // SOURCE_H