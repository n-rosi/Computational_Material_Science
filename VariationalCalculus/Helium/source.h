#ifndef SOURCE_H
#define SOURCE_H

#define PI 3.14159265358979323846

#define dim1 4
#define dim2 4
#define dim3 4
#define dim4 4

typedef double VectorType[dim1];
typedef double MatrixType[dim1][dim2]; 
typedef double TensorType[dim1][dim3][dim2][dim4];

void Q_prqs(const VectorType *alpha, TensorType *Q);
void S_pq(const VectorType *alpha, MatrixType *S);
void T_pq(const VectorType *alpha, MatrixType *T);
void A_pq(const VectorType *alpha, MatrixType *A);
void h_pq(const VectorType *alpha, MatrixType *h, MatrixType *T, MatrixType *A);
void F_pq(VectorType *c, TensorType *Q, MatrixType *h, MatrixType *F);
double GroundStateEnergy(VectorType *c, TensorType *Q, MatrixType *h);
void normalize_eigenvector(VectorType *c, MatrixType *S);

#endif // SOURCE_H