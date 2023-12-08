// function for matrix vector product
// function for eigenvalues-eigenstates computation

#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <complex.h>

#include "gsl/gsl_eigen.h"
#include "gsl/gsl_vector_complex_double.h"
#include "gsl/gsl_vector_double.h"

#include "source.h"
#include "linalg.h"

void copyVector(VectorType *destination, const VectorType *source, size_t size) {
    for (size_t i = 0; i < size; ++i) {
        (*destination)[i] = (*source)[i];
    }
}

void FillGSLMatrix(gsl_matrix *M_gsl, MatrixType *M){

    //gsl_complex value;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            gsl_matrix_set(M_gsl, i, j, (*M)[i][j]); 
        }
    }

}

void FillGSLMatrixComplex (gsl_matrix_complex *M_gsl, MatrixType *M){

    gsl_complex complexValue;

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            complexValue = (*M)[i][j] + 0*I;
            gsl_matrix_complex_set(M_gsl, i, j, complexValue); 
        }
    }

}

void complexToReal(gsl_vector_complex *complexVec, double *realVec) {
    size_t size = complexVec->size;

    for (size_t i = 0; i < size; ++i) {
        gsl_complex z = gsl_vector_complex_get(complexVec, i);
        realVec[i] = GSL_REAL(z);
    }
}

void fillC(gsl_vector *fillingVector, double *newVec){
    size_t size = fillingVector->size;
    for (size_t i = 0; i < size; ++i) {
        newVec[i] = gsl_vector_get(fillingVector, i);
    }
}

void RealSymmetricDefinite_GeneralizedEigenvalueProblem(VectorType *c, MatrixType *S, MatrixType *F, VectorType *new_c, double *E){

    /* Get problem dimensions */
    sizeVectorType size_v;
    sizeMatrixType size_M;

    DiscoverVectorSize(c, &size_v);
    DiscoverMatrixSize(S, &size_M);

    /* Define gsl matrices and vectors */
    gsl_vector_view evec_g;
    gsl_vector *eval = gsl_vector_alloc(size_v);
    gsl_matrix *evec = gsl_matrix_alloc(size_M[0], size_M[1]);

    gsl_matrix *S_gsl = gsl_matrix_alloc(size_M[0], size_M[1]);
    gsl_matrix *F_gsl = gsl_matrix_alloc(size_M[0], size_M[1]);

    /* Turn Schrodinger operators into gsl matrices */
    FillGSLMatrix(S_gsl,S);
    FillGSLMatrix(F_gsl,F);

    /* Create a workspace for the generalized eigenvalue problem  and solve it */
    gsl_eigen_gensymmv_workspace *workspace = gsl_eigen_gensymmv_alloc(size_v);
    gsl_eigen_gensymmv(F_gsl, S_gsl, eval, evec, workspace);

    /* Sort eigenvalues and eigenvectors to find ground state */
    gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_ASC);

    /* Modify coefficient vector (new_c) and Ground State Energy (Eg) */
    *E = gsl_vector_get(eval,0);
    evec_g = gsl_matrix_column(evec,0);
    fillC(&evec_g.vector, *new_c);
}

void ComplexHermitianDefinite_GeneralizedEigenvalueProblem(VectorType *c, MatrixType *S, MatrixType *F, VectorType *new_c, double *E){

    /* Get problem dimensions */
    sizeVectorType size_v;
    sizeMatrixType size_M;

    DiscoverVectorSize(c, &size_v);
    DiscoverMatrixSize(S, &size_M);

    /* Define gsl matrices and vectors */
    gsl_vector_complex_view evec_g;
    gsl_vector *eval = gsl_vector_alloc(size_v);
    gsl_matrix_complex *evec = gsl_matrix_complex_alloc(size_M[0], size_M[1]);

    gsl_matrix_complex *S_gsl = gsl_matrix_complex_alloc(size_M[0], size_M[1]);
    gsl_matrix_complex *F_gsl = gsl_matrix_complex_alloc(size_M[0], size_M[1]);

    /* Turn Schrodinger operators into gsl matrices */
    FillGSLMatrixComplex(S_gsl,S);
    FillGSLMatrixComplex(F_gsl,F);

    /* Create a workspace for the generalized eigenvalue problem  and solve it */
    gsl_eigen_genhermv_workspace *workspace = gsl_eigen_genhermv_alloc(size_v);
    gsl_eigen_genhermv(F_gsl, S_gsl, eval, evec, workspace);

    /* Sort eigenvalues and eigenvectors to find ground state */
    gsl_eigen_genhermv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC);

    /* Modify coefficient vector (new_c) and Ground State Energy (Eg) */
    *E = gsl_vector_get(eval,0);
    evec_g = gsl_matrix_complex_column(evec,0);
    complexToReal(&evec_g.vector, *new_c);
}

void DiscoverVectorSize(VectorType *v, sizeVectorType *sizeV){
    *sizeV = (sizeVectorType) (sizeof(*v)/sizeof((*v)[0]));
}

void DiscoverMatrixSize(MatrixType *M, sizeMatrixType *sizeM){
    (*sizeM)[0] = (sizeVectorType) (sizeof(*M) / sizeof((*M)[0]));
    (*sizeM)[1] = (sizeVectorType) (sizeof((*M)[0]) / sizeof((*M)[0][0]));
}







































// double* MatrixVectorProduct(VectorType *v, MatrixType *M){

//     sizeVectorType size_v;
//     sizeMatrixType size_M;

//     double *b = calloc(size_M[0], sizeof(double));

//     DiscoverVectorSize(v, &size_v);
//     DiscoverMatrixSize(M, &size_M);

//     for (int i=0; i<size_M[0]; i++){
//         for(int j=0; j<size_M[1]; j++){
//             b[i] += (*M)[i][j] * (*v)[j];
//         }
//     }

//     return b;
// }


