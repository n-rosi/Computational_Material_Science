// function for matrix vector product
// function for eigenvalues-eigenstates computation

#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>

#include "gsl/gsl_eigen.h"
#include "gsl/gsl_vector_double.h"

#include "source.h"
#include "linalg.h"

void copyVector(VectorType *destination, const VectorType *source, size_t size) {

    // Copy an existing vector in a destination vector.

    for (size_t i = 0; i < size; ++i) {
        (*destination)[i] = (*source)[i];
    }
}

void FillGSLMatrix(gsl_matrix *M_gsl, MatrixType *M){

    // Fill a GSL matrix with the values contained in a MatrixType.

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            gsl_matrix_set(M_gsl, i, j, (*M)[i][j]); 
        }
    }

}

void fillC(gsl_vector *fillingVector, double *newVec){

    // Fill a vector of double with the values contained in a GSL vector.

    size_t size = fillingVector->size;
    for (size_t i = 0; i < size; ++i) {
        newVec[i] = gsl_vector_get(fillingVector, i);
    }
}

void DiscoverVectorSize(VectorType *v, sizeVectorType *sizeV){
    *sizeV = (sizeVectorType) (sizeof(*v)/sizeof((*v)[0]));
}
void DiscoverMatrixSize(MatrixType *M, sizeMatrixType *sizeM){
    (*sizeM)[0] = (sizeVectorType) (sizeof(*M) / sizeof((*M)[0]));
    (*sizeM)[1] = (sizeVectorType) (sizeof((*M)[0]) / sizeof((*M)[0][0]));
}

void RealSymmetricDefinite_GeneralizedEigenvalueProblem(MatrixType *S, MatrixType *F, VectorType *c, double *E){

    // Solve eigenvalue problem.

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
    fillC(&evec_g.vector, *c);
}










































