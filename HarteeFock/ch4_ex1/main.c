#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "gsl/gsl_matrix.h"
#include "gsl/gsl_eigen.h"
#include "gsl/gsl_math.h"

#include "source.h"
#include "linalg.h"

void main(void){

    /* Declare problem parameters */
    const VectorType alpha_values = {0.298073, 1.242567, 5.782948, 38.474970};
    const VectorType *alpha = &alpha_values;

    VectorType c_values = {1., 1., 1., 1.};
    VectorType *c = &c_values;
    VectorType *new_c = calloc(sizeof(c_values),0.);

    double groundStateEnergy = 1000.0;
    double new_groundStateEnergy = 1000.0;
    double *Eg = &groundStateEnergy;
    double *new_Eg = &new_groundStateEnergy;
    double delta = 1e-6;

    /* Declare hamiltonian quantities */
    TensorType *Q = malloc(sizeof(TensorType));
    MatrixType *h = malloc(sizeof(MatrixType));
    MatrixType *T = malloc(sizeof(MatrixType));
    MatrixType *A = malloc(sizeof(MatrixType));
    MatrixType *S = malloc(sizeof(MatrixType));
    MatrixType *F = malloc(sizeof(MatrixType));

    /* Compute hamiltonian quantities */
    Q_prqs(alpha,Q);
    S_pq(alpha,S); 
    h_pq(alpha,h,T,A);

    /*  */
    while(fabs(Eg-new_Eg)>delta){

        groundStateEnergy = new_groundStateEnergy;
        F_pq(c,Q,h,F); 
        GeneralizedEigenvalueProblem(c,S,F,new_c,new_Eg);
        copyVector(c, new_c, sizeof(VectorType) / sizeof(double));
        printf("%f\n",new_Eg);

    }
}