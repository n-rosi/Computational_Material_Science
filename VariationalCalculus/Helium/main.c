#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "source.h"
#include "linalg.h"

void main(void){

    /* Declare problem parameters */
    VectorType *alpha = malloc(sizeof(VectorType));
    VectorType *c = malloc(sizeof(VectorType));
    VectorType *new_c = calloc(sizeof(VectorType),0.);

    /* Declare hamiltonian quantities */
    TensorType *Q = malloc(sizeof(TensorType));
    MatrixType *h = malloc(sizeof(MatrixType));
    MatrixType *T = malloc(sizeof(MatrixType));
    MatrixType *A = malloc(sizeof(MatrixType));
    MatrixType *S = malloc(sizeof(MatrixType));
    MatrixType *F = malloc(sizeof(MatrixType));

    /* Declare other variables */
    double groundStateEnergy_value = 0.0;
    double Eg_value = 1000.0;
    double new_Eg_value = 0.0;
    double *groundStateEnergy = &groundStateEnergy_value;
    double *Eg = &Eg_value;
    double *new_Eg = &new_Eg_value;
    double delta = 1e-6;

    /* Assign problem parameters */
    (*alpha)[0] = 0.298073;
    (*alpha)[1] = 1.242567; 
    (*alpha)[2] = 5.782948;
    (*alpha)[3] = 38.474970;

    (*c)[0] = 0.25;
    (*c)[1] = 0.25; 
    (*c)[2] = 0.25;
    (*c)[3] = 0.25;

    /* Compute hamiltonian quantities */
    Q_prqs(alpha,Q); 
    S_pq(alpha,S); 
    h_pq(alpha,h,T,A);
    normalize_eigenvector(c,S);
    F_pq(c,Q,h,F);

    /* Self-consistency loop */
    while(fabs(*Eg-*new_Eg)>delta){

        *Eg = *new_Eg;
        RealSymmetricDefinite_GeneralizedEigenvalueProblem(S,F,new_c,new_Eg);
        normalize_eigenvector(new_c,S);
        copyVector(c, new_c, sizeof(VectorType) / sizeof(double));
        F_pq(c,Q,h,F);
    }

    groundStateEnergy_value = GroundStateEnergy(new_c, Q, h);
    printf("Helium Ground State Energy: %f Hartee.\n", groundStateEnergy_value);
}