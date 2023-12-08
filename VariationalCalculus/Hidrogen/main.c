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
    VectorType *alpha = malloc(sizeof(VectorType));
    VectorType *c = malloc(sizeof(VectorType));

    (*alpha)[0] = 13.00773;
    (*alpha)[1] = 1.962079; 
    (*alpha)[2] = 0.444529;
    (*alpha)[3] = 0.1219492;

    double Eg_value = 0.0;
    double *Eg = &Eg_value;
    double groundStateEnergy_value = 0.0;
    double *groundStateEnergy = &groundStateEnergy_value;

    /* Declare hamiltonian quantities */
    MatrixType *h = malloc(sizeof(MatrixType));
    MatrixType *T = malloc(sizeof(MatrixType));
    MatrixType *A = malloc(sizeof(MatrixType));
    MatrixType *S = malloc(sizeof(MatrixType));

    /* Compute hamiltonian quantities */
    S_pq(alpha,S); 
    h_pq(alpha,h,T,A);

    /* Solve generalized eigenvalue problem hC = ESC */
    RealSymmetricDefinite_GeneralizedEigenvalueProblem(S,h,c,Eg);
    printf("\nHidrogen ground state energy in Hartee: %f\n\n", (*Eg));
}