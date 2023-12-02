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
    VectorType *new_c = calloc(sizeof(VectorType),0.);

    (*alpha)[0] = 0.298073;
    (*alpha)[1] = 1.242567; 
    (*alpha)[2] = 5.782948;
    (*alpha)[3] = 38.474970;

    (*c)[0] = 10.;
    (*c)[1] = 5.; 
    (*c)[2] = 7.;
    (*c)[3] = 0.;

    double groundStateEnergy_value = 0.0;
    double Eg_value = 1000.0;
    double new_Eg_value = 0.0;
    double *groundStateEnergy = &groundStateEnergy_value;
    double *Eg = &Eg_value;
    double *new_Eg = &new_Eg_value;
    double delta = 1e-10;

    /* Declare hamiltonian quantities */
    TensorType *Q = malloc(sizeof(TensorType));
    MatrixType *h = malloc(sizeof(MatrixType));
    MatrixType *T = malloc(sizeof(MatrixType));
    MatrixType *A = malloc(sizeof(MatrixType));
    MatrixType *S = malloc(sizeof(MatrixType));
    MatrixType *F = malloc(sizeof(MatrixType));

    /* Compute hamiltonian quantities */
    Q_prqs(alpha,Q); // CHECKED, LOOKS OK
    S_pq(alpha,S); // CHECKED, LOOKS OK
    h_pq(alpha,h,T,A); // CHECKED, LOOKS OK
    Cnormalization(c,S);

    /*  */  
    while(fabs(*Eg - *new_Eg)>delta){
        *Eg = *new_Eg;
        F_pq(c,Q,h,F); 
        GeneralizedEigenvalueProblem(c,S,F,new_c,new_Eg);
        for(int i=0; i<4;i++){
            printf("%f ", (*new_c)[i]);
        }
        printf("\n%f\n", (*new_Eg));
        GroundStateEnergy(new_c, Q, h, groundStateEnergy);
        Cnormalization(new_c,S);
        copyVector(c, new_c, sizeof(VectorType) / sizeof(double));
        printf("%f\n", *groundStateEnergy);
    }
}