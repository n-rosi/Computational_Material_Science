#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "source.h"

void T_pq(const VectorType *alpha, MatrixType *T){

    // Build kinetic energy

    for(int p=0; p<4; p++){
        for(int q=0; q<4; q++){
            (*T)[p][q] = (3 * (*alpha)[p] * (*alpha)[q] * pow(PI, (1.5) )) / pow( ((*alpha)[p]+(*alpha)[q]) , (2.5) );
        }
    }
}

void A_pq(const VectorType *alpha, MatrixType *A){

    // Build n-e Coulomb energy

    for(int p=0; p<4; p++){
        for(int q=0; q<4; q++){
            (*A)[p][q] = - (2 * PI) / ( (*alpha)[p]+(*alpha)[q] );
        }
    }
}

void h_pq(const VectorType *alpha, MatrixType *h, MatrixType *T, MatrixType *A){

    // Build independent particle term

    T_pq(alpha,T);
    A_pq(alpha,A);

    for(int p=0; p<4; p++){
        for(int q=0; q<4; q++){
            (*h)[p][q] = (*T)[p][q] + (*A)[p][q];
        }
    }

}

void S_pq(const VectorType *alpha, MatrixType *S){

    // Build <Chi_p|Chi_q> product
    
    for(int p=0; p<4; p++){
        for(int q=0; q<4; q++){
            (*S)[p][q] = pow( ( PI / ( (*alpha)[p] + (*alpha)[q] ) ) , (1.5) );
        }
    }
}