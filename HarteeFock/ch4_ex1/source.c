#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "source.h"

void Q_prqs(const VectorType *alpha, TensorType *Q){

    // Build e-e Coulomb Energy

    for(int p=1; p<=4; p++){
        for(int r=1; r<=4; r++){
            for(int q=1; q<=4; q++){
                for(int s=1; s<=4; s++){
                    (*Q)[p][r][q][s] = (2 * pow(PI, 2.5)) / ( ((*alpha)[p]+(*alpha)[q])*((*alpha)[r]+(*alpha)[s]) * 
                                                            pow(((*alpha)[p]+(*alpha)[q]+(*alpha)[r]+(*alpha)[s]), 0.5) );
                }
            }
        }
    }

}

void S_pq(const VectorType *alpha, MatrixType *S){

    // Build <Chi_p|Chi_q> product

    for(int p=1; p<=4; p++){
        for(int q=1; q<=4; q++){
            (*S)[p][q] = pow((PI/((*alpha)[p]+(*alpha)[q])), 1.5);
        }
    }
}

void T_pq(const VectorType *alpha, MatrixType *T){

    // Build kinetic energy

    for(int p=1; p<=4; p++){
        for(int q=1; q<=4; q++){
            (*T)[p][q] = (3 * (*alpha)[p] * (*alpha)[q] * pow(PI, 1.5)) / pow(((*alpha)[p]+(*alpha)[q]), 2.5);
        }
    }
}

void A_pq(const VectorType *alpha, MatrixType *A){

    // Build n-e Coulomb energy

    for(int p=1; p<=4; p++){
        for(int q=1; q<=4; q++){
            (*A)[p][q] = - (2 * PI) / ((*alpha)[p]+(*alpha)[q]);
        }
    }
}

void h_pq(const VectorType *alpha, MatrixType *h, MatrixType *T, MatrixType *A){

    // Build independent particle term

    T_pq(alpha,T);
    A_pq(alpha,A);

    for(int p=1; p<=4; p++){
        for(int q=1; q<=4; q++){
            (*h)[p][q] = (*T)[p][q] + (*A)[p][q];
        }
    }

}

void F_pq(VectorType *c, TensorType *Q, MatrixType *h, MatrixType *F){

    // Combine all terms together
    MatrixType *suppVar = calloc(0, sizeof(MatrixType));

    for(int p=1; p<=4; p++){
        for(int r=1; r<=4; r++){
            for(int q=1; q<=4; q++){
                for(int s=1; s<=4; s++){
                    (*suppVar)[r][s] = (*suppVar)[r][s] + (*Q)[p][q][r][s]*(*c)[r]*(*c)[s];
                }
            }
        }
    }

    for(int p=1; p<=4; p++){
        for(int q=1; q<=4; q++){
            (*F)[p][q] = (*h)[p][q] + (*suppVar)[p][q];
        }
    }

}