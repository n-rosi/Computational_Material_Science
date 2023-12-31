#include <stdlib.h>
#include <math.h>

#include "source.h"

void normalize_eigenvector(VectorType *c, MatrixType *S){

    /* Smallest eigenvector normalization */

    double normalizationFactor2 = 0;

    for(int p=0; p<4; p++){
        for(int q=0; q<4; q++){
            normalizationFactor2 = normalizationFactor2 + (*c)[p]*(*S)[p][q]*(*c)[q];
        }
    }

    for(int p=0; p<4; p++){
        (*c)[p] = (*c)[p] / pow(normalizationFactor2, 0.5);
    }

}

double GroundStateEnergy(VectorType *c, TensorType *Q, MatrixType *h){

    /* Compute the ground state energy */

    double term1=0, term2=0;

    for(int p=0; p<4; p++){
        for(int q=0; q<4; q++){
            term1 = term1 + (*c)[p]*(*h)[p][q]*(*c)[q];
        }
    }

    for(int p=0; p<4; p++){
        for(int q=0; q<4; q++){
            for(int r=0; r<4; r++){
                for(int s=0; s<4; s++){
                    term2 = term2 + (*Q)[p][r][q][s]*(*c)[p]*(*c)[q]*(*c)[r]*(*c)[s];
                }
            }
        }
    }
    return 2*term1+term2;
}

void F_pq(VectorType *c, TensorType *Q, MatrixType *h, MatrixType *F){

    /* Hamiltonian operator in matrix form */

    MatrixType *suppVar = malloc(sizeof(MatrixType));

    for(int p=0; p<4; p++){
        for(int q=0; q<4; q++){
            (*suppVar)[p][q] = 0.;
        }
    }

    for(int p=0; p<4; p++){
        for(int r=0; r<4; r++){
            for(int q=0; q<4; q++){
                for(int s=0; s<4; s++){
                    (*suppVar)[q][p] = (*suppVar)[q][p] + (*Q)[p][r][q][s]*(*c)[r]*(*c)[s];
                }
            }
        }
    }

    for(int p=0; p<4; p++){
        for(int q=0; q<4; q++){
            (*F)[p][q] = (*h)[p][q] + (*suppVar)[p][q];
        }
    }

}

void T_pq(const VectorType *alpha, MatrixType *T){

    /* Kinetic energy */

    for(int p=0; p<4; p++){
        for(int q=0; q<4; q++){
            (*T)[p][q] = (3 * (*alpha)[p] * (*alpha)[q] * pow(PI, (1.5) )) / pow( ((*alpha)[p]+(*alpha)[q]) , (2.5) );
        }
    }
}

void A_pq(const VectorType *alpha, MatrixType *A){

    /* n-e Coulomb energy */

    for(int p=0; p<4; p++){
        for(int q=0; q<4; q++){
            (*A)[p][q] = - (2 * PI) / ( (*alpha)[p]+(*alpha)[q] );
        }
    }
}

void h_pq(const VectorType *alpha, MatrixType *h, MatrixType *T, MatrixType *A){

    /* single particle hamiltonian (kinetic energy and n-e Coulomb energy) */

    T_pq(alpha,T);
    A_pq(alpha,A);

    for(int p=0; p<4; p++){
        for(int q=0; q<4; q++){
            (*h)[p][q] = (*T)[p][q] + 2*(*A)[p][q];
        }
    }

}

void S_pq(const VectorType *alpha, MatrixType *S){

    /* <Chi_p|Chi_q> product */

    for(int p=0; p<4; p++){
        for(int q=0; q<4; q++){
            (*S)[p][q] = pow( ( PI / ( (*alpha)[p] + (*alpha)[q] ) ) , (1.5) );
        }
    }
}

void Q_prqs(const VectorType *alpha, TensorType *Q){

    /* e-e Coulomb Energy */

    for(int p=0; p<4; p++){
        for(int r=0; r<4; r++){
            for(int q=0; q<4; q++){
                for(int s=0; s<4; s++){
                    (*Q)[p][r][q][s] = (2 * pow(PI, 2.5)) / ( ( (*alpha)[p]+(*alpha)[q] ) * ( (*alpha)[r]+(*alpha)[s] ) * pow(((*alpha)[p]+(*alpha)[q]+(*alpha)[r]+(*alpha)[s]), 0.5) );
                }
            }
        }
    }

}