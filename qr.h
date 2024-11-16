#ifndef __ISOMAP_QR
#define __ISOMAP_QR
#include "matrix.h"
typedef struct {
  double *betas;
  matrix *qrT;
} houseHolderFactor;
void naive_gram_schmidt(matrix *a, matrix **q, matrix **r);
void house(matrix *v, double *beta);
houseHolderFactor *houseHolderQR(matrix *a);
void getExplicitQRFromHouseholder(houseHolderFactor *hhf, matrix **q,
                                  matrix **r);
void implilcitQx(houseHolderFactor *hhf, matrix *x);
void implilcitQTx(houseHolderFactor *hhf, matrix *x);

#endif
