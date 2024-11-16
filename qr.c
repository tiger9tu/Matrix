#include "qr.h"
#include "matrix.h"
#include <math.h>
#include <stdio.h>

void house(matrix *v, double *beta) {
  /* v is initialized as x */
  int m = v->height;
  matrix *subx = subVectorRef(v, 1, m);
  double sigma = innerProductVector(subx, subx);
  if (sigma == 0 && v->data[0] >= 0) {
    *beta = 0;
  } else if (sigma == 0 && v->data[0] < 0) {
    *beta = -2;
  } else {
    // calculate x1 - ||x|| while reduce cancelation error
    double normx = sqrt(sigma + (v->data[0]) * (v->data[0]));
    if (v->data[0] <= 0) {
      v->data[0] = v->data[0] - normx;
    } else {
      v->data[0] = -sigma / (v->data[0] + normx);
    }
    *beta = 2 * v->data[0] * v->data[0] / (sigma + v->data[0] * v->data[0]);
    rescaleMatrix(v, 1 / v->data[0]);
  }
}

houseHolderFactor *makeHouseHolderFactor(int n) {
  houseHolderFactor *out =
      (houseHolderFactor *)malloc(sizeof(houseHolderFactor));
  out->betas = (double *)malloc(sizeof(double) * (n));
  return out;
}

/*===========================================================================
 * Implicit Q.c, where Q is the one Householder factor Q, and c is a vector
 *=========================================================================*/
void implicitQix(matrix *vi, double beta, matrix *x) {
  double projLength = innerProductVector(vi, x);
  rescaleMatrixAdd(x, vi, 1, -beta * projLength);
}

void implilcitQx(houseHolderFactor *hhf, matrix *x) {
  for (int i = hhf->qrT->height - 1; i >= 0; i--) {
    matrix *vi = makeMatrix(1, hhf->qrT->width);
    matrix *qrTi =
        subVectorRef(hhf->qrT, i * hhf->qrT->width, (i + 1) * hhf->qrT->width);
    vi->data[i] = 1;
    for (size_t j = i + 1; j < vi->height; j++) {
      vi->data[j] = qrTi->data[j];
    }

    implicitQix(vi, hhf->betas[i], x);
    freeMatrix(vi);
  }
}

void implilcitQTx(houseHolderFactor *hhf, matrix *x) {
  for (int i = 0; i < hhf->qrT->height; i++) {
    matrix *vi = makeMatrix(1, hhf->qrT->width);
    matrix *qrTi =
        subVectorRef(hhf->qrT, i * hhf->qrT->width, (i + 1) * hhf->qrT->width);
    vi->data[i] = 1;
    for (size_t j = i + 1; j < vi->height; j++) {
      vi->data[j] = qrTi->data[j];
    }

    implicitQix(vi, hhf->betas[i], x);
    freeMatrix(vi);
  }
}

void printTranspose(matrix *a) {
  matrix *tmp = transposeMatrix(a);
  printMatrix(tmp);
  freeMatrix(tmp);
}

/*===========================================================================
 * Implicit
 *=========================================================================*/
houseHolderFactor *houseHolderQR(matrix *a) {
  // allocate factor
  houseHolderFactor *hhf = makeHouseHolderFactor(a->width);
  hhf->qrT = transposeMatrix(copyMatrix(a));

  matrix *xbuffer = makeMatrix(1, a->height);
  for (int i = 0; i < a->width && i < (a->height - 1); i++) {
    matrix *vi = subVectorRef(hhf->qrT, i * a->height + i, (i + 1) * a->height);
    for (int l = 0; l < a->height - i; l++) {
      xbuffer->data[l] = vi->data[l];
    }

    house(vi, &(hhf->betas[i]));

    // update the submatrix
    for (int j = i + 1; j < a->width; j++) {
      matrix *vj =
          subVectorRef(hhf->qrT, j * a->height + i, (j + 1) * a->height);
      implicitQix(vi, hhf->betas[i], vj);
    }

    // update the current column's top entry R(i, i)
    double vTc = 0;
    for (int l = 0; l < a->height - i; l++) {
      vTc += xbuffer->data[l] * vi->data[l];
    }

    vi->data[0] = xbuffer->data[0] - hhf->betas[i] * vTc;
  }

  return hhf;
}

void getExplicitQRFromHouseholder(houseHolderFactor *hhf, matrix **q,
                                  matrix **r) {

  assert(*q == NULL && *r == NULL,
         "The q and r matrices must be initalized to NULL.");

  *q = makeMatrix(hhf->qrT->width, hhf->qrT->width);
  *r = makeMatrix(hhf->qrT->height, hhf->qrT->width);
  for (int i = 0; i < (*r)->height; i++) {
    for (int j = i; j < (*r)->width; j++) {
      (*r)->data[i * (*r)->width + j] = hhf->qrT->data[j * hhf->qrT->width + i];
    }
  }

  matrix *I = eyeMatrix((*q)->width);
  // initialize Q as Im
  plusMatrix((*q), I);
  for (size_t i = 0; i < (*q)->width; i++) {
    // this part can be replaced by implicit opeartions,
    // which would be more memory efficient
    matrix *vi = makeMatrix(1, hhf->qrT->width);

    // subVectorRef(hhf->qrT, i * hhf->qrT->width, (i + 1) * hhf->qrT->width);

    for (size_t k = i + 1; k < hhf->qrT->width; k++) {
      vi->data[k] = hhf->qrT->data[i * hhf->qrT->width + k];
    }

    vi->data[i] = 1;

    matrix *viT = transposeMatrix(vi);
    matrix *betaviviT = multiplyMatrix(vi, viT);
    rescaleMatrix(betaviviT, -hhf->betas[i]);

    matrix *Qi = addMatrix(I, betaviviT);
    // timesMatrix(Q, I);
    matrix *QQi = multiplyMatrix((*q), Qi);
    copyData(QQi, (*q), 0, QQi->height, 0, QQi->width);

    freeMatrix(viT);
    freeMatrix(betaviviT);
    freeMatrix(Qi);
    freeMatrix(QQi);
    freeMatrix(vi);
  }

  freeMatrix(I);
}

void naive_gram_schmidt(matrix *a, matrix **q, matrix **r) {
  int i, j;
  double norm;
  double *qPtr;
  double *aPtr;
  double *rPtr;

  // assert that q and r matrices are initialized to NULL
  // assert(*q == NULL && *r == NULL);

  // Allocate Q and R matrices
  *q = makeMatrix(a->width, a->height);
  *r = makeMatrix(a->width, a->height);

  // Perform Gram-Schmidt Process
  for (j = 0; j < a->width; j++) {
    // Step 1: Set q_j = a_j (copy column j of a into q)
    qPtr = (*q)->data + j;
    aPtr = a->data + j;
    for (i = 0; i < a->height; i++) {
      *qPtr = *aPtr;
      qPtr += (*q)->width;
      aPtr += a->width;
    }

    // Step 2: Orthogonalize column j against all previous columns
    for (i = 0; i < j; i++) {
      // Compute R(i, j) = dot product of q_i and a_j
      double r_ij = dotProduct((*q)->data + i, a->data + j, a->height,
                               (*q)->width, a->width);
      rPtr = (*r)->data + (i * (*r)->width) + j;
      *rPtr = r_ij;

      // Subtract projection from q_j
      qPtr = (*q)->data + j;
      double *q_iPtr = (*q)->data + i;
      for (int k = 0; k < a->height; k++) {
        *qPtr -= (*q_iPtr) * r_ij;
        qPtr += (*q)->width;
        q_iPtr += (*q)->width;
      }
    }

    // Step 3: Normalize q_j and set R(j, j)
    norm = sqrt(dotProduct((*q)->data + j, (*q)->data + j, a->height,
                           (*q)->width, (*q)->width));
    rPtr = (*r)->data + (j * (*r)->width) + j;
    *rPtr = norm;

    // Divide q_j by its norm to normalize it
    qPtr = (*q)->data + j;
    for (i = 0; i < a->height; i++) {
      *qPtr = (norm == 0 ? 0 : *qPtr / norm);
      qPtr += (*q)->width;
    }
  }
}
