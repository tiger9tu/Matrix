#include "matrix.h"
#include "qr.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <math.h>
#include <stdio.h>
#include <time.h>

// Function to generate a random upper triangular matrix of size n x n
void generateRandomUpperTriangularSquareMatrix(int n, matrix **r,
                                               gsl_matrix **gsl_r) {
  *r = makeMatrix(n, n);
  *gsl_r = gsl_matrix_alloc(n, n);
  printf("gslr size = %d\n", (*gsl_r)->size2);
  // Use the current time as the seed for the random number generator
  srand(time(NULL));

  // Generate the upper triangular matrix
  for (size_t i = 0; i < n; i++) {
    for (size_t j = i; j < n; j++) {
      // Assign a random number to the upper triangular part (including the
      // diagonal)
      int randn = rand() % 100;
      (*r)->data[i * n + j] =
          randn; // Generate a random integer between 0 and 99
      gsl_matrix_set(*gsl_r, i, j, (double)randn);
    }
  }
}

// Function to generate a random upper triangular matrix of size n x n
void generateRandomSquareMatrix(int n, matrix **r, gsl_matrix **gsl_r) {
  *r = makeMatrix(n, n);
  *gsl_r = gsl_matrix_alloc(n, n);

  // Use the current time as the seed for the random number generator
  srand(time(NULL));

  // Generate the upper triangular matrix
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      // Assign a random number to the upper triangular part (including the
      // diagonal)
      int randn = rand() % 100;
      (*r)->data[i * n + j] =
          randn; // Generate a random integer between 0 and 99
      gsl_matrix_set(*gsl_r, i, j, (double)randn);
    }
  }
}

// Function to print a gsl_matrix
void printGslMatrix(const gsl_matrix *matrix) {
  size_t rows = matrix->size1;
  size_t cols = matrix->size2;

  for (size_t i = 0; i < rows; i++) {
    for (size_t j = 0; j < cols; j++) {
      // Access the matrix element using gsl_matrix_get
      printf("%6.2f ", gsl_matrix_get(matrix, i, j));
    }
    printf("\n");
  }
}

void getNormalMatrixFromGsl(matrix *a, gsl_matrix *b) {
  for (size_t i = 0; i < a->height; i++) {
    for (size_t j = 0; j < a->width; j++) {
      // Access the matrix element using gsl_matrix_get
      a->data[i * a->width + j] = gsl_matrix_get(b, i, j);
    }
  }
}

double calculateFrobeniusNorm(matrix *a) {
  double norm = 0.0;
  for (int i = 0; i < a->height; i++) {
    for (int j = 0; j < a->width; j++) {
      norm += a->data[i * a->width + j] * a->data[i * a->width + j];
    }
  }
  return sqrt(norm);
}

double calculateGslFrobeniusNorm(gsl_matrix *a) {
  double norm = 0.0;
  for (int i = 0; i < a->size1; i++) {
    for (int j = 0; j < a->size2; j++) {
      norm += gsl_matrix_get(a, i, j) * gsl_matrix_get(a, i, j);
    }
  }
  return sqrt(norm);
}

double giveError(matrix *a, matrix *aHat) {
  double normA = calculateFrobeniusNorm(a);
  matrix *aMinusAHat = makeMatrix(a->width, a->height);
  matrixAdd(aMinusAHat, a);
  matrixMinus(aMinusAHat, aHat);
  double error = calculateFrobeniusNorm(aMinusAHat);
  double relativeError = error / normA;

  return relativeError;
}

matrix *generateVandMatrix(int m, int n) {
  matrix *A = makeMatrix(n, m);
  for (int i = 0; i < m; i++) {
    double alpha_i = (double)i / (m - 1);
    for (int j = 0; j < n; j++) {
      set_element(A, i, j, pow(alpha_i, j));
    }
  }
  return A;
}

// Function to generate vector b
matrix *generateVandVec(int m) {
  matrix *b = makeMatrix(1, m);
  for (int i = 0; i < m; i++) {
    double alpha_i = (double)i / (m - 1);
    set_element(b, i, 0, exp(sin(4 * alpha_i)));
  }
  return b;
}

matrix *solveSystemR(matrix *R, matrix *b) {
  int n = R->width;
  matrix *x = makeMatrix(1, n);

  for (int i = n - 1; i >= 0; i--) {
    double sum = 0.0;
    for (int j = i + 1; j < n; j++) {
      sum += get_element(R, i, j) * get_element(x, j, 0);
    }
    double value;
    double pivot = get_element(R, i, i);
    if (pivot == 0) {
      value = 0;
    } else {
      value = (get_element(b, i, 0) - sum) / pivot;
    }
    set_element(x, i, 0, value);
  }

  return x;
}

gsl_matrix *getGslFromNormal(matrix *a) {
  gsl_matrix *gsla = gsl_matrix_alloc(a->height, a->width);
  for (size_t i = 0; i < a->height; i++) {
    for (size_t j = 0; j < a->width; j++) {
      gsl_matrix_set(gsla, i, j, a->data[i * a->width + j]);
    }
  }
  return gsla;
}

void solve_linear_system(const gsl_matrix *A, const gsl_vector *b,
                         gsl_vector *x) {
  // A: Coefficient matrix (input)
  // b: Right-hand side vector (input)
  // x: Solution vector (output)

  // Make a copy of the input matrix A because LU decomposition modifies the
  // matrix
  gsl_matrix *LU = gsl_matrix_alloc(A->size1, A->size2);
  gsl_matrix_memcpy(LU, A);

  // Create a permutation object for row exchanges during LU decomposition
  gsl_permutation *perm = gsl_permutation_alloc(A->size1);
  int signum;
  // Perform LU decomposition of A
  gsl_linalg_LU_decomp(LU, perm, &signum);

  // Solve the system Ax = b using the LU decomposition
  gsl_linalg_LU_solve(LU, perm, b, x);

  // Free allocated memory
  gsl_matrix_free(LU);
  gsl_permutation_free(perm);
}

void usinggslSolveNormalEqu(matrix *A, matrix *b) {
  gsl_matrix *gslA = getGslFromNormal(A);
  gsl_vector *gsl_b = gsl_vector_alloc(b->height);
  for (size_t i = 0; i < b->height; i++) {
    gsl_vector_set(gsl_b, i, b->data[i]);
  }

  gsl_matrix *gslATA = gsl_matrix_alloc(A->width, A->width);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, gslA, gslA, 0.0, gslATA);
  gsl_matrix *gslATb = gsl_vector_alloc(A->width);
  gsl_blas_dgemv(CblasTrans, 1.0, gslA, gsl_b, 0.0, gslATb);

  gsl_vector *x = gsl_vector_alloc(gslA->size2);

  // solve ATAx = ATb:
  solve_linear_system(gslATA, gslATb, x);

  for (size_t i = 0; i < x->size; ++i) {
    printf("v[%d] = %lf\n", i, gsl_vector_get(x, i));
  }
}

matrix *generateTestMatrix() {
  matrix *a = makeMatrix(3, 4);
  for (int i = 0; i < a->height; i++)
    for (int j = 0; j < a->width; j++)
      a->data[i * 3 + j] = i * 3 + j;
  return a;
}

void a_b_() {
  matrix *a = generateTestMatrix();
  printf("The original matrix A:\n");
  printMatrix(a);
  // gram-schimit decomposition
  matrix *gramSchimidtQ = NULL;
  matrix *gramSchimidtR = NULL;

  naive_gram_schmidt(a, &gramSchimidtQ, &gramSchimidtR);
  printf("The gram schimidt decomposition:");
  printf("Q: \n");
  printMatrix(gramSchimidtQ);
  printf("R: \n");
  printMatrix(gramSchimidtR);

  houseHolderFactor *hhf = houseHolderQR(a);
  matrix *houseHolderQ = NULL;
  matrix *houseHolderR = NULL;
  getExplicitQRFromHouseholder(hhf, &houseHolderQ, &houseHolderR);
  printf("The Householder decomposition:");
  printf("Q: \n");
  printMatrix(houseHolderQ);
  printf("R: \n");
  printMatrix(houseHolderR);
  printf("\n\n\n");
}

void c_() {

  matrix *a = generateTestMatrix();
  houseHolderFactor *hhf = houseHolderQR(a);
  matrix *x = makeMatrix(1, 4);
  for (size_t i = 0; i < x->height; i++)
    x->data[i] = i + 1;

  matrix *y = copyMatrix(x);
  matrix *r = NULL;
  matrix *q = NULL;

  getExplicitQRFromHouseholder(hhf, &q, &r);
  printf("Q: \n");
  printMatrix(q);

  printf("x: \n");
  printMatrix(x);
  printf("explicit computation of QX: \n");
  printMatrix(multiplyMatrix(q, x));
  printf("Implicit computation of Qx: \n");
  implilcitQx(hhf, x);
  printMatrix(x);

  printf("\ny: \n");
  printMatrix(y);
  printf("explicit computation of QTy: \n");
  printMatrix(multiplyMatrix(transposeMatrix(q), y));
  printf("implicit computation of QTy: \n");
  implilcitQTx(hhf, y);
  printMatrix(y);
}

void d_() {

  size_t tries = 1;
  FILE *file = fopen("output.txt", "w");
  for (int n = 1; n <= 15; n++) {
    double hhER = 0;
    double hhEQ = 0;
    double hhEQR = 0;
    double hhEQQ = 0;
    double gsER = 0;
    double gsEQ = 0;
    double gsEQR = 0;
    double gsEQQ = 0;
    for (size_t try = 0; try < tries; try++) {
      // setting up the random sample
      matrix *r = NULL;
      matrix *b = NULL;

      gsl_matrix *gsl_r = NULL;
      gsl_matrix *gsl_b = NULL;

      generateRandomUpperTriangularSquareMatrix(n, &r, &gsl_r);
      generateRandomSquareMatrix(n, &b, &gsl_b);

      gsl_matrix *gsl_t = gsl_matrix_alloc(n, n);
      gsl_linalg_QR_decomp_r(gsl_b, gsl_t);

      gsl_matrix *gsl_bq = gsl_matrix_alloc(n, n);
      gsl_matrix *gsl_br = gsl_matrix_alloc(n, n);

      gsl_linalg_QR_unpack_r(gsl_b, gsl_t, gsl_bq, gsl_br);

      matrix *q = makeMatrix(n, n);
      getNormalMatrixFromGsl(q, gsl_bq);

      gsl_matrix *gsl_a = gsl_matrix_alloc(n, n);
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, gsl_bq, gsl_r, 0.0,
                     gsl_a);
      matrix *a = makeMatrix(n, n);
      getNormalMatrixFromGsl(a, gsl_a);

      // my QR decomposition
      houseHolderFactor *hhf = houseHolderQR(a);
      matrix *myq = NULL;
      matrix *myr = NULL;
      getExplicitQRFromHouseholder(hhf, &myq, &myr);

      // printf("R error of the Householder:\n");
      hhER += giveError(r, myr);
      if (hhER > 0.5) {
        printf("hhER too large, r = \n");
        printMatrix(r);
        printf("my r : \n");
        printMatrix(myr);

        printf("a : \n");
        printMatrix(a);
      }

      // error type 2
      hhEQ += giveError(q, myq);

      // error type 3
      //  printf("QR error of the Householder:\n");
      hhEQR += giveError(a, multiplyMatrix(myq, myr));

      // error type 4
      //  printf("QTQ error of the Householder:\n");
      hhEQQ +=
          giveError(eyeMatrix(n), multiplyMatrix(transposeMatrix(myq), myq));

      // printf("\n\n\n");

      // gramSchmidt
      matrix *gramschmidtq = NULL;
      matrix *gramschmidtr = NULL;
      naive_gram_schmidt(a, &gramschmidtq, &gramschmidtr);
      // double normR = calculateFrobeniusNorm(r);
      // printf("R error of the Gram Schimidt:\n");
      gsER += giveError(r, gramschmidtr);

      // error type 2
      // printf("Q error of the Gram Schimidt.\n");
      gsEQ += giveError(q, gramschmidtq);

      // error type 3
      //  printf("QR error of the Gram Schimidt:\n");
      gsEQR += giveError(a, multiplyMatrix(gramschmidtq, gramschmidtr));

      // error type 4
      //  printf("QTQ error of the Gram Schimidt:\n");
      gsEQQ +=
          giveError(eyeMatrix(n), multiplyMatrix(transposeMatrix(gramschmidtq),
                                                 gramschmidtq));
      // fprintf(file, "\n");
    }
    fprintf(file, "%lf %lf %lf %lf ", hhER / tries, hhEQ / tries, hhEQR / tries,
            hhEQQ / tries);
    fprintf(file, "%lf %lf %lf %lf\n", gsER / tries, gsEQ / tries,
            gsEQR / tries, gsEQQ / tries);
  }

  // gsl_linalg_QR_decomp_r()
}

void e_() {
  int m = 100;
  int n = 15;
  matrix *A = generateVandMatrix(m, n);
  matrix *b = generateVandVec(m);

  // (i)
  houseHolderFactor *hhf = houseHolderQR(A);
  matrix *b1 = copyMatrix(b);
  implilcitQTx(hhf, b1); // now b is Qtb
  matrix *x = solveSystemR(transposeMatrix(hhf->qrT) /* same effect as R*/, b1);
  //   printMatrix(x);
  printf("The (i) gives result c14 = %.17f\n\n", x->data[14]);

  // (ii)
  matrix *augAb = makeMatrix(A->width + 1, A->height);
  copyData(A, augAb, 0, A->height, 0, A->width);
  for (size_t i = 0; i < b->height; i++)
    augAb->data[(i + 1) * augAb->width - 1] = b->data[i];

  houseHolderFactor *aughhf = houseHolderQR(augAb);
  matrix *qr = transposeMatrix(aughhf->qrT);
  //   printf("qr->width = %d, qr->height = %d\n", qr->width, qr->height);
  matrix *Rnn = makeMatrix(n, n);
  copyData(qr, Rnn, 0, n, 0, n);
  matrix *Rnplus1 = makeMatrix(1, n);
  for (size_t i = 0; i < n; i++)
    Rnplus1->data[i] = qr->data[(i + 1) * (n + 1) - 1];

  matrix *x2 = solveSystemR(Rnn, Rnplus1);
  printf("The (ii) gives result c14 = %.17f\n\n", x->data[14]);
}

int main(int argc, char **argv) {

  a_b_();

  c_();

  d_();
  e_();
  return 0;
}
