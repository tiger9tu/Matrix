#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*===========================================================================
 * assert
 * If the assertion is non-zero (i.e. true), then it returns.
 * If the assertion is zero (i.e. false), then it display the string and
 * aborts the program.
 * This is ment to act like Python's assert keyword.
 *=========================================================================*/
void assert(int assertion, char *message) {
  if (assertion == 0) {
    fprintf(stderr, "%s\n", message);
    exit(1);
  }
}

/*===========================================================================
 * readMatrix
 * Reads a file containing a Matrix
 *=========================================================================*/
matrix *readMatrix(char *filename) {
  FILE *fh;
  matrix *data;
  float myValue;
  int width, height, i, elements;
  int scan_test;
  double *ptr;

  if ((fh = fopen(filename, "r")) == NULL) {
    fprintf(stderr, "Error: Cannot open %s\n", filename);
    exit(1);
  }

  scan_test = fscanf(fh, "%d", &width);
  assert(scan_test != EOF, "Failed to read from file.");

  scan_test = fscanf(fh, "%d", &height);
  assert(scan_test != EOF, "Failed to read from file.");

  data = makeMatrix(width, height);
  elements = width * height;

  ptr = data->data;
  for (i = 0; i < elements; i++) {
    scan_test = fscanf(fh, "%f", &myValue);
    assert(scan_test != EOF, "Failed to read from file. File is incomplete.");
    *(ptr++) = myValue;
  }

  fclose(fh);
  return data;
}

/*===========================================================================
 * makeMatrix
 * Makes a matrix with a width and height parameters.
 *=========================================================================*/
matrix *makeMatrix(int width, int height) {
  matrix *out;
  assert(width > 0 && height > 0, "New matrix must be at least a 1 by 1");
  out = (matrix *)malloc(sizeof(matrix));

  assert(out != NULL, "Out of memory.");

  out->width = width;
  out->height = height;
  out->data = (double *)malloc(sizeof(double) * width * height);

  assert(out->data != NULL, "Out of memory.");

  memset(out->data, 0.0, width * height * sizeof(double));

  return out;
}

/*===========================================================================
 * copyMatrix
 * Copies a matrix. This function uses scaleMatrix, because scaling matrix
 * by 1 is the same as a copy.
 *=========================================================================*/
matrix *copyMatrix(matrix *m) { return scaleMatrix(m, 1); }

void copyData(matrix *from, matrix *to, int startRow, int endRow, int startCol,
              int endCol) {
  //   assert(from->height == to->height && from->width == to->width,
  //          "copyData: matrices must have the same size");
  for (size_t i = startRow; i < endRow; i++)
    for (size_t j = startCol; j < endCol; j++)
      to->data[i * to->width + j] = from->data[i * from->width + j];
}

/*===========================================================================
 * freeMatrix
 * Frees the resources of a matrix
 *=========================================================================*/
void freeMatrix(matrix *m) {
  if (m != NULL) {
    if (m->data != NULL) {
      free(m->data);
      m->data = NULL;
    }

    free(m);
    m = NULL;
  }
  return;
}

/*===========================================================================
 * writeMatrix
 * Write a matrix to a file.
 *=========================================================================*/
void writeMatrix(matrix *m, char *filename) {
  FILE *fh;
  int i, j;
  double *ptr = m->data;

  if ((fh = fopen(filename, "w")) == NULL) {
    fprintf(stderr, "Error: Cannot open %s\n", filename);
    exit(1);
  }

  fprintf(fh, "%d %d\n", m->width, m->height);

  for (i = 0; i < m->height; i++) {
    for (j = 0; j < m->width; j++) {
      fprintf(fh, " %2.5f", *(ptr++));
    }
    fprintf(fh, "\n");
  }

  close(fh);
  return;
}

/*===========================================================================
 * printMatrix
 * Prints a matrix. Great for debugging.
 *=========================================================================*/
void printMatrix(matrix *m) {
  int i, j;
  // double* ptr = m->data;
  printf("%d %d\n", m->width, m->height);
  for (i = 0; i < m->height; i++) {
    for (j = 0; j < m->width; j++) {
      // printf(" %9.6f", *(ptr++));
      printf(" %.6f", m->data[i * m->width + j]);
    }
    printf("\n");
  }
  return;
}

/*===========================================================================
 * eyeMatrix
 * Returns an identity matrix of size n by n, where n is the input
 * parameter.
 *=========================================================================*/
matrix *eyeMatrix(int n) {
  int i;
  matrix *out;
  double *ptr;

  assert(n > 0, "Identity matrix must have value greater than zero.");

  out = makeMatrix(n, n);
  ptr = out->data;
  for (i = 0; i < n; i++) {
    *ptr = 1.0;
    ptr += n + 1;
  }

  return out;
}

/*===========================================================================
 * traceMatrix
 * Given an "m rows by n columns" matrix, it returns the sum of the elements
 * along the diagonal. This is know as the matrix 'trace'.
 *=========================================================================*/
double traceMatrix(matrix *m) {
  int i;
  int size;
  double *ptr = m->data;
  double sum = 0.0;

  if (m->height < m->width) {
    size = m->height;
  } else {
    size = m->width;
  }

  for (i = 0; i < size; i++) {
    sum += *ptr;
    ptr += m->width + 1;
  }

  return sum;
}

/*===========================================================================
 * meanMatrix
 * Given an "m rows by n columns" matrix, it returns a matrix with 1 row and
 * n columns, where each element represents the mean of that full column.
 *=========================================================================*/
matrix *meanMatrix(matrix *m) {
  int i, j;
  double *ptr;
  matrix *out;

  assert(m->height > 0, "Height of matrix cannot be zero.");

  out = makeMatrix(m->width, 1);

  for (i = 0; i < m->width; i++) {
    out->data[i] = 0.0;
    ptr = &m->data[i];
    for (j = 0; j < m->height; j++) {
      out->data[i] += *ptr;
      ptr += out->width;
    }
    out->data[i] /= (double)m->height;
  }
  return out;
}

/*===========================================================================
 * covarianceMatrix
 * Given an "m rows by n columns" matrix, it returns a matrix with n row and
 * n columns, where each element represents covariance of 2 columns.
 *=========================================================================*/
matrix *covarianceMatrix(matrix *m) {
  int i, j, k, L = 0;
  matrix *out;
  matrix *mean;
  double *ptrA;
  double *ptrB;
  double *ptrOut;

  assert(m->height > 1, "Height of matrix cannot be zero or one.");

  mean = meanMatrix(m);
  out = makeMatrix(m->width, m->width);
  ptrOut = out->data;

  for (i = 0; i < m->width; i++) {
    for (j = 0; j < m->width; j++) {
      ptrA = &m->data[i];
      ptrB = &m->data[j];
      *ptrOut = 0.0;
      for (k = 0; k < m->height; k++) {
        *ptrOut += (*ptrA - mean->data[i]) * (*ptrB - mean->data[j]);
        ptrA += m->width;
        ptrB += m->width;
      }
      *ptrOut /= m->height - 1;
      ptrOut++;
    }
  }

  freeMatrix(mean);
  return out;
}

/*===========================================================================
 * transposeMatrix
 * Given an matrix, returns the transpose.
 *=========================================================================*/
matrix *transposeMatrix(matrix *m) {
  matrix *out = makeMatrix(m->height, m->width);
  double *ptrOut;
  double *ptrM = m->data;
  int i, j;

  for (i = 0; i < m->height; i++) {
    ptrOut = &out->data[i];
    for (j = 0; j < m->width; j++) {
      *ptrOut = *ptrM;
      ptrM++;
      ptrOut += out->width;
    }
  }

  return out;
}

/*===========================================================================
 * multiplyMatrix
 * Given a two matrices, returns the multiplication of the two.
 *=========================================================================*/
matrix *multiplyMatrix(matrix *a, matrix *b) {
  int i, j, k;
  matrix *out;
  double *ptrOut;
  double *ptrA;
  double *ptrB;

  assert(a->width == b->height,
         "Matrices have incorrect dimensions. a->width != b->height");

  out = makeMatrix(b->width, a->height);
  ptrOut = out->data;

  for (i = 0; i < a->height; i++) {

    for (j = 0; j < b->width; j++) {
      ptrA = &a->data[i * a->width];
      ptrB = &b->data[j];

      *ptrOut = 0;
      for (k = 0; k < a->width; k++) {
        *ptrOut += *ptrA * *ptrB;
        ptrA++;
        ptrB += b->width;
      }
      ptrOut++;
    }
  }

  return out;
}

matrix *addMatrix(matrix *a, matrix *b) {

  matrix *c;
  assert(a->width == b->width, "Matrices width different");
  assert(a->height == b->height, "Matrices height.");
  c = makeMatrix(a->width, a->height);
  for (int i = 0; i < a->height; i++)
    for (int j = 0; j < a->width; j++) {
      int idx = i * a->width + j;
      c->data[idx] = a->data[idx] + b->data[idx];
    }
  return c;
}

void plusMatrix(matrix *a, matrix *b) {
  assert(a->width == b->width, "Matrices width different");
  assert(a->height == b->height, "Matrices height.");
  for (int i = 0; i < a->height; i++)
    for (int j = 0; j < a->width; j++) {
      int idx = i * a->width + j;
      a->data[idx] += b->data[idx];
    }
}

/*===========================================================================
 * scaleMatrix
 * Given a matrix and a double value, this returns a new matrix where each
 * element in the input matrix is multiplied by the double value
 *=========================================================================*/
matrix *scaleMatrix(matrix *m, double value) {
  int i, elements = m->width * m->height;
  matrix *out = makeMatrix(m->width, m->height);
  double *ptrM = m->data;
  double *ptrOut = out->data;

  for (i = 0; i < elements; i++) {
    *(ptrOut++) = *(ptrM++) * value;
  }

  return out;
}

void rescaleMatrix(matrix *m, double scale) {
  for (int i = 0; i < m->height; i++)
    for (int j = 0; j < m->width; j++)
      m->data[i * m->width + j] *= scale;
}

/*===========================================================================
 * rowSwap
 * Given a matrix, this algorithm will swap rows p and q, provided
 * that p and q are less than or equal to the height of matrix A and p
 * and q are different values.
 *=========================================================================*/
void rowSwap(matrix *a, int p, int q) {
  int i;
  double temp;
  double *pRow;
  double *qRow;

  assert(a->height > 2, "Matrix must have at least two rows to swap.");
  assert(p < a->height && q < a->height,
         "Values p and q must be less than the height of the matrix.");

  // If p and q are equal, do nothing.
  if (p == q) {
    return;
  }

  pRow = a->data + (p * a->width);
  qRow = a->data + (q * a->width);

  // Swap!
  for (i = 0; i < a->width; i++) {
    temp = *pRow;
    *pRow = *qRow;
    *qRow = temp;
    pRow++;
    qRow++;
  }

  return;
}

double innerProductVector(matrix *x, matrix *y) {
  double l2norm = 0;
  for (int j = 0; j < x->height; j++) {
    l2norm += x->data[j] * y->data[j];
  }
  return l2norm;
}

/*===========================================================================
 * subVectorRef
 * s is the start, e is the index after the last element.
 *=========================================================================*/
matrix *subVectorRef(matrix *a, int s, int e) {
  matrix *subVec = (matrix *)malloc(sizeof(matrix));
  int length = e - s;
  subVec->height = length;
  subVec->width = 1;
  subVec->data = a->data + s;
  return subVec;
}

/*===========================================================================
 * dotProductMatrix
 * Given a two matrices (or the same matrix twice) with identical widths and
 * different heights, this method returns a a->height by b->height matrix of
 * the cross product of each matrix.
 *
 * Dot product is essentially the sum-of-squares of two vectors.
 *
 * Also, if the second paramter is NULL, it is assumed that we
 * are performing a cross product with itself.
 *=========================================================================*/
matrix *dotProductMatrix(matrix *a, matrix *b) {
  matrix *out;
  double *ptrOut;
  double *ptrA;
  double *ptrB;
  int i, j, k;

  if (b != NULL) {
    assert(a->width == b->width,
           "Matrices must be of the same dimensionality.");
  }

  // Are we computing the sum of squares of the same matrix?
  if (a == b || b == NULL) {
    b = a; // May not appear safe, but we can do this without risk of losing b.
  }

  out = makeMatrix(b->height, a->height);
  ptrOut = out->data;

  for (i = 0; i < a->height; i++) {
    ptrB = b->data;

    for (j = 0; j < b->height; j++) {
      ptrA = &a->data[i * a->width];

      *ptrOut = 0;
      for (k = 0; k < a->width; k++) {
        *ptrOut += *ptrA * *ptrB;
        ptrA++;
        ptrB++;
      }
      ptrOut++;
    }
  }

  return out;
}

/*===========================================================================
 * matrixDotDiagonal
 * Given a two matrices (or the same matrix twice) with identical widths and
 * heights, this method returns a 1 by a->height matrix of the cross
 * product of each matrix along the diagonal.
 *
 * Dot product is essentially the sum-of-squares of two vectors.
 *
 * If the second paramter is NULL, it is assumed that we
 * are performing a cross product with itself.
 *=========================================================================*/
matrix *dotDiagonalMatrix(matrix *a, matrix *b) {
  matrix *out;
  double *ptrOut;
  double *ptrA;
  double *ptrB;
  int i, j;

  if (b != NULL) {
    assert(a->width == b->width && a->height == b->height,
           "Matrices must be of the same dimensionality.");
  }

  // Are we computing the sum of squares of the same matrix?
  if (a == b || b == NULL) {
    b = a; // May not appear safe, but we can do this without risk of losing b.
  }

  out = makeMatrix(1, a->height);
  ptrOut = out->data;
  ptrA = a->data;
  ptrB = b->data;

  for (i = 0; i < a->height; i++) {
    *ptrOut = 0;
    for (j = 0; j < a->width; j++) {
      *ptrOut += *ptrA * *ptrB;
      ptrA++;
      ptrB++;
    }
    ptrOut++;
  }

  return out;
}

/* x = xscale*x + yscale* y*/
void rescaleMatrixAdd(matrix *x, matrix *y, double xScale, double yScale) {
  // int length = x->height > x->width ? x->height : x->width;
  for (int i = 0; i < x->height; i++) {
    for (size_t j = 0; j < x->width; j++) {
      x->data[i * x->width + j] = xScale * x->data[i * x->width + j] +
                                  yScale * y->data[i * x->width + j];
    }
  }
}

void matrixAdd(matrix *x, matrix *y) { rescaleMatrixAdd(x, y, 1, 1); }

void matrixMinus(matrix *x, matrix *y) { rescaleMatrixAdd(x, y, 1, -1); }

double dotProduct(const double *a, const double *b, int length, int strideA,
                  int strideB) {
  double sum = 0.0;
  for (int i = 0; i < length; i++) {
    sum += (*a) * (*b);
    a += strideA;
    b += strideB;
  }
  return sum;
}

// Function to set an element in the matrix
void set_element(matrix *mat, int row, int col, double value) {
  mat->data[row * mat->width + col] = value;
}

// Function to get an element from the matrix
double get_element(matrix *mat, int row, int col) {
  return mat->data[row * mat->width + col];
}
