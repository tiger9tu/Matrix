#ifndef __ISOMAP_MATRIX
#define __ISOMAP_MATRIX

/*
 * This is my personal matrix C library.
 * Now I understand why the NR library was created.
 * I didn't realize how quickly I would run into
 * memory problems when working with matricies.
 *
 * None of the matrix operations below change the
 * input matrices if an input is required. They all
 * return a new matrix with the new changes.
 *
 * Because memory issues are so prevelant, be sure
 * to use the "freeMatrix" function to reclaim some
 * of that memory.
 *
 * Most of the algorithms come from the following books:
 *
 * Numerical Recipes in C (2nd Edition)
 * 1992
 * By Press, Teukolsky, Vetterling, and Flannery
 *
 * Matric Computations (3rd Edition)
 * 1996
 * By Golub and Loan
 */

typedef struct _matrix {
  int height;
  int width;
  double *data;
} matrix;

//============================
// Stuff I wish C had.
//============================
void assert(int assertion, char *message);

//============================
// Catch and release functions
//============================
matrix *readMatrix(char *filename);
matrix *makeMatrix(int width, int height);
matrix *copyMatrix(matrix *m);
void copyData(matrix *from, matrix *to, int startRow, int endRow, int startCol,
              int endCol);
void freeMatrix(matrix *m);
void writeMatrix(matrix *m, char *filename);
void printMatrix(matrix *m);

//============================
// Basic Matrix operations
//============================
matrix *eyeMatrix(int n);
double traceMatrix(matrix *m);
matrix *transposeMatrix(matrix *m);
matrix *meanMatrix(matrix *m);
matrix *multiplyMatrix(matrix *a, matrix *b);
matrix *addMatrix(matrix *a, matrix *b);
void plusMatrix(matrix *a, matrix *b);
matrix *scaleMatrix(matrix *m, double value);
matrix *covarianceMatrix(matrix *m);
void rowSwap(matrix *a, int p, int q); // This method changes the input matrix.
double innerProductVector(matrix *x, matrix *y);
matrix *dotProductMatrix(matrix *a, matrix *b);
matrix *dotDiagonalMatrix(matrix *a, matrix *b);
matrix *L2_distance(matrix *a, matrix *b);
matrix *subVectorRef(matrix *a, int s, int e);
void rescaleMatrix(matrix *m, double scale);
void rescaleMatrixAdd(matrix *x, matrix *y, double xScale, double yScale);
void matrixAdd(matrix *x, matrix *y);
void matrixMinus(matrix *x, matrix *y);

double dotProduct(const double *a, const double *b, int length, int strideA,
                  int strideB);
void set_element(matrix *mat, int row, int col, double value);
double get_element(matrix *mat, int row, int col);
#endif
