#ifndef __ISOMAP_MATRIX
#define __ISOMAP_MATRIX

typedef struct _matrix {
  int height;
  int width;
  double *data;
} matrix;

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
