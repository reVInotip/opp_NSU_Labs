#pragma once

#include <stdbool.h>

typedef struct matrix {
    double* data;
    unsigned width;
    unsigned height;
    unsigned size;
}* Matrix;

enum ERR_CODE_MATRIX {
    ALLOCATE_MEMORY_ERR = 10,
    DATA_LESS_THEN_MATRIX_SIZE_ERR = 11
};

Matrix createMatrix(unsigned height, unsigned width, int *errCode);
void fillMatrix(Matrix matrix, int seed);
void resizeMatrix(Matrix matrix, const unsigned width, const unsigned height, int *errCode);
void destroyMatrix(Matrix matrix);
void multTransposeMatrixOnMatrix(Matrix result, const Matrix first, const Matrix second);
void multMatrixOnMatrix(Matrix result, const Matrix first, const Matrix second);
void printMatrix(const Matrix matrix);
void copy(Matrix dest, const Matrix matrix);
bool isCorrectCalcualtion(const Matrix first, const Matrix second);