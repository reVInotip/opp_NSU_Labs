#pragma once

typedef struct matrix {
    double* data;
    unsigned width;
    unsigned height;
    unsigned size;
}* Matrix;

enum ERR_CODE_MATRIX {
    FILE_NOT_FOUND_ERR = 10,
    ALLOCATE_MEMORY_ERR = 11,
    NOT_FOUND_MATRIX_WIDTH_ERR = 12,
    NOT_FOUND_MATRIX_HEIGHT_ERR = 13,
    DATA_LESS_THEN_MATRIX_SIZE_ERR = 14
};

Matrix createMatrixFromFile(const char *fileName, int *errCode);
Matrix createMatrix(const unsigned width, const unsigned height, int *errCode);
void destroyMatrix(Matrix matrix);
void multMatrixOnMatrix(Matrix result, Matrix first, Matrix second);
void printMatrix(const Matrix matrix);