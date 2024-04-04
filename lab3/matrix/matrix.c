#include "matrix.h"
#include <stdio.h>
#include <errno.h>
#include <assert.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>

Matrix createMatrix(unsigned height, unsigned width, int *errCode) {
    assert(errCode);

    Matrix matrix = (Matrix)malloc(sizeof(struct matrix));
    if (!matrix) {
        fprintf(stderr, "Allocate memory error!\n");
        *errCode = ALLOCATE_MEMORY_ERR;
        return NULL;
    }

    matrix->height = height;
    matrix->width = width;
    matrix->size = height * width;

    if (matrix->size) {
        double *buffer = (double*)malloc(matrix->size * sizeof(double));
        if (!buffer) {
            free(matrix);
            fprintf(stderr, "Allocate memory error!\n");
            *errCode = ALLOCATE_MEMORY_ERR;
            return NULL;
        }
        matrix->data = buffer;
    } else {
        matrix->data = NULL;
    }

    return matrix;
}

void resizeMatrix(Matrix matrix, const unsigned width, const unsigned height, int *errCode) {
    assert(errCode);

    matrix->width = width;
    matrix->height = height;
    matrix->size = width * height;

    double *buffer = (double*)realloc(matrix->data, matrix->size * sizeof(double));
    if (!buffer) {
        free(matrix);
        fprintf(stderr, "Allocate memory error!\n");
        *errCode = ALLOCATE_MEMORY_ERR;
        return;
    }

    matrix->data = buffer;
}

void fillMatrix(Matrix matrix, int seed) {
    assert(matrix);

    srand(seed);
    for (int i = 0; i < matrix->size; ++i) {
        matrix->data[i] = (double)((rand() % 50) - 25);
    }
}

void destroyMatrix(Matrix matrix) {
    if (matrix == NULL) {
        return;
    }

    free(matrix->data);
    free(matrix);
    matrix = NULL;
}

void printMatrix(const Matrix matrix) {
    assert(matrix);
    for (unsigned i = 0; i < matrix->size; ++i) {
        if (i % (matrix->width) == 0) {
            printf("\n");
        }
        printf("%lf ", matrix->data[i]);
    }
    printf("\n");
}

void copy(Matrix dest, const Matrix matrix) {
    for (int i = 0; i < matrix->height; ++i) {
        for (int j = 0; j < matrix->width; ++j) {
            dest->data[i * dest->width + j] = matrix->data[i * matrix->width + j];
        }
    }
}

void multTransposeMatrixOnMatrix(Matrix result, const Matrix first, const Matrix second) {
    unsigned width = (first->width > second->width) * second->width + (first->width <= second->width) * first->width;
    for (int col = 0; col < result->height; ++col) {
        for (int row = 0; row < result->width; ++row) {
            result->data[col * result->width + row] = 0;
            for (int k = 0; k < width; ++k) {
                result->data[col * result->width + row] += first->data[col * first->width + k] * second->data[row * second->width + k];
            }
        }
    }
}

void multMatrixOnMatrix(Matrix result, const Matrix first, const Matrix second) {
    unsigned width = (first->width > second->width) * first->width + (first->width <= second->width) * second->width;
    for (int col = 0; col < result->height; ++col) {
        for (int row = 0; row < result->width; ++row) {
            result->data[col * result->width + row] = 0;
            for (int k = 0; k < width; ++k) {
                result->data[col * result->width + row] += first->data[col * first->width + k] * second->data[k * second->width + row];
            }
        }
    }
}

bool isCorrectCalcualtion(const Matrix first, const Matrix second) {
    bool isCorrect = false;
    unsigned width = (first->width == second->width) * first->width;
    unsigned height = (first->height == second->height) * first->height;
    if (!width || !height) {
        return false;
    }

    for (int col = 0; col < height; ++col) {
        for (int row = 0; row < width; ++row) {
            if (first->data[col * width + row] != second->data[col * width + row]) {
                return false;
            }
        }
    }

    return true;
}