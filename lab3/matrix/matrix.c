#include "matrix.h"
#include <stdio.h>
#include <errno.h>
#include <assert.h>
#include <stdlib.h>

// 0 - width, 1 - height
void getMatrSizeParametrs(FILE *dataFile, int parametrs[2]) {
    char character;
    int isEOF = 0;
    while ((character = fgetc(dataFile)) != EOF) {
        if (character == 'w') {
            isEOF = fscanf(dataFile, " %d", &parametrs[0]);
        } else if (character == 'h') {
            isEOF = fscanf(dataFile, " %d", &parametrs[1]);
        }

        if (isEOF == EOF) {
            break;
        }

        if ((parametrs[0] != -1) && (parametrs[1] != -1)) {
            break;
        }
    }
}

Matrix createMatrixFromFile(const char *fileName, int *errCode) {
    assert(fileName);
    assert(errCode);

    FILE *dataFile = fopen(fileName, "r");
    if (!dataFile) {
        fprintf(stderr, "File %s not found!\n", fileName);
        *errCode = FILE_NOT_FOUND_ERR;
        return NULL;
    }

    int matrParametrs[2] = {-1, -1};
    getMatrSizeParametrs(dataFile, matrParametrs);

    if (matrParametrs[0] == -1) {
        fclose(dataFile);
        fprintf(stderr, "Not found matrix width error!\n");
        *errCode = NOT_FOUND_MATRIX_WIDTH_ERR;
        return NULL;
    } else if (matrParametrs[0] == -1) {
        fclose(dataFile);
        fprintf(stderr, "Not found matrix height error!\n");
        *errCode = NOT_FOUND_MATRIX_HEIGHT_ERR;
        return NULL;
    }

    Matrix matrix = (Matrix)malloc(sizeof(struct matrix));
    if (!matrix) {
        fclose(dataFile);
        fprintf(stderr, "Allocate memory error!\n");
        *errCode = ALLOCATE_MEMORY_ERR;
        return NULL;
    }

    matrix->width = matrParametrs[0];
    matrix->height = matrParametrs[1];
    matrix->size = matrParametrs[0] * matrParametrs[1];

    double *buffer = (double*)malloc(matrix->size * sizeof(double));
    if (!buffer) {
        free(matrix);
        fclose(dataFile);
        fprintf(stderr, "Allocate memory error!\n");
        *errCode = ALLOCATE_MEMORY_ERR;
        return NULL;
    }

    for (unsigned i = 0; i < matrix->size; ++i) {
        if (fscanf(dataFile, "%lf", &buffer[i]) == EOF) {
            free(matrix);
            free(buffer);
            fclose(dataFile);
            fprintf(stderr, "Count elemets is not equal to matrix size!\n");
            *errCode = DATA_LESS_THEN_MATRIX_SIZE_ERR;
            return NULL;
        }
    }

    fclose(dataFile);

    matrix->data = buffer;
    return matrix;
}

Matrix createMatrix(const unsigned width, const unsigned height, int *errCode) {
    assert(errCode);

    Matrix matrix = (Matrix)malloc(sizeof(struct matrix));
    if (!matrix) {
        fprintf(stderr, "Allocate memory error!\n");
        *errCode = ALLOCATE_MEMORY_ERR;
        return NULL;
    }

    matrix->width = width;
    matrix->height = height;
    matrix->size = width * height;

    double *buffer = (double*)malloc(matrix->size * sizeof(double));
    if (!buffer) {
        free(matrix);
        fprintf(stderr, "Allocate memory error!\n");
        *errCode = ALLOCATE_MEMORY_ERR;
        return NULL;
    }

    matrix->data = buffer;
    return matrix;
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