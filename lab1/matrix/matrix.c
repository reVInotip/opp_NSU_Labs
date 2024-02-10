#include "matrix.h"
#include <assert.h>
#include <stdlib.h>

void CreateMatrix(Matrix* matr, unsigned width, unsigned length) {
    assert(matr);
    int** matrix = malloc(width * length * sizeof(int));
    matr->matrix = matrix;
    matr->length = length;
    matr->width = width;
}

void DestroyMatrix(Matrix* matr) {
    assert(matr);
    free(matr->matrix);
    matr->length = 0;
    matr->width = 0;
}

void FillMatrix(Matrix* matr) {

}