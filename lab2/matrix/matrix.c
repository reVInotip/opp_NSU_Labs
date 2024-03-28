#include "matrix.h"
#include "../vector/vector.h"
#include "../settings.h"
#include <assert.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>

void InitMatrix(Matrix* matr, const unsigned blockLength, unsigned countBlocksInRow) {
    assert(matr);
    matr->size = (unsigned long)blockLength * blockLength * countBlocksInRow * countBlocksInRow;
    matr->matrix = calloc(
        matr->size,
        sizeof(int));
    matr->blockLength = blockLength;
    matr->countBlocksInRow = countBlocksInRow;

    unsigned offset = blockLength * countBlocksInRow;
    int currentOffset = 0;
    unsigned long fakeMatrixOneSize =
        (unsigned long)blockLength * blockLength *
        (countBlocksInRow - 1) * countBlocksInRow;
    long fakeMatrixTwoSize =
        (long)blockLength * blockLength * countBlocksInRow;
    for (int i = 0; i < matr->size; i += (int)offset) {
        if (i + currentOffset < matr->size) {
            matr->matrix[i + currentOffset] = -4;
        }
        if (
            (i + 1 + currentOffset < matr->size) &&
            (currentOffset % blockLength != blockLength - 1)
        ) {
            matr->matrix[i + 1 + currentOffset] = 1;
        }
        if (
            (i - 1 + currentOffset > 0) &&
            (currentOffset % blockLength != 0)
        ) {
            matr->matrix[i - 1 + currentOffset] = 1;
        }
        if (i + blockLength + currentOffset < fakeMatrixOneSize) {
            matr->matrix[i + blockLength + currentOffset] = 1;
        }
        if (i - (int)blockLength + currentOffset >= fakeMatrixTwoSize) {
            matr->matrix[i - blockLength + currentOffset] = 1;
        }
        ++currentOffset;
    }
}

void DestroyMatrix(Matrix* matr) {
    assert(matr);
    assert(matr->matrix);
    free(matr->matrix);
    matr->matrix = NULL;
    matr->blockLength = 0;
    matr->countBlocksInRow = 0;
}

void MultMatrixOnVector(Vector* result, const Vector* vec, const Matrix* matr) {
    assert(result);
    assert(vec);
    assert(matr);

    unsigned matrRowLength = matr->blockLength * matr->countBlocksInRow;
    if (
        (matrRowLength != vec->length) ||
        (result->length != vec->length)
    ) {
        printf("bad args\n");
        return;
    }

    #ifdef OpenMP_V1
    #pragma omp parallel for schedule(static) collapse(2)
    #elif defined(OpenMP_V2)
    #pragma omp for schedule(static) collapse(2)
    #endif
    for (int i = 0; i < matrRowLength; ++i) {
        for (int j = 0; j < matrRowLength; ++j) {
            result->data[i] += vec->data[j] * (double)matr->matrix[matrRowLength * i + j];
        }
    }
}

void PrintMatrix(const Matrix* matr) {
    assert(matr);
    for (unsigned i = 0; i < matr->size; ++i) {
        if (i % (matr->blockLength * matr->countBlocksInRow) == 0) {
            printf("\n");
        }
        printf("%d ", matr->matrix[i]);
    }
    printf("\n");
}

#ifdef STATIC
void SumRows(int* sumRows, const Matrix* matr) {
    const unsigned matrRowLength = matr->blockLength * matr->countBlocksInRow;

    for (unsigned i = 0; i < matr->blockLength * matr->countBlocksInRow; ++i) {
        for (unsigned j = 0; j < matrRowLength; ++j) {
            sumRows[i] += matr->matrix[matrRowLength * i + j];
        }
    }
}
#endif
