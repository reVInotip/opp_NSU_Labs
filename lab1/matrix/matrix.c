#include "matrix.h"
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
