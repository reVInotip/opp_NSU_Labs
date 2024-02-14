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
}

void DestroyMatrix(Matrix* matr) {
    assert(matr);
    free(matr->matrix);
    matr->blockLength = 0;
    matr->countBlocksInRow = 0;
}

void FillMatrix(Matrix* matr) {
    assert(matr);
    unsigned offset = matr->blockLength * matr->countBlocksInRow;
    int currentOffset = 0;
    unsigned long fakeMatrixOneSize =
        (unsigned long)matr->blockLength * matr->blockLength *
        (matr->countBlocksInRow - 1) * matr->countBlocksInRow;
    long fakeMatrixTwoSize =
        (long)matr->blockLength * matr->blockLength *
        matr->countBlocksInRow;
    for (int i = 0; i < matr->size; i += (int)offset) {
        if (i + currentOffset < matr->size) {
            matr->matrix[i + currentOffset] = -4;
        }
        if (
            (i + 1 + currentOffset < matr->size) &&
            (currentOffset % matr->blockLength != matr->blockLength - 1)
        ) {
            matr->matrix[i + 1 + currentOffset] = 1;
        }
        if (
            (i - 1 + currentOffset > 0) &&
            (currentOffset % matr->blockLength != 0)
        ) {
            matr->matrix[i - 1 + currentOffset] = 1;
        }
        if (i + matr->blockLength + currentOffset < fakeMatrixOneSize) {
            matr->matrix[i + matr->blockLength + currentOffset] = 1;
        }
        if (i - (int)matr->blockLength + currentOffset >= fakeMatrixTwoSize) {
            matr->matrix[i - matr->blockLength + currentOffset] = 1;
        }
        ++currentOffset;
    }
}

void PrintMatrix(Matrix* matr) {
    assert(matr);
    for (unsigned i = 0; i < matr->size; ++i) {
        if (i % (matr->blockLength * matr->countBlocksInRow) == 0) {
            printf("\n");
        }
        printf("%d ", matr->matrix[i]);
    }
    printf("\n");
}
