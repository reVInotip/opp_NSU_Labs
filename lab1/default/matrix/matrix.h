#pragma once

#include "../vector/vectror.h"

typedef struct {
    int* matrix;
    unsigned blockLength;
    unsigned countBlocksInRow;
    unsigned long size;
} Matrix;

void InitMatrix(Matrix* matr, const unsigned blockLength, unsigned countBlocksInRow);
void DestroyMatrix(Matrix* matr);
void PrintMatrix(const Matrix* matr);
void MultMatrixOnVector(Vector* result, const Vector* vec, const Matrix* matr);
