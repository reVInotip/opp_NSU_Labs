#pragma once

enum ERR_CODE_MATR {
    OK
};

typedef struct {
    int** matrix;
    unsigned width;
    unsigned length;
} Matrix;

void CreateMatrix(Matrix* matr, unsigned width, unsigned length);
void DestroyMatrix(Matrix* matr);
void FillMatrix(Matrix* matr);