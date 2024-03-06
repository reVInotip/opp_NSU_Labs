#pragma once

#include "../settings.h"
#include "../vector/vector.h"

#ifdef DEFAULT

    typedef struct {
        int* matrix;
        unsigned blockLength;
        unsigned countBlocksInRow;
        unsigned long size;
    } Matrix;

    void InitMatrix(Matrix* matr, const unsigned blockLength, unsigned countBlocksInRow);
#elif defined(MPI_V1) || defined(MPI_V2)
    #include "../shift_table/shift_table.h"

    typedef struct {
        int* matrix;
        const ShiftTable* table;
        unsigned blockLength;
        unsigned blockWidth;
        unsigned countBlocksInRow;
        unsigned long size;
    } Matrix;

    void InitMatrix(Matrix* matr, const ShiftTable* table, const unsigned blockLength, const unsigned blockWidth, const unsigned countBlocksInRow);
#endif

void DestroyMatrix(Matrix* matr);
void PrintMatrix(const Matrix* matr);

#ifdef DEFAULT
    void MultMatrixOnVector(Vector* result, const Vector* vec, const Matrix* matr);
#elif defined(MPI_V1) || defined(MPI_V2)
    void MultMatrixOnVector(Vector* result, Vector* vec, const Matrix* matr);
#endif

#ifdef STATIC
    void SumRows(int* sumRows, const Matrix* matr);
#endif
