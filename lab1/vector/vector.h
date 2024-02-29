#pragma once

#include "../settings.h"

#if defined (MPI_V1) || defined (MPI_V2)
    #include "../shift_table/shift_table.h"
#endif

typedef struct {
    double* data;
    unsigned length;
    #ifdef MPI_V1
        unsigned realLength;
    #endif
    #if defined (MPI_V1) || defined (MPI_V2)
        const ShiftTable* table;
    #endif
} Vector;

#ifdef DEFAULT
    void InitVector(Vector* vec, const unsigned length);
    double ScalarProduct(const Vector* vec1, const Vector* vec2);
    double Norm(const Vector* vec);
    void Copy(Vector* result, const Vector* vec);
#elifdef MPI_V1
    void InitVector(Vector* vec, const ShiftTable* table, const unsigned length);
    double ScalarProduct(Vector* vec1, Vector* vec2);
    double Norm(const Vector* vec);
    void Copy(Vector* result, const Vector* vec);
#elifdef MPI_V2
    void InitVector(Vector* vec, const ShiftTable* table, const unsigned length);
    double ScalarProduct(Vector* vec1, Vector* vec2);
    double Norm(Vector* vec);
    void Copy(Vector* result, Vector* vec);
#endif
void DestoryVector(Vector* vec);
#ifdef RANDOM
    void FillVector(Vector* vec);
#elifdef STATIC
    void FillVector(Vector* vec, const int* rowSums);
#endif
void PrintVector(const Vector* vec);
void Subtraction(Vector* result, const Vector* vec1, const Vector* vec2);
void Addition(Vector* result, const Vector* vec1, const Vector* vec2);
void MultOnConst(Vector* result, const double value);
