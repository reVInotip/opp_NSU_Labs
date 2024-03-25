#pragma once

#include "../settings.h"


typedef struct {
    double* data;
    unsigned length;
} Vector;

void InitVector(Vector* vec, const unsigned length);
double ScalarProduct(const Vector* vec1, const Vector* vec2);
double Norm(const Vector* vec);
void Copy(Vector* result, const Vector* vec);
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
