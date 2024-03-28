#pragma once

#include "../settings.h"

typedef struct {
    double* data;
    unsigned length;
} Vector;

void InitVector(Vector* vec, const unsigned length);
double ScalarProduct(const Vector* vec1, const Vector* vec2);
void Copy(Vector* result, const Vector* vec);
void DestoryVector(Vector* vec);

#ifdef RANDOM
void FillVector(Vector* vec);
#elif defined(STATIC)
void FillVector(Vector* vec, const int* rowSums);
#endif

void PrintVector(const Vector* vec);
void SubtractionWithMultOnConst(Vector* result, Vector* vec1, Vector* vec2,
    const double constantForVec1, const double constantForVec2);
void AdditionWithMultOnConst(Vector* result, Vector* vec1, Vector* vec2,
    const double constantForVec1, const double constantForVec2);
