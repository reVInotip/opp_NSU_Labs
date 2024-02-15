#pragma once

typedef struct {
    double* data;
    unsigned length;
} Vector;

void InitVector(Vector* vec, const unsigned length);
void DestoryVector(Vector* vec);
void FillVector(Vector* vec);
void PrintVector(const Vector* vec);
double ScalarProduct(const Vector* vec1, const Vector* vec2);
void Subtraction(Vector* result, const Vector* vec1, const Vector* vec2);
void Addition(Vector* result, const Vector* vec1, const Vector* vec2);
double Norm(const Vector* vec);
void Copy(Vector* result, const Vector* vec);
void MultOnConst(Vector* result, const double value);
