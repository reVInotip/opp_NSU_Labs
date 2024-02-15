#pragma once

typedef struct {
    int* data;
    unsigned length;
} Vector;

void InitVector(Vector* vec, const unsigned length);
void DestoryVector(Vector* vec);
void FillVector(Vector* vec);
void PrintVector(const Vector* vec);
