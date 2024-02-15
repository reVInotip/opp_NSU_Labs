#include "vectror.h"
#include <assert.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>

void InitVector(Vector* vec, const unsigned length) {
    assert(vec);
    vec->data = calloc((unsigned long)length, sizeof(int));
    vec->length = length;
}

void DestoryVector(Vector* vec) {
    assert(vec);
    assert(vec->data);
    free(vec->data);
    vec->data = NULL;
    vec->length = 0;
}

void FillVector(Vector *vec) {
    assert(vec);
    assert(vec->data);
    srand(time(NULL));

    const unsigned countNumbers = (unsigned int)(vec->length * 0.2);
    unsigned index = 0;
    int value = 0;
    for (unsigned i = 0; i < countNumbers; ++i) {
        index = rand() % (vec->length);
        value = (rand() % 100) - 50;
        vec->data[index] = value;
    }
}

void PrintVector(const Vector *vec) {
    assert(vec);
    assert(vec->data);

    for (unsigned i = 0; i < vec->length; ++i) {
        printf("%lf ", vec->data[i]);
    }
    printf("\n");
}

double ScalarProduct(const Vector* vec1, const Vector* vec2) {
    assert(vec1);
    assert(vec1->data);
    assert(vec2);
    assert(vec2->data);

    const unsigned minLength = vec1->length > vec2->length ? vec2->length : vec1->length;
    double result = 0;
    for (int i = 0; i < minLength; ++i) {
        result += vec1->data[i] * vec2->data[i];
    }

    return result;
}


void Subtraction(Vector *result, const Vector *vec1, const Vector *vec2) {
    assert(result);
    assert(vec1);
    assert(vec2);

    unsigned minLength = vec1->length > vec2->length ? vec2->length : vec1->length;
    if (minLength > result->length) {
        minLength = result->length;
    }

    for (int i = 0; i < minLength; ++i) {
        result->data[i] = vec1->data[i] - vec2->data[i];
    }
}


void Addition(Vector* result, const Vector* vec1, const Vector* vec2) {
    assert(result);
    assert(vec1);
    assert(vec2);

    unsigned minLength = vec1->length > vec2->length ? vec2->length : vec1->length;
    if (minLength > result->length) {
        minLength = result->length;
    }

    for (int i = 0; i < minLength; ++i) {
        result->data[i] = vec1->data[i] + vec2->data[i];
    }
}

double Norm(const Vector* vec) {
    assert(vec);
    assert(vec->data);

    double result = 0;
    for (int i = 0; i < vec->length; ++i) {
        result += vec->data[i] * vec->data[i];
    }

    result = sqrt(result);
    return result;
}

void Copy(Vector* result, const Vector* vec) {
    assert(result);
    assert(vec);

    if (vec->length != result->length) {
        double* newData = realloc(result->data, vec->length * sizeof(int));
        if (newData) {
            result->data = newData;
        } else {
            return;
        }
    }
    result->length = vec->length;

    for (int i = 0; i < vec->length; ++i) {
        result->data[i] = vec->data[i];
    }
}

void MultOnConst(Vector* result, const double value) {
    assert(result);

    for (int i = 0; i < result->length; ++i) {
        result->data[i] *= value;
    }
}
