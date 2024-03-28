#include "vector.h"
#include "../settings.h"
#include <assert.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>

#define FAKE_BORDER 100
#define REAL_BORDER 50

void InitVector(Vector* vec, const unsigned length) {
    assert(vec);

    vec->data = calloc((unsigned long)length, sizeof(double));
    vec->length = length;
}

void DestoryVector(Vector* vec) {
    assert(vec);
    assert(vec->data);

    free(vec->data);
    vec->data = NULL;
    vec->length = 0;
}

#ifdef RANDOM
void FillVector(Vector *vec) {
    assert(vec);
    assert(vec->data);
    srand(time(NULL));

    const unsigned countNumbers = (unsigned int)(vec->length * 0.2) + 1;
    unsigned index = 0;
    int value = 0;
    for (unsigned i = 0; i < countNumbers; ++i) {
        index = rand() % (vec->length);
        value = (rand() % FAKE_BORDER) - REAL_BORDER;
        vec->data[index] = value;
    }
}
#elif defined(STATIC)
void FillVector(Vector *vec, const int* rowSums) {
    assert(vec);
    assert(vec->data);

    for (unsigned i = 0; i < vec->length; ++i) {
        vec->data[i] = rowSums[i];
    }
}
#endif

void SubtractionWithMultOnConst(
    Vector* result, Vector* vec1, Vector* vec2,
    const double constantForVec1, const double constantForVec2
) {
    assert(result);
    assert(vec1);
    assert(vec2);

    unsigned minLength = vec1->length > vec2->length ? vec2->length : vec1->length;
    if (minLength > result->length) {
        minLength = result->length;
    }

    #ifdef OpenMP_V1
    #pragma omp parallel for schedule(static)
    #elif defined(OpenMP_V2)
    #pragma omp for schedule(static)
    #endif
    for (unsigned i = 0; i < minLength; ++i) {
        vec1->data[i] *= constantForVec1;
        vec2->data[i] *= constantForVec2;
        result->data[i] = vec1->data[i] - vec2->data[i];
    }
}

void AdditionWithMultOnConst(
    Vector* result, Vector* vec1, Vector* vec2,
    const double constantForVec1, const double constantForVec2
) {
    assert(result);
    assert(vec1);
    assert(vec2);

    unsigned minLength = vec1->length > vec2->length ? vec2->length : vec1->length;
    if (minLength > result->length) {
        minLength = result->length;
    }

    #ifdef OpenMP_V1
    #pragma omp parallel for schedule(static)
    #elif defined(OpenMP_V2)
    #pragma omp for schedule(static)
    #endif
    for (unsigned i = 0; i < minLength; ++i) {
        vec1->data[i] *= constantForVec1;
        vec2->data[i] *= constantForVec2;
        result->data[i] = vec1->data[i] + vec2->data[i];
    }
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

    #ifdef OpenMP_V1
    #pragma omp parallel for schedule(static)
    #elif defined(OpenMP_V2)
    #pragma omp for schedule(static)
    #endif
    for (unsigned i = 0; i < vec->length; ++i) {
        result->data[i] = vec->data[i];
    }
}

double ScalarProduct(const Vector* vec1, const Vector* vec2) {
    assert(vec1);
    assert(vec1->data);
    assert(vec2);
    assert(vec2->data);

    const unsigned minLength = vec1->length > vec2->length ? vec2->length : vec1->length;
    double result = 0;

    #pragma omp parallel for schedule(static) reduction(+:result)
    for (int i = 0; i < minLength; ++i) {
        result += vec1->data[i] * vec2->data[i];
    }

    return result;
}

void MultOnConst(Vector* result, const double value) {
    assert(result);

    #ifdef OpenMP_V1
    #pragma omp parallel for schedule(static)
    #elif defined(OpenMP_V2)
    #pragma omp for schedule(static)
    #endif
    for (int i = 0; i < result->length; ++i) {
        result->data[i] *= value;
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
