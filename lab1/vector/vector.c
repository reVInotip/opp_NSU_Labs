#include "vectror.h"
#include <assert.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>

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
        index = 1 + (rand() % (vec->length - 1));
        value = (rand() % 100) - 50;
        vec->data[index] = value;
    }
}

void PrintVector(const Vector *vec) {
    assert(vec);
    assert(vec->data);

    for (unsigned i = 0; i < vec->length; ++i) {
        printf("%d\n", vec->data[i]);
    }
}
