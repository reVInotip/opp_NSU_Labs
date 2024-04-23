#include "function.h"
#include <stdlib.h>
#include <assert.h>

enum ERR_CODE_FUNC {
    MemoryAllocateError = 1,
    OK = 0
};

Function initFunction(const unsigned size, const unsigned lengthByCoord[3], int *status) {
    if (size == 0) {
        *status = OK;
        return NULL;
    }

    Function func = (Function)malloc(sizeof(struct function));
    if (!func) {
        *status = MemoryAllocateError;
        return NULL;
    }

    double *buf1 = (double*)calloc(size, sizeof(double));
    double *buf2 = (double*)calloc(size, sizeof(double));
    if (!buf1 || !buf2) {
        *status = MemoryAllocateError;
        if (buf1) {
            free(buf1);
        } else if (buf2) {
            free(buf2);
        }
        
        free(func);
        return NULL;
    }

    func->data[0] = buf1;
    func->data[1] = buf2;
    func->size = size;
    func->currentBufferIndex = 0;
    func->lengthByCoord[0] = lengthByCoord[0];
    func->lengthByCoord[1] = lengthByCoord[1];
    func->lengthByCoord[2] = lengthByCoord[2];

    return func;
}

void destroyFunction(Function func) {
    assert(func);

    free(func->data[0]);
    free(func->data[1]);
    free(func);
}

double getPrevious(const Function func, const int coords[3]) {
    int index = coords[2] * func->lengthByCoord[0] * func->lengthByCoord[1] + coords[1] * func->lengthByCoord[1] + coords[0];
    if (coords[0] < 0 || coords[1] < 0 || coords[2] < 0 || index >= func->size) {
        return 0;
    }

    return func->data[(func->currentBufferIndex + 1) % 2][index];
}

double* getCurrentBuffer(const Function func) {
    return func->data[func->currentBufferIndex];
}

void put(const Function func, const int coords[3], const double value) {
    int index = coords[2] * func->lengthByCoord[0] * func->lengthByCoord[1] + coords[1] * func->lengthByCoord[1] + coords[0];
    if (coords[0] < 0 || coords[1] < 0 || coords[2] < 0 || index >= func->size) {
        return;
    }

    func->data[func->currentBufferIndex][index] = value;
}