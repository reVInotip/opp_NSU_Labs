#include "../../include/node_function.h"
#include "../../include/configuration.h"
#include "../../include/function.h"
#include <stdlib.h>
#include <assert.h>

#include <stdio.h>

enum ERR_CODE_FUNC {
    MemoryAllocateError = 1,
    OK = 0
};

NodeFunction* initNodeFunction(const unsigned size, const int totalProcessesNumber, const int rank, int *status) {
    if (size == 0) {
        *status = OK;
        return NULL;
    }

    NodeFunction *func = (NodeFunction*)malloc(sizeof(NodeFunction));
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
    func->lengthByZCoord = N_z / totalProcessesNumber + (N_z % totalProcessesNumber > rank);

    return func;
}

void destroyNodeFunction(NodeFunction *func) {
    assert(func);

    free(func->data[0]);
    free(func->data[1]);
    free(func);
}

double getPreviousValue(const NodeFunction *func, const int x, const int y, const int z, const double *h) {
    int index = z * N_x * N_y + y * N_x + x;
    if ((x < 0 || y < 0 || z < 0) || (x * h[0] > D_x || y * h[1] > D_y || z * h[2] > D_z)) {
        return getFunctionValue(x, y, z, h);
    }

    if (index >= func->size) {
        return 0;
    }

    return func->data[(func->currentBufferIndex + 1) % 2][index];
}

double getCurrentValue(const NodeFunction *func, const int x, const int y, const int z, const double *h) {
    int index = z * N_x * N_y + y * N_x + x;
    if ((x < 0 || y < 0 || z < 0) || (x * h[0] > D_x || y * h[1] > D_y || z * h[2] > D_z)) {
        return getFunctionValue(x, y, z, h);
    }

    if (index >= func->size) {
        return 0;
    }

    return func->data[func->currentBufferIndex][index];
}

double* getCurrentBuffer(const NodeFunction *func) {
    return func->data[func->currentBufferIndex];
}

void put(const NodeFunction *func, const int x, const int y, const int z, const double value) {
    int index = z * N_x * N_y + y * N_x + x;
    if (x < 0 || y < 0 || z < 0 || index >= func->size) {
        return;
    }

    func->data[func->currentBufferIndex][index] = value;
}