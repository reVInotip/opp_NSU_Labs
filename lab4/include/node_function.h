#pragma once

typedef struct function {
    double *data[2];
    unsigned size;
    unsigned lengthByZCoord;
    unsigned char currentBufferIndex;
} NodeFunction;

NodeFunction* initNodeFunction(const unsigned size, const int totalProcessesNumber, const int rank, int *status);
void destroyNodeFunction(NodeFunction *func);
double getPreviousValue(const NodeFunction *func, const int x, const int y, const int z);
double* getCurrentBuffer(const NodeFunction *func);
void put(const NodeFunction *func, const int x, const int y, const int z, const double value);
double getCurrentValue(const NodeFunction *func, const int x, const int y, const int z);

inline void swap(NodeFunction *func) {
    func->currentBufferIndex = (func->currentBufferIndex + 1) % 2;
}
