#pragma once

typedef struct function {
    double *data[2];
    unsigned size;
    unsigned lengthByCoord[3];
    unsigned char currentBufferIndex;
} Function;

Function* initFunction(const unsigned size, const unsigned lengthByCoord[3], int *status);
void destroyFunction(Function *func);
double getPrevious(const Function *func, const int coords[3]);
double* getCurrentBuffer(const Function *func);
void put(const Function *func, const int coords[3], const double value);
double getCurrent(const Function *func, const int coords[3]);

void inline swap(Function *func) {
    func->currentBufferIndex = (func->currentBufferIndex + 1) % 2;
}
