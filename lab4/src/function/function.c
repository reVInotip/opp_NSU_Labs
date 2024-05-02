#include "../../include/function.h"

void getNodeCoords(const double h[3], const int nodeNumber[3], double nodeCoords[3]) {
    for (int i = 0; i < 3; ++i) {
        nodeCoords[i] = nodeNumber[i] * h[i] - 1; // start coords: (0, 0, 0)
    }
}

double getFunctionValue(const int x, const int y, const int z, const double h[3]){
    double realCoords[3] = {0, 0, 0};
    int coords[3] = {x, y, z};
    getNodeCoords(h, coords, realCoords);
    return realCoords[0] * realCoords[0] + realCoords[1] * realCoords[1] + realCoords[2] * realCoords[2];
}