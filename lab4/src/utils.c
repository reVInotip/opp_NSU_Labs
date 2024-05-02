#include "../include/utils.h"
#include "../include/configuration.h"
#include "../include/node_function.h"
#include "../include/function.h"
#include <mpi.h>
#include <stdio.h>
#include <math.h>

void isCorrectCalculation(const double h[3], const NodeFunction *phi, const int rank, const int startZCoord) {
    double delta = 0, maxDelta = 0;

    for (int k = 0, z = startZCoord; k < (int)phi->lengthByZCoord; ++k, ++z) {
        for (int j = 0; j < N_y; ++j) {
            for (int i = 0; i < N_x; ++i) {

                delta = fabs(getFunctionValue(i, j, z, h) - getCurrentValue(phi, i, j, k));

                if (delta > maxDelta) {
                    maxDelta = delta;
                }

                if (delta > DELTA_MAX) {
                    printf("Function value: %lf, result: %lf, node (%f, %f, %f)\n",
                        getFunctionValue(i, j, z, h), getCurrentValue(phi, i, j, k), i * h[0] - 1, j * h[1] - 1, z * h[2] - 1);
                }
            }
        }
    }

    MPI_Allreduce(MPI_IN_PLACE, &maxDelta, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    if (!rank) {
        if (maxDelta > DELTA_MAX) {
            printf("INCORRECT! Max delta: %lf\n", maxDelta);
        } else {
            printf("CORRECT! Max delta: %lf\n", maxDelta);
        }
    }
}