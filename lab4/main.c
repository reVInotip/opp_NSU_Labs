/*
    OPP Laboratory work number 4: Parallel implemetation of the Jacobi method in three-dimensional space.
    Copyright by Grigoriy Novikov.

    ⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿
    ⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⠿⠯⠉⠉⠉⠉⠙⠛⢿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿
    ⣿⣿⣿⣿⣿⣿⣿⣟⣓⣒⡂⠄⠄⠄⠄⠄⠄⠄⠄⠄⠈⠙⢿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿
    ⣿⣿⣿⣿⣿⠯⠭⠭⠉⠁⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⠻⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿
    ⣿⣿⣿⣿⠿⠁⠄⠄⠤⠁⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⠹⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿
    ⣿⣿⣿⣯⣭⣍⣉⣁⡀⠄⡀⠐⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⠹⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿
    ⣿⣿⣟⣒⣒⣒⣒⣀⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿
    ⣿⣿⣷⠤⠶⠶⠶⠦⠤⠄⠄⣀⣀⣠⣤⣤⣤⣤⠄⠄⢀⣀⣀⣀⣀⣀⢸⣿⣿⣿⣿⣿⠟⠛⠛⠻⢿⣿⣿⣿
    ⣿⣿⡭⠭⠭⠭⠉⠁⠄⠄⢀⣠⣤⣤⣤⣤⣤⡀⠄⠄⢀⣭⠥⠤⢤⣄⠸⣿⣿⣿⡿⠁⠄⠄⠄⠄⠄⠹⣿⣿
    ⣿⣿⣯⣭⠍⠁⠄⠄⠄⠄⠸⣇⠄⠄⠘⢁⡼⠃⠄⣤⢼⡀⠈⠁⣠⡟⠄⣿⣿⣿⠁⠄⠄⠄⠄⠄⠄⠄⢹⣿
    ⣿⣭⣭⣽⣧⣤⣀⣀⠄⠐⠂⠈⠙⠓⠚⠋⠄⠄⠄⢿⠄⠙⠓⠛⠋⠄⠄⣿⣿⡇⠄⠄⠄⠄⠄⠄⠄⠄⠘⣿
    ⡟⠛⠛⣷⡒⠒⠒⠒⠂⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⢸⡄⠄⠄⠄⠄⠄⠄⣿⣿⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⣿
    ⣷⣭⣽⣿⣍⡉⠉⠁⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⣧⠄⠄⠄⠄⠄⠄⣿⣿⠄⠄⣿⠛⣿⠄⣿⣠⣿⠄⣿
    ⣿⣗⢿⣾⡷⠶⠤⠤⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⠁⠄⠄⠄⠄⠄⢠⣿⣿⠄⠄⣿⢁⡿⠄⣿⣿⠄⠄⣿
    ⣿⣿⣯⣉⣉⣭⠭⠭⠥⠤⠤⠄⠄⠄⠄⠄⠄⠄⣀⣀⡀⠄⠄⠄⠄⠄⣸⣿⣏⠄⠄⣿⣤⣿⠄⣿⠄⣷⠄⣿
    ⣿⣿⣿⣿⣿⣿⡀⠒⠒⠒⠒⠄⠄⠄⠄⠄⠄⠄⠉⠉⠁⠄⠄⠄⠄⢠⣿⣿⣿⡄⠄⠄⠄⠄⠄⠄⠄⠄⢠⣿
    ⣿⣿⣿⣿⣿⣿⣿⣭⣭⣄⣐⡒⠂⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⢀⣾⣿⣿⣿⡇⠄⠄⠄⠄⠄⠄⠄⠄⢸⣿
    ⣿⣿⣿⣿⣿⣿⣿⣿⣖⠒⠒⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⣠⣾⣿⣿⣿⣿⣷⠄⠄⠄⠄⠄⠄⠄⠄⣿⣿
    ⣿⣿⣿⣿⣿⣿⣿⣿⣿⣷⢦⣀⠤⠄⠄⠄⠄⠄⠄⠄⠄⣠⣶⣿⣿⣿⣿⣿⣿⣿⣆⠄⠄⠄⠄⠄⢀⣼⣿⣿
    ⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣭⣭⣻⡷⠦⣤⣤⣤⠤⠶⢿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣦⣤⣤⣶⣿⣿⣿⣿
    ⣿⣿⣿⣿⣿⣿⣿⣿⢿⣿⣷⠒⠒⠄⠄⠄⠄⠄⠄⠄⢸⡏⠛⢿⡛⢻⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿
    ⣿⣿⣿⣿⣿⠿⢿⣿⠿⣯⠄⠄⠄⠄⠄⠄⠄⠄⠄⠄⢸⣷⡀⠈⢻⡄⢻⡟⢿⡟⠻⣿⣿⣿⣿⣿⣿⣿⣿⣿
    ⣿⣿⣟⠛⣛⡻⢿⣿⡍⢿⣧⣀⠄⠄⠄⠄⠄⠄⠄⠄⠄⣫⡿⣆⠄⢻⣿⣿⣶⣧⣀⡬⠛⢿⣿⣿⣿⣿⣿⣿
    ⠭⠿⠿⣿⣿⣿⣿⣿⣯⠘⣿⠙⠓⠶⣤⣤⣄⣠⣤⣴⠾⠋⠄⠙⣧⡀⢻⣿⣿⣿⡿⠄⣠⣶⣿⣿⣿⣿⣿⣿
    ⣭⡭⠭⠽⢿⣿⣿⣿⠟⣰⣿⡀⠄⠄⠄⠄⣿⢿⣻⠿⣇⠄⠄⣠⡾⠛⠷⠯⣭⡥⠶⠟⠋⠄⣿⣿⣿⣿⣿⣿
    ⣿⡛⠻⢷⣤⣴⡶⠞⠛⠉⠙⠻⢦⣄⣀⣀⣿⡌⠉⠰⣿⣰⡟⠋⠄⠄⠄⠄⠄⠄⠄⠄⡀⠄⢹⣿⣿⣿⣿⣿
    ⠂⠄⠄⠄⠐⣧⢲⣖⠒⠂⠄⠄⠄⠄⠉⢹⣿⠄⠲⢀⣼⢿⡇⠄⠄⠄⠄⠄⠄⠄⠄⣸⠁⠄⢸⣿⣿⣿⣿⣿
    ⠒⠄⠄⠄⠄⣿⣿⡭⠭⠄⠄⠄⠄⠄⠄⢸⡏⠻⣿⣿⠁⢸⡇⠄⠄⠄⠄⠄⠄⠄⢀⣿⡞⠄⠄⣿⣿⣿⣿⣿
    ⠄⠄⠄⠄⠄⠛⠛⠓⠒⠄⠄⠄⠄⠄⠄⠘⠃⠄⠛⠛⠂⠘⠃⠄⠄⠄⠄⠄⠄⠄⠘⠛⠄⠄⠄⠛⠛⠛⠛⠛
*/

#include "parser/parser.h"
#include "function/function.h"
#include <stdlib.h>
#include <assert.h>
#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <float.h>

typedef struct {
    int D[3];
    int N[3];
    double h[3];
    double H[3];
    int a;
    double epsilon;
    double precision;
} Configuration;

/*typedef struct {
    int dims[3];
    int periods[3];
    int coords[3];
    MPI_Comm comm3D;
    MPI_Comm comm2D[3];
} CartSpace;*/

typedef struct buffer {
    double *send;
    double *recieve;
    unsigned size;
} Buffer;

void initBuffer(Buffer *buf, const unsigned size) {
    buf->send = (double*)calloc(size, sizeof(double));
    assert(buf->send);
    buf->recieve = (double*)calloc(size, sizeof(double));
    assert(buf->recieve);

    buf->size = size;
}

void destroyBuffer(Buffer *buf) {
    assert(buf);

    free(buf->send);
    free(buf->recieve);
    buf->size = 0;
}

void initConfig(Configuration *conf, const unsigned counts[3]) {
    assert(conf);
    conf->a = 1e5;
    conf->epsilon = 1e-8;
    conf->precision = 0.5;

    for (int i = 0; i < 3; ++i) {
        conf->D[i] = 2;
        conf->N[i] = counts[i];
        conf->h[i] = (double)2 / (counts[i] - 1);
        conf->H[i] = conf->h[i] * conf->h[i];
    }
}

void isCorrectCalculation(const double *funcValues, const Configuration conf) {
    bool isCorrect = true;
    int coords[3] = {0};
    double delta = 0, middleDelta = 0;

    for (int k = 0; k < conf.N[2]; ++k) {
        coords[2] = k;

        for (int j = 0; j < conf.N[1]; ++j) {
            coords[1] = j;

            for (int i = 0; i < conf.N[0]; ++i) {
                coords[0] = i;
                delta = (i * i + j * j + k * k) - (funcValues[k * conf.N[0] * conf.N[1] + j * conf.N[1] + i]);

                if (delta >= conf.precision) {
                    printf("Incorrect value in (%d, %d, %d) node!\nCorrect: %d, result: %lf\nPrecision: %lf\n",
                        i, j, k, i * i + j * j + k * k, funcValues[k * conf.N[0] * conf.N[1] + j * conf.N[1] + i], delta);
                    isCorrect = false;
                }

                middleDelta += delta;
            }
        }
    }

    if (isCorrect) {
        middleDelta = middleDelta / (conf.N[0] * conf.N[1] * conf.N[2]);
        printf("CORRECT! Middle precision: %lf\n", middleDelta);
    }
}

/*void initCartSpace(CartSpace *space, const int totalProcessesNumber, const int rank) {
    assert(space);
    space->periods[0] = 0; space->periods[1] = 0; space->periods[2] = 0;

    MPI_Dims_create(totalProcessesNumber, 3, space->dims);
    MPI_Cart_create(MPI_COMM_WORLD, 3, space->dims, space->periods, 1, &(space->comm3D));

    int remainDimsForXYSpace[3] = {1, 1, 0};
    int remainDimsForXZSpace[3] = {1, 0, 1};
    int remainDimsForZYSpace[3] = {0, 1, 1};

    MPI_Cart_sub(space->comm3D, remainDimsForXYSpace, &(space->comm2D[0]));
    MPI_Cart_sub(space->comm3D, remainDimsForXZSpace, &(space->comm2D[1]));
    MPI_Cart_sub(space->comm3D, remainDimsForZYSpace, &(space->comm2D[2]));

    MPI_Cart_coords(space->comm3D, rank, 3, space->coords);
}*/

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, totalProcessesNumber;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &totalProcessesNumber);

    unsigned countVirtualNodes[3] = {0, 0, 0};
    parseAgrs(argc, argv, &countVirtualNodes[0], &countVirtualNodes[1], &countVirtualNodes[2]);

    Configuration conf;
    //CartSpace space;
    initConfig(&conf, countVirtualNodes);
    //initCartSpace(&space, totalProcessesNumber, rank);

    const double calculationConstant = 2.0 / conf.H[0] + 2.0 / conf.H[1] + 2.0 / conf.H[2] + conf.a;
    
    //double phi[conf.N[0] * conf.N[1] * conf.N[2]] = {0};
    int status = 0;
    Function phi = initFunction((conf.N[0] * conf.N[1] * conf.N[2]) / totalProcessesNumber +
        ((conf.N[0] * conf.N[1] * conf.N[2]) % totalProcessesNumber > rank), (unsigned*)conf.N, &status);
    
    if (status) {
        return EXIT_FAILURE;
    }

    int boundaryK[totalProcessesNumber];
    boundaryK[0] = conf.N[2] / totalProcessesNumber + (conf.N[2] % totalProcessesNumber > 0);
    for (int i = 1; i < totalProcessesNumber; ++i) {
        boundaryK[i] = boundaryK[i - 1] + conf.N[2] / totalProcessesNumber + (conf.N[2] % totalProcessesNumber > 0);
    }
    
    Buffer buffer[2];
    if (rank != 0) {
        initBuffer(&buffer[0], conf.N[0] * conf.N[1]);
    } else if (rank != totalProcessesNumber - 1) {
        initBuffer(&buffer[1], conf.N[0] * conf.N[1]);
    }
    
    MPI_Request desc[4] = {MPI_REQUEST_NULL};
    const int subSize = conf.N[2] / totalProcessesNumber + (conf.N[2] % totalProcessesNumber > rank);
    int virtualNodeCoords[3] = {0};
    double values[2] = {0};
    double result = 0;
    double maxDiff = DBL_MAX, localMaxDiff = DBL_MAX;

    int z = 0;

    while (maxDiff > conf.epsilon) {
        localMaxDiff = conf.epsilon;

        if (rank != 0) {
            //printf("porn %d\n", rank);
            //printf("%d\n", buffer[0].size);
            MPI_Isend(buffer[0].send, buffer[0].size, MPI_DOUBLE, rank - 1, 101, MPI_COMM_WORLD, &desc[0]);
            MPI_Irecv(buffer[0].recieve, buffer[0].size, MPI_DOUBLE, rank - 1, MPI_ANY_TAG, MPI_COMM_WORLD, &desc[1]);
            //printf("%d\n", buffer[0].size);
        } else if (rank != totalProcessesNumber - 1) {
            //printf("porn %d\n", rank);
            MPI_Isend(buffer[1].send, buffer[1].size, MPI_DOUBLE, rank + 1, 103, MPI_COMM_WORLD, &desc[2]);
            MPI_Irecv(buffer[1].recieve, buffer[1].size, MPI_DOUBLE, rank + 1, MPI_ANY_TAG, MPI_COMM_WORLD, &desc[3]);
            //printf("%d\n", buffer[1].size);
        }

        // calculating independent values in inner nodes
        int n = rank > 0 ? boundaryK[rank - 1] : 0;
        int startIndex = (rank != 0);
        int endIndex = (rank != totalProcessesNumber - 1);
        for (int k = startIndex; k < subSize - endIndex; ++k) {
            virtualNodeCoords[2] = k;

            for (int j = 0; j < conf.N[1]; ++j) {
                virtualNodeCoords[1] = j;

                for (int i = 0; i < conf.N[0]; ++i) {
                    virtualNodeCoords[0] = i;
                    
                    for (int m = 0; m < 3; ++m) {
                        ++virtualNodeCoords[m];
                        values[0] = getPrevious(phi, virtualNodeCoords);
                        virtualNodeCoords[m] -= 2;
                        values[1] = getPrevious(phi, virtualNodeCoords);
                        ++virtualNodeCoords[m];

                        result += (values[0] + values[1]) / conf.H[m];
                    }
                    
                    // subtract right part
                    result -= 6 - conf.a * (i * i + j * j + n * n);
                    result /= calculationConstant;

                    if (fabs(result - getPrevious(phi, virtualNodeCoords)) > localMaxDiff) {
                        localMaxDiff = fabs(result - getPrevious(phi, virtualNodeCoords));
                    }

                    put(phi, virtualNodeCoords, result);

                    result = 0;
                }
            }

            ++n;
        }
        
        // waiting for recording buffer & calculating dependent on outer nodes inner nodes
        if (rank != 0) {
            //printf("%d\n", desc[0]);
            MPI_Wait(&desc[0], MPI_STATUS_IGNORE);
            //printf("porno\n");
            MPI_Wait(&desc[1], MPI_STATUS_IGNORE);
            //printf("porno\n");

            int n = rank > 0 ? boundaryK[rank - 1] : 0;
            virtualNodeCoords[2] = n;
            for (int j = 0; j < conf.N[1]; ++j) {
                virtualNodeCoords[1] = j;

                for (int i = 0; i < conf.N[0]; ++i) {
                    virtualNodeCoords[0] = i;
                    
                    for (int m = 0; m < 2; ++m) {
                        ++virtualNodeCoords[m];
                        values[0] = getPrevious(phi, virtualNodeCoords);
                        virtualNodeCoords[m] -= 2;
                        values[1] = getPrevious(phi, virtualNodeCoords);
                        ++virtualNodeCoords[m];

                        result += (values[0] + values[1]) / conf.H[m];
                    }

                    values[0] = buffer[0].recieve[j * conf.N[1] + i];
                    values[1] = getPrevious(phi, virtualNodeCoords);

                    result += (values[0] + values[1]) / conf.H[2];
                    
                    // subtract right part
                    result -= 6 - conf.a * (i * i + j * j + n * n);
                    result /= calculationConstant;

                    buffer[0].send[j * conf.N[1] + i] = result;

                    printf("Res_%d_%d: %lf %lf val: %lf %lf, node: (%d %d %d)\n", z, rank, result, getPrevious(phi, virtualNodeCoords), values[0], values[1], i, j, n);

                    if (fabs(result - getPrevious(phi, virtualNodeCoords)) > localMaxDiff) {
                        localMaxDiff = fabs(result - getPrevious(phi, virtualNodeCoords));
                    }

                    put(phi, virtualNodeCoords, result);

                    result = 0;
                }
            }
        } else if (rank != totalProcessesNumber - 1) {
            //printf("%d\n", desc[2]);
            MPI_Wait(&desc[2], MPI_STATUS_IGNORE);
            MPI_Wait(&desc[3], MPI_STATUS_IGNORE);
            //printf("porno\n");

            virtualNodeCoords[2] = 0;
            for (int j = 0; j < conf.N[1]; ++j) {
                virtualNodeCoords[1] = j;

                for (int i = 0; i < conf.N[0]; ++i) {
                    virtualNodeCoords[0] = i;
                    
                    for (int m = 0; m < 2; ++m) {
                        ++virtualNodeCoords[m];
                        values[0] = getPrevious(phi, virtualNodeCoords);
                        virtualNodeCoords[m] -= 2;
                        values[1] = getPrevious(phi, virtualNodeCoords);
                        ++virtualNodeCoords[m];

                        result += (values[0] + values[1]) / conf.H[m];
                    }

                    values[0] = getPrevious(phi, virtualNodeCoords);
                    values[1] = buffer[1].recieve[j * conf.N[1] + i];

                    result += (values[0] + values[1]) / conf.H[2];
                    
                    // subtract right part
                    result -= 6 - conf.a * (i * i + j * j);
                    result /= calculationConstant;

                    printf("Res_%d_%d: %lf %lf val: %lf %lf, node: (%d %d %d)\n", z, rank, result, getPrevious(phi, virtualNodeCoords), values[0], values[1], i, j, boundaryK[rank]);

                    buffer[1].send[j * conf.N[1] + i] = result;

                    if (fabs(result - getPrevious(phi, virtualNodeCoords)) > localMaxDiff) {
                        localMaxDiff = fabs(result - getPrevious(phi, virtualNodeCoords));
                    }

                    put(phi, virtualNodeCoords, result);

                    result = 0;
                }
            }
        }

        MPI_Allreduce(&localMaxDiff, &maxDiff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        if (!rank) {
            //printf("%.9f\n", maxDiff);
        }

        if (z == 2) {
            break;
        }
        ++z;

        swap(phi);
    }

    int recvcounts[totalProcessesNumber];
    int displs[totalProcessesNumber];
    for (int i = 0; i < totalProcessesNumber; ++i) {
        recvcounts[i] = (conf.N[0] * conf.N[1] * conf.N[2]) / totalProcessesNumber +
            ((conf.N[0] * conf.N[1] * conf.N[2]) % totalProcessesNumber > i);
        displs[i] = i > 0 ? displs[i - 1] + recvcounts[i] : 0;
    }

    double funcValues[conf.N[0] * conf.N[1] * conf.N[2]];
    MPI_Gatherv(getCurrentBuffer(phi), phi->size, MPI_DOUBLE, funcValues, recvcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (!rank) {
        isCorrectCalculation(funcValues, conf);
    }

    MPI_Finalize();
    return EXIT_SUCCESS;
}