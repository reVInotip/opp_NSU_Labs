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

void isCorrectCalculation(const Function *phi, const double precision, const int rank, int startZCoord) {
    bool isCorrect = true;
    int coords[3] = {0};
    double delta = 0, maxDelta = 0;

    for (int k = 0; k < phi->lengthByCoord[2]; ++k, ++startZCoord) {
        coords[2] = k;

        for (int j = 0; j < phi->lengthByCoord[1]; ++j) {
            coords[1] = j;

            for (int i = 0; i < phi->lengthByCoord[0]; ++i) {
                coords[0] = i;
                delta = fabs((i * i + j * j + startZCoord * startZCoord) - getCurrent(phi, coords));

                if (delta >= precision && delta > 240) {
                    printf("Incorrect value in (%d, %d, %d)_%d node!\nCorrect: %d, result: %lf\nPrecision: %lf\n",
                        i, j, startZCoord, rank, i * i + j * j + startZCoord * startZCoord, getCurrent(phi, coords), delta);
                    isCorrect = false;
                }

                if (delta > maxDelta) {
                    maxDelta = delta;
                };
            }
        }
    }
    
    MPI_Allreduce(MPI_IN_PLACE, &isCorrect, 1, MPI_C_BOOL, MPI_LAND, MPI_COMM_WORLD);

    if (isCorrect) {
        MPI_Allreduce(MPI_IN_PLACE, &maxDelta, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        if (!rank) {
            printf("CORRECT! Max precision: %lf\n", maxDelta);
        }
    }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, totalProcessesNumber;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &totalProcessesNumber);

    unsigned countVirtualNodes[3] = {0, 0, 0};
    parseAgrs(argc, argv, &countVirtualNodes[0], &countVirtualNodes[1], &countVirtualNodes[2]);

    Configuration conf;
    initConfig(&conf, countVirtualNodes);

    const double calculationConstant = 2.0 / conf.H[0] + 2.0 / conf.H[1] + 2.0 / conf.H[2] + conf.a;
    
    int status = 0;
    unsigned lengthByCoord[3] = {conf.N[0], conf.N[1], conf.N[2] / totalProcessesNumber + (conf.N[2] % totalProcessesNumber > rank)};
    Function *phi = initFunction((conf.N[2] / totalProcessesNumber + (conf.N[2] % totalProcessesNumber > rank)) *
        conf.N[0] * conf.N[1], (unsigned*)lengthByCoord, &status);
    
    if (status) {
        return EXIT_FAILURE;
    }

    int boundaryK[totalProcessesNumber];
    boundaryK[0] = conf.N[2] / totalProcessesNumber + (conf.N[2] % totalProcessesNumber > 0);
    for (int i = 1; i < totalProcessesNumber; ++i) {
        boundaryK[i] = boundaryK[i - 1] + conf.N[2] / totalProcessesNumber + (conf.N[2] % totalProcessesNumber > i);
    }
    
    Buffer buffer[2];
    if (rank != 0) {
        initBuffer(&buffer[0], conf.N[0] * conf.N[1]);
    }
    
    if (rank != totalProcessesNumber - 1) {
        initBuffer(&buffer[1], conf.N[0] * conf.N[1]);
    }
    
    MPI_Request desc[4] = {MPI_REQUEST_NULL};
    const int subSize = conf.N[2] / totalProcessesNumber + (conf.N[2] % totalProcessesNumber > rank);
    int virtualNodeCoords[3] = {0};
    double values[2] = {0};
    double result = 0;
    double maxDiff = DBL_MAX, localMaxDiff = DBL_MAX, intermediateLocalMaxDiff = 0;
    const int startZCoord = rank > 0 ? boundaryK[rank - 1] + 1 : 0;

    while (maxDiff > conf.epsilon) {
        localMaxDiff = conf.epsilon;

        if (rank != 0) {
            MPI_Isend(buffer[0].send, buffer[0].size, MPI_DOUBLE, rank - 1, 101, MPI_COMM_WORLD, &desc[0]);
            MPI_Irecv(buffer[0].recieve, buffer[0].size, MPI_DOUBLE, rank - 1, MPI_ANY_TAG, MPI_COMM_WORLD, &desc[1]);
        }
        
        if (rank != totalProcessesNumber - 1) {
            MPI_Isend(buffer[1].send, buffer[1].size, MPI_DOUBLE, rank + 1, 103, MPI_COMM_WORLD, &desc[2]);
            MPI_Irecv(buffer[1].recieve, buffer[1].size, MPI_DOUBLE, rank + 1, MPI_ANY_TAG, MPI_COMM_WORLD, &desc[3]);
        }

        // calculating independent values in inner nodes
        int realZCoord = startZCoord;
        int startIndex = (rank != 0);
        int endIndex = (rank != totalProcessesNumber - 1);
        for (int k = startIndex; k < subSize - endIndex; ++k) {
            virtualNodeCoords[2] = k;

            for (int j = 0; j < conf.N[1]; ++j) {
                virtualNodeCoords[1] = j;

                for (int i = 0; i < conf.N[0]; ++i) {
                    virtualNodeCoords[0] = i;
                    
                    for (int coordNumber = 0; coordNumber < 3; ++coordNumber) {
                        ++virtualNodeCoords[coordNumber];
                        values[0] = getPrevious(phi, virtualNodeCoords);
                        virtualNodeCoords[coordNumber] -= 2;
                        values[1] = getPrevious(phi, virtualNodeCoords);
                        ++virtualNodeCoords[coordNumber];

                        result += (values[0] + values[1]) / conf.H[coordNumber];
                    }
                    
                    // subtract right part
                    result -= 6 - conf.a * (i * i + j * j + realZCoord * realZCoord);
                    result /= calculationConstant;

                    if (fabs(result - getPrevious(phi, virtualNodeCoords)) > localMaxDiff) {
                        localMaxDiff = fabs(result - getPrevious(phi, virtualNodeCoords));
                    }

                    put(phi, virtualNodeCoords, result);

                    result = 0;
                }
            }

            ++realZCoord;
        }
        
        intermediateLocalMaxDiff = localMaxDiff;

        // waiting for recording buffer & calculating dependent on outer nodes inner nodes
        if (rank != 0) {
            MPI_Wait(&desc[0], MPI_STATUS_IGNORE);
            MPI_Wait(&desc[1], MPI_STATUS_IGNORE);

            realZCoord = boundaryK[rank - 1];
            virtualNodeCoords[2] = 0;
            for (int j = 0; j < conf.N[1]; ++j) {
                virtualNodeCoords[1] = j;

                for (int i = 0; i < conf.N[0]; ++i) {
                    virtualNodeCoords[0] = i;
                    
                    for (int coordNumber = 0; coordNumber < 2; ++coordNumber) {
                        ++virtualNodeCoords[coordNumber];
                        values[0] = getPrevious(phi, virtualNodeCoords);
                        virtualNodeCoords[coordNumber] -= 2;
                        values[1] = getPrevious(phi, virtualNodeCoords);
                        ++virtualNodeCoords[coordNumber];

                        result += (values[0] + values[1]) / conf.H[coordNumber];
                    }

                    values[0] = buffer[0].recieve[j * conf.N[0] + i];
                    ++virtualNodeCoords[2];
                    values[1] = getPrevious(phi, virtualNodeCoords);
                    --virtualNodeCoords[2];

                    result += (values[0] + values[1]) / conf.H[2];
                    
                    // subtract right part
                    result -= 6 - conf.a * (i * i + j * j + realZCoord * realZCoord);
                    result /= calculationConstant;

                    buffer[0].send[j * conf.N[0] + i] = result;

                    if (fabs(result - getPrevious(phi, virtualNodeCoords)) > localMaxDiff) {
                        localMaxDiff = fabs(result - getPrevious(phi, virtualNodeCoords));
                    }

                    put(phi, virtualNodeCoords, result);

                    result = 0;
                }
            }
        }

        if (rank != totalProcessesNumber - 1) {
            MPI_Wait(&desc[2], MPI_STATUS_IGNORE);
            MPI_Wait(&desc[3], MPI_STATUS_IGNORE);

            localMaxDiff = intermediateLocalMaxDiff;

            realZCoord = boundaryK[rank] - 1;
            virtualNodeCoords[2] = phi->lengthByCoord[2] - 1;
            for (int j = 0; j < conf.N[1]; ++j) {
                virtualNodeCoords[1] = j;

                for (int i = 0; i < conf.N[0]; ++i) {
                    virtualNodeCoords[0] = i;
                    
                    for (int coordNumber = 0; coordNumber < 2; ++coordNumber) {
                        ++virtualNodeCoords[coordNumber];
                        values[0] = getPrevious(phi, virtualNodeCoords);
                        virtualNodeCoords[coordNumber] -= 2;
                        values[1] = getPrevious(phi, virtualNodeCoords);
                        ++virtualNodeCoords[coordNumber];

                        result += (values[0] + values[1]) / conf.H[coordNumber];
                    }

                    --virtualNodeCoords[2];
                    values[0] = getPrevious(phi, virtualNodeCoords);
                    ++virtualNodeCoords[2];
                    values[1] = buffer[1].recieve[j * conf.N[0] + i];

                    result += (values[0] + values[1]) / conf.H[2];
                    
                    // subtract right part
                    result -= 6 - conf.a * (i * i + j * j + realZCoord * realZCoord);
                    result /= calculationConstant;

                    buffer[1].send[j * conf.N[0] + i] = result;

                    if (fabs(result - getPrevious(phi, virtualNodeCoords)) > localMaxDiff) {
                        localMaxDiff = fabs(result - getPrevious(phi, virtualNodeCoords));
                    }

                    put(phi, virtualNodeCoords, result);

                    result = 0;
                }
            }
        }

        MPI_Allreduce(&localMaxDiff, &maxDiff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        swap(phi);
    }

    isCorrectCalculation(phi, conf.precision, rank, rank > 0 ? boundaryK[rank - 1] : 0);

    if (rank != 0) {
        destroyBuffer(&buffer[0]);
    }
    
    if (rank != totalProcessesNumber - 1) {
        destroyBuffer(&buffer[1]);
    }

    destroyFunction(phi);

    MPI_Finalize();
    return EXIT_SUCCESS;
}