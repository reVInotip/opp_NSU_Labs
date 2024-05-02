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

#include "../include/node_function.h"
#include "../include/function.h"
#include "../include/configuration.h"
#include "../include/utils.h"
#include "../include/buffer.h"
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <stdbool.h>
#include <float.h>
#include <stdio.h>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, totalProcessesNumber;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &totalProcessesNumber);

    int boundaryKForPrevRank = 0;
    for (int i = 0; i < rank; ++i) {
        boundaryKForPrevRank = boundaryKForPrevRank + N_z / totalProcessesNumber + (N_z % totalProcessesNumber > i);
    }

    int boundaryKForCurrRank = boundaryKForPrevRank + N_z / totalProcessesNumber + (N_z % totalProcessesNumber > rank);

    double h[3] = {0, 0, 0};
    double squareH[3] = {0, 0, 0};
    h[0] = (double)D_x / (N_x - 1);
    h[1] = (double)D_y / (N_y - 1);
    h[2] = (double)D_z / (N_z - 1);
    for (int i = 0; i < 3; ++i) {
        squareH[i] = h[i] * h[i];
    }

    double result = 0;
    double maxDiff = DBL_MAX, localMaxDiff = DBL_MAX, intermediateLocalMaxDiff = 0;
    MPI_Request desc[4] = {MPI_REQUEST_NULL};

    const int startZCoord = (boundaryKForPrevRank + 1) * (rank != 0);
    const int rightBorderByZ = (N_z / totalProcessesNumber + (N_z % totalProcessesNumber > rank)) - (rank != totalProcessesNumber - 1);
    const int leftBorderByZ = (rank != 0);

    const double calculationConstant = 2.0 / squareH[0] + 2.0 / squareH[1] + 2.0 / squareH[2] + A;
    
    int status = 0;
    NodeFunction *phi = initNodeFunction((N_z / totalProcessesNumber + (N_z % totalProcessesNumber > rank)) *
        N_x * N_y, totalProcessesNumber, rank, &status);
    if (status) {
        return EXIT_FAILURE;
    }
    
    Buffer buffer[2];
    if (rank != 0) {
        initBuffer(&buffer[0], N_x * N_y);
    }
    if (rank != totalProcessesNumber - 1) {
        initBuffer(&buffer[1], N_x * N_y);
    }

    double start = MPI_Wtime();

    while (maxDiff > EPSILON) {
        localMaxDiff = EPSILON;

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
        for (int k = leftBorderByZ; k < rightBorderByZ; ++k) {
            for (int j = 0; j < N_y; ++j) {
                for (int i = 0; i < N_x; ++i) {

                    result += (getPreviousValue(phi, i, j, k - 1) + getPreviousValue(phi, i, j, k + 1)) / squareH[2] +
                        (getPreviousValue(phi, i, j - 1, k) + getPreviousValue(phi, i, j + 1, k)) / squareH[1] +
                        (getPreviousValue(phi, i - 1, j, k) + getPreviousValue(phi, i + 1, j, k)) / squareH[0];

                    // subtract right part
                    result -= 6.0 - A * getFunctionValue(i, j, realZCoord, h);
                    result /= calculationConstant;

                    if (fabs(result - getPreviousValue(phi, i, j, k)) > localMaxDiff) {
                        localMaxDiff = fabs(result - getPreviousValue(phi, i, j, k));
                    }

                    put(phi, i, j, k, result);

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

            realZCoord = boundaryKForPrevRank;
            for (int j = 0; j < N_y; ++j) {
                for (int i = 0; i < N_x; ++i) {

                    result += (buffer[0].recieve[j * N_x + i] + getPreviousValue(phi, i, j, 1)) / squareH[2] +
                        (getPreviousValue(phi, i, j - 1, 0) + getPreviousValue(phi, i, j + 1, 0)) / squareH[1] +
                        (getPreviousValue(phi, i - 1, j, 0) + getPreviousValue(phi, i + 1, j, 0)) / squareH[0];
                    
                    // subtract right part
                    result -= 6.0 - A * getFunctionValue(i, j, realZCoord, h);
                    result /= calculationConstant;

                    buffer[0].send[j * N_x + i] = result;

                    if (fabs(result - getPreviousValue(phi, i, j, 0)) > localMaxDiff) {
                        localMaxDiff = fabs(result - getPreviousValue(phi, i, j, 0));
                    }

                    put(phi, i, j, 0, result);

                    result = 0;
                }
            }
        }

        if (rank != totalProcessesNumber - 1) {
            MPI_Wait(&desc[2], MPI_STATUS_IGNORE);
            MPI_Wait(&desc[3], MPI_STATUS_IGNORE);

            localMaxDiff = intermediateLocalMaxDiff;

            realZCoord = boundaryKForCurrRank - 1;
            int z = phi->lengthByZCoord - 1;
            for (int j = 0; j < N_y; ++j) {
                for (int i = 0; i < N_x; ++i) {

                    result += (buffer[1].recieve[j * N_x + i] + getPreviousValue(phi, i, j, z - 1)) / squareH[2] +
                        (getPreviousValue(phi, i, j - 1, z) + getPreviousValue(phi, i, j + 1, z)) / squareH[0] +
                        (getPreviousValue(phi, i - 1, j, z) + getPreviousValue(phi, i + 1, j, z)) / squareH[1];
                    
                    // subtract right part
                    result -= 6.0 - A * getFunctionValue(i, j, realZCoord, h);
                    result /= calculationConstant;

                    buffer[1].send[j * N_x + i] = result;

                    if (fabs(result - getPreviousValue(phi, i, j, z)) > localMaxDiff) {
                        localMaxDiff = fabs(result - getPreviousValue(phi, i, j, z));
                    }

                    put(phi, i, j, z, result);

                    result = 0;
                }
            }
        }

        MPI_Allreduce(&localMaxDiff, &maxDiff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        swap(phi);
    }

    if (!rank) {
        printf("Time taken: %lf sec.\n", MPI_Wtime() - start);
    }

    isCorrectCalculation(h, phi, rank, boundaryKForPrevRank);

    if (rank != 0) {
        destroyBuffer(&buffer[0]);

    }
    if (rank != totalProcessesNumber - 1) {
        destroyBuffer(&buffer[1]);
    }

    destroyNodeFunction(phi);

    MPI_Finalize();
    return EXIT_SUCCESS;
}