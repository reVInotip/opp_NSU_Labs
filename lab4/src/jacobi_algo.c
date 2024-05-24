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

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int rank, totalProcessesNumber;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &totalProcessesNumber);

    int boundaryKForPrevRank = 0;
    for (int i = 0; i < rank; ++i)
    {
        boundaryKForPrevRank = boundaryKForPrevRank + N_z / totalProcessesNumber + (N_z % totalProcessesNumber > i);
    }

    int boundaryKForCurrRank = boundaryKForPrevRank + N_z / totalProcessesNumber + (N_z % totalProcessesNumber > rank);

    double h[3] = {0, 0, 0};
    double squareH[3] = {0, 0, 0};
    h[0] = (double)D_x / (N_x - 1);
    h[1] = (double)D_y / (N_y - 1);
    h[2] = (double)D_z / (N_z - 1);
    for (int i = 0; i < 3; ++i)
    {
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
    NodeFunction *phi = initNodeFunction((N_z / totalProcessesNumber + (N_z % totalProcessesNumber > rank)) * N_x * N_y,
        totalProcessesNumber, rank, &status);
    if (status)
    {
        return EXIT_FAILURE;
    }

    Buffer buffer[2];
    if (rank != 0)
    {
        initBuffer(&buffer[0], N_x * N_y);
    }
    if (rank != totalProcessesNumber - 1)
    {
        initBuffer(&buffer[1], N_x * N_y);
    }

    double start = MPI_Wtime();

    while (maxDiff > EPSILON)
    {
        localMaxDiff = EPSILON;

        if (rank != 0)
        {
            MPI_Isend(buffer[0].send, buffer[0].size, MPI_DOUBLE, rank - 1, 101, MPI_COMM_WORLD, &desc[0]);
            MPI_Irecv(buffer[0].recieve, buffer[0].size, MPI_DOUBLE, rank - 1, 103, MPI_COMM_WORLD, &desc[1]);
        }

        if (rank != totalProcessesNumber - 1)
        {
            MPI_Isend(buffer[1].send, buffer[1].size, MPI_DOUBLE, rank + 1, 103, MPI_COMM_WORLD, &desc[2]);
            MPI_Irecv(buffer[1].recieve, buffer[1].size, MPI_DOUBLE, rank + 1, 101, MPI_COMM_WORLD, &desc[3]);
        }

        // calculating independent values in inner nodes
        int realZCoord = startZCoord;
        //printf("%d , %d, %d\n", leftBorderByZ, rightBorderByZ, phi->lengthByZCoord);
        for (int k = leftBorderByZ; k < rightBorderByZ; ++k)
        {
            for (int j = 0; j < N_y; ++j)
            {
                for (int i = 0; i < N_x; ++i)
                {
                    int m = ((k == rightBorderByZ - 1) && (rank == totalProcessesNumber - 1)) ? realZCoord + 1 : k + 1;
                    int m1 = ((j + 1) * h[1] > D_y) ? realZCoord : k;
                    int m11 = (j - 1 < 0) ? realZCoord : k;
                    int m2 = ((i + 1) * h[0] > D_x) ? realZCoord : k;
                    int m22 = (i - 1 < 0) ? realZCoord : k;
                    result = (getPreviousValue(phi, i, j, k - 1, h) + getPreviousValue(phi, i, j, m, h)) / squareH[2] +
                              (getPreviousValue(phi, i, j - 1, m11, h) + getPreviousValue(phi, i, j + 1, m1, h)) / squareH[1] +
                              (getPreviousValue(phi, i - 1, j, m22, h) + getPreviousValue(phi, i + 1, j, m2, h)) / squareH[0];

                    // subtract right part
                    result -= 6.0 - A * getFunctionValue(i, j, realZCoord, h);
                    result /= calculationConstant;

                    if (fabs(result - getPreviousValue(phi, i, j, k, h)) > localMaxDiff)
                    {
                        localMaxDiff = fabs(result - getPreviousValue(phi, i, j, k, h));
                    }

                    if (i * h[0] - 1 == 1 && j * h[1] - 1 == 1) {
                        //printf("LMD_i: %lf (%lf %lf %lf), result: %lf, prev: %lf, k: (%d, %d, %d)\n",
                        //  localMaxDiff, i * h[0] - 1, j * h[1] - 1, realZCoord * h[2] - 1, result, getPreviousValue(phi, i, j - 1, m1, h), i, j-1, m1);
                    }

                    put(phi, i, j, k, result);
                }
            }

            ++realZCoord;
        }

        intermediateLocalMaxDiff = localMaxDiff;

        // waiting for recording buffer & calculating dependent on outer nodes inner nodes
        if (rank != 0)
        {
            MPI_Wait(&desc[0], MPI_STATUS_IGNORE);
            MPI_Wait(&desc[1], MPI_STATUS_IGNORE);

            realZCoord = boundaryKForPrevRank;
            for (int j = 0; j < N_y; ++j)
            {
                for (int i = 0; i < N_x; ++i)
                {
                    int m = (phi->lengthByZCoord == 1 ? realZCoord : 0) + 1;
                    int m1 = ((j + 1) * h[1] > D_y) ? realZCoord : 0;
                    int m11 = (j - 1 < 0) ? realZCoord : 0;
                    int m2 = ((i + 1) * h[0] > D_x) ? realZCoord : 0;
                    int m22 = (i - 1 < 0) ? realZCoord : 0;
                    result = (buffer[0].recieve[j * N_x + i] + getPreviousValue(phi, i, j, m, h)) / squareH[2] +
                              (getPreviousValue(phi, i, j - 1, m11, h) + getPreviousValue(phi, i, j + 1, m1, h)) / squareH[1] +
                              (getPreviousValue(phi, i - 1, j, m22, h) + getPreviousValue(phi, i + 1, j, m2, h)) / squareH[0];

                    // subtract right part
                    result -= 6.0 - A * getFunctionValue(i, j, realZCoord, h);
                    result /= calculationConstant;

                    buffer[0].send[j * N_x + i] = result;

                    if (fabs(result - getPreviousValue(phi, i, j, 0, h)) > localMaxDiff)
                    {
                        localMaxDiff = fabs(result - getPreviousValue(phi, i, j, 0, h));
                    }

                    if (i * h[0] - 1 == 1 && j * h[1] - 1 == 1) {
                        //printf("LMD_o1: %lf (%lf %lf %lf), result: %lf, prev: %lf, k: (%d, %d, %d)\n",
                        //    localMaxDiff, i * h[0] - 1, j * h[1] - 1, realZCoord * h[2] - 1, result, getPreviousValue(phi, i, j - 1, m11, h), i, j-1, m11);
                    }

                    put(phi, i, j, 0, result);
                }
            }
        }

        if (rank != totalProcessesNumber - 1)
        {
            MPI_Wait(&desc[2], MPI_STATUS_IGNORE);
            MPI_Wait(&desc[3], MPI_STATUS_IGNORE);

            localMaxDiff = intermediateLocalMaxDiff;

            realZCoord = boundaryKForCurrRank - 1;
            int z = phi->lengthByZCoord - 1;
            for (int j = 0; j < N_y; ++j)
            {
                for (int i = 0; i < N_x; ++i)
                {
                    int m1 = ((j + 1) * h[1] > D_y) ? realZCoord : z;
                    int m11 = (j - 1 < 0) ? realZCoord : z;
                    int m2 = ((i + 1) * h[0] > D_x) ? realZCoord : z;
                    int m22 = (i - 1 < 0) ? realZCoord : z;
                    result = (buffer[1].recieve[j * N_x + i] + getPreviousValue(phi, i, j, z - 1, h)) / squareH[2] +
                              (getPreviousValue(phi, i, j - 1, m11, h) + getPreviousValue(phi, i, j + 1, m1, h)) / squareH[1] +
                              (getPreviousValue(phi, i - 1, j, m22, h) + getPreviousValue(phi, i + 1, j, m2, h)) / squareH[0];

                    // subtract right part
                    result -= 6.0 - A * getFunctionValue(i, j, realZCoord, h);
                    result /= calculationConstant;

                    buffer[1].send[j * N_x + i] = result;

                    if (fabs(result - getPreviousValue(phi, i, j, z, h)) > localMaxDiff)
                    {
                        localMaxDiff = fabs(result - getPreviousValue(phi, i, j, z, h));
                    }

                    if (i * h[0] - 1 == 1 && j * h[1] - 1 == 1) {
                        //printf("LMD_o2: %lf (%lf %lf %lf), result: %lf, prev: %lf, k: %d\n",
                        //    localMaxDiff, i * h[0] - 1, j * h[1] - 1, realZCoord * h[2] - 1, result, getPreviousValue(phi, i, j - 1, z, h), z);
                    }

                    put(phi, i, j, z, result);
                }
            }
        }

        MPI_Allreduce(&localMaxDiff, &maxDiff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        swap(phi);
    }

    if (!rank)
    {
        printf("Time taken: %lf sec.\n", MPI_Wtime() - start);
    }

    isCorrectCalculation(h, phi, rank, boundaryKForPrevRank);

    if (rank != 0)
    {
        destroyBuffer(&buffer[0]);
    }
    if (rank != totalProcessesNumber - 1)
    {
        destroyBuffer(&buffer[1]);
    }

    destroyNodeFunction(phi);

    MPI_Finalize();
    return EXIT_SUCCESS;
}