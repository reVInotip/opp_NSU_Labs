#include "matrix/matrix.h"
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

enum ERR_CODE_MAIN {
    WRONG_MARTIX_SIZE_ERR = 20
};

int main(int argc, char **argv) {
    int errorCode = 0;
    unsigned matrixSizes[3] = {0, 0, 0};
    int dims[2]={0, 0};
    int rank;
    
    MPI_Init(argc, argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // create new communicator
    MPI_Comm comm2D, comm1D[2];

    int size;
    int periods[2]={0,0};
    int coords[2];
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Dims_create(size, 2, dims);
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &comm2D);

    // create subcommunicatros, 0 - x, 1 - y
    int remain_dims[2] = {0, 0};
    for (int i  = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            remain_dims[j] = (i == j);
        }
        MPI_Cart_sub(comm2D, remain_dims, &comm1D[i]);
    }

    int sendcountsForA[size], sendcountsForB[size];
    int displsForA[size], displsForB[size];

    // create full matrix on 0 process
    if (!rank) {
        Matrix A = createMatrixFromFile("./A", &errorCode);
        Matrix B = createMatrixFromFile("./B", &errorCode);
        Matrix C = createMatrixFromFile("./C", &errorCode);
        if (errorCode) {
            destroyMatrix(A);
            destroyMatrix(B);
            destroyMatrix(C);
            return errorCode;
        }

        if ((A->width != B->height) || (B->width != C->width) || (A->height != C->height)) {
            fprintf(stderr, "Can not multiply this matrices");
            errorCode = WRONG_MARTIX_SIZE_ERR;
            return errorCode;
        }

        matrixSizes[0] = A->width;
        matrixSizes[1] = B->width;
        matrixSizes[2] = A->height;

        MPI_Bcast(matrixSizes, 3, MPI_INT, 0, MPI_COMM_WORLD);

        // create displs and counts for matrix A
        int i = 0;
        for (int process = 0; process < size; ++process) {
            if (process % dims[1] == 0 && process != 0) {
                ++i;
            }
            sendcountsForA[process] = A->height / dims[0] + (A->height % dims[0] > i);
        }
    }

    // create submatrices
    unsigned subAHeight = matrixSizes[2] / dims[0] + (matrixSizes[2] % dims[0] > dims[0]);
    unsigned subBWidth = matrixSizes[1] / dims[1] + (matrixSizes[1] / dims[1] > dims[1]);
    Matrix subA = createMatrix(matrixSizes[0], subAHeight, &errorCode);
    Matrix subB = createMatrix(subBWidth, matrixSizes[0], &errorCode);
    Matrix subC = createMatrix(subBWidth, subAHeight, &errorCode);
    if (errorCode) {
        destroyMatrix(subA);
        destroyMatrix(subB);
        destroyMatrix(subC);
        return errorCode;
    }

    MPI_Datatype vectorTypeForB;

    if (!rank) {
        MPI_Datatype types[2];
        // create MPI type for matrix B
        int vectorSize;
        MPI_Type_vector(1, matrixSizes[0], matrixSizes[1], MPI_DOUBLE, &types[0]);
        MPI_Type_extent(types[0], &vectorSize);
        MPI_Type_create_resized(types[0], 0, vectorSize, &vectorTypeForB);
        MPI_Type_commit(&vectorTypeForB);
    }

    return EXIT_SUCCESS;
}