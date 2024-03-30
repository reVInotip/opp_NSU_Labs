#include "matrix/matrix.h"
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

enum ERR_CODE_MAIN {
    WRONG_MARTIX_SIZE_ERR = 20,
    CANT_DIVIDE_MATRIX_ON_SUBMATRIX_ERR = 21
};

int main(int argc, char **argv) {
    int errorCode = 0;
    unsigned matrixSizes[3] = {0, 0, 0};
    int dims[2]={0, 0};
    int rank;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // create new communicator
    MPI_Comm comm2D, commRows, commColumns;

    int size;
    int periods[2] = {0,0};
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Dims_create(size, 2, dims);
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &comm2D);
    
    MPI_Comm_rank(comm2D, &rank);

    // create subcommunicatros, 0 - x, 1 - y
    int remainDimsForRows[2] = {1, 0};
    int remainDimsForColumns[2] = {0, 1};
    MPI_Cart_sub(comm2D, remainDimsForRows, &commRows);
    MPI_Cart_sub(comm2D, remainDimsForColumns, &commColumns);

    int defaultRank, comm2DRank, commRowsRank, commColumnsRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &defaultRank);
    MPI_Comm_rank(comm2D, &comm2DRank);
    MPI_Comm_rank(commRows, &commRowsRank);
    MPI_Comm_rank(commColumns, &commColumnsRank);

    int coords[2] = {0, 0};
    MPI_Cart_coords(comm2D, rank, 2, coords);

    printf("Ranks: default - %d, comm2D - %d, commRows - %d, commColumns - %d, x - %d, y - %d\n",
        defaultRank, comm2DRank, commRowsRank, commColumnsRank, coords[0], coords[1]);

    Matrix A;
    Matrix B;
    Matrix C;

    // create full matrix on 0 process
    if (!rank) {
        A = createMatrixFromFile("./A", &errorCode);
        B = createMatrixFromFile("./B", &errorCode);
        C = createMatrix(B->width, A->height, &errorCode);
        if (errorCode) {
            destroyMatrix(A);
            destroyMatrix(B);
            destroyMatrix(C);
            return errorCode;
        }

        if ((A->width != B->height) || (B->width != C->width) || (A->height != C->height)) {
            fprintf(stderr, "Can not multiply this matrices");
            errorCode = WRONG_MARTIX_SIZE_ERR;
            destroyMatrix(A);
            destroyMatrix(B);
            destroyMatrix(C);
            return errorCode;
        } else if (A->height % dims[0] != 0 || B->width % dims[1] != 0) {
            fprintf(stderr, "Matrix sizes should be multiplies of count processes!");
            errorCode = CANT_DIVIDE_MATRIX_ON_SUBMATRIX_ERR;
            destroyMatrix(A);
            destroyMatrix(B);
            destroyMatrix(C);
            return errorCode;
        }

        matrixSizes[0] = A->width;
        matrixSizes[1] = B->width;
        matrixSizes[2] = A->height;
    }

    MPI_Bcast(matrixSizes, 3, MPI_INT, 0, MPI_COMM_WORLD);

    // create submatrices
    unsigned subAHeight = matrixSizes[2] / dims[0];
    unsigned subBWidth = matrixSizes[1] / dims[1];
    Matrix subA = createMatrix(matrixSizes[0], subAHeight, &errorCode);
    Matrix subB = createMatrix(subBWidth, matrixSizes[0], &errorCode);
    Matrix subC = createMatrix(subBWidth, subAHeight, &errorCode);
    if (errorCode) {
        destroyMatrix(subA);
        destroyMatrix(subB);
        destroyMatrix(subC);
        return errorCode;
    }

    MPI_Datatype vectorTypeForB, vectorTypeForC;

    int sendcountsForC[size], sendcountsForB[dims[1]];
    int displsForC[size], displsForB[dims[1]];

    /*if (!rank) {
        MPI_Datatype types[2];
        // create MPI type for matrix B
        MPI_Aint vectorSize;
        MPI_Type_vector(matrixSizes[0], subBWidth, matrixSizes[1], MPI_DOUBLE, &types[0]);
        MPI_Type_extent(types[0], &vectorSize);
        MPI_Type_create_resized(types[0], 0, vectorSize, &vectorTypeForB);
        MPI_Type_commit(&vectorTypeForB);

        // calc displs and send counts for B
        for (int i = 0; i < dims[1]; ++i) {
            sendcountsForB[i] = 1;
            displsForB[i] = i;
        }

        // create MPI type for matrix C
        MPI_Type_vector(matrixSizes[0] / dims[0], matrixSizes[1] / dims[1], matrixSizes[1], MPI_DOUBLE, &types[1]);
        MPI_Type_extent(types[1], &vectorSize);
        MPI_Type_create_resized(types[1], 0, vectorSize, &vectorTypeForC);
        MPI_Type_commit(&vectorTypeForC);

        // calc displs and send counts for C
        for (int i = 0; i < dims[0]; ++i) {
            for (int j = 0; j < dims[1]; ++j) {
                sendcountsForC[i * dims[1] + j] = 1;
                displsForC[i * dims[1] + j] = i * subBWidth + j;
            }
        }
    }*/

    // send B by Y coordinate
    //MPI_Scatterv(B, sendcountsForB, displsForB, vectorTypeForB, subB, subB->size, MPI_DOUBLE, 0, comm1D[1]);

    //MPI_Comm_rank(commRows, &rank);

    // send A by X coordinate
    if (coords[0] == 0) {
        MPI_Scatter(A, subA->size , MPI_DOUBLE, subA, subA->size, MPI_DOUBLE, 0, commColumns);
    }

    // MPI_Comm_rank(commColumns, &rank);
    MPI_Bcast(subA, subA->size, MPI_DOUBLE, 0, commRows);

    printf("Rank: %d\n", rank);
    printMatrix(subA);

    for (int k = 0; k < subC->size; ++k) {
        for (int i = 0; i < subA->width; ++i) {
            for (int j = 0; j < subB->height; ++j) {
                subC->data[k] = subA->data[i * subA->width + j] + subB->data[j * subA->width + i];
            }
        }
    }

    MPI_Finalize();

    return EXIT_SUCCESS;
}