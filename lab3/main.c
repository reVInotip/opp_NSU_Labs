#include "matrix/matrix.h"
#include "parser/parser.h"
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <mpi.h>

#define A_SEED 1239032
#define B_SEED 5982847
#define SIZE_CAP 200

enum ERR_CODE_MAIN {
    WRONG_MARTIX_SIZE_ERR = 20,
    INCORRECT_CALCULATION = 21,
    INCORRECT_PROCESS_COUNT = 22,
    OK = 0
};

typedef struct {
    // A height
    unsigned aHeight;
    // B width
    unsigned bWidth;
    // length of other sides of both matrices
    unsigned side;
} MatrixSizes;


int main(int argc, char **argv) {
    int errorCode = 0;
    int dims[2]={0, 0};
    int totalProcessesNumber;
    int rank;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &totalProcessesNumber);

    // create new communicators
    MPI_Comm comm2D, commRows, commColumns;

    
    int periods[2] = {0,0};
    MPI_Dims_create(totalProcessesNumber, 2, dims);
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &comm2D);

    // create subcommunicatros, 0 - x, 1 - y
    int remainDimsForRows[2] = {1, 0};
    int remainDimsForColumns[2] = {0, 1};
    MPI_Cart_sub(comm2D, remainDimsForRows, &commRows);
    MPI_Cart_sub(comm2D, remainDimsForColumns, &commColumns);

    int coords[2] = {0, 0};
    MPI_Cart_coords(comm2D, rank, 2, coords);

    // initialize matrices
    Matrix A = createMatrix(0, 0, &errorCode);
    Matrix B = createMatrix(0, 0, &errorCode);
    Matrix C = createMatrix(0, 0, &errorCode);
    if (errorCode) {
        return errorCode;
    }

    MatrixSizes matrixSizes = {0, 0, 0};

    // create type for struct matrixSizes
    MPI_Datatype MATRIX_SIZES;
    
    int blocklengths[3] = {1, 1, 1};
    MPI_Aint displacements[3] = {offsetof(MatrixSizes, aHeight), offsetof(MatrixSizes, bWidth), offsetof(MatrixSizes, side)};
    MPI_Datatype types[3] = {MPI_UNSIGNED, MPI_UNSIGNED, MPI_UNSIGNED};
    MPI_Type_create_struct(3, blocklengths, displacements, types, &MATRIX_SIZES);
    MPI_Type_commit(&MATRIX_SIZES);

    int heightsForSubA[totalProcessesNumber];
    int widthsForSubB[totalProcessesNumber];
    // create full matrix on 0 process
    if (!rank) {
        // parse matrices sizes
        parseAgrs(argc, argv, &matrixSizes.bWidth, &matrixSizes.aHeight, &matrixSizes.side, &errorCode);
        if (errorCode) {
            return errorCode;
        }

        resizeMatrix(A, matrixSizes.side, matrixSizes.aHeight, &errorCode);
        resizeMatrix(B, matrixSizes.bWidth, matrixSizes.side, &errorCode);
        resizeMatrix(C, matrixSizes.bWidth, matrixSizes.aHeight, &errorCode);
        if (errorCode) {
            return errorCode;
        }

        fillMatrix(A, A_SEED);
        fillMatrix(B, B_SEED);

        int commRowsRank = 0;
        int commColumnsRank = 0;
        int coords2D[2] = {0, 0};
        for (int i = 0; i < totalProcessesNumber; ++i) {
            MPI_Cart_coords(comm2D, i, 2, coords2D);
            MPI_Cart_rank(commRows, &coords2D[0], &commRowsRank);
            MPI_Cart_rank(commColumns, &coords2D[1], &commColumnsRank);

            heightsForSubA[i] = A->height / dims[0] + (A->height % dims[0] > commRowsRank);
            widthsForSubB[i] = B->width / dims[1] + (B->width % dims[1] > commColumnsRank);
        }
    }

    MPI_Bcast(&matrixSizes, 1, MATRIX_SIZES, 0, MPI_COMM_WORLD);
    MPI_Bcast(heightsForSubA, totalProcessesNumber, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(widthsForSubB, totalProcessesNumber, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Type_free(&MATRIX_SIZES);

    // create submatrices
    unsigned subAHeight = heightsForSubA[rank];
    unsigned subBWidth = widthsForSubB[rank];
    Matrix subA = createMatrix(subAHeight, matrixSizes.side, &errorCode);
    Matrix subB = createMatrix(subBWidth, matrixSizes.side, &errorCode); // transpose
    Matrix subC = createMatrix(subAHeight, subBWidth, &errorCode);
    if (errorCode) {
        return errorCode;
    }

    // send B by Y coordinate
    if (coords[0] == 0) {
        MPI_Datatype B_VECTOR_TYPE, type;
        // create MPI type for matrix B
        MPI_Aint lb, vectorSize;
        MPI_Type_vector(matrixSizes.side, 1, matrixSizes.bWidth, MPI_DOUBLE, &type);
        MPI_Type_create_resized(type, 0, sizeof(double), &B_VECTOR_TYPE);
        MPI_Type_commit(&B_VECTOR_TYPE);

        int commColumnsProcessesNumber;
        MPI_Comm_size(commColumns, &commColumnsProcessesNumber);

        int displacementsForB[commColumnsProcessesNumber];
        int countVectorsForB[commColumnsProcessesNumber];
        int displ = 0;
        for (int i = 0; i < commColumnsProcessesNumber; ++i) {
            displacementsForB[i] = displ;
            countVectorsForB[i] = matrixSizes.bWidth / dims[1] + (matrixSizes.bWidth % dims[1] > i);
            displ += countVectorsForB[i];
        }

        MPI_Scatterv(B->data, countVectorsForB, displacementsForB, B_VECTOR_TYPE, subB->data, subB->size, MPI_DOUBLE, 0, commColumns);

        MPI_Type_free(&type);
        MPI_Type_free(&B_VECTOR_TYPE);
    }
    MPI_Bcast(subB->data, subB->size, MPI_DOUBLE, 0, commRows);

    // send A by X coordinate
    if (coords[1] == 0) {
        int commRowsProcessesNumber;
        MPI_Comm_size(commRows, &commRowsProcessesNumber);

        int displacementsForA[commRowsProcessesNumber];
        int sizesForA[commRowsProcessesNumber];
        int displ = 0;
        for (int i = 0; i < commRowsProcessesNumber; ++i) {
            displacementsForA[i] = displ;
            sizesForA[i] = (matrixSizes.aHeight / dims[0] + (matrixSizes.aHeight % dims[0] > i)) * matrixSizes.side;
            displ += sizesForA[i];
        }

        MPI_Scatterv(A->data, sizesForA, displacementsForA, MPI_DOUBLE, subA->data, subA->size, MPI_DOUBLE, 0, commRows);
    }
    MPI_Bcast(subA->data, subA->size, MPI_DOUBLE, 0, commColumns);

    multTransposeMatrixOnMatrix(subC, subA, subB);

    // gather matrix C
    if (rank == 0) {
        MPI_Datatype C_SUBARRAY_TYPE[totalProcessesNumber];
        int senderCoords[2] = {0, 0};
        for (int i = 1; i < totalProcessesNumber; ++i) {
            int sizes[2] = {matrixSizes.aHeight, matrixSizes.bWidth};
            int subsizes[2] = {heightsForSubA[i], widthsForSubB[i]};

            MPI_Type_create_subarray(2, sizes, subsizes, senderCoords, MPI_ORDER_C, MPI_DOUBLE, &C_SUBARRAY_TYPE[i]);
            MPI_Type_commit(&C_SUBARRAY_TYPE[i]);
        }

        copy(C, subC);
        int displ = 0;
        for (int i = 1; i < totalProcessesNumber; ++i) {
            MPI_Cart_coords(comm2D, i, 2, senderCoords);

            displ = senderCoords[0] * matrixSizes.bWidth * heightsForSubA[0] + senderCoords[1] * widthsForSubB[0];
            MPI_Recv(&(C->data[displ]), 1, C_SUBARRAY_TYPE[i], i, 1, comm2D, MPI_STATUS_IGNORE);
        }

        for (int i = 1; i < totalProcessesNumber; ++i) {
            MPI_Type_free(&C_SUBARRAY_TYPE[i]);
        }
    } else {
        MPI_Send(subC->data, subC->size, MPI_DOUBLE, 0, 1, comm2D);
    }
    
    int returnedValue = OK;
    if (rank == 0) {
        Matrix trueC = createMatrix(C->height, C->width, &errorCode);
        if (errorCode) {
            return errorCode;
        }
        
        multMatrixOnMatrix(trueC, A, B);

        if (trueC->size < SIZE_CAP) {
            printf("Matrix A\n==================");
            printMatrix(A);
            printf("Matrix B\n==================");
            printMatrix(B);
            printf("Matrix C\n==================");
            printMatrix(C);
            printf("True matrix C\n==================");
            printMatrix(trueC);
        }

        if (isCorrectCalcualtion(C, trueC)) {
            printf("*********CORRECT*********\n");
        } else {
            printf("********INCORRECT********\n");
            returnedValue = INCORRECT_CALCULATION;
        }
        
        destroyMatrix(trueC);
        destroyMatrix(A);
        destroyMatrix(B);
        destroyMatrix(C);
    }
    destroyMatrix(subA);
    destroyMatrix(subB);
    destroyMatrix(subC);
    MPI_Comm_free(&comm2D);
    MPI_Comm_free(&commRows);
    MPI_Comm_free(&commColumns);
    MPI_Finalize();

    return returnedValue;
}