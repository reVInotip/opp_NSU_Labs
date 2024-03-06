#include "parser/parser.h"
#include "matrix/matrix.h"
#include "vector/vector.h"
#include "settings.h"
#include <stdio.h>
#include <stdlib.h>

#ifdef  DEFAULT
    #include <time.h>
#elif defined(MPI_V1) || defined(MPI_V2)
    #include <mpi.h>
    #include "shift_table/shift_table.h"
#endif

const double E = 1e-5;

#ifdef DEFAULT
    void IterationAlgo(Vector* x, const Vector* b, const Matrix* A) {
        struct timespec start;
        struct timespec end;
        Vector r;
        Vector z;
        double alpha = 0;
        double beta = 0;

        InitVector(&r, b->length);
        InitVector(&z, b->length);

        Vector res;
        Vector res1;
        InitVector(&res, b->length);
        InitVector(&res1, b->length);

        clock_gettime(CLOCK_MONOTONIC_RAW, &start);

        MultMatrixOnVector(&res, x, A);
        Subtraction(&r, b, &res);
        Copy(&z, &r);

        while (Norm(&r)/Norm(b) > E) {
            const double mult = ScalarProduct(&r, &r);

            MultMatrixOnVector(&res, &z, A);
            alpha = mult / ScalarProduct(&res, &z);

            Copy(&res1, &z);
            MultOnConst(&res1, alpha);
            Addition(x, x, &res1);

            MultOnConst(&res, alpha);
            Subtraction(&r, &r, &res);

            beta = ScalarProduct(&r, &r) / mult;

            MultOnConst(&z, beta);
            Addition(&z, &r, &z);
        }

        clock_gettime(CLOCK_MONOTONIC_RAW, &end);
        printf("==============\n");
        printf("Time taken: %lf sec.\n",
            end.tv_sec - start.tv_sec + 0.000000001*(end.tv_nsec - start.tv_nsec));
        printf("==============\n");

        DestoryVector(&r);
        DestoryVector(&z);
        DestoryVector(&res);
        DestoryVector(&res1);
    }
#elif defined(MPI_V1) || defined(MPI_V2)
    void IterationAlgo(Vector* x, Vector* b, const Matrix* A, const ShiftTable* table) {
        double start = 0;
        double end = 0;
        Vector r;
        Vector z;
        double alpha = 0;
        double beta = 0;

        InitVector(&r, table, b->length);
        InitVector(&z, table, b->length);

        Vector res;
        Vector res1;
        InitVector(&res, table, b->length);
        InitVector(&res1, table, b->length);

        start = MPI_Wtime();

        MultMatrixOnVector(&res, x, A);
        int rank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        Subtraction(&r, b, &res);
        Copy(&z, &r);

        while (Norm(&r)/Norm(b) > E) {
            const double mult = ScalarProduct(&r, &r);

            MultMatrixOnVector(&res, &z, A);
            alpha = mult / ScalarProduct(&res, &z);

            Copy(&res1, &z);
            MultOnConst(&res1, alpha);
            Addition(x, x, &res1);

            MultOnConst(&res, alpha);
            Subtraction(&r, &r, &res);

            beta = ScalarProduct(&r, &r) / mult;

            MultOnConst(&z, beta);
            Addition(&z, &r, &z);
        }

        end = MPI_Wtime();

        rank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (!rank) {
            printf("==============\n");
            printf("Time taken: %lf sec.\n", end - start);
            printf("==============\n");
        }

        DestoryVector(&r);
        DestoryVector(&z);
        DestoryVector(&res);
        DestoryVector(&res1);
    }
#endif

int main(int argc, char** argv) {
    unsigned plateWidth = 0;
    unsigned plateLength = 0;

    int code = ParseAgrs(argc, argv, &plateWidth, &plateLength);
    if (code) {
        return code;
    }

    #ifdef DEFAULT
        Matrix A;
        InitMatrix(&A, plateLength, plateWidth);
        printf("Matrix A\n");
        PrintMatrix(&A);
    #elif defined(MPI_V1) || defined(MPI_V2)
        MPI_Init(&argc, &argv);

        int size = 0;
        int rank = 0;
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        const unsigned countRowsInMatrix = plateLength * plateWidth;

        if (size > countRowsInMatrix) {
            if (!rank) {
                printf("bad count process");
            }
            return EXIT_FAILURE;
        }

        ShiftTable* table = calloc(1, sizeof(ShiftTable));
        InitShiftTable(table, countRowsInMatrix);

        Matrix A;
        const unsigned countRowsForProcess = countRowsInMatrix / size + ((countRowsInMatrix % size) > rank);
        InitMatrix(&A, table, plateLength, countRowsForProcess, plateWidth);

        printf("Matrix A_%d\n", rank);
        PrintMatrix(&A);
    #endif

    Vector x;
    Vector b;
    #ifdef DEFAULT
        InitVector(&x, plateLength * plateWidth);
        InitVector(&b, plateLength * plateWidth);
        #ifdef RANDOM
            FillVector(&b);
        #elifdef STATIC
            int* sumRows = (int*)calloc((unsigned long)plateLength * plateWidth, sizeof(int));
            SumRows(sumRows, &A);

            FillVector(&b, sumRows);
            free(sumRows);
        #endif

        IterationAlgo(&x, &b, &A);
        printf("\n vector b\n");
        PrintVector(&b);
        printf("\nvector x (answer)\n");
        PrintVector(&x);

        DestoryVector(&x);
        DestoryVector(&b);
        DestroyMatrix(&A);
    #elifdef MPI_V1
        InitVector(&x, table, plateLength * plateWidth);
        InitVector(&b, table, plateLength * plateWidth);
        #ifdef RANDOM
            if (!rank) {
                FillVector(&b);
            }
        #elifdef STATIC
            int* sumRows = (int*)calloc((unsigned long)plateLength * plateWidth, sizeof(int));
            SumRows(sumRows, &A);

            MPI_Allgatherv(&(sumRows)[A.table->shifts[rank]], (int)A.blockWidth, MPI_INT, //send
            &(sumRows)[0], (const int*)A.table->countElementsForProcess,
            (const int*)A.table->shifts, MPI_INT, MPI_COMM_WORLD); //recvest

            if (!rank) {
                FillVector(&b, sumRows);
            }
            free(sumRows);
        #endif


        MPI_Bcast(b.data, (int)b.length, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        IterationAlgo(&x, &b, &A, table);

        if (!rank) {
            printf("\n vector b\n");
            PrintVector(&b);
            printf("\nvector x (answer)\n");
            PrintVector(&x);
        }

        DestoryVector(&x);
        DestoryVector(&b);
        DestroyMatrix(&A);
        DestroyShiftTable(table);
        free(table);

        MPI_Finalize();
    #elifdef MPI_V2
        InitVector(&x, table, countRowsForProcess);
        InitVector(&b, table, countRowsForProcess);
        #ifdef RANDOM
            FillVector(&b);
        #elifdef STATIC
            int* sumRows = (int*)calloc((unsigned long)countRowsForProcess, sizeof(int));
            SumRows(sumRows, &A);

            FillVector(&b, sumRows);
            free(sumRows);
        #endif

        IterationAlgo(&x, &b, &A, table);
        printf("\n vector b_%d\n", rank);
        PrintVector(&b);
        printf("\nvector x_%d (answer)\n", rank);
        PrintVector(&x);

        DestoryVector(&x);
        DestoryVector(&b);
        DestroyMatrix(&A);
        DestroyShiftTable(table);
        free(table);

        MPI_Finalize();
    #endif
    return EXIT_SUCCESS;
}
