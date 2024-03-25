#include "parser/parser.h"
#include "matrix/matrix.h"
#include "vector/vector.h"
#include "settings.h"
#include <stdio.h>
#include <stdlib.h>

#ifdef DEFAULT
    #include <time.h>
#endif

#if defined(OpenMP_V1) || defined (OpenMP_V2)
    #include <omp.h>
#endif

const double E = 1e-5;

double IterationAlgo(Vector* x, const Vector* b, const Matrix* A) {
    #ifdef DEFAULT
        struct timespec start;
        struct timespec end;
    #elif defined(OpenMP_V1) || defined (OpenMP_V2)
        double start = 0;
        double end = 0;
    #endif
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

    #ifdef DEFAULT
        clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    #elif defined(OpenMP_V1) || defined (OpenMP_V2)
        start = omp_get_wtime();
    #endif

    #ifdef OpenMP_V2
        #pragma omp parallel firstprivate(alpha) firstprivate(beta)
        {
    #endif

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

    #ifdef OpenMP_V2
        }
    #endif

    #ifdef DEFAULT
        clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    #elif defined(OpenMP_V1) || defined (OpenMP_V2)
        end = omp_get_wtime();
    #endif

    DestoryVector(&r);
    DestoryVector(&z);
    DestoryVector(&res);
    DestoryVector(&res1);

    #ifdef DEFAULT
        return end.tv_sec - start.tv_sec + 0.000000001*(end.tv_nsec - start.tv_nsec);
    #elif defined(OpenMP_V1) || defined (OpenMP_V2)
        return end - start;
    #endif
}

int main(int argc, char** argv) {
    unsigned plateWidth = 0;
    unsigned plateLength = 0;

    int code = ParseAgrs(argc, argv, &plateWidth, &plateLength);
    if (code) {
        return code;
    }
    Matrix A;
    InitMatrix(&A, plateLength, plateWidth);
    printf("Matrix A\n");
    PrintMatrix(&A);

    Vector x;
    Vector b;
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

    double workTime = IterationAlgo(&x, &b, &A);

    printf("\n vector b\n");
    PrintVector(&b);
    printf("\nvector x (answer)\n");
    PrintVector(&x);

    printf("==============\n");
    printf("Time taken: %lf sec.\n", workTime);
    printf("==============\n");

    DestoryVector(&x);
    DestoryVector(&b);
    DestroyMatrix(&A);

    return EXIT_SUCCESS;
}
