#include "parser/parser.h"
#include "matrix/matrix.h"
#include "vector/vector.h"
#include "settings.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

#ifdef DEFAULT
#include <time.h>
#endif

#if defined(OpenMP_V1) || defined (OpenMP_V2)
#include <omp.h>
#endif

const double E = 1e-5;

double IterationAlgo(Vector* x, Vector* b, const Matrix* A) {
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

    double normB = 0;
    double normR = 0;
    double productR = 0;

    #if defined(OpenMP_V1) || defined(DEFAULT)
    MultMatrixOnVector(&res, x, A);
    SubtractionWithMultOnConst(&r, b, &res, 1, 1);
    Copy(&z, &r);

    #pragma omp parallel for reduction(+:normB) reduction(+:normR)
    for (int i = 0; i < b->length; ++i) {
        normB += b->data[i] * b->data[i];
        normR += r.data[i] * r.data[i];
    }

    productR = normR;
    normB = sqrt(normB);
    normR = sqrt(normR);

    while (normR/normB > E) {
        const double mult = productR;

        MultMatrixOnVector(&res, &z, A);
        alpha = mult / ScalarProduct(&res, &z);

        Copy(&res1, &z);
        AdditionWithMultOnConst(x, x, &res1, 1, alpha);

        SubtractionWithMultOnConst(&r, &r, &res, 1, alpha);

        productR = ScalarProduct(&r, &r);

        beta = productR / mult;

        AdditionWithMultOnConst(&z, &r, &z, 1, beta);

        normR = sqrt(productR);
    }
    #elif defined(OpenMP_V2)
    double result = 0;
    #pragma omp parallel firstprivate(alpha) firstprivate(beta)
    {
        MultMatrixOnVector(&res, x, A);
        SubtractionWithMultOnConst(&r, b, &res, 1, 1);
        Copy(&z, &r);

        #pragma omp for reduction(+:normB) reduction(+:normR)
        for (int i = 0; i < b->length; ++i) {
            normB += b->data[i] * b->data[i];
            normR += r.data[i] * r.data[i];
        }

        #pragma omp single
        {
            productR = normR;
            normB = sqrt(normB);
            normR = sqrt(normR);
        }

        while (normR/normB > E) {
            const double mult = productR;

            MultMatrixOnVector(&res, &z, A);

            #pragma omp single
            result = 0;

            #pragma omp for reduction(+:result)
            for (int i = 0; i < res.length; ++i) {
                result += res.data[i] * z.data[i];
            }

            alpha = mult / result;

            Copy(&res1, &z);
            AdditionWithMultOnConst(x, x, &res1, 1, alpha);

            SubtractionWithMultOnConst(&r, &r, &res, 1, alpha);

            #pragma omp single
            productR = 0;

            #pragma omp for reduction(+:productR)
            for (int i = 0; i < r.length; ++i) {
                productR += r.data[i] * r.data[i];
            }

            beta = productR / mult;

            AdditionWithMultOnConst(&z, &r, &z, 1, beta);

            #pragma omp single
            normR = sqrt(productR);
        }
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

    Vector x;
    Vector b;
    InitVector(&b, plateLength * plateWidth);
    #ifdef RANDOM
    FillVector(&b);
    #elif defined(STATIC)
    int* sumRows = (int*)calloc((unsigned long)plateLength * plateWidth, sizeof(int));
    SumRows(sumRows, &A);

    FillVector(&b, sumRows);
    free(sumRows);
    #endif

    int countThreads[5] = {1, 2, 4, 8, 12};
    double minTimes[5] = {INT_MAX, INT_MAX, INT_MAX, INT_MAX, INT_MAX};
    double workTime = 0;

    printf("==============\n");
    for (int j = 0; j < 5; ++j) {
        omp_set_num_threads(countThreads[j]);
        for (int i = 0; i < 5; ++i) {
            InitVector(&x, plateLength * plateWidth);
            workTime = IterationAlgo(&x, &b, &A);
            if (minTimes[j] > workTime) {
                minTimes[j] = workTime;
            }
            DestoryVector(&x);
        }
        printf("Time taken: %lf sec. Count threads: %d\n", minTimes[j], countThreads[j]);
        printf("==============\n");
    }

   
    DestoryVector(&b);
    DestroyMatrix(&A);

    return EXIT_SUCCESS;
}