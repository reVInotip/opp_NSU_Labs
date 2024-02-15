#include "parser/parser.h"
#include "matrix/matrix.h"
#include "vector/vectror.h"
#include <stdio.h>
#include <stdlib.h>

#define E 0.00001

void IterationAlgo(Vector* x, const Vector* b, const Matrix* A) {
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

    MultMatrixOnVector(&res, x, A);
    Subtraction(&r, b, &res);
    Copy(&z, &r);

    while (Norm(&r)/Norm(b) < E) {
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

    DestoryVector(&r);
    DestoryVector(&z);
    DestoryVector(&res);
    DestoryVector(&res1);
}

int main(int argc, char** argv) {
    unsigned blockWidth = 0;
    unsigned blockLength = 0;

    int code = ParseAgrs(argc, argv, &blockWidth, &blockLength);
    if (code) {
        return code;
    }

    Matrix A;
    InitMatrix(&A, blockLength, blockWidth);

    Vector x;
    Vector b;
    InitVector(&x, blockLength * blockWidth);
    InitVector(&b, blockLength * blockWidth);
    FillVector(&b);

    IterationAlgo(&x, &b, &A);
    PrintVector(&x);

    DestoryVector(&x);
    DestoryVector(&b);
    DestroyMatrix(&A);

    return EXIT_SUCCESS;
}
