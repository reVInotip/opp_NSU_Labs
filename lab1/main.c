#include "parser/parser.h"
#include "matrix/matrix.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv) {
    unsigned blockWidth = 0;
    unsigned blockLength = 0;

    int code = ParseAgrs(argc, argv, &blockWidth, &blockLength);
    if (code) {
        return code;
    }

    Matrix matr;
    printf("%u\n", blockLength * blockWidth);
    InitMatrix(&matr, blockLength, blockWidth);
    PrintMatrix(&matr);
    DestroyMatrix(&matr);

    return EXIT_SUCCESS;
}
