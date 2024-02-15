#include "parser/parser.h"
#include "matrix/matrix.h"
#include "vector/vectror.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv) {
    unsigned blockWidth = 0;
    unsigned blockLength = 0;

    int code = ParseAgrs(argc, argv, &blockWidth, &blockLength);
    if (code) {
        return code;
    }

    Vector vec;
    InitVector(&vec, blockLength * blockWidth);
    PrintVector(&vec);
    FillVector(&vec);
    printf("------\n");
    PrintVector(&vec);
    DestoryVector(&vec);

    return EXIT_SUCCESS;
}
