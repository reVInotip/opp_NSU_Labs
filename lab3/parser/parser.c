#include "parser.h"
#include <stdbool.h>
#include <stdio.h>

int toNumber(const char* strNumber) {
    int size = 1;
    while (strNumber[size - 1] != '\0') {
        ++size;
    }

    int deg = 1;
    int result = 0;
    for (int i = size - 2; i >= 0; --i) {
        result += (strNumber[i] - '0') * deg;
        deg *= 10;
    }

    return result;
}

void parseAgrs(int argc, char **argv, unsigned *bWidth, unsigned *aHeight, unsigned *side, int *errCode) {
    if (argc != 7) {
        *errCode = BAD_COUNT_OF_ARGS;
        return;
    }

    // a, b, s
    bool isAHeight = false, isBWidth = false, isABSide = false;
    for (int i = 0; i < argc; ++i) {
        if (argv[i][0] == '-') {
            int j = 1;
            while (argv[i][j] != '\0' && j <= 2) {
                if (argv[i][j] == 'a') {
                    isAHeight = true;
                } else if (argv[i][j] == 'b') {
                    isBWidth = true;
                } else if (argv[i][j] == 's') {
                    isABSide = true;
                }
                ++j;
            }
            if (j > 2) {
                *errCode = BAD_ARGUMENT;
                return;
            }
        } else if (isAHeight) {
            *aHeight = toNumber(argv[i]);
            isAHeight = false;
        } else if (isBWidth) {
            *bWidth = toNumber(argv[i]);
            isBWidth = false;
        } else if (isABSide) {
            *side = toNumber(argv[i]);
            isABSide = false;
        }
    }

    if (*bWidth == 0) {
        *errCode = MISS_WIDTH;
    } else if (*aHeight == 0) {
        *errCode = MISS_HEIGHT;
    } else if (*side == 0) {
        *errCode = MISS_SIDE;
    }

    return;
}
