#include "parser.h"
#include <stdbool.h>
#include <stdio.h>

int ToNumber(const char* strNumber) {
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

int ParseAgrs(int argc, char** argv, unsigned* width, unsigned* length) {
    if (argc != 5) {
        return BadCountOfArgs;
    }

    bool isLength = false, isWidth = false;
    for (int i = 0; i < argc; ++i) {
        if (argv[i][0] == '-') {
            int j = 0;
            while (argv[i][j] != '\0' && j <= 2) {
                if (argv[i][j] == 'l') {
                    isLength = true;
                } else if (argv[i][j] == 'w') {
                    isWidth = true;
                }
                ++j;
            }
            if (j > 2) {
                return BadArgument;
            }
        } else if (isLength) {
            *length = ToNumber(argv[i]);
            isLength = false;
        } else if (isWidth) {
            *width = ToNumber(argv[i]);
            isWidth = false;
        }
    }

    if (*width == 0) {
        return MissWidth;
    } else if (*length == 0) {
        return MissLength;
    }

    return OK;
}
