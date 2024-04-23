#include "parser.h"
#include <stdbool.h>
#include <stdio.h>

enum ERR_CODE_PARS {
    MissN_x = 5,
    MissN_y = 4,
    MissN_z = 3,
    BadArgument = 2,
    BadCountOfArgs = 1,
    OK = 0,
};

int toNumber(const char *strNumber) {
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

int parseAgrs(int argc, char **argv, unsigned *N_x, unsigned *N_y, unsigned *N_z) {
    if (argc != 7) {
        return BadCountOfArgs;
    }

    bool isN_x = false, isN_y = false, isN_z = false;
    for (int i = 0; i < argc; ++i) {
        if (argv[i][0] == '-') {
            int j = 1;
            while (argv[i][j] != '\0' && j < 4) {
                if (argv[i][j] == 'x') {
                    isN_x = true;
                    break;
                } else if (argv[i][j] == 'y') {
                    isN_y = true;
                    break;
                } else if (argv[i][j] == 'z') {
                    isN_z = true;
                    break;
                } else if (argv[i][j] == 'N' || argv[i][j] == '_') {
                    ++j;
                } else {
                    return BadArgument;
                }
            }
            if (j >= 4) {
                return BadArgument;
            }
        } else if (isN_x) {
            *N_x = toNumber(argv[i]);
            isN_x = false;
        } else if (isN_y) {
            *N_y = toNumber(argv[i]);
            isN_y = false;
        } else if (isN_z) {
            *N_z = toNumber(argv[i]);
            isN_z = false;
        }
    }

    if (*N_x == 0) {
        return MissN_x;
    } else if (*N_y == 0) {
        return MissN_y;
    } else if (*N_z == 0) {
        return MissN_z;
    }

    return OK;
}
