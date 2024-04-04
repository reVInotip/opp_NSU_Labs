#pragma once

enum ERR_CODE_PARS {
    MISS_HEIGHT = 5,
    MISS_WIDTH = 4,
    MISS_SIDE = 3,
    BAD_ARGUMENT = 2,
    BAD_COUNT_OF_ARGS = 1,
};

void parseAgrs(int argc, char **argv, unsigned *bWidth, unsigned *aHeight, unsigned *side, int *errCode);