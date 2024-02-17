#pragma once

enum ERR_CODE_PARS {
    MissLength = 4,
    MissWidth = 3,
    BadArgument = 2,
    BadCountOfArgs = 1,
    OK = 0,
};

int ParseAgrs(int argc, char** argv, unsigned* width, unsigned* length);