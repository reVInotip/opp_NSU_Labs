#include "parser/parser.h"
#include <stdio.h>
#include <stdlib.h>

void FillMatrix() {

}

int main(int argc, char** argv) {
    unsigned width = 0, length = 0;

    int code = ParseAgrs(argc, argv, &width, &length);
    if (code) {
        return code; 
    }

    

    return EXIT_SUCCESS;
}