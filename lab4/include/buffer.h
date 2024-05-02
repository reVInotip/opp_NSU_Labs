#pragma once

typedef struct buffer {
    double *send;
    double *recieve;
    unsigned size;
} Buffer;

void initBuffer(Buffer *buf, const unsigned size);
void destroyBuffer(Buffer *buf);