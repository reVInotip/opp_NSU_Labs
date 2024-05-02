#include "../include/buffer.h"
#include <stdlib.h>
#include <assert.h>

void initBuffer(Buffer *buf, const unsigned size) {
    buf->send = (double*)calloc(size, sizeof(double));
    assert(buf->send);
    buf->recieve = (double*)calloc(size, sizeof(double));
    assert(buf->recieve);

    buf->size = size;
}

void destroyBuffer(Buffer *buf) {
    assert(buf);

    free(buf->send);
    free(buf->recieve);
    buf->size = 0;
}