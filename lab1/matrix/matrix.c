#include "matrix.h"
#include "../vector/vector.h"
#include "../settings.h"
#include <assert.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>

#ifdef MPI_V1
    #include <mpi.h>
#elifdef MPI_V2
    #include <mpi.h>
    #define TAG 101
#endif

#ifdef DEFAULT
    void InitMatrix(Matrix* matr, const unsigned blockLength, unsigned countBlocksInRow) {
        assert(matr);
        matr->size = (unsigned long)blockLength * blockLength * countBlocksInRow * countBlocksInRow;
        matr->matrix = calloc(
            matr->size,
            sizeof(int));
        matr->blockLength = blockLength;
        matr->countBlocksInRow = countBlocksInRow;

        unsigned offset = blockLength * countBlocksInRow;
        int currentOffset = 0;
        unsigned long fakeMatrixOneSize =
            (unsigned long)blockLength * blockLength *
            (countBlocksInRow - 1) * countBlocksInRow;
        long fakeMatrixTwoSize =
            (long)blockLength * blockLength * countBlocksInRow;
        for (int i = 0; i < matr->size; i += (int)offset) {
            if (i + currentOffset < matr->size) {
                matr->matrix[i + currentOffset] = -4;
            }
            if (
                (i + 1 + currentOffset < matr->size) &&
                (currentOffset % blockLength != blockLength - 1)
            ) {
                matr->matrix[i + 1 + currentOffset] = 1;
            }
            if (
                (i - 1 + currentOffset > 0) &&
                (currentOffset % blockLength != 0)
            ) {
                matr->matrix[i - 1 + currentOffset] = 1;
            }
            if (i + blockLength + currentOffset < fakeMatrixOneSize) {
                matr->matrix[i + blockLength + currentOffset] = 1;
            }
            if (i - (int)blockLength + currentOffset >= fakeMatrixTwoSize) {
                matr->matrix[i - blockLength + currentOffset] = 1;
            }
            ++currentOffset;
        }
    }

    void DestroyMatrix(Matrix* matr) {
        assert(matr);
        assert(matr->matrix);
        free(matr->matrix);
        matr->matrix = NULL;
        matr->blockLength = 0;
        matr->countBlocksInRow = 0;
    }

    void MultMatrixOnVector(Vector* result, const Vector* vec, const Matrix* matr) {
        assert(result);
        assert(vec);
        assert(matr);

        unsigned matrRowLength = matr->blockLength * matr->countBlocksInRow;
        if (
            (matrRowLength != vec->length) ||
            (result->length != vec->length)
        ) {
            printf("bad args\n");
            return;
        }

        for (int i = 0; i < matrRowLength; ++i) {
            for (int j = 0; j < matrRowLength; ++j) {
                result->data[i] += vec->data[j] * (double)matr->matrix[matrRowLength * i + j];
            }
        }
    }
#elif defined(MPI_V1) || defined(MPI_V2)
    void InitMatrix(Matrix* matr, const ShiftTable* table, const unsigned blockLength, const unsigned blockWidth, const unsigned countBlocksInRow) {
        assert(matr);

        matr->size = (unsigned long)blockLength * blockWidth * countBlocksInRow;
        matr->matrix = calloc(
            matr->size,
            sizeof(int));
        matr->blockLength = blockLength;
        matr->blockWidth = blockWidth;
        matr->countBlocksInRow = countBlocksInRow;
        matr->table = table;

        int rank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        const unsigned offset = blockLength * countBlocksInRow;
        int currentOffset = table->shifts[rank];

        int countRows = currentOffset;
        for (int i = 0; i < matr->size; i += (int)offset) {
            if (i + currentOffset < matr->size) {
                matr->matrix[i + currentOffset] = -4;
            }
            if (
                (i + 1 + currentOffset < matr->size) &&
                (currentOffset % blockLength != blockLength - 1)
            ) {
                matr->matrix[i + 1 + currentOffset] = 1;
            }
            if (
                (i - 1 + currentOffset > 0) &&
                (currentOffset % blockLength != 0)
            ) {
                matr->matrix[i - 1 + currentOffset] = 1;
            }
            if (
                (i / offset + countRows < blockLength * (countBlocksInRow - 1)) &&
                (i + blockLength + currentOffset < matr->size)
            ) {
                matr->matrix[i + blockLength + currentOffset] = 1;
            }
            if (
                (i / offset + countRows >= blockLength) &&
                (i - (int)blockLength + currentOffset < matr->size)
            ) {
                matr->matrix[i - blockLength + currentOffset] = 1;
            }
            ++currentOffset;
        }
    }

    void DestroyMatrix(Matrix* matr) {
        assert(matr);
        assert(matr->matrix);

        free(matr->matrix);
        matr->matrix = NULL;
        matr->blockLength = 0;
        matr->countBlocksInRow = 0;
        matr->blockWidth = 0;
        matr->size = 0;
        matr->table = NULL;
    }
#endif

#ifdef MPI_V1
    void MultMatrixOnVector(Vector* result, Vector* vec, const Matrix* matr) {
        assert(result);
        assert(vec);
        assert(matr);

        const unsigned matrRowLength = matr->blockLength * matr->countBlocksInRow;
        if (
            (matrRowLength != vec->realLength) ||
            (result->realLength != vec->realLength)
        ) {
            printf("bad args\n");
            return;
        }

        int rank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        unsigned startIndex = matr->table->shifts[rank];

        if ((startIndex + matr->blockWidth > result->realLength) || (startIndex > result->realLength)) {
            printf("bad args\n");
            return;
        }

        if (matr->blockWidth == vec->length) {
            MPI_Allgatherv(&(vec->data)[startIndex], (int)matr->blockWidth, MPI_DOUBLE, //send
            &(vec->data)[0], (const int*)matr->table->countElementsForProcess,
            (const int*)matr->table->shifts, MPI_DOUBLE, MPI_COMM_WORLD); //recvest
            vec->length = vec->realLength;
        }
        result->length = matr->blockWidth;

        for (unsigned i = 0; i < matr->blockWidth; ++i) {
            for (unsigned j = 0; j < matrRowLength; ++j) {
                result->data[startIndex] += vec->data[j] * (double)matr->matrix[matrRowLength * i + j];
            }
            ++startIndex;
        }
    }
#elifdef MPI_V2
    void MultMatrixOnVector(Vector* result, Vector* vec, const Matrix* matr) {
        assert(result);
        assert(vec);
        assert(matr);

        const unsigned matrRowLength = matr->blockLength * matr->countBlocksInRow;
        if (
            (matr->blockWidth != vec->length)
        ) {
            printf("bad args1\n");
            return;
        }
        int rank = 0;
        int size = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        int currRank = rank;
        for (int process = 0; process < size; ++process) {
            const unsigned startIndex = vec->table->shifts[currRank];
            const unsigned endIndex = startIndex + vec->table->countElementsForProcess[currRank];
            vec->length = vec->table->countElementsForProcess[currRank];

            for (unsigned i = 0; i < matr->blockWidth; ++i) {
                for (unsigned j = startIndex, k = 0; (j < endIndex) && (k < vec->length); ++j, ++k) {
                    result->data[i] += vec->data[k] * (double)matr->matrix[matrRowLength * i + j];
                }
            }

            MPI_Sendrecv_replace(vec->data, (int)vec->table->maxCountElementsForProcess, MPI_DOUBLE, (rank + size - 1) % size, TAG,
                (rank + 1) % size, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            currRank = (currRank + 1) % size;
        }

        vec->length = vec->table->countElementsForProcess[currRank];
    }
#endif

void PrintMatrix(const Matrix* matr) {
    assert(matr);
    for (unsigned i = 0; i < matr->size; ++i) {
        if (i % (matr->blockLength * matr->countBlocksInRow) == 0) {
            printf("\n");
        }
        printf("%d ", matr->matrix[i]);
    }
    printf("\n");
}

#ifdef STATIC
    #ifdef DEFAULT
        void SumRows(int* sumRows, const Matrix* matr) {
            const unsigned matrRowLength = matr->blockLength * matr->countBlocksInRow;

            for (unsigned i = 0; i < matr->blockLength * matr->countBlocksInRow; ++i) {
                for (unsigned j = 0; j < matrRowLength; ++j) {
                    sumRows[i] += matr->matrix[matrRowLength * i + j];
                }
            }
        }
    #elifdef MPI_V2
        void SumRows(int* sumRows, const Matrix* matr) {
            const unsigned matrRowLength = matr->blockLength * matr->countBlocksInRow;

            for (unsigned i = 0; i < matr->blockWidth; ++i) {
                for (unsigned j = 0; j < matrRowLength; ++j) {
                    sumRows[i] += matr->matrix[matrRowLength * i + j];
                }
            }
        }
    #elifdef MPI_V1
        void SumRows(int* sumRows, const Matrix* matr) {
            const unsigned matrRowLength = matr->blockLength * matr->countBlocksInRow;

            int rank = 0;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            unsigned startIndex = matr->table->shifts[rank];

            for (unsigned i = 0; i < matr->blockWidth; ++i) {
                for (unsigned j = 0; j < matrRowLength; ++j) {
                    sumRows[startIndex] += matr->matrix[matrRowLength * i + j];
                }
                ++startIndex;
            }
        }
    #endif
#endif
