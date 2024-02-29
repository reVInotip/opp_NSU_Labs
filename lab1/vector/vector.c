#include "vector.h"
#include "../settings.h"
#include <assert.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>

#if defined(MPI_V1) || defined (MPI_V2)
    #include <mpi.h>
    #include "../shift_table/shift_table.h"

    #define TAG 102
#endif

#define FAKE_BORDER 100
#define REAL_BORDER 50

#if defined(MPI_V1) || defined (MPI_V2)
void InitVector(Vector* vec, const ShiftTable* table, const unsigned length) {
    assert(vec);

    vec->data = calloc((unsigned long)length, sizeof(double));
    vec->length = length;
    #ifdef MPI_V1
        vec->realLength = length;
    #endif
    vec->table = table;
}
#elifdef DEFAULT
    void InitVector(Vector* vec, const unsigned length) {
        assert(vec);

        vec->data = calloc((unsigned long)length, sizeof(double));
        vec->length = length;
    }
#endif

void DestoryVector(Vector* vec) {
    assert(vec);
    assert(vec->data);

    free(vec->data);
    vec->data = NULL;
    vec->length = 0;
    #ifdef MPI_V1
        vec->realLength = 0;
    #endif
    #if defined(MPI_V1) || defined (MPI_V2)
        vec->table = NULL;
    #endif
}

#ifdef DEFAULT
    #ifdef RANDOM
        void FillVector(Vector *vec) {
            assert(vec);
            assert(vec->data);
            srand(time(NULL));

            const unsigned countNumbers = (unsigned int)(vec->length * 0.2) + 1;
            unsigned index = 0;
            int value = 0;
            for (unsigned i = 0; i < countNumbers; ++i) {
                index = rand() % (vec->length);
                value = (rand() % FAKE_BORDER) - REAL_BORDER;
                vec->data[index] = value;
            }
        }
    #elifdef STATIC
        void FillVector(Vector *vec, const int* rowSums) {
            assert(vec);
            assert(vec->data);

            for (unsigned i = 0; i < vec->length; ++i) {
                vec->data[i] = rowSums[i];
            }
        }
    #endif

    void Subtraction(Vector *result, const Vector *vec1, const Vector *vec2) {
        assert(result);
        assert(vec1);
        assert(vec2);

        unsigned minLength = vec1->length > vec2->length ? vec2->length : vec1->length;
        if (minLength > result->length) {
            minLength = result->length;
        }

        for (unsigned i = 0; i < minLength; ++i) {
            result->data[i] = vec1->data[i] - vec2->data[i];
        }
    }

    void Addition(Vector* result, const Vector* vec1, const Vector* vec2) {
        assert(result);
        assert(vec1);
        assert(vec2);

        unsigned minLength = vec1->length > vec2->length ? vec2->length : vec1->length;
        if (minLength > result->length) {
            minLength = result->length;
        }

        for (unsigned i = 0; i < minLength; ++i) {
            result->data[i] = vec1->data[i] + vec2->data[i];
        }
    }

    void Copy(Vector* result, const Vector* vec) {
        assert(result);
        assert(vec);

        if (vec->length != result->length) {
            double* newData = realloc(result->data, vec->length * sizeof(int));
            if (newData) {
                result->data = newData;
            } else {
                return;
            }
        }
        result->length = vec->length;

        for (unsigned i = 0; i < vec->length; ++i) {
            result->data[i] = vec->data[i];
        }
    }

    double ScalarProduct(const Vector* vec1, const Vector* vec2) {
        assert(vec1);
        assert(vec1->data);
        assert(vec2);
        assert(vec2->data);

        const unsigned minLength = vec1->length > vec2->length ? vec2->length : vec1->length;
        double result = 0;
        for (int i = 0; i < minLength; ++i) {
            result += vec1->data[i] * vec2->data[i];
        }

        return result;
    }

    double Norm(const Vector* vec) {
        assert(vec);
        assert(vec->data);

        double result = 0;
        for (int i = 0; i < vec->length; ++i) {
            result += vec->data[i] * vec->data[i];
        }

        result = sqrt(result);
        return result;
    }

    void MultOnConst(Vector* result, const double value) {
        assert(result);

        for (int i = 0; i < result->length; ++i) {
            result->data[i] *= value;
        }
    }

    void PrintVector(const Vector *vec) {
        assert(vec);
        assert(vec->data);

        for (unsigned i = 0; i < vec->length; ++i) {
            printf("%lf ", vec->data[i]);
        }
        printf("\n");
    }
#elifdef MPI_V1
    #ifdef RANDOM
        void FillVector(Vector *vec) {
            assert(vec);
            assert(vec->data);
            srand(time(NULL));

            const unsigned countNumbers = (unsigned int)(vec->length * 0.2) + 1;
            unsigned index = 0;
            int value = 0;
            for (unsigned i = 0; i < countNumbers; ++i) {
                index = rand() % (vec->length);
                value = (rand() % FAKE_BORDER) - REAL_BORDER;
                vec->data[index] = value;
            }
        }
    #elifdef STATIC
        void FillVector(Vector *vec, const int* rowSums) {
            assert(vec);
            assert(vec->data);

            for (unsigned i = 0; i < vec->length; ++i) {
                vec->data[i] = rowSums[i];
            }
        }
    #endif

    void Subtraction(Vector *result, const Vector *vec1, const Vector *vec2) {
        assert(result);
        assert(vec1);
        assert(vec2);

        int rank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        unsigned startIndex = result->table->shifts[rank];
        unsigned endIndex = startIndex + result->table->countElementsForProcess[rank];
        if (result->length == result->realLength) {
            startIndex = 0;
            endIndex = result->realLength;
        }

        if ((endIndex > vec1->realLength) || (endIndex > vec2->realLength) || (endIndex > result->realLength)) {
            printf("bad args5\n");
            return;
        }

        for (unsigned i = startIndex; i < endIndex; ++i) {
            result->data[i] = vec1->data[i] - vec2->data[i];
        }

        unsigned minLength = vec1->length > vec2->length ? vec2->length : vec1->length;
        if (minLength > result->length) {
            minLength = result->length;
        }

        result->length = minLength;
    }

    void Addition(Vector* result, const Vector* vec1, const Vector* vec2) {
        assert(result);
        assert(vec1);
        assert(vec2);

        int rank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        unsigned startIndex = result->table->shifts[rank];
        unsigned endIndex = startIndex + result->table->countElementsForProcess[rank];
        if (result->length == result->realLength) {
            startIndex = 0;
            endIndex = result->realLength;
        }

        if ((endIndex > vec1->realLength) || (endIndex > vec2->realLength) || (endIndex > result->realLength)) {
            printf("bad args5\n");
            return;
        }

        for (unsigned i = startIndex; i < endIndex; ++i) {
            result->data[i] = vec1->data[i] + vec2->data[i];
        }

        unsigned minLength = vec1->length > vec2->length ? vec2->length : vec1->length;
        if (minLength > result->length) {
            minLength = result->length;
        }

        result->length = minLength;
    }

    void Copy(Vector* result, const Vector* vec) {
        assert(result);
        assert(vec);

        int rank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        unsigned startIndex = vec->table->shifts[rank];
        unsigned endIndex = startIndex + vec->table->countElementsForProcess[rank];
        if (vec->length == vec->realLength) {
            startIndex = 0;
            endIndex = vec->realLength;
        }

        if ((endIndex > vec->realLength) || (endIndex > vec->realLength) || (endIndex > result->realLength)) {
            printf("bad args5\n");
            return;
        }

        for (unsigned i = startIndex; i < endIndex; ++i) {
            result->data[i] = vec->data[i];
        }

        const unsigned minLength = vec->length > result->length ? result->length : vec->length;
        result->length = minLength;
    }

    double ScalarProduct(Vector* vec1, Vector* vec2) {
        assert(vec1);
        assert(vec1->data);
        assert(vec2);
        assert(vec2->data);

        double result = 0;

        int rank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        const unsigned startIndex = vec1->table->shifts[rank];
        const unsigned endIndex = startIndex + vec1->table->countElementsForProcess[rank];
        if ((endIndex > vec1->realLength) || (endIndex > vec1->realLength)) {
            printf("bad args5\n");
            return 0;
        }


        for (unsigned i = startIndex; i < endIndex; ++i) {
            result += vec1->data[i] * vec2->data[i];
        }

        MPI_Allreduce(&result, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        return result;
    }

    double Norm(const Vector* vec) {
        assert(vec);
        assert(vec->data);

        int rank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        const unsigned startIndex = vec->table->shifts[rank];
        const unsigned endIndex = startIndex + vec->table->countElementsForProcess[rank];
        if (endIndex > vec->realLength) {
            printf("bad args5\n");
            return 0;
        }

        double result = 0;
        for (unsigned i = startIndex; i < endIndex; ++i) {
            result += vec->data[i] * vec->data[i];
        }

        MPI_Allreduce(&result, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        result = sqrt(result);
        return result;
    }

    void MultOnConst(Vector* result, const double value) {
        assert(result);

        int rank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        unsigned startIndex = result->table->shifts[rank];
        unsigned endIndex = startIndex + result->table->countElementsForProcess[rank];
        if (result->length == result->realLength) {
            startIndex = 0;
            endIndex = result->realLength;
        }

        if (endIndex > result->realLength) {
            printf("bad args5\n");
            return;
        }

        for (unsigned i = startIndex; i < endIndex; ++i) {
            result->data[i] *= value;
        }
    }

    void PrintVector(const Vector *vec) {
        assert(vec);
        assert(vec->data);

        int rank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        unsigned startIndex = vec->table->shifts[rank];

        if (vec->length == vec->realLength) {
            startIndex = 0;
        }

        for (unsigned i = 0; i < vec->length; ++i) {
            printf("%lf ", vec->data[startIndex + i]);
        }
        printf("\n");
    }
#elifdef MPI_V2
    double ScalarProduct(Vector* vec1, Vector* vec2) {
        assert(vec1);
        assert(vec1->data);
        assert(vec2);
        assert(vec2->data);

        if (vec1->length != vec2->length) {
            printf("bad args\n");
            return 0;
        }

        double result = 0;

        for (unsigned i = 0; i < vec1->length; ++i) {
            result += vec1->data[i] * vec2->data[i];
        }

        MPI_Allreduce(&result, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        return result;
    }

    #ifdef RANDOM
        void FillVector(Vector *vec) {
            assert(vec);
            assert(vec->data);
            srand(time(NULL));

            int rank = 0;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);

            const unsigned countNumbers = (unsigned int)(vec->length * 0.2) + 1;
            unsigned index = 0;
            int value = 0;
            for (unsigned i = 0; i < countNumbers; ++i) {
                index = (rand() + (rank * 11)) % (vec->length);
                value = ((rand() + (rank * 11)) % FAKE_BORDER) - REAL_BORDER;
                vec->data[index] = value;
            }
        }
    #elifdef STATIC
        void FillVector(Vector *vec, const int* rowSums) {
            assert(vec);
            assert(vec->data);

            for (unsigned i = 0; i < vec->length; ++i) {
                vec->data[i] = rowSums[i];
            }
        }
    #endif

    void Subtraction(Vector *result, const Vector *vec1, const Vector *vec2) {
        assert(result);
        assert(vec1);
        assert(vec2);

        if (
            (vec1->length != vec2->length) ||
            (vec2->length != result->length)
        ) {
            printf("bad args2\n");
            return;
        }

        for (int i = 0; i < result->length; ++i) {
            result->data[i] = vec1->data[i] - vec2->data[i];
        }
    }

    void Addition(Vector* result, const Vector* vec1, const Vector* vec2) {
        assert(result);
        assert(vec1);
        assert(vec2);

        for (int i = 0; i < result->length; ++i) {
            result->data[i] = vec1->data[i] + vec2->data[i];
        }
    }

    void Copy(Vector* result, Vector* vec) {
        assert(result);
        assert(vec);

        if (vec->length != result->length) {
            printf("bad args\n");
            return;
        }

        for (int i = 0; i < vec->length; ++i) {
            result->data[i] = vec->data[i];
        }
    }

    double Norm(Vector* vec) {
        assert(vec);
        assert(vec->data);

        int rank = 0;
        int size = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        int currRank = rank;
        double result = 0;
        for (int process = 0; process < size; ++process) {
            vec->length = vec->table->countElementsForProcess[currRank];

            for (int i = 0; i < vec->length; ++i) {
                result += vec->data[i] * vec->data[i];
            }

            MPI_Sendrecv_replace(vec->data, (int)vec->table->maxCountElementsForProcess, MPI_DOUBLE, (rank + size - 1) % size,
            TAG, (rank + 1) % size, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            currRank = (currRank + 1) % size;
        }

        vec->length = vec->table->countElementsForProcess[currRank];

        MPI_Allreduce(&result, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        result = sqrt(result);
        return result;
    }

    void MultOnConst(Vector* result, const double value) {
        assert(result);

        for (int i = 0; i < result->length; ++i) {
            result->data[i] *= value;
        }
    }

    void PrintVector(const Vector *vec) {
        assert(vec);
        assert(vec->data);

        for (unsigned i = 0; i < vec->length; ++i) {
            printf("%lf ", vec->data[i]);
        }
        printf("\n");
    }
#endif
