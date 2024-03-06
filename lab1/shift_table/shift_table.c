#include "shift_table.h"
#include "../settings.h"
#include <assert.h>
#include <stdlib.h>
#include <mpi.h>

#if defined(MPI_V1) || defined(MPI_V2)

void InitShiftTable(ShiftTable* table, const unsigned countRowsInMatrix) {
    assert(table);

    int size = 0;
    int rank = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    table->shifts = calloc((unsigned long)size, sizeof(int));
    table->countElementsForProcess = calloc((unsigned long)size, sizeof(int));
    table->length = size;

    unsigned currShift = 0;
    for (unsigned i = 0; i < size; ++i) {
        table->shifts[i] = currShift;
        currShift += countRowsInMatrix / size + ((countRowsInMatrix % size) > i);
        table->countElementsForProcess[i] = countRowsInMatrix / size + ((countRowsInMatrix % size) > i);
    }

#ifdef MPI_V2
    table->maxCountElementsForProcess = (countRowsInMatrix / size) + ((countRowsInMatrix % size) > 0);
#endif
}

void DestroyShiftTable(ShiftTable* table) {
    assert(table);
    assert(table->shifts);

    free(table->shifts);
    free(table->countElementsForProcess);
    table->length = 0;
}

#endif
