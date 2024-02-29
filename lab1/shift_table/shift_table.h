#pragma once

#include "../settings.h"

#if defined(MPI_V1) || defined(MPI_V2)

typedef struct {
    unsigned* shifts;
    unsigned* countElementsForProcess;
    unsigned length;
#ifdef MPI_V2
    unsigned maxCountElementsForProcess;
#endif
} ShiftTable;

void InitShiftTable(ShiftTable* table, const unsigned countRowsInMatrix);
void DestroyShiftTable(ShiftTable* table);
unsigned GetShiftForCurrRank(const ShiftTable* table);

#endif
