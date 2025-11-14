#include "gaussian_elimination.h"
#include <math.h>

int findPartialPivot(double** rows, int r1, int r2, int col) {
    if (fabs(rows[r1][col]) > fabs(rows[r2][col])) {
        return r1;
    } else if (fabs(rows[r2][col]) > fabs(rows[r1][col])) {
        return r2;
    }
    return -1;
}

void swapRows(double** rows, int r1, int r2) {
    double* tmp = rows[r1];
    rows[r1] = rows[r2];
    rows[r2] = tmp;
}

void eliminateRow(double **rows, int r1, int r2, int col, int rLen){
    double factor = rows[r2][col] / rows[r1][col];
    for (int i = 0; i < rLen; i++) {
        rows[r2][i] -= factor * rows[r1][i];
    }
}
