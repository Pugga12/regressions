#pragma once
int findPartialPivot(double** rows, int r1, int r2, int col);
void swapRows(double** rows, int r1, int r2);
void eliminateRow(double** rows, int r1, int r2, int col, int rLen);