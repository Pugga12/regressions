#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include "gaussian_elimination.h"

static inline double square(double x) {
	return x * x;
}

static inline double pow3(double x) {
	return x * x * x;
}

static inline double pow4(double x) {
	return x * x * x * x;
}

bool linReg(const int n, const double *x, const double *y, double* m, double* b, double* r) {
	double sumX = 0.0;
	double sumY = 0.0;
	double sumX2 = 0.0;
	double sumXY = 0.0;
	if (m == NULL || b == NULL) return false;

	for (int i = 0;i<n;i++) {
		sumX += x[i];
		sumY += y[i];
		sumX2 += square(x[i]);
		sumXY += x[i] * y[i];
	}

	double slopeDenom = (n * sumX2 - square(sumX));
	if (slopeDenom == 0) {
		*m = 0.0;
		*b = 0.0;
		*r = 0.0;
		return false;
	}
	*m = (n * sumXY - sumX * sumY) / slopeDenom;
	*b = (sumY - *m * sumX) / n;
	return true;
}

bool expReg(const int n, const double *x, const double *y, double* a, double* cr) {
	if (x == NULL || y == NULL || a == NULL || cr == NULL) return false;
	double *y_log = malloc(sizeof(double)*n);
	if (y_log == NULL) {
		printf("malloc of size %d bytes failed", sizeof(double)*n);
		free(y_log);
		return false;
	}
	double m = 0;
	double b = 0;
	double r = 0;

	for (int i = 0; i < n; i++) {
		if (y[i] <= 0) {
			// regression fail (ln is not defined when its input <= 0)
			free(y_log);
			return false;
		}
		y_log[i] = log(y[i]);
	}

	bool result = linReg(n, x, y_log, &m, &b, &r);
	if (!result) {
		// linreg fail
		free(y_log);
		return false;
	} 
	*a = exp(b);
	*cr = exp(m);
	free(y_log);
	return true;
}

bool quadReg(const int n, const double x[], const double y[], double* a, double* b, double* c, double* r) {
	double sumX = 0.0;
	double sumY = 0.0;
	double sumX2 = 0.0;
	double sumX3 = 0.0;
	double sumX4 = 0.0;
	double sumXY = 0.0;
	double sumX2Y = 0.0;
	if (a == NULL || b == NULL || c == NULL) return false;
	for (int i = 0; i<n; i++) {
		sumX += x[i];
		sumY += y[i];
		sumX2 += square(x[i]);
		sumX3 += pow3(x[i]);
		sumX4 += pow4(x[i]);
		sumXY += x[i] * y[i];
		sumX2Y += square(x[i]) * y[i];
	}
	double matrix[3][4] = {
		{sumX4, sumX3, sumX2, sumX2Y},
		{sumX3, sumX2, sumX, sumXY},
		{sumX2, sumX, n, sumY}
	};
	double* rows[3] = {matrix[0], matrix[1], matrix[2]};
	
	if (findPartialPivot(rows, 0, 1, 0) != 0) swapRows(rows, 0, 1);
	eliminateRow(rows, 0, 1, 0, 4);
	eliminateRow(rows, 0, 2, 0, 4);
	if (findPartialPivot(rows, 1, 2, 1) != 0) swapRows(rows, 1, 2);
	eliminateRow(rows, 1, 2, 1, 4);

	*c = rows[2][3] / rows[2][2];
	*b = (rows[1][3] - rows[1][2] * (*c)) / rows[1][1];
	*a = (rows[0][3] - rows[0][2] * (*c) - rows[0][1] * (*b)) / rows[0][0];
	return true;
}

int main() {
	double *x_array;
	double *y_array;
	int numPoints;
	int choice;

	printf("Select a regression\n(1) LinReg(ax+b)\n(2) QuadReg(ax^2+bx+c)\n(3) ExpReg (A*r^x)\n");
	scanf("%d", &choice);
	if (choice < 1 || choice > 3) {
		printf("Invalid selection");
		exit(1);
	}
	
	printf("Enter the number of points to use: \n");
	scanf("%d", &numPoints);
	if (numPoints <= 0) {
		printf("Invalid input entered: %d\n", numPoints);
		exit(1);
	}

	x_array = (double *)malloc(sizeof(double)*numPoints);
	y_array = (double *)malloc(sizeof(double)*numPoints);

	if (x_array == NULL || y_array == NULL) {
		printf("malloc of size %d failed\n", numPoints);
		exit(1);
	}

	for (int i = 0; i < numPoints; i++) {
		printf("x[%d]= ?\n", i);
		scanf("%lf", &x_array[i]);
		printf("y[%d]= ?\n", i);
		scanf("%lf", &y_array[i]);
	}

	double a, b, c, r;
	bool result = false;
	switch (choice) {
		case 1:
			result = linReg(numPoints, x_array, y_array, &a, &b, &r);
			if (result) {
				printf("y = %fx + %f\n", a, b);
			} else {
				printf("Not successful!");
			}
			break;
		case 2:
			result = quadReg(numPoints, x_array, y_array, &a, &b, &c, &r);
			if (result) {
				printf("y = %fx^2 + %fx + %f", a, b, c);
			} else {
				printf("Not successful!");
			}
			break;
		case 3:
			result = expReg(numPoints, x_array, y_array, &a, &b);
			if (result) {
				printf("y = %f * %f^x", a, b);
			} else {
				printf("Not successful!");
			}
			break;
	}

	free(x_array);
	free(y_array);
	return !result;
}
