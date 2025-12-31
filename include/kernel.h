#ifndef KERNEL_H
#define KERNEL_H

typedef double real;

int ker(double *matent, real *matesc, real *base, int nlin, int ncol);
void print_matrix(real *mat, int rows, int cols);

#endif