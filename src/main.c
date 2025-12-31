#include <stdio.h>
#include <stdlib.h>
#include "kernel.h"

int main() {
    int nlin, ncol;

    printf("Dimensões (linhas colunas): ");
    if (scanf("%d %d", &nlin, &ncol) != 2) return 1;

    double *matent = malloc(nlin * ncol * sizeof(double));
    real *matesc = malloc(nlin * ncol * sizeof(real));
    real *base = malloc(ncol * ncol * sizeof(real));

    if (!matent || !matesc || !base) return 1;

    printf("Matriz (%d elementos): ", nlin * ncol);
    for (int i = 0; i < nlin; i++) {
        for (int j = 0; j < ncol; j++) {
            scanf("%lf", &matent[i * ncol + j]);
        }
    }

    printf("\n[Entrada]\n");
    print_matrix(matent, nlin, ncol);

    int dim = ker(matent, matesc, base, nlin, ncol);

    printf("[RREF]\n");
    print_matrix(matesc, nlin, ncol);

    printf("Dimensão Kernel: %d\n", dim);

    if (dim > 0) {
        printf("[Base]\n");
        print_matrix(base, dim, ncol);
    } else {
        printf("Kernel trivial.\n");
    }

    free(matent);
    free(matesc);
    free(base);

    return 0;
}