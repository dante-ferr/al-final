#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define AT(row, col, num_cols) ((row) * (num_cols) + (col))
#define EPSILON 1e-9

typedef double real;

int ker(double *matent, real *matesc, real *base, int nlin, int ncol);
void print_matrix(real *mat, int rows, int cols);

int ker(double *matent, real *matesc, real *base, int nlin, int ncol) {
    for (int i = 0; i < nlin * ncol; i++) {
        matesc[i] = matent[i];
    }

    int pivot_row = 0;
    int *pivot_cols = (int*)calloc(ncol, sizeof(int));
    int rank = 0;

    for (int col = 0; col < ncol && pivot_row < nlin; col++) {
        int sel = pivot_row;
        for (int i = pivot_row + 1; i < nlin; i++) {
            if (fabs(matesc[AT(i, col, ncol)]) > fabs(matesc[AT(sel, col, ncol)])) {
                sel = i;
            }
        }

        if (fabs(matesc[AT(sel, col, ncol)]) < EPSILON) {
            continue;
        }

        if (sel != pivot_row) {
            for (int k = 0; k < ncol; k++) {
                real temp = matesc[AT(sel, k, ncol)];
                matesc[AT(sel, k, ncol)] = matesc[AT(pivot_row, k, ncol)];
                matesc[AT(pivot_row, k, ncol)] = temp;
            }
        }

        real pivot_val = matesc[AT(pivot_row, col, ncol)];
        for (int k = 0; k < ncol; k++) {
            matesc[AT(pivot_row, k, ncol)] /= pivot_val;
        }

        for (int i = 0; i < nlin; i++) {
            if (i != pivot_row) {
                real factor = matesc[AT(i, col, ncol)];
                for (int k = 0; k < ncol; k++) {
                    matesc[AT(i, k, ncol)] -= factor * matesc[AT(pivot_row, k, ncol)];
                }
            }
        }

        pivot_cols[col] = 1;
        pivot_row++;
        rank++;
    }

    int nullity = ncol - rank;

    if (nullity == 0) {
        free(pivot_cols);
        return 0;
    }

    int current_base_vec = 0;
    
    for (int col = 0; col < ncol; col++) {
        if (!pivot_cols[col]) {
            for (int row = 0; row < ncol; row++) {
                int base_idx = AT(current_base_vec, row, ncol);
                
                if (row == col) {
                    base[base_idx] = 1.0;
                } else if (!pivot_cols[row]) {
                    base[base_idx] = 0.0;
                } else {
                    int r_pivot = -1;
                    for(int r = 0; r < nlin; r++) {
                        if(fabs(matesc[AT(r, row, ncol)] - 1.0) < EPSILON) {
                             r_pivot = r; 
                             break;
                        }
                    }
                    
                    if (r_pivot != -1) {
                         base[base_idx] = -matesc[AT(r_pivot, col, ncol)];
                    } else {
                         base[base_idx] = 0.0;
                    }
                }
            }
            current_base_vec++;
        }
    }

    free(pivot_cols);
    return nullity;
}

void print_matrix(real *mat, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        printf("| ");
        for (int j = 0; j < cols; j++) {
            double val = mat[AT(i, j, cols)];
            if(fabs(val) < EPSILON) val = 0.0;
            printf("%8.3f ", val);
        }
        printf("|\n");
    }
    printf("\n");
}

int main() {
    int nlin, ncol;

    printf("Dimensões (linhas colunas): ");
    if (scanf("%d %d", &nlin, &ncol) != 2) return 1;

    double *matent = (double*)malloc(nlin * ncol * sizeof(double));
    real *matesc = (real*)malloc(nlin * ncol * sizeof(real));
    real *base = (real*)malloc(ncol * ncol * sizeof(real));

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