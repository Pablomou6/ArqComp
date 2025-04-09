#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "counter.h"

/*
    Entradas:
        a: Matriz de coeficientes del sistema (float[n x n])
        b: Vector de términos independientes (float[n])
        x: Vector solución (float[n])
        tol: Tolerancia para la convergencia (float)
        max_iter: Número máximo de iteraciones (int)

    Variables auxiliares:
        x_new: Vector nueva solución (float[n])
*/

int n;

void originalJacobi(float** a, float* b, float* x, float tol, int max_iter) {
    int iter;
    float *x_new = (float*)aligned_alloc(64, n*sizeof(float));
    for(iter = 0; iter < max_iter; iter++) {
        float norm2 = 0;
        for(int i = 0; i < n; i++) {
            float sigma = 0.0;
            for(int j = 0; j < n; j++) {
                if(i != j) {
                    sigma += a[i][j] * x[j];
                }
            }
            x_new[i] = (b[i] - sigma) / a[i][i];
            norm2 += (x_new[i] - x[i]) * (x_new[i] - x[i]);
        }
        // Copiar los valores de x_new a x
        for (int i = 0; i < n; i++) {
            x[i] = x_new[i];
        }
        if(sqrt(norm2) < tol) {
            break;
        }
    }
    free(x_new);
    printf("Iteraciones: %d\n", iter);
}

int main(int argc, char *argv[]) {
    srand(time(NULL));

    n = atoi(argv[1]);

    double ck = 0;
    float** a = (float**)aligned_alloc(64, n * sizeof(float*));
    for (int i = 0; i < n; i++) {
        a[i] = (float*)aligned_alloc(64, n * sizeof(float));
    }
    float* x = (float*)aligned_alloc(64, n*sizeof(float));
    float* b = (float*)aligned_alloc(64, n*sizeof(float));
    float tol = 1e-8;
    int max_iter = 20000;

    for(int i = 0; i < n; i++) {
        float row_sum = 0.0;
        for(int j = 0; j < n; j++) {
            a[i][j] = (float)rand() / RAND_MAX;
            row_sum += a[i][j];
        }
        a[i][i] += row_sum;
        b[i] = (float)rand() / RAND_MAX;
        x[i] = 0.0;
    }
    
    start_counter();
    originalJacobi(a, b, x, tol, max_iter);
    ck = get_counter();

    printf("Ciclos: %.0f\n", ck);
    return 0;   
}