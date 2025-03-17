#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "counter.h"
#include <time.h>

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

int n = 200;

void originalJacobi(float a[n][n], float b[n], float x[n], float tol, int max_iter) {
    float *x_new = (float*)malloc(n*sizeof(float));
    for(int iter = 0; iter < max_iter; iter++) {
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
        x = x_new;
        if(sqrt(norm2) < tol) {
            break;
        }
    }
}

int main(int argc, char *argv[]) {
    srand(time(NULL));

    double ck = 0;
    float a[n][n];
    float b[n];
    float x[n];
    float tol = 1e-6;
    int max_iter = 20000;

    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            a[i][j] = (float)rand() / RAND_MAX;
        }
        b[i] = (float)rand() / RAND_MAX;
        x[i] = 0;
    }
    
    start_counter();
    originalJacobi(a, b, x, tol, max_iter);
    ck = get_counter();

    printf("Ciclos: %.0f\n", ck);
    return 0;   
}