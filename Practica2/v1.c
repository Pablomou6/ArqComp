/*
    Entradas:
        a: Matriz de coeficientes del sistema (float[n x n])
        b: Vector de términos independientes (float[n])
        x: Vector solución (float[n])
        tol: Tolerancia para la convergencia (float)
        max_iter: Número máximo de iteraciones (int)

    Variables auxiliares:
        x_new: Vector nueva solución (float[n])

    Cómputo:
        Para iter (int) desde 0 hasta iter:
            norm2 = 0: norma del vector al cuadrado (float)
            Para i (int) desde 0 hasta n:
                sigma = 0.0 (float)
                Para j desde 0 hasta n:
                    Si i ≠ j:
                        sigma += a[i][j] * x[j]
                x_new[i] = (b[i] - sigma) / a[i][i]
                norm2 += (x_new[i] - x[i]) ^ 2
            x = x_new
            Si srqt(norm2) < tol:
                Terminar

    Salida:
        Imprimir valor de norm2
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int n = 0;

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