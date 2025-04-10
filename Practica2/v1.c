#include <stdio.h>
#include <stdlib.h>
#include <math.h>
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

//Definimos ciertas constantes
#define TOL 1e-8
#define MAX_ITER 20000

//Definimos variables globales
int n = 0;
int iter = 0;
float norm2 = 0;
float* x_new = NULL;

void originalJacobi(float** a, float* b, float* x, float tol, int max_iter) {
    
    for(iter = 0; iter < max_iter; iter++) {
        norm2 = 0;

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
        
        for (int i = 0; i < n; i++) {
            x[i] = x_new[i];
        }

        if(sqrtf(norm2) < tol) {
            break;
        }
    }
}



int main(int argc, char* argv[]) {

    if(argc != 2) {
        printf("Error de entrada, se debe especificar el tamaño de la matriz como argumento.\n");
        return EXIT_FAILURE;
    }

    //Recuperamos el valor de n introducido por parámetros
    n = atoi(argv[1]);
    
    //Reservamos memoria dinámica para las estructuras de datos que dependen de n
    float** a = (float**)aligned_alloc(64, n * sizeof(float*));
    if (a == NULL) {
        printf("Error al reservar memoria para la matriz.\n");
        return EXIT_FAILURE;
    }
    for (int i = 0; i < n; i++) {
        a[i] = (float*)aligned_alloc(64, n * sizeof(float));
        if (a[i] == NULL) {
            printf("Error al reservar memoria para la fila %d de la matriz.\n", i);
            for (int j = 0; j < i; j++) {
                free(a[j]);
            }
            free(a);
            return EXIT_FAILURE;
        }
    }

    float* b = (float*)aligned_alloc(64, n * sizeof(float));
    if(b == NULL) {
        printf("Error al reservar memoria para el vector de términos independientes.\n");
        for (int i = 0; i < n; i++) {
            free(a[i]);
        }
        free(a);
        return EXIT_FAILURE;
    }

    float* x = (float*)aligned_alloc(64, n * sizeof(float));
    if(x == NULL) {
        printf("Error al reservar memoria para el vector solución.\n");
        for (int i = 0; i < n; i++) {
            free(a[i]);
        }
        free(a);
        free(b);
        return EXIT_FAILURE;
    }
    
    x_new = (float*)aligned_alloc(64, n * sizeof(float));
    if(x_new == NULL) {
        printf("Error al reservar memoria para el vector nueva solución.\n");
        for (int i = 0; i < n; i++) {
            free(a[i]);
        }
        free(a);
        free(b);
        free(x);
        return EXIT_FAILURE;
    }

    //Inicializamos la matriz y los vectores.
    //Declaramos la semilla para  rand() con N.
    srand(n);

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

    //Declaramos la variable que almacenará el número de ciclos
    double ck = 0;

    start_counter();
    originalJacobi(a, b, x, TOL, MAX_ITER);
    ck = get_counter();

    printf("Ciclos: %.0f\n", ck);
    if(iter == MAX_ITER) {
        printf("No se ha alcanzado la convergencia en %d iteraciones.\n", MAX_ITER);
    } else {
        printf("Convergencia alcanzada en %d iteraciones.\n", iter);
    }
    printf("norma2: %.6f\n", sqrt(norm2));

    //Ahora liberamos la memoria reservada
    for (int i = 0; i < n; i++) {
        free(a[i]);
    }
    free(a);
    free(b);
    free(x);
    free(x_new);

    return EXIT_SUCCESS;
}

