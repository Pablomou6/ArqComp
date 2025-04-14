#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "counter.h"
#include <immintrin.h>

//Declaramos constantes que usaremos
#define TOL 1e-8
#define MAX_ITER 20000
#define ALIGNMENT 32

//Definimos como variable global el tamaño de la matriz
int n = 0;

//Función que implementa el método de Jacobi
void v3Jacobi(float** a, float* b, float* x, float tol, int max_iter) {
    double ck = 0.0;
    float* x_new = (float*)aligned_alloc(ALIGNMENT, n * sizeof(float));
    int iter = 0;
    float norm2 = 0.0;

    start_counter();

    for (iter = 0; iter < max_iter; iter++) {
        norm2 = 0.0;

        for (int i = 0; i < n; i++) {
            float sigma = 0.0f;

            int j;
            for (j = 0; j <= n - 8; j += 8) {
                __m256 va = _mm256_loadu_ps(&a[i][j]);
                __m256 vx = _mm256_loadu_ps(&x[j]);

                // Cargar máscara manualmente con if en bucle
                float a_vals[8], x_vals[8];
                _mm256_storeu_ps(a_vals, va);
                _mm256_storeu_ps(x_vals, vx);

                for (int k = 0; k < 8; k++) {
                    if (j + k != i) sigma += a_vals[k] * x_vals[k];
                }
            }

            for (; j < n; j++) {
                if (j != i) sigma += a[i][j] * x[j];
            }

            x_new[i] = (b[i] - sigma) / a[i][i];
            float diff = x_new[i] - x[i];
            norm2 += diff * diff;
        }

        for (int i = 0; i < n; i++) {
            x[i] = x_new[i];
        }

        if (sqrtf(norm2) < tol) break;
    }

    ck = get_counter();
    printf("Ciclos: %.2lf\n", ck);
    printf("Iteraciones: %d\n", iter);
    printf("Norma: %f\n", sqrtf(norm2));

    free(x_new);
}

int main(int argc, char* argv[]) {
    //Declaramos las variables que usaremos
    float** a = NULL;
    float* b = NULL;
    float* x = NULL;

    //Comprobamos que se ha introducido el tamaño de la matriz
    if(argc != 2) {
        printf("Error: se debe introducir el tamaño de la matriz como argumento.\n");
        printf("Uso: %s <tamaño de la matriz>\n", argv[0]);
        return EXIT_FAILURE;
    }

    //Recuperamos el tamaño de la matriz
    n = atoi(argv[1]);

    //Reservamos memoria para la matriz (la almacenamos en un vector plano)
    a = (float**)aligned_alloc(ALIGNMENT, n * sizeof(float*));
    if(a == NULL) {
        printf("Error: no se ha podido reservar memoria para la matriz de coeficientes.\n");
        return EXIT_FAILURE;
    }
    //Se reserva memoria para cada fila de la matriz
    for(int i = 0; i < n; i++) {
        a[i] = (float*)aligned_alloc(ALIGNMENT, n * sizeof(float));
        if(a[i] == NULL) {
            printf("Error: no se ha podido reservar memoria para la fila %d de la matriz de coeficientes.\n", i);
            //Liberamos la memoria reservada hasta el momento
            for(int j = 0; j < i; j++) {
                free(a[j]);
            }
            free(a);
            return EXIT_FAILURE;
        }
    }

    //Se reserva memoria para el vector de términos independientes
    b = (float*)aligned_alloc(ALIGNMENT, n * sizeof(float));
    if(b == NULL) {
        printf("Error: no se ha podido reservar memoria para el vector de términos independientes.\n");
        //Liberamos la memoria reservada para la matriz
        for(int i = 0; i < n; i++) {
            free(a[i]);
        }
        free(a);
        return EXIT_FAILURE;
    }

    //Se reserva memoria para el vector solución
    x = (float*)aligned_alloc(ALIGNMENT, n * sizeof(float));
    if(x == NULL) {
        printf("Error: no se ha podido reservar memoria para el vector solución.\n");
        //Liberamos la memoria reservada para el vector de términos independientes
        for(int i = 0; i < n; i++) {
            free(a[i]);
        }
        free(a);
        free(b);
        return EXIT_FAILURE;
    }

    //Inicializamos la semilla para la generación de números aleatorios
    srand(n);

    //Inicializamos la matriz
    for(int i = 0; i < n; i++) {
        float row_sum = 0.0;
        for(int j = 0; j < n; j++) {
            a[i][j] = (float)rand() / RAND_MAX;
            row_sum += a[i][j];
        }
        a[i][i] += row_sum;
    }

    //Inicializamos el vector de términos independientes
    for(int i = 0; i < n; i++) {
        b[i] = (float)rand() / RAND_MAX;
    }

    //Inicializamos el vector solución
    for(int i = 0; i < n; i++) {
        x[i] = 0.0;
    }

    //Llamamos a la función que implementa el método de Jacobi
    v3Jacobi(a, b, x, TOL, MAX_ITER);

    //Liberamos la memoria reservada para la matriz y los vectores
    for(int i = 0; i < n; i++) {
        free(a[i]);
    }
    free(a);
    free(b);
    free(x);

    return EXIT_SUCCESS;
}