#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "counter.h"
#include <immintrin.h>

//Declaramos constantes que usaremos
#define TOL 1e-8
#define MAX_ITER 20000
#define ALIGNMENT 64

//Definimos como variable global el tamaño de la matriz
int n = 0;

//Función que implementa el método de Jacobi
void v3Jacobi(float** a, float* b, float* x, float tol, int max_iter) {
    double ck = 0.0;
    float* x_new = (float*)_mm_malloc(n * sizeof(float), ALIGNMENT);
    int iter = 0;
    float norm2 = 0.0;

    start_counter();

    for (iter = 0; iter < max_iter; iter++) {
        norm2 = 0.0;

        for (int i = 0; i < n; i += 2) { // Procesamos dos filas por iteración
            float sigma1 = 0.0f, sigma2 = 0.0f;

            int j;
            for (j = 0; j <= n - 8; j += 8) {
                // Cargamos 8 elementos de la fila i
                __m256 va1 = _mm256_loadu_ps(&a[i][j]);
                __m256 vx = _mm256_loadu_ps(&x[j]);

                // Cargamos 8 elementos de la fila i+1 (si existe)
                __m256 va2 = (i + 1 < n) ? _mm256_loadu_ps(&a[i + 1][j]) : _mm256_setzero_ps();

                // Convertimos a arreglos para manejar la máscara
                float a_vals1[8], a_vals2[8], x_vals[8];
                _mm256_storeu_ps(a_vals1, va1);
                _mm256_storeu_ps(a_vals2, va2);
                _mm256_storeu_ps(x_vals, vx);

                for (int k = 0; k < 8; k++) {
                    if (j + k != i) sigma1 += a_vals1[k] * x_vals[k];
                    if (i + 1 < n && j + k != i + 1) sigma2 += a_vals2[k] * x_vals[k];
                }
            }

            // Procesamos los elementos restantes
            for (; j < n; j++) {
                if (j != i) sigma1 += a[i][j] * x[j];
                if (i + 1 < n && j != i + 1) sigma2 += a[i + 1][j] * x[j];
            }

            // Calculamos los nuevos valores de x[i] y x[i+1]
            x_new[i] = (b[i] - sigma1) / a[i][i];
            if (i + 1 < n) {
                x_new[i + 1] = (b[i + 1] - sigma2) / a[i + 1][i + 1];
            }

            // Calculamos las diferencias para la norma
            float diff1 = x_new[i] - x[i];
            norm2 += diff1 * diff1;

            if (i + 1 < n) {
                float diff2 = x_new[i + 1] - x[i + 1];
                norm2 += diff2 * diff2;
            }
        }

        // Actualizamos el vector x
        float *temp = x;
        x = x_new;
        x_new = temp;

        // Verificamos la convergencia
        if (sqrtf(norm2) < tol) break;
    }

    ck = get_counter();
    printf("Ciclos: %.2lf\n", ck);
    printf("Iteraciones: %d\n", iter);
    printf("Norma: %f\n", sqrtf(norm2));

    _mm_free(x_new);
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
    a = (float**)_mm_malloc(n * sizeof(float*), ALIGNMENT);
    if(a == NULL) {
        printf("Error: no se ha podido reservar memoria para la matriz de coeficientes.\n");
        return EXIT_FAILURE;
    }
    //Se reserva memoria para cada fila de la matriz
    for(int i = 0; i < n; i++) {
        a[i] = (float*)_mm_malloc(n * sizeof(float), ALIGNMENT);
        if(a[i] == NULL) {
            printf("Error: no se ha podido reservar memoria para la fila %d de la matriz de coeficientes.\n", i);
            //Liberamos la memoria reservada hasta el momento
            for(int j = 0; j < i; j++) {
                _mm_free(a[j]);
            }
            _mm_free(a);
            return EXIT_FAILURE;
        }
    }

    //Se reserva memoria para el vector de términos independientes
    b = (float*)_mm_malloc(n * sizeof(float), ALIGNMENT);
    if(b == NULL) {
        printf("Error: no se ha podido reservar memoria para el vector de términos independientes.\n");
        //Liberamos la memoria reservada para la matriz
        for(int i = 0; i < n; i++) {
            _mm_free(a[i]);
        }
        _mm_free(a);
        return EXIT_FAILURE;
    }

    //Se reserva memoria para el vector solución
    x = (float*)_mm_malloc(n * sizeof(float), ALIGNMENT);
    if(x == NULL) {
        printf("Error: no se ha podido reservar memoria para el vector solución.\n");
        //Liberamos la memoria reservada para el vector de términos independientes
        for(int i = 0; i < n; i++) {
            _mm_free(a[i]);
        }
        _mm_free(a);
        _mm_free(b);
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
    _mm_free(a);
    _mm_free(b);
    _mm_free(x);

    return EXIT_SUCCESS;
}