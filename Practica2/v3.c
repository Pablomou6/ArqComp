#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "counter.h"
#include <immintrin.h>

//Declaramos constantes que usaremos
#define TOL 1e-8
#define MAX_ITER 20000

//Definimos como variable global el tamaño de la matriz
int n = 0;

//Función que implementa el método de Jacobi
void v3Jacobi(float** a, float* b, float* x, float tol, int max_iter) {
    // Declaramos las variables necesarias
    double ck = 0.0;
    float* x_new = (float*)aligned_alloc(64, n * sizeof(float));
    int iter = 0;
    float norm2 = 0.0;

    // Implementamos el pseudocódigo del método de Jacobi e iniciamos el contador
    start_counter();
    for(iter = 0; iter < max_iter; iter++) {
        norm2 = 0.0;

        // Iteramos sobre cada fila de la matriz
        for(int i = 0; i < n; i++) {
            float sigma = 0.0;
            
            // Preparamos registros vectoriales para acumular resultados parciales
            __m512 sum_vec_512 = _mm512_setzero_ps(); // Para AVX-512 (16 floats)
            __m256 sum_vec_256 = _mm256_setzero_ps(); // Para AVX2 (8 floats)
            __m128 sum_vec_128 = _mm_setzero_ps();    // Para SSE (4 floats)
            
            // Cálculo de sigma usando AVX-512 para j < i
            int j = 0;
            for(; j <= i - 16; j += 16) {
                __m512 a_vec = _mm512_loadu_ps(&a[i][j]);
                __m512 x_vec = _mm512_loadu_ps(&x[j]);
                sum_vec_512 = _mm512_fmadd_ps(a_vec, x_vec, sum_vec_512);
            }
            
            // Procesar elementos restantes con AVX2 (8 floats a la vez)
            for(; j <= i - 8; j += 8) {
                __m256 a_vec = _mm256_loadu_ps(&a[i][j]);
                __m256 x_vec = _mm256_loadu_ps(&x[j]);
                sum_vec_256 = _mm256_fmadd_ps(a_vec, x_vec, sum_vec_256);
            }
            
            // Procesar elementos restantes con SSE (4 floats a la vez)
            for(; j <= i - 4; j += 4) {
                __m128 a_vec = _mm_loadu_ps(&a[i][j]);
                __m128 x_vec = _mm_loadu_ps(&x[j]);
                sum_vec_128 = _mm_add_ps(sum_vec_128, _mm_mul_ps(a_vec, x_vec));
            }
            
            // Procesar elementos restantes de forma escalar
            for(; j < i; j++) {
                sigma += a[i][j] * x[j];
            }
            
            // Cálculo de sigma usando AVX-512 para j > i
            j = i + 1;
            for(; j <= n - 16; j += 16) {
                __m512 a_vec = _mm512_loadu_ps(&a[i][j]);
                __m512 x_vec = _mm512_loadu_ps(&x[j]);
                sum_vec_512 = _mm512_fmadd_ps(a_vec, x_vec, sum_vec_512);
            }
            
            // Procesar elementos restantes con AVX2 (8 floats a la vez)
            for(; j <= n - 8; j += 8) {
                __m256 a_vec = _mm256_loadu_ps(&a[i][j]);
                __m256 x_vec = _mm256_loadu_ps(&x[j]);
                sum_vec_256 = _mm256_fmadd_ps(a_vec, x_vec, sum_vec_256);
            }
            
            // Procesar elementos restantes con SSE (4 floats a la vez)
            for(; j <= n - 4; j += 4) {
                __m128 a_vec = _mm_loadu_ps(&a[i][j]);
                __m128 x_vec = _mm_loadu_ps(&x[j]);
                sum_vec_128 = _mm_add_ps(sum_vec_128, _mm_mul_ps(a_vec, x_vec));
            }
            
            // Procesar elementos restantes de forma escalar
            for(; j < n; j++) {
                sigma += a[i][j] * x[j];
            }
            
            // Reducir los resultados vectoriales a un único valor escalar
            // Reducción AVX-512
            float sum_array_512[16] __attribute__((aligned(64)));
            _mm512_store_ps(sum_array_512, sum_vec_512);
            for(int k = 0; k < 16; k++) {
                sigma += sum_array_512[k];
            }
            
            // Reducción AVX2
            float sum_array_256[8] __attribute__((aligned(32)));
            _mm256_store_ps(sum_array_256, sum_vec_256);
            for(int k = 0; k < 8; k++) {
                sigma += sum_array_256[k];
            }
            
            // Reducción SSE
            float sum_array_128[4] __attribute__((aligned(16)));
            _mm_store_ps(sum_array_128, sum_vec_128);
            for(int k = 0; k < 4; k++) {
                sigma += sum_array_128[k];
            }
            
            // Calculamos el nuevo valor del elemento i del vector solución
            x_new[i] = (b[i] - sigma) / a[i][i];
        }
        
        // Calculamos la norma y actualizamos el vector solución con vectorización
        __m512 norm_vec_512 = _mm512_setzero_ps();
        int i = 0;
        
        // Procesar con AVX-512 (16 floats a la vez)
        for(; i <= n - 16; i += 16) {
            __m512 x_new_vec = _mm512_loadu_ps(&x_new[i]);
            __m512 x_vec = _mm512_loadu_ps(&x[i]);
            __m512 diff_vec = _mm512_sub_ps(x_new_vec, x_vec);
            
            // Actualizar x con x_new
            _mm512_storeu_ps(&x[i], x_new_vec);
            
            // Calcular la contribución a la norma al cuadrado
            norm_vec_512 = _mm512_fmadd_ps(diff_vec, diff_vec, norm_vec_512);
        }
        
        // Procesar con AVX2 (8 floats a la vez)
        __m256 norm_vec_256 = _mm256_setzero_ps();
        for(; i <= n - 8; i += 8) {
            __m256 x_new_vec = _mm256_loadu_ps(&x_new[i]);
            __m256 x_vec = _mm256_loadu_ps(&x[i]);
            __m256 diff_vec = _mm256_sub_ps(x_new_vec, x_vec);
            
            // Actualizar x con x_new
            _mm256_storeu_ps(&x[i], x_new_vec);
            
            // Calcular la contribución a la norma al cuadrado
            norm_vec_256 = _mm256_fmadd_ps(diff_vec, diff_vec, norm_vec_256);
        }
        
        // Procesar con SSE (4 floats a la vez)
        __m128 norm_vec_128 = _mm_setzero_ps();
        for(; i <= n - 4; i += 4) {
            __m128 x_new_vec = _mm_loadu_ps(&x_new[i]);
            __m128 x_vec = _mm_loadu_ps(&x[i]);
            __m128 diff_vec = _mm_sub_ps(x_new_vec, x_vec);
            
            // Actualizar x con x_new
            _mm_storeu_ps(&x[i], x_new_vec);
            
            // Calcular la contribución a la norma al cuadrado
            norm_vec_128 = _mm_add_ps(_mm_mul_ps(diff_vec, diff_vec), norm_vec_128);
        }
        
        // Procesar elementos restantes de forma escalar
        for(; i < n; i++) {
            float diff = x_new[i] - x[i];
            x[i] = x_new[i];
            norm2 += diff * diff;
        }
        
        // Reducir los resultados vectoriales de la norma
        // Reducción AVX-512
        float norm_array_512[16] __attribute__((aligned(64)));
        _mm512_store_ps(norm_array_512, norm_vec_512);
        for(int k = 0; k < 16; k++) {
            norm2 += norm_array_512[k];
        }
        
        // Reducción AVX2
        float norm_array_256[8] __attribute__((aligned(32)));
        _mm256_store_ps(norm_array_256, norm_vec_256);
        for(int k = 0; k < 8; k++) {
            norm2 += norm_array_256[k];
        }
        
        // Reducción SSE
        float norm_array_128[4] __attribute__((aligned(16)));
        _mm_store_ps(norm_array_128, norm_vec_128);
        for(int k = 0; k < 4; k++) {
            norm2 += norm_array_128[k];
        }
        
        // Comprobamos si la norma es menor que la tolerancia
        if(sqrtf(norm2) < tol) {
            break;
        }
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
    a = (float**)aligned_alloc(64, n * sizeof(float*));
    if(a == NULL) {
        printf("Error: no se ha podido reservar memoria para la matriz de coeficientes.\n");
        return EXIT_FAILURE;
    }
    //Se reserva memoria para cada fila de la matriz
    for(int i = 0; i < n; i++) {
        a[i] = (float*)aligned_alloc(64, n * sizeof(float));
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
    b = (float*)aligned_alloc(64, n * sizeof(float));
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
    x = (float*)aligned_alloc(64, n * sizeof(float));
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