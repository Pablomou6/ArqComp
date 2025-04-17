#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "counter.h"
#include <immintrin.h>
#include <stdint.h>

//Declaramos constantes que usaremos
#define TOL 1e-8
#define MAX_ITER 20000
#define ALIGNMENT 64

//Definimos como variable global el tamaño de la matriz
int n = 0;

//Función que implementa el método de Jacobi
void v3Jacobi(float** a, float* b, float* x, float tol, int max_iter) {
    double ck = 0.0;
    
    //Reservamos memoria para el vector solución. Además, comprobamos que se alinease correctamente
    float* x_new = (float*)_mm_malloc(n * sizeof(float), ALIGNMENT);
    if (((uintptr_t)x_new % ALIGNMENT) != 0) {
        printf("Vector b NO alineado: %p\n", (void*)b);
    } 

    int iter = 0;
    float norm2 = 0.0;

    //Iniciamos el contador de ciclos
    start_counter();

    /*
        Reservamos memoria para un vector que almacene los elementos de la diagonal. Además, reservamos memoria para dos vectores temporales (uno por fila
        que se procese simultáneamente) para almacenar los resultados de las multiplicaciones. Comprobamos que se alinean correctamente.
        Como se dice en el guión, se deben incluír en la medida de ciclos lo que implique una sobrecarga al algoritmo secuencial optimizado.
    */
    float* aII = (float*)_mm_malloc(n * sizeof(float), ALIGNMENT);
    if (((uintptr_t)aII % ALIGNMENT) != 0) {
        printf("Vector aII NO alineado: %p\n", (void*)b);
    }

    float* temp1 = (float*)_mm_malloc(8 * sizeof(float), ALIGNMENT);
    if (((uintptr_t)temp1 % ALIGNMENT) != 0) {
        printf("Vector temp1 NO alineado: %p\n", (void*)b);
    }

    float* temp2 = (float*)_mm_malloc(8 * sizeof(float), ALIGNMENT);
    if (((uintptr_t)temp2 % ALIGNMENT) != 0) {
        printf("Vector temp2 NO alineado: %p\n", (void*)b);
    }

    /*
        Como SIMD realiza operaciones vectoriales, haciendo que una sola operación se aplique a varios elementos, controlar la diagonal es más complejo.
        Nuestra solución se encuentra en este bucle, en el cual hacemos una copia de la diagonal de la matriz y la sustituimos por ceros. De esta forma,
        podemos realizar las operacionse de forma vectorial sin preocuparnos de la diagonal (su valor de la operación será 0, ya que es una multiplicación)
        y no afectará al resultado. Cuando necesitemos los valores de la diagonal para el cálculo de la nueva solución, simplemente los recuperamos de aII.
    */
    for (int i = 0; i < n; i++) {
        aII[i] = a[i][i];
        a[i][i] = 0.0f;
    }

    for (iter = 0; iter < max_iter; iter++) {
        norm2 = 0.0;

        //Procesamos dos filas simultáneamente
        for (int i = 0; i < n; i += 2) { 
            float sigma1 = 0.0f, sigma2 = 0.0f;

            int j;
            for (j = 0; j <= n - 8; j += 8) {
                //Cargamos 8 elementos de la fila i
                __m256 va1 = _mm256_load_ps(&a[i][j]);
                __m256 vx = _mm256_load_ps(&x[j]);

                //Cargamos 8 elementos de la fila i+1 en caso de que exista
                __m256 va2 = (i + 1 < n) ? _mm256_load_ps(&a[i + 1][j]) : _mm256_setzero_ps();

                //Realizamos las multiplicaciones 
                __m256 mul1 = _mm256_mul_ps(va1, vx);
                __m256 mul2 = _mm256_mul_ps(va2, vx);

                //Almacenamos en los vectores temporales los resultados de las multiplicaciones
                _mm256_store_ps(temp1, mul1);
                _mm256_store_ps(temp2, mul2);

                //Sumamos los resultados de las multiplicaciones
                for (int k = 0; k < 8; k++) {
                    sigma1 += temp1[k];
                    if (i + 1 < n) sigma2 += temp2[k];
                }
            }

            //Debemos tener en cuenta que puede no ser exacto el bucle anterior. Procesamos los elementos restantes
            for (; j < n; j++) {
                sigma1 += a[i][j] * x[j];
                if (i + 1 < n) sigma2 += a[i + 1][j] * x[j];
            }

            //Calculamos los nuevos valores de x[i] y x[i+1] usando la diagonal almacenada en aII
            x_new[i] = (b[i] - sigma1) / aII[i];
            if (i + 1 < n) {
                x_new[i + 1] = (b[i + 1] - sigma2) / aII[i + 1];
            }

            //Calculamos las diferencias para la norma
            float diff1 = x_new[i] - x[i];
            norm2 += diff1 * diff1;
            if (i + 1 < n) {
                float diff2 = x_new[i + 1] - x[i + 1];
                norm2 += diff2 * diff2;
            }
        }

        //Actualizamos el vector x con intercambio de punteros
        float* temp = x;
        x = x_new;
        x_new = temp;

        //Verificamos la convergencia
        if (sqrtf(norm2) < tol) {
            break;
        }
    }
    ck = get_counter();

    printf("Ciclos: %.2lf\n", ck);
    printf("Iteraciones: %d\n", iter);
    printf("Norma: %lf\n", sqrtf(norm2));

    _mm_free(x_new);
    _mm_free(aII);
    _mm_free(temp1);
    _mm_free(temp2);
}

int main(int argc, char* argv[]) {
    //Declaramos las variables que usaremos
    float** a = NULL;
    float* b = NULL;
    float* x = NULL;

    //Comprobamos que se ha introducido el tamaño de la matriz
    if(argc != 2 && argc != 3) {
        printf("Error: se debe introducir el tamaño de la matriz como argumento.\n");
        printf("Uso: %s <tamaño de la matriz>\n", argv[0]);
        return EXIT_FAILURE;
    }

    //Recuperamos el tamaño de la matriz
    n=atoi(argv[1]);
    if(n<=0){
        printf("El tamaño de la matriz debe ser mayor que 0.\n");
        return EXIT_FAILURE;
    }

    //Reservamos memoria para la matriz
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
    //Comprobamos que la memoria se ha alineado correctamente
    for (int i = 0; i < n; i++) {
        if (((uintptr_t)a[i] % ALIGNMENT) != 0) {
            printf("Fila %d NO alineada: %p\n", i, (void*)a[i]);
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
    //Comprobamos que la memoria se ha alineado correctamente
    if (((uintptr_t)b % ALIGNMENT) != 0) {
        printf("Vector b NO alineado: %p\n", (void*)b);
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
    //Comprobamos que la memoria se ha alineado correctamente
    if (((uintptr_t)x % ALIGNMENT) != 0) {
        printf("Vector b NO alineado: %p\n", (void*)b);
    } 

    //Inicializamos la semilla para la generación de números aleatorios
    srand(n);

    //Inicializamos la matriz. Sumamos a la diagonal el sumatorio de cada fila, haciendo que la matriz sea diagonal dominante
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