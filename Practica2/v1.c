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

//Definimos como variable global el tamaño de la matriz
int n = 0;

//Función que implementa el método de Jacobi
void v1Jacobi(float** a, float* b, float* x, float tol, int max_iter) {
    //Declaramos las variables necesarias
    double ck = 0.0;
    float* x_new = (float*)aligned_alloc(64, n * sizeof(float));
    int iter = 0;
    float norm2 = 0.0;

    //Implementamos el pseudocódigo del método de Jacobi e iniciamos el contador
    start_counter();
    for(iter = 0; iter < max_iter; iter++) {
        norm2 = 0.0;

        //Iteramos sobre cada fila de la matriz
        for(int i = 0; i < n; i++) {
            float sigma = 0.0;

            //Calculamos la suma de los productos de los elementos de la fila i con los elementos del vector solución
            for(int j = 0; j < n; j++) {
                if(i != j) {
                    sigma += a[i][j] * x[j];
                }
            }

            //Calculamos el nuevo valor del elemento i del vector solución
            x_new[i] = (b[i] - sigma) / a[i][i];
            //Calculamos la norma del vector solución
            norm2 += (x_new[i] - x[i]) * (x_new[i] - x[i]);
        }

        //Actualizamos el vector solución
        for(int i = 0; i < n; i++) {
            x[i] = x_new[i];
        }

        //Comprobamos si la norma es menor que la tolerancia
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

    //Comprobamos que se ha introducido el tamaño de la matriz. También se pudo introducir el número de hilos, pero no lo necesitamos.
    if(argc != 2 && argc != 3) {
        printf("Error: se debe introducir el tamaño de la matriz como argumento.\n");
        printf("Uso: %s <tamaño de la matriz>\n", argv[0]);
        return EXIT_FAILURE;
    }

    //Recuperamos el tamaño de la matriz
    n = atoi(argv[1]);
    if(n <= 0) {
        printf("El tamaño de la matriz debe ser mayor que 0.\n");
        return EXIT_FAILURE;
    }

    //Reservamos memoria para la matriz
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
        //Liberamos la memoria reservada para la matriz y el vector de términos independientes
        for(int i = 0; i < n; i++) {
            free(a[i]);
        }
        free(a);
        free(b);
        return EXIT_FAILURE;
    }

    //Establecemos la semilla para la generación de números aleatorios
    srand(n);

    //Inicializamos la matriz de coeficientes, sumando a la diagonal el sumatorio de cada fila, haciendo que la matriz sea diagonal dominante
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
    v1Jacobi(a, b, x, TOL, MAX_ITER);

    //Liberamos la memoria reservada de las filas de la matriz
    for(int i = 0; i < n; i++) {
        free(a[i]);
    }
    //Liberamos la memoria reservada de la matriz, el vector de términos independientes y el vector solución
    free(a);
    free(b);
    free(x);

    return EXIT_SUCCESS;
}   

