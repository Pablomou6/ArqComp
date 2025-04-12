#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "counter.h"

/*
    Cambios respecto v1:
    - Se calcula el inverso de la diagonal (a[i][i]) al principio de cada iteración. Evitando una operación compleja.
    - Se calcula sigma en un bucle desenrollado, teniendo en cuenta que la diagonal no se tiene que sumar. Reduce iteraciones y es más eficiente.
    - Se calcula la normal en el mismo bucle que se copia x_new a x; mejorando la localidad de la caché.
    - Almacenamos la difrencia en una variable. De esta forma, la calculamos 1 vez y accedemos a ella 2 veces, en vez de calcularla 2 veces.  
    - La matriz se almacena en un vector plano (tamaño n*n), de forma que es más eficiente acceder a ella. Además, este cambio nos permite un 
    acceso secuencial a la memoria, mejorando la localidad de la caché.
    NOTA: Dado que la matriz la almacenamos como un vector plano, las filas se almacenan de forma contigua. Para acceder a [i][j] se accede a a[i*n+j].
    Por ejemplo, par acceder a la posición real [5][5] para n = 7, se accede realmente a [4][4] en la matriz almacenada. Por lo que, el elemento 
    será el que está tras 4 filas enteras (desde 0 a 3) (4 filas de 7 elementos, i * n) y 4 elementos (4 elementos de la quinta fila, + j). 
    Por lo que se accede a a[4*7+4] = a[32]. Además, como hemos mencionado, nos permite un acceso secuencia a la memoria, ya que dada una fila, por 
    ejemplo i = 1, se accede a n + j, siendo j = {0,1,2,3,...,n-1}
     
*/

//Declaramos constantes que usaremos
#define TOL 1e-8
#define MAX_ITER 20000
#define BLOCK_SIZE 64

//Definimos como variable global el tamaño de la matriz
int n = 0;

//Función que implementa el método de Jacobi
//! Función de v1, sin optimizar. Probar diferentes cousas.
void v2Jacobi(float** a, float* b, float* x, float tol, int max_iter) {
    //Declaramos las variables necesarias
    double ck = 0.0;
    float* x_new = (float*)aligned_alloc(64, n * sizeof(float));
    int iter = 0;
    float norm2 = 0.0, sigma = 0.0;

    //Implementamos el pseudocódigo del método de Jacobi e iniciamos el contador
    start_counter();
    for(iter = 0; iter < max_iter; iter++) {
        norm2 = 0.0;

        //Iteramos sobre cada fila de la matriz
        for(int i = 0; i < n; i++) {
            sigma = 0.0;

            //! Calculamos la suma de los productos de los elementos de la fila i con los elementos del vector solución
            // Suma de los elementos antes de la diagonal (desenrollado de 4 en 4)
            int j;
            for (j = 0; j <= i - 4; j += 4) {
                sigma += a[i][j] * x[j];
                sigma += a[i][j + 1] * x[j + 1];
                sigma += a[i][j + 2] * x[j + 2];
                sigma += a[i][j + 3] * x[j + 3];
            }
            // Procesar los elementos restantes
            for (; j < i; j++) {
                sigma += a[i][j] * x[j];
            }

            // Suma de los elementos después de la diagonal (desenrollado de 4 en 4)
            for (j = i + 1; j <= n - 4; j += 4) {
                sigma += a[i][j] * x[j];
                sigma += a[i][j + 1] * x[j + 1];
                sigma += a[i][j + 2] * x[j + 2];
                sigma += a[i][j + 3] * x[j + 3];
            }
            // Procesar los elementos restantes
            for (; j < n; j++) {
                sigma += a[i][j] * x[j];
            }

            //Calculamos el nuevo valor del elemento i del vector solución
            x_new[i] = (b[i] - sigma) / a[i][i];

            //! Calculamos la norma del vector solución
            float diff = x_new[i] - x[i];
            norm2 += diff * diff;
        }

        //! Actualizamos el vector solución con los nuevos valores
        float* temp = x;
        x = x_new;
        x_new = temp;

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
    v2Jacobi(a, b, x, TOL, MAX_ITER);

    //Liberamos la memoria reservada para la matriz y los vectores
    for(int i = 0; i < n; i++) {
        free(a[i]);
    }
    free(a);
    free(b);
    free(x);

    return EXIT_SUCCESS;
}