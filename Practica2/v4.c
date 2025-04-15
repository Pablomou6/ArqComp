#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "counter.h"
#include <omp.h>
#include <stdint.h>

#define ALIGNMENT 64

int n=0;
int num_threads=0; //a partir de 16, el rendimiento

void Jacobi(float **a, float *b, float *x, float tol, int max_iter){
    /*
        Variables auxiliares:
        x_new: Vector nueva solución (float[n])
     */
    double ck = 0.0;
    int iter;
    float *x_new = (float *)aligned_alloc(ALIGNMENT, n * sizeof(float));
    if (x_new == NULL) {
        printf("Error al reservar memoria para x_new\n");
        return;
    }
    float norm2 = 0.0;

    //! Iniciamos el contador de ciclos
    start_counter();

    for (iter = 0; iter < max_iter; iter++) {
        norm2 = 0.0;

        // Paralelizamos el bucle externo con OpenMP
        #pragma omp parallel
        {
            float local_norm2 = 0.0; // Variable local para la reducción

            #pragma omp for schedule(dynamic) reduction(+:norm2)
            for (int i = 0; i < n; i++) {
                float sigma = 0.0;

                // Paralelizamos el cálculo de sigma con reducción
                #pragma omp parallel for reduction(+:sigma)
                for (int j = 0; j < n; j++) {
                    if (i != j) {
                        sigma += a[i][j] * x[j];
                    }
                }

                // Calculamos el nuevo valor de x[i]
                x_new[i] = (b[i] - sigma) / a[i][i];

                // Calculamos la diferencia para la norma
                float diff = x_new[i] - x[i];
                local_norm2 += diff * diff;
            }

            // Actualizamos la norma global
            #pragma omp atomic
            norm2 += local_norm2;
        }

        // Actualizamos el vector x
        #pragma omp parallel for
        for (int i = 0; i < n; i++) {
            x[i] = x_new[i];
        }

        // Verificamos la convergencia
        if (sqrtf(norm2) < tol) {
            break;
        }
    }

    //! Detenemos el contador de ciclos
    ck = get_counter();

    printf("Iteraciones: %d\nNorma: %.5lf\n", iter, sqrtf(norm2));
    printf("Ciclos totales = %.2lf\n", ck);

    free(x_new);
}

int main(int argc, char *argv[]){
    /*
    Entradas:
    a: Matriz de coeficientes del sistema (float[n x n])
    b: Vector de términos independientes (float[n])
    x: Vector solución (float[n])
    tol: Tolerancia para la convergencia (float)
    max_iter: Número máximo de iteraciones (int)
    */
    float **a=NULL; // Matriz de coeficientes
    float *b=NULL; // vector de terminos independientes
    float *x=NULL; // vector solución 
    float tol = 1e-8; // Tolerancia de 1e-8
    int max_iter = 20000; // 20000 iteracciones máximas

    if(argc != 3) {
        printf("Error: se necesitan 2 argumentos: n y num_threads\n");
        return EXIT_FAILURE;
    }

    n = atoi(argv[1]);
    if (n <= 0) {
        printf("Error: el tamaño de la matriz debe ser un número entero positivo.\n");
        return EXIT_FAILURE;
    }

    num_threads = atoi(argv[2]);
    if (num_threads <= 0) {
        printf("Error: el número de hilos debe ser un número entero positivo.\n");
        return EXIT_FAILURE;
    }
    omp_set_num_threads(num_threads);

    //!reservamos memoria para la matriz a de esta forma
    a=(float **)aligned_alloc(ALIGNMENT, n*sizeof(float *));
    if(a==NULL){
        printf("Error al reservar memoria para la matriz de coeficientes a\n");
        return EXIT_FAILURE;
    }
    //!reservamos memoria para cada una de las filas de la matriz a
    for(int i=0; i<n ; i++){
        a[i]=(float *)aligned_alloc(ALIGNMENT, n*sizeof(float));
        if(a[i]==NULL){
            printf("Error al reservar memoria para las filas de la matriz a\n");
            //!para que en caso de que se reservase memoria hasta el momento, se libere debido al fallo
            for(int j=0;j<i;j++){
                free(a[j]);
            }
            free(a);
            return EXIT_FAILURE;
        }
    }

    //!reservamos memoria para el vector b
    b=(float*)aligned_alloc(ALIGNMENT, n*sizeof(float));
    if(b==NULL){
        printf("Error al reservar memoria para el vector de terminos independientes b\n");
        for(int i=0;i<n;i++){
            free(a[i]);
        }
        free(a);
        return EXIT_FAILURE;
    }

    //!reservamos memoria para el vector solución x
    x=(float*)aligned_alloc(ALIGNMENT, n*sizeof(float));
    if(x==NULL){
        printf("Error al reservar memoria para el vector solución x\n");
        for(int i=0;i<n;i++){
            free(a[i]);
        }
        free(a);
        free(b);
        return EXIT_FAILURE;
    }

    srand(n); //! Inicializamos la semilla para la generación de números aleatorios

    //!inicializamos la variable a como se indica en el pdf
    for(int i=0;i<n;i++){
        float sum_acum = 0.0;    
        for(int j=0;j<n;j++){
            a[i][j]=(float) rand() / RAND_MAX;
            //!hallamos la suma de los elementos de la fila
            sum_acum+=a[i][j];
        }
        //!sumamos a cada valor de la diagonal el valor de la suma de los elementos de esa fila
        a[i][i]+=sum_acum;
    }

    //!inicializamos el vector b con valores aleatorios
    for(int i=0;i<n;i++){
        b[i]=(float) rand() / RAND_MAX;
    }

    //!inicializamos el vector solución x a 0
    for(int i=0;i<n;i++){
        x[i]=0.0;
    }

    //!llamamos al método Jacobi
    Jacobi(a,b,x,tol,max_iter);

    //!liberamos la memoria reservada y finalizamos el programa
    for(int i=0;i<n;i++){
        free(a[i]);
    }
    free(a);
    free(b);
    free(x);
    printf("Fin del programa\n");

    return 0;
}
            
