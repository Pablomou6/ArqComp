#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "counter.h"

#define ALIGNMENT 64

//Declaramos el tamaño de la matriz
int n = 0;

/*
    Entradas:
    a: Matriz de coeficientes del sistema (float[n x n])
    b: Vector de términos independientes (float[n])
    x: Vector solución (float[n])
    tol: Tolerancia para la convergencia (float)
    max_iter: Número máximo de iteraciones (int)
*/

void Jacobi(float **a, float *b, float *x, float tol, int max_iter){
        double ck = 0.0;
        int iter = 0;
        float* x_new = (float*)aligned_alloc(ALIGNMENT, n * sizeof(float));
        float norm2 = 0.0;

        //Iniciamos el contador de ciclos
        start_counter();
        for(iter=0;iter<max_iter;iter++){
            norm2=0.0; 

            //Respecto a v1, desarollamos el bucle de las filas de la matriz. Nos añade eficiencia.
            for (int i = 0; i < n; i += 2) {

                float sigma1 = 0.0, sigma2 = 0.0;
    
                /*
                    A su vez, el bucle de las columnas lo desarrollamos de 4 en 4, permitiendo procesar 4 elementos en una iteración. Como se hace en
                    dos filas simultáneamente, procesamos 8 elementos en una iteración.
                    Decidimos comprobar en todo momento que no sumásemos la diagonal ya que, al procesar dos filas a la vez, la diagonal varía y esto
                    supondría un aumento en el número de bucles si dividiésemos el bucle en dos.
                */
                int j = 0;
                for (; j + 3 < n; j += 4) {
                    if (i != j) {
                        sigma1 += a[i][j] * x[j];
                    }
                    if (i + 1 < n && i + 1 != j) {
                        sigma2 += a[i + 1][j] * x[j];
                    }
                    if (i != j + 1) {
                        sigma1 += a[i][j + 1] * x[j + 1];
                    }
                    if (i + 1 < n && i + 1 != j + 1) {
                        sigma2 += a[i + 1][j + 1] * x[j + 1];
                    }
                    if (i != j + 2) {
                        sigma1 += a[i][j + 2] * x[j + 2];
                    }
                    if (i + 1 < n && i + 1 != j + 2) {
                        sigma2 += a[i + 1][j + 2] * x[j + 2];
                    }
                    if (i != j + 3) {
                        sigma1 += a[i][j + 3] * x[j + 3];
                    }
                    if (i + 1 < n && i + 1 != j + 3) {
                        sigma2 += a[i + 1][j + 3] * x[j + 3];
                    }
                }
                //Dado que el programa debe funcionar para todo n, debemos tener en cuenta el caso de que n no sea múltiplo de 4. Procesamos los restantes
                for (; j < n; j++) {
                    if (i != j) {
                        sigma1 += a[i][j] * x[j];
                    }
                    if (i + 1 < n && i + 1 != j) {
                        sigma2 += a[i + 1][j] * x[j];
                    }
                }

                /*
                    Antes de llegar a la conclusión de que desenrollar el bucle por filas nos reducía ciclos, probamos también los bloques
                    con un blocking_size de 64, tamaño el cual calculamos que debería ser el mejor a priori. Es importante esta afirmación
                    previa ya que el código comentado del blocking tiene el bucle de las columnas desenrollado y dividido en 2 para 
                    evitar comprobar la diagonal, pero como ya mencionamos antes, al desenrollar el bucle por filas se nos presentaba el problema de
                    tener 2 diagonales diferentes, siendo mejor desenrollar el bucle en menos elementos por fila (que al tener dos filas simultáneamente
                    conseguíamos el mismo resultado) y comprobar con if las diagonales. El código de los bloques es el siguiente:

                    Iteramos sobre bloques de filas
                    for (int bi = 0; bi < n; bi += BLOCK_SIZE) {
                        int bi_end = (bi + BLOCK_SIZE > n) ? n : bi + BLOCK_SIZE;

                        Iteramos sobre cada fila dentro del bloque
                        for (int i = bi; i < bi_end; i++) {
                            float sigma = 0.0;

                            Suma de los elementos antes de la diagonal (bloques de columnas)
                            for (int bj = 0; bj < i; bj += BLOCK_SIZE) {
                                int bj_end = (bj + BLOCK_SIZE > i) ? i : bj + BLOCK_SIZE;

                                for (int j = bj; j < bj_end; j++) {
                                    sigma += a[i][j] * x[j];
                                }
                            }

                            Suma de los elementos después de la diagonal (bloques de columnas)
                            for (int bj = i + 1; bj < n; bj += BLOCK_SIZE) {
                                int bj_end = (bj + BLOCK_SIZE > n) ? n : bj + BLOCK_SIZE;

                                for (int j = bj; j < bj_end; j++) {
                                    sigma += a[i][j] * x[j];
                                }
                            }
                            (calculos del vector solución y diferencia)
                        }
                    }
                
                */
                
                //Calculamos el nuevo valor del elemento i del vector solución
                x_new[i] = (b[i] - sigma1) / a[i][i];
                if (i + 1 < n) {
                    x_new[i + 1] = (b[i + 1] - sigma2) / a[i + 1][i + 1];
                }
    
                //Sumamos el cuadrado de las diferencias. Respecto a v1, decidimos almacenar la diferencia en una variable para así solo calcularla una
                //vez y no dos.
                float diff1 = x_new[i] - x[i];
                norm2 += diff1 * diff1;
                if (i + 1 < n) {
                    float diff2 = x_new[i + 1] - x[i + 1];
                    norm2 += diff2 * diff2;
                }
            }
    
            //Actualizamos el vector solución. Respecto a v1, lo hacemos con intercambio de punteros.
            float* temp = x;
            x = x_new;
            x_new = temp;
    
            //Verificamos la convergencia
            if (sqrtf(norm2) < tol) {
                break;
            }
        }    
        //Detenemos el contador de ciclos
        ck = get_counter();
    
        printf("Iteracion %d\nNorma= %.5lf\n", iter, sqrtf(norm2));
        printf("Ciclos totales = %.2lf\n", ck);
    
        free(x_new);
}

int main(int argc, char *argv[]){
    //Declaramos las variables que usaremos
    float **a=NULL; //Matriz de coeficientes
    float *b=NULL; //Vector de terminos independientes
    float *x=NULL; //Vector solución 
    float tol = 1e-8; //Tolerancia de 1e-8
    int max_iter = 20000; //20000 iteracciones máximas

    //Comprobamos que se ha introducido el tamaño de la matriz. También se pudo introducir el número de hilos, pero no lo necesitamos.
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

    //Reservamos memoria para la matriz a
    a=(float **)aligned_alloc(ALIGNMENT, n*sizeof(float *));
    if(a==NULL){
        printf("Error al reservar memoria para la matriz de coeficientes a\n");
        return EXIT_FAILURE;
    }
    //Reservamos memoria para cada una de las filas de la matriz a
    for(int i=0; i<n ; i++){
        a[i]=(float *)aligned_alloc(ALIGNMENT, n*sizeof(float));
        if(a[i]==NULL){
            printf("Error al reservar memoria para las filas de la matriz a\n");
            //Para que en caso de que se reservase memoria hasta el momento, se libere debido al fallo
            for(int j=0;j<i;j++){
                free(a[j]);
            }
            free(a);
            return EXIT_FAILURE;
        }
    }

    //Reservamos memoria para el vector b
    b=(float*)aligned_alloc(ALIGNMENT, n*sizeof(float));
    if(b==NULL){
        printf("Error al reservar memoria para el vector de terminos independientes b\n");
        for(int i=0;i<n;i++){
            free(a[i]);
        }
        free(a);
        return EXIT_FAILURE;
    }

    //Reservamos memoria para el vector solución x
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

    //Inicializamos la semilla para la generación de números aleatorios
    srand(n); 

    //Inicializamos la variable a como se indica en el pdf
    for(int i=0;i<n;i++){
        float sum_acum = 0.0;    
        for(int j=0;j<n;j++){
            a[i][j]=(float) rand() / RAND_MAX;
            //Acumulamos el valor de la fila para luego sumarlo a la diagonal
            sum_acum+=a[i][j];
        }
        //Sumamos a cada valor de la diagonal el valor de la suma de los elementos de esa fila
        a[i][i]+=sum_acum;
    }

    //Inicializamos el vector b con valores aleatorios
    for(int i=0;i<n;i++){
        b[i]=(float) rand() / RAND_MAX;
    }

    //Inicializamos el vector solución x a 0
    for(int i=0;i<n;i++){
        x[i]=0.0;
    }

    //Llamamos al método Jacobi
    Jacobi(a,b,x,tol,max_iter);

    //Liberamos la memoria reservada y finalizamos el programa
    for(int i=0;i<n;i++){
        free(a[i]);
    }
    free(a);
    free(b);
    free(x);
    printf("Fin del programa\n");

    return 0;
}
            
