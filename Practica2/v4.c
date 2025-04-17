#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <omp.h>
#include "counter.h"

#define ALIGNMENT 64

int n=0;

void Jacobi(float **a, float *b, float *x, float tol, int max_iter){
    double ck=0.0;
    int iter;
    float *x_new=(float *)aligned_alloc(ALIGNMENT,n*sizeof(float));
    float norm2=0.0;
    
    //Iniciamos el contador de ciclos
    start_counter();

    for(iter=0;iter<max_iter;iter++){
        //Declaramos variables
        int i, j;
        float sigma1, sigma2, diff1, diff2;

        norm2=0.0;

        /*
            Utilizamos un solo pragma, ya que de esta forma paralelizamos todo el bucle, haciendo que cada hilo procese dos filas a la vez.
            Este pragma nos paraleliza el bucle for de las filas, haciendo que cada bucle procese dos filas a la vez con sus respectvas columnas. Tiene
            un scheduling static, ya que nos proporciona mejor resultado que un scheduling dynamic. La diferencia es que en el scheduling static, cada hilo
            procesa un número fijo de iteraciones, mientras que en el scheduling dynamic, cada hilo procesa un número variable de iteraciones. En este caso,
            el scheduling static es más eficiente, ya que cada hilo procesa un número fijo de iteraciones y no hay sobrecarga de gestión de hilos.
            Además, le declaramos explícitamente las privacidades de las variables.
            Por último, la clausula reduction(+:norm2) nos proporciona un resultado más eficiente, ya que cada hilo calcula su propio valor de norm2 y
            al final se suman todos los valores de norm2 de cada hilo, evitando conflictos entre hilos.
        */
        #pragma omp parallel for schedule(static) shared(a,b,x,x_new,n) private(i,j,diff1,diff2,sigma1,sigma2) reduction(+:norm2)
        for(i=0; i<n; i+=2){
            sigma1=0, sigma2=0;

            //Aprovechamos las mejoras de la versión v2, como el desenrollo de bucles
            for (j = 0; j <= n - 4; j += 4) {
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
            
            //Calculamos, en caso de ser necesario, las iteraciones restantes
            int k = 0;
            for (k = j; k < n; k++) {
                if (i != k) {
                    sigma1 += a[i][k] * x[k];
                }
                if (i + 1 < n && i + 1 != k) {
                    sigma2 += a[i + 1][k] * x[k];
                }
            }

            //Calculamos los nuevos valores de x[i] y x[i+1]
            x_new[i]=(b[i]-sigma1)/a[i][i];
            if(i+1<n){
                x_new[i+1]=(b[i+1]-sigma2)/a[i+1][i+1];
            }
            
            //Calculamos las diferencias para la norma
            diff1=(x_new[i]-x[i]);
            norm2+= diff1 * diff1;
            if(i+1<n){
                diff2=(x_new[i+1]-x[i+1]);
                norm2+= diff2 * diff2;
            }
        }
        
        //Actualizamos el vector x con intercambio de punteros
        float *temp=x;
        x=x_new;
        x_new=temp;

        //Calculamos la condición de parada
        if(sqrtf(norm2)<tol){
            break;
        }
    }    
    //Recuperamos el contador de ciclos
    ck = get_counter();

    printf("Iteracion %d\nNorma= %lf\n", iter, sqrtf(norm2));
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

    //Comprobamos que se han introducido los argumentos necesarios
    if(argc != 3){
        printf("Uso: %s <tamaño de la matriz> <número de hilos>\n", argv[0]);
        return EXIT_FAILURE;
    }

    n=atoi(argv[1]);
    if(n<=0){
        printf("El tamaño de la matriz debe ser mayor que 0.\n");
        return EXIT_FAILURE;
    }

    int num_threads = atoi(argv[2]);
    if (num_threads <= 0) {
        printf("El número de hilos debe ser mayor que 0.\n");
        return EXIT_FAILURE;
    }

    //Configuramos el número de hilos para usar en OpenMP
    omp_set_num_threads(num_threads);

    //Reservamos memoria para la matriz a de esta forma
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
            
