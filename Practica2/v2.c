#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "counter.h"

#define ALIGNMENT 64

int n=0;

void Jacobi(float **a, float *b, float *x, float tol, int max_iter){
    /*
        Variables auxiliares:
        x_new: Vector nueva solución (float[n])
     */
        double ck=0.0;
        int iter;
        float *x_new=(float *)aligned_alloc(ALIGNMENT,n*sizeof(float));
        float norm2=0.0;
        //!iniciamos el contador de ciclos
        start_counter();
        for(iter=0;iter<max_iter;iter++){
            norm2=0.0; //!norma del vector al cuadrado (float) 
            for (int i = 0; i < n; i += 2) {

                float sigma1 = 0.0, sigma2 = 0.0;
    
                // Fusión de bucles internos
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
                for (; j < n; j++) {
                    if (i != j) {
                        sigma1 += a[i][j] * x[j];
                    }
                    if (i + 1 < n && i + 1 != j) {
                        sigma2 += a[i + 1][j] * x[j];
                    }
                }
    
                x_new[i] = (b[i] - sigma1) / a[i][i];
                if (i + 1 < n) {
                    x_new[i + 1] = (b[i + 1] - sigma2) / a[i + 1][i + 1];
                }
    
                // Sumar el cuadrado de las diferencias
                float diff1 = x_new[i] - x[i];
                norm2 += diff1 * diff1;
                if (i + 1 < n) {
                    float diff2 = x_new[i + 1] - x[i + 1];
                    norm2 += diff2 * diff2;
                }
            }
    
            // Actualizar el vector x
            float* temp = x;
            x = x_new;
            x_new = temp;
    
            // Verificar convergencia
            if (sqrtf(norm2) < tol) {
                break;
            }
        }    
        //! Detenemos el contador de ciclos
        ck = get_counter();
    
        printf("Iteracion %d\nNorma= %.5lf\n", iter, sqrtf(norm2));
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

    n=atoi(argv[1]);

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
            
