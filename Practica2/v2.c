#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int n = 0;

/*
    Cambios respecto v1:
    - Se calcula el inverso de la diagonal (a[i][i]) al principio de cada iteración. Evitando una operación compleja.
    - Se calcula sigma en un bucle desenrrollado, teniendo en cuenta que la diagonal no se tiene que sumar. Reduce iteraciones y es más eficiente.
    - Se calcula la normal en el mismo bucle que se copia x_new a x; mejorando la localidad de la caché.
    - Almacenamos la difrencia en una variable. De esta forma, la calculamos 1 vez y accedemos a ella 2 veces, en vez de calcularla 2 veces. 
    - Cambiamos la forma en que se almacenan las variables. Usamos memoria dinámica de forma que se alineen a la caché. 
    NOTA: Dado que la matriz la almacenamos como un vector plano, las filas se almacenan de forma contigua. Para acceder a [i][j] se accede a a[i*n+j].
    Por ejemplo, par acceder a la posición real [5][5] para n = 7, se accede realmente a [4][4] en la matriz almacenada. Por lo que, el elemento 
    será el que está tras 4 filas enteras (desde 0 a 3) (4 filas de 7 elementos, i * n) y 4 elementos (4 elementos de la quinta fila, + j). 
    Por lo que se accede a a[4*7+4] = a[32]. 
     
*/

void v2Jacobi(float* a, float* b, float* x, float tol, int max_iter) {
    float *x_new = (float*)malloc(n*sizeof(float));
    int iter = 0;

    for(iter = 0; iter < max_iter; iter++) {
        float norm2 = 0;

        //Iteraciones sobre las filas
        for(int i = 0; i < n; i++) {
            float sigma = 0.0;
            //Precalculamos el inverso de la diagonal
            float aInv = 1.0 / a[i * n + i];
            
            int j = 0;  
            for(j = 0; j <= n - 4; j += 4) {
                if(i != j) {
                    sigma += a[i * n + j] * x[j];
                }
                if(i != j + 1) {
                    sigma += a[i * n + j + 1] * x[j + 1];
                }
                if(i != j + 2) {
                    sigma += a[i * n + j + 2] * x[j + 2];
                }
                if(i != j + 3) {
                    sigma += a[i * n + j + 3] * x[j + 3];
                }
            }

            //Ahora, vamos a procesar los últimos elementos
            for(; j < n; j++) {
                if(i != j) {
                    sigma += a[i * n + j] * x[j];
                }
            }
            x_new[i] = (b[i] - sigma) * aInv;
        }
        
        //Fusiono en un bucle la copia de x_new a x con el cálculo de la norma
        for(int i = 0; i < n; i++) {
            //No calculo 2 veces la diferencia. En cambio, la guardo en una variable y accedo a ella dos veces.
            float diff = (x_new[i] - x[i]);
            norm2 += diff * diff;

            x[i] = x_new[i];
        }

        if(sqrt(norm2) < tol) {
            break;
        }

    }

    free(x_new);    
    printf("Iteraciones: %d\n", iter);
}

int main(int argc, char *argv[]) {
    srand(time(NULL));

    double ck = 0;
    float* a = (float*)aligned_alloc(64, n*n*sizeof(float));
    float* b = (float*)aligned_alloc(64, n*sizeof(float));
    float* x = (float*)aligned_alloc(64, n*sizeof(float));
    float tol = 1e-8;
    int max_iter = 20000;

    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            a[i * n + j] = (float)rand() / RAND_MAX;
        }
        b[i] = (float)rand() / RAND_MAX;
        x[i] = 0;
    }
    
    start_counter();
    v2Jacobi(a, b, x, tol, max_iter);
    ck = get_counter();

    printf("Ciclos: %.0f\n", ck);

    free(a);
    free(b);
    free(x);

    return 0;   
}