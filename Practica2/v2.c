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
*/

void v2Jacobi(float a[n][n], float b[n], float x[n], float tol, int max_iter) {
    float *x_new = (float*)malloc(n*sizeof(float));
    
    for(int iter = 0; iter < max_iter; iter++) {
        float norm2 = 0;

        //Iteraciones sobre las filas
        for(int i = 0; i < n; i++) {
            float sigma = 0.0;
            //Precalculamos el inverso de la diagonal
            float aInv = 1.0 / a[i][i];
            
            int j = 0;  
            for(j = 0; j <= n - 4; j += 4) {
                if(i != j) {
                    sigma += a[i][j] * x[j];
                }
                if(i != j + 1) {
                    sigma += a[i][j + 1] * x[j + 1];
                }
                if(i != j + 2) {
                    sigma += a[i][j + 2] * x[j + 2];
                }
                if(i != j + 3) {
                    sigma += a[i][j + 3] * x[j + 3];
                }
            }

            //Ahora, vamos a procesar los últimos elementos
            for(; j < n; j++) {
                if(i != j) {
                    sigma += a[i][j] * x[j];
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

        free(x_new);
        printf("Iteraciones: %d\n", iter);
    }
}

int main(int argc, char *argv[]) {
    srand(time(NULL));

    double ck = 0;
    float a[n][n];
    float b[n];
    float x[n];
    float tol = 1e-8;
    int max_iter = 20000;

    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            a[i][j] = (float)rand() / RAND_MAX;
        }
        b[i] = (float)rand() / RAND_MAX;
        x[i] = 0;
    }
    
    start_counter();
    v2Jacobi(a, b, x, tol, max_iter);
    ck = get_counter();

    printf("Ciclos: %.0f\n", ck);
    return 0;   
}