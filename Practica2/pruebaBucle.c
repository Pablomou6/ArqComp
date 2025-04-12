#include <stdio.h>
#include <stdlib.h>

int main() {
    int n = 0;
    float sum = 0.0;

    printf("Introduce un tamaño para la matriz.\n");
    scanf(" %d", &n);

    float** matriz = (float**)aligned_alloc(64, n*sizeof(float*));
    if(matriz == NULL) {
        printf("Error: no se ha podido reservar memoria para la matriz de coeficientes.\n");
        return EXIT_FAILURE;
    }
    //Se reserva memoria para cada fila de la matriz
    for(int i = 0; i < n; i++) {
        matriz[i] = (float*)aligned_alloc(64, n * sizeof(float));
        if(matriz[i] == NULL) {
            printf("Error: no se ha podido reservar memoria para la fila %d de la matriz de coeficientes.\n", i);
            //Liberamos la memoria reservada hasta el momento
            for(int j = 0; j < i; j++) {
                free(matriz[j]);
            }
            free(matriz);
            return EXIT_FAILURE;
        }
    }

    // Inicializamos la matriz 
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) {
                matriz[i][j] = -n+1; // Diagonal principal con -n+1
            } else {
                matriz[i][j] = 1; // Resto de columnas con 1
            }
        }
    }

    /* Imprimimos la matriz
    printf("La matriz es:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%6.2f ", matriz[i][j]);
        }
        printf("\n");
    }*/

    /*
        Ahora mismo, en una matriz n = 10, tendremos ue la diagonal es -9 y el resto suma 9, por lo que, para probar el bucl haremos lo siguiente:
        Recorremos las filas con el bucle dividido como en la práctica, para evitar la diagonal, de forma ue haremos un sumatorio. Si algún sumatorio
        da resultado diferente de n-1 (si cuenta la diagonal, deberá de dar 0), se deben añadir las comprobaciones.
    */

    //Hacemos un primer bucle para recorrer todas las filas
    for(int i = 0; i < n; i++) {
        sum = 0.0;
        //Ahora hacemos un bucle que recorrerá, de 4 en 4, todas las columnas previas a la diagonal (i == j) y sumará los valores
        int j;
        for(j = 0; j <= i - 4; j += 4) {
            sum += matriz[i][j];
            sum += matriz[i][j+1];
            sum += matriz[i][j+2];
            sum += matriz[i][j+3];
        }
        //El siguiente bucle procesará el resto de elementos hasta la diagonal (excluyéndola)
        for (; j < i; j++) {
            sum += matriz[i][j];
        }

        //Ahora haremos lo mismo, pero para las columnas posteriores a la diagonal
        for(j = i + 1; j <= n - 4; j += 4) {
            sum += matriz[i][j];
            sum += matriz[i][j+1];
            sum += matriz[i][j+2];
            sum += matriz[i][j+3];
        }
        //El siguiente bucle procesará el resto de elementos hasta la diagonal (excluyéndola)
        for (; j < n; j++) {
            sum += matriz[i][j];
        }

        printf("Sumatorio de la fila %d: %f\n", i, sum);

    }

    /*
        Otra alternativa será no hacer un sumatorio, sino que, como ahora no nos importa el rendimiento, hacer los 4 if's ue comprueben el criterio.
        En caso de que coincida, imprimirá un mensaje.
    */
    //Hacemos un primer bucle para recorrer todas las filas
    for(int i = 0; i < n; i++) {
        
        int j;
        for(j = 0; j <= i - 4; j += 4) {
            if (j == i) printf("Ha coincidido.\n");
            if (j + 1 == i) printf("Ha coincidido.\n");
            if (j + 2 == i) printf("Ha coincidido.\n");
            if (j + 3 == i) printf("Ha coincidido.\n");
        }
        //El siguiente bucle procesará el resto de elementos hasta la diagonal (excluyéndola)
        for (; j < i; j++) {
            if(j == i) printf("Ha coincidido.\n");
        }

        //Ahora haremos lo mismo, pero para las columnas posteriores a la diagonal
        for(j = i + 1; j <= n - 4; j += 4) {
            if (j == i) printf("Ha coincidido.\n");
            if (j + 1 == i) printf("Ha coincidido.\n");
            if (j + 2 == i) printf("Ha coincidido.\n");
            if (j + 3 == i) printf("Ha coincidido.\n");
        }
        //El siguiente bucle procesará el resto de elementos hasta la diagonal (excluyéndola)
        for (; j < n; j++) {
            if(j == i) printf("Ha coincidido.\n");
        }

    }

    return EXIT_SUCCESS;
}