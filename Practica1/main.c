#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

#define CACHE_LINE_SIZE 64


////////////////// FUNCIONES DE MEDICIÓN DE TIEMPOS //////////////////////
void start_counter();
double get_counter();
double mhz();

/* Initialize the cycle counter */
 static unsigned cyc_hi = 0;
 static unsigned cyc_lo = 0;


 /* Set *hi and *lo to the high and low order bits of the cycle counter.
 Implementation requires assembly code to use the rdtsc instruction. */
 void access_counter(unsigned *hi, unsigned *lo)
 {
 asm("rdtsc; movl %%edx,%0; movl %%eax,%1" /* Read cycle counter */
 : "=r" (*hi), "=r" (*lo) /* and move results to */
 : /* No input */ /* the two outputs */
 : "%edx", "%eax");
 }

 /* Record the current value of the cycle counter. */
 void start_counter()
 {
 access_counter(&cyc_hi, &cyc_lo);
 }

 /* Return the number of cycles since the last call to start_counter. */
 double get_counter()
 {
 unsigned ncyc_hi, ncyc_lo;
 unsigned hi, lo, borrow;
 double result;

 /* Get cycle counter */
 access_counter(&ncyc_hi, &ncyc_lo);

 /* Do double precision subtraction */
 lo = ncyc_lo - cyc_lo;
 borrow = lo > ncyc_lo;
 hi = ncyc_hi - cyc_hi - borrow;
 result = (double) hi * (1 << 30) * 4 + lo;
 if (result < 0) {
 fprintf(stderr, "Error: counter returns neg value: %.0f\n", result);
 }
 return result;
 }

double mhz(int verbose, int sleeptime)
 {
 double rate;

 start_counter();
 sleep(sleeptime);
 rate = get_counter() / (1e6*sleeptime);
 if (verbose)
 printf("\n Processor clock rate = %.1f MHz\n", rate);
 return rate;
 }
 /////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[]) {
    //Declaramos la semilla para los número aleatorios
    srand(time(NULL));

    int R = -1; //Valor inda desconocido (tamaño do vector)
    int D = -1; //Valor inda desconocido (parámetro localidad, espaciado entre elemntos)
    double* A = NULL; //Tamaño de A inda desconocido
    double sum = 0.0;
    double* S = NULL; //Tamaño de S; 10

    //Variable para la medida de tiempos
    double ck = 0.0;

    //Reservamos y comprobamos que la reseva fuese exitosa
    S = malloc(10 * sizeof(double));
    if(!S) {
        printf("No se ha reservado memoria para el vector S.\n");
        return EXIT_FAILURE;
    }
    
    //Recuperamos los valores de D y L
    if(argc != 3) {
        printf("No se especificaron los valores de D y R\n");
        return EXIT_FAILURE;
    }
    D = atoi(argv[1]);
    R = atoi(argv[2]);

    //Hacemos el vector para el acceso indirecto
    int ind[R];
    for(int i = 0; i < R; i++) {
        ind[i] = i * D;
    }

    //En aligned_alloc, el primer parámetro es el alineamiento (tamño de la línea caché, debe ser potencia de 2)
    //El segundo parámetro es el tamaño del vector
    //Debemos asegurarnos que el tamaño reservado, sea múltiplo de, en este caso, el tamaño de la línea caché
    size_t size = R * sizeof(double);
    if (size % CACHE_LINE_SIZE != 0) {
        //Nos aseguramos que sea múltiplo de 64, ya que la dvisión se encarga de truncar
        size = ((size / CACHE_LINE_SIZE) + 1) * CACHE_LINE_SIZE; 
    }

    A = (double*)aligned_alloc(CACHE_LINE_SIZE, size);
    //Debemos comprobar que se haya reservado memoria correctamente
    if (!A) {
        printf("Error al asignar memoria para A\n");
        return EXIT_FAILURE;
    }
    
    //inicializamos el vector
    for(int i = 0; i < R; i++) {
        A[i] = ((double)rand() / RAND_MAX) * 3.8 - 1.9;
    }

    //declaramos una variable para contabilizar los accesos
    int accesos = 0;

    //Hacemos la operación de reducción 10 veces. Inicializamos el contador, el código siguiente será medido
    start_counter();
    for(int i = 0; i < 10; i++) {
        for(int j = 0; j < R; j++) {
            sum += A[ind[j]];
            accesos++;
        }
        S[i] = sum;
        sum = 0.0;
    }
    ck=get_counter();
    printf("Ciclos = %1.10lf \n",ck);
    printf("Accesos medios por cada ciclo: %lf\n", (double)accesos/ck);
    //Las dos líneas superiores terminan la medición

    printf("Los resultados son:\n");
    for(int i = 0; i < 10; i++) {
        printf("Experimento número %d: %lf\n", i, S[i]);
    }

    double aux = 0.0;
    for(int i = 0; i < 10; i++) {
        aux += S[i];
    }
    
    printf("La media de S es: %lf\n", aux/10);
    printf("\n");

    return 0;
}