#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <math.h>

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

    //Número de elementos de A
    unsigned long R = -1; 
    //Valor do parámetro de localidad, espaciado entre elemntos
    unsigned long D = -1;
    //Declaramos o vector A 
    double* A = NULL; 
    //Declaramos o vector S, que almacenará os resultados. O tamaño deste será 10
    double* S = NULL;
    //Variable para contabilizar os ciclos
    double ck = 0.0; 

    //Reservamos memoria para S e comprobamos que a reseva foi exitosa
    S = malloc(10 * sizeof(double));
    if(!S) {
        printf("No se ha reservado memoria para el vector S.\n");
        return EXIT_FAILURE;
    }
    
    //Comprobamos que hai dous argumentos. En caso de non ser así, especificámolo por pantalla e acabamos a ejecución
    if(argc != 3) {
        printf("No se especificaron los valores de D y R\n");
        return EXIT_FAILURE;
    }
    //Recuperamos o valor de D. R calculase como o valor de L introducido polo script, multiplicado por 8/D; por moi grande que sea D, ceil(8/D) será 1 mínimo
    D = atoi(argv[1]);
    if(D < 8) {
        R = atoi(argv[2]) * ceil(8.0/D);
    }
    else {
        R = atoi(argv[2]);
    }

    //Calculamos o tamaño do vector A  
    unsigned long size = R * D;  

    /*  
        En aligned_alloc, o primeiro parámetro é o alineamiento (tamño da línea caché, debe ser potencia de 2)
        O segundo parámetro é o tamaño do vector
    */
    //Reservamos memoria para A. Debemos especificar o tamaño e o alineamiento
    A = (double*)aligned_alloc(CACHE_LINE_SIZE, size*sizeof(double));
    //Comprobamos que se reservase memoria correctamente
    if (!A) {
        printf("Error al asignar memoria para A\n");
        return EXIT_FAILURE;
    }
    
    //Inicializamos o vector A con valores aleatorios no intervalo [-1, 2)
    for(int i = 0; i < size; i++) {
        //Nos aseguramos que nunca genere un 1 aleatoriamente al sumar el +1.0
        A[i] = ((double)rand() / (RAND_MAX + 1.0)) * 3 - 1;
    }

    //Facemos a operación de reducción as 10 veces. Cada unha, almacena o resultado no vector S.
    //Ademais, o contador é inicializado, mediremos os ciclos.
    start_counter();
    for(int i = 0; i < 10; i++) {
        S[i] = 0.0;
        for(int j = 0; j < R; j++) { //Como debemos ter en conta o D, accedemos a A con saltos de D. Debido a esto, o bucle vai de 0 a R.
            S[i] += A[j*D];
        }
    }
    ck = get_counter();
    
    //Printeamos os resultados. Serán os mismos, xa que o vector non cambia
    printf("Los resultados son:\n");
    for(int i = 0; i < 10; i++) {
        printf("Experimento número %d: %lf\n", i, S[i]);
    }

    //Calculamos a media do vector S, que será igual ao valor de cada elemento, xa que son iguais
    double avgS = 0.0;
    for(int i = 0; i < 10; i++) {
        avgS += S[i];
    }
    avgS = avgS / 10;
    printf("Media del vector S: %lf\n", avgS);

    //Imprimimos os ciclos totales das 10 iteracións
    printf("Ciclos totales (10 iteraciones) = %1.10lf \n",ck);

    //Calculamos a media de ciclos en cada iteración (Dividimos o total entre as 10 iteracións)
    double avgck = ck/10;
    printf("Media de ciclos: %lf\n", avgck);

    //Calculamos a media de ciclos por acceso a memoria
    /*
        Esto poderíase facer de 2 formas. Poderíanse contabilizar os accesos manualmente ou, como se fai neste caso,
        temos en conta que o bucle interno realiza R accesos a memoria. Dado que hai 10 iteracións, o total de accesos
        a memoria son os R accesos 10 veces.
    */
    double avgAccesosCiclo = ck/(10*R);
    printf("Media de ciclos por acceso a memoria: %lf\n", avgAccesosCiclo);
    
    printf("\n");

    FILE* doc = fopen("accesoDirectoResultadosDouble.txt", "a+");
    fprintf(doc, "Resultados para D = %ld, L = %d\n", D, atoi(argv[2]));
    fprintf(doc, "Media de ciclos: %lf\n", avgck);
    fprintf(doc, "Media de ciclos por acceso a memoria: %lf\n", avgAccesosCiclo);
    fprintf(doc, "\n");

    //liberamos la memoria
    free(A);
    free(S);

    return 0;
}