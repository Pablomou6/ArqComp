#!/bin/bash

#Solicitamos un nodo con 64 cores y 64 GB de memoria durante 2 minutos
#SBATCH -n 1 -c 64 -t 00:02:00 --mem=64G
# Ponemos nombre a nuestro trabajo para poder identificarlo.
# ATENCIÓN - Debes sustituir el NN por el número de equipo.
#SBATCH --job-name p2acg03

#Es de especial interés representar la ganancia en velocidad (también llamada aceleración o
#speedup) de la versión secuencial optimizada con respecto a la versión inicial compilada con -
#O0, la ganancia en velocidad de los códigos de los apartados iii) y iv) con respecto a la versión
#secuencial optimizada, y finalmente las versiones desarrolladas en los apartados ii), iii) y iv)
#respecto de la versión inicial compilada con -O3).
#En el caso del apartado iv), además de las gráficas que se consideren adecuadas, es
#imprescindible representar en una gráfica para cada valor de N la ganancia en velocidad conseguida.

#El texto anterior nos lleva a las siguientes comparaciones:
    #- secuencial optimizada (v2) vs inicial con -O0
    #- ⁠v3 vs v2
    #- ⁠v4 vs v2
    #- ⁠v2 vs v1 -O3
    #- ⁠v3 vs v1 -O3 
    #- ⁠v4 vs v1 -O3
    #- ⁠v4 depende de N

doc="v1_250.txt"

gcc -Wall -o v1 v1.c -O0 -lm
#gcc -Wall -o v1 v1.c -O3 -lm
#gcc -Wall -o v2 v2.c -lm
#gcc -Wall -o v3 v3.c -mavx2 -mfma -lm
#gcc -Wall -o v4 v4.c -fopenmp -lm

# Eliminamos el contenido del fichero si existe
> $doc

for i in {1..15}
do
    # Ejecutamos 15 veces el programa correspondiente y guardamos los resultados en un fichero
    ./v1 250 >> $doc
done

#Ahora ordenamos el fichero para que el primer valor sea el más pequeño y el último el más grande
sort -n $doc -o $doc
