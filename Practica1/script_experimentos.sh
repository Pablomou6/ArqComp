#!/bin/bash
# Solicitamos un nodo con 64 cores y 256 GB de memoria durante 2 horas
#SBATCH -n 1 -c 64 -t 02:00:00 --mem=256G
# Ponemos nombre a nuestro trabajo para poder identificarlo.
# ATENCIÓN - Debes sustituir el NN por el número de equipo.
#SBATCH --job-name p1acgNN

# Sustituir los valores de Di y Li por los calculados para la realización de la práctica.
#
gcc main.c -o acp1 -msse3 -O0

doc1="experimento1.txt"
doc2="experimento2.txt"
doc3="experimento3.txt"

for D in 1 #,4,64,256,1024}
do
	for R in {3072,9216,81920,122880,327680,655360,1310720}
	do
		for i in {1..10}; do
			./acp1 $D $R #>> $doc1 #Redireccionamos el output al documento que queremos, el doble ">" es para que no sobreescriba
		done
	done
done


