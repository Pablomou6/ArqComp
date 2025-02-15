#!/bin/bash
# Solicitamos un nodo con 64 cores y 256 GB de memoria durante 2 horas
#SBATCH -n 1 -c 64 -t 02:00:00 --mem=256G
# Ponemos nombre a nuestro trabajo para poder identificarlo.
# ATENCIÓN - Debes sustituir el NN por el número de equipo.
#SBATCH --job-name p1acgNN

# Sustituir los valores de Di y Li por los calculados para la realización de la práctica.
#
gcc main.c -o acp1 -msse3 -O0

#TODO
#for i in {1..10} Comentamos esto pq xa se ejecutan 10 veces no código, preguntar se é necesario 
#manter este ou o outro, mais que nada polo vector S
#do
	for D in {1} #,4,64,256,1024}
	do
		for L in {3072,9216,81920,122880,327680,655360,1310720}
		do
			./acp1 $D $L
		done
	done
#done


