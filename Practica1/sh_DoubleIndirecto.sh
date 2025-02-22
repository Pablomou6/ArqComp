#!/bin/bash
# Solicitamos un nodo con 64 cores y 64 GB de memoria durante 10 minutos
#SBATCH -n 1 -c 64 -t 00:10:00 --mem=64G
# Ponemos nombre a nuestro trabajo para poder identificarlo.
# ATENCIÓN - Debes sustituir el NN por el número de equipo.
#SBATCH --job-name p1acg03

# Sustituir los valores de Di y Li por los calculados para la realización de la práctica.

#gcc accesoIndirecto.c -o acp1 -msse3 -O0 -lm
gcc accesoDirecto.c -o acp2 -msse3 -O0 -lm
#gcc accesoIndirectoInt.c -o acp3 -msse3 -O0 -lm


for D in {1,4,32,64,512}
do
	for L in {384,1152,10240,15360,40960,81920,163840} 
	do
		for i in {1..10}
		do
			echo "Valores: D = $D, L = $L; iteración = $i"
			#./acp1 $D $L
			./acp2 $D $L
			#./acp3 $D $L 
		done
		echo "-----------------------------------------"
		echo ""
	done
done
