#!/bin/bash
# Solicitamos un nodo con 64 cores y 64 GB de memoria durante 10 minutos
#SBATCH -n 1 -c 64 -t 00:10:00 --mem=64G
# Ponemos nombre a nuestro trabajo para poder identificarlo.
# ATENCIÓN - Debes sustituir el NN por el número de equipo.
#SBATCH --job-name p1acg03

# Sustituir los valores de Di y Li por los calculados para la realización de la práctica.


#Utilizamos el mismo script para los programas con tipo de dato double, ya que ambos tienen los mismos valores de D y R, solo cambia el .c que se ejecuta

gcc accesoIndirectoInt.c -o acp3 -msse3 -O0

#Facemos un primeiro bucle para D = 1, xa que ten uns valores de R concretos.
for D in 1
do
	for R in {6144,18432,163840,245760,655360,1310720,2621440}
	do
		for i in {1..10}
		do
			echo "Valores: D = $D, R = $R; iteración = $i"
			./acp3 $D $R 
		done
		echo "-----------------------------------------"
		echo ""
	done
done

#Da mesma forma, para D = 4 temos uns valores de R concretos.
for D in 4
do
    for R in {1536,4608,40960,61440,163840,327680,655360}
	do
		for i in {1..10}
		do
			echo "Valores: D = $D, R = $R; iteración = $i"
			./acp3 $D $R 
		done
		echo "-----------------------------------------"
		echo ""
	done
done

#Como R = L*(16/D), unha vez que D = 16, temos que R = L. Desta forma, con D's máis grandes; 16/D será menor que 1. Polo tanto, con D >= 16; R = L. 
for D in {64,256,1024}
do
	for R in {384,1152,10240,15360,40960,81920,163840}
	do
		for i in {1..10}
		do
			echo "Valores: D = $D, R = $R; iteración = $i"
			./acp3 $D $R 
		done
		echo "-----------------------------------------"
		echo ""
	done
done
