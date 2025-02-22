#!/bin/bash
# Solicitamos un nodo con 64 cores y 256 GB de memoria durante 2 horas
#SBATCH -n 1 -c 64 -t 02:00:00 --mem=256G
# Ponemos nombre a nuestro trabajo para poder identificarlo.
# ATENCIÓN - Debes sustituir el NN por el número de equipo.
#SBATCH --job-name p1acgNN

# Sustituir los valores de Di y Li por los calculados para la realización de la práctica.


#Utilizamos el mismo script para los programas con tipo de dato double, ya que ambos tienen los mismos valores de D y R, solo cambia el .c que se ejecuta

gcc accesoIndirecto.c -o acp1 -msse3 -O0
#gcc accesoDirecto.c -o acp3 -msse3 -O0

doc1="accesoIndirectoResultadosDouble.txt"
#doc2="accesoDirectoResultadosDouble.txt"

#Facemos un primeiro bucle para D = 1, xa que ten uns valores de R concretos.
for D in 1
do
	for R in {3072,9216,81920,122880,327680,655360,1310720}
	do
		for i in {1..10}
		do
			echo "Valores: D = $D, R = $R; iteración = $i"
			./acp1 $D $R 
		done
		echo "-----------------------------------------"
		echo ""
	done
done

#Da mesma forma, para D = 4 temos uns valores de R concretos.
for D in 4
do
    for R in {768,2304,20480,30720,81920,163840,327680}
	do
		for i in {1..10}
		do
			echo "Valores: D = $D, R = $R; iteración = $i"
			./acp1 $D $R 
		done
		echo "-----------------------------------------"
		echo ""
	done
done

#Como R = L*(8/D), unha vez que D = 8, temos que R = L. Desta forma, con D's máis grandes; 8/D será menor que 1. Polo tanto, con D >= 8; R = L. 
for D in {64,256,1024}
do
	for R in {384,1152,10240,15360,40960,81920,163840}
	do
		for i in {1..10}
		do
			echo "Valores: D = $D, R = $R; iteración = $i"
			./acp1 $D $R 
		done
		echo "-----------------------------------------"
		echo ""
	done
done
