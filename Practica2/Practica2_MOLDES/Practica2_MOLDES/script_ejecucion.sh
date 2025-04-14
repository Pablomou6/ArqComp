#!/bin/bash

#Solicitamos un nodo con 64 cores y 64 GB de memoria durante 2 minutos
#SBATCH -n 1 -c 64 -t 00:02:00 --mem=64G
# Ponemos nombre a nuestro trabajo para poder identificarlo.
# ATENCIÓN - Debes sustituir el NN por el número de equipo.
#SBATCH --job-name p2acg03

make

make run