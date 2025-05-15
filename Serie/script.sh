#!/bin/bash

icc -o heat_s heat_s.c diffusion_s.c faux_s.c

echo "" >../Resultados/resultados-serie.txt

echo -e "Resultados serie\n" >>../Resultados/resultados-serie.txt
./heat_s ../Cards/card >>../Resultados/resultados-serie.txt
