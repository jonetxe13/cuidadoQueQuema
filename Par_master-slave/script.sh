#!/bin/bash

mpicc -o heat_pMS heat_pMS.c diffusion_pMS.c faux_pMS.c

echo "" >../Resultados/resultados-totales-ms.txt

echo "9 nodo" >>../Resultados/resultados-totales-ms.txt
mpirun -mca btl openib,vader,self,tcp -mca mpi_cuda_support 0 -hostfile /export/home/nodos/hostfinal -map-by node -bind-to none -n 9 heat_pMS ../Cards/card >>../Resultados/resultados-totales-ms.txt

echo "17 nodo" >>../Resultados/resultados-totales-ms.txt
mpirun -mca btl openib,vader,self,tcp -mca mpi_cuda_support 0 -hostfile /export/home/nodos/hostfinal -map-by node -bind-to none -n 17 heat_pMS ../Cards/card >>../Resultados/resultados-totales-ms.txt

echo "25 nodo" >>../Resultados/resultados-totales-ms.txt
mpirun -mca btl openib,vader,self,tcp -mca mpi_cuda_support 0 -hostfile /export/home/nodos/hostfinal -map-by node -bind-to none -n 25 heat_pMS ../Cards/card >>../Resultados/resultados-totales-ms.txt

echo "33 nodo" >>../Resultados/resultados-totales-ms.txt
mpirun -mca btl openib,vader,self,tcp -mca mpi_cuda_support 0 -hostfile /export/home/nodos/hostfinal -map-by node -bind-to none -n 33 heat_pMS ../Cards/card >>../Resultados/resultados-totales-ms.txt

echo "41 nodo" >>../Resultados/resultados-totales-ms.txt
mpirun -mca btl openib,vader,self,tcp -mca mpi_cuda_support 0 -hostfile /export/home/nodos/hostfinal -map-by node -bind-to none -n 41 heat_pMS ../Cards/card >>../Resultados/resultados-totales-ms.txt

echo "49 nodo" >>../Resultados/resultados-totales-ms.txt
mpirun -mca btl openib,vader,self,tcp -mca mpi_cuda_support 0 -hostfile /export/home/nodos/hostfinal -map-by node -bind-to none -n 49 heat_pMS ../Cards/card >>../Resultados/resultados-totales-ms.txt
