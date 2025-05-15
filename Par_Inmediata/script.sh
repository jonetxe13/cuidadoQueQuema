#!/bin/bash

mpicc -o heat_pAsinc heat_pAsinc.c diffusion_pAsinc.c faux_pAsinc.c

touch ../Resultados/resultados-totales-asinc.txt

echo "2 nodo" >>../Resultados/resultados-totales-asinc.txt
mpirun -mca btl openib,vader,self,tcp -mca mpi_cuda_support 0 -hostfile /export/home/nodos/hostfinal -map-by node -bind-to none -n 2 heat_pAsinc ../Cards/card >>../Resultados/resultados-totales-asinc.txt

echo "4 nodo" >>../Resultados/resultados-totales-asinc.txt
mpirun -mca btl openib,vader,self,tcp -mca mpi_cuda_support 0 -hostfile /export/home/nodos/hostfinal -map-by node -bind-to none -n 4 heat_pAsinc ../Cards/card >>../Resultados/resultados-totales-asinc.txt

echo "8 nodo" >>../Resultados/resultados-totales-asinc.txt
mpirun -mca btl openib,vader,self,tcp -mca mpi_cuda_support 0 -hostfile /export/home/nodos/hostfinal -map-by node -bind-to none -n 8 heat_pAsinc ../Cards/card >>../Resultados/resultados-totales-asinc.txt

echo "12 nodo" >>../Resultados/resultados-totales-asinc.txt
mpirun -mca btl openib,vader,self,tcp -mca mpi_cuda_support 0 -hostfile /export/home/nodos/hostfinal -map-by node -bind-to none -n 12 heat_pAsinc ../Cards/card >>../Resultados/resultados-totales-asinc.txt

echo "16 nodo" >>../Resultados/resultados-totales-asinc.txt
mpirun -mca btl openib,vader,self,tcp -mca mpi_cuda_support 0 -hostfile /export/home/nodos/hostfinal -map-by node -bind-to none -n 16 heat_pAsinc ../Cards/card >>../Resultados/resultados-totales-asinc.txt

echo "20 nodo" >>../Resultados/resultados-totales-asinc.txt
mpirun -mca btl openib,vader,self,tcp -mca mpi_cuda_support 0 -hostfile /export/home/nodos/hostfinal -map-by node -bind-to none -n 20 heat_pAsinc ../Cards/card >>../Resultados/resultados-totales-asinc.txt

echo "24 nodo" >>../Resultados/resultados-totales-asinc.txt
mpirun -mca btl openib,vader,self,tcp -mca mpi_cuda_support 0 -hostfile /export/home/nodos/hostfinal -map-by node -bind-to none -n 24 heat_pAsinc ../Cards/card >>../Resultados/resultados-totales-asinc.txt

echo "32 nodo" >>../Resultados/resultados-totales-asinc.txt
mpirun -mca btl openib,vader,self,tcp -mca mpi_cuda_support 0 -hostfile /export/home/nodos/hostfinal -map-by node -bind-to none -n 32 heat_pAsinc ../Cards/card >>../Resultados/resultados-totales-asinc.txt

echo "48 nodo" >>../Resultados/resultados-totales-asinc.txt
mpirun -mca btl openib,vader,self,tcp -mca mpi_cuda_support 0 -hostfile /export/home/nodos/hostfinal -map-by node -bind-to none -n 48 heat_pAsinc ../Cards/card >>../Resultados/resultados-totales-asinc.txt
