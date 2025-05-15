#!/bin/bash

mpicc -o heat_pSinc heat_pSinc.c diffusion_pSinc.c faux_pSinc.c

echo "" >../Resultados/resultados-totales-sinc.txt

echo "2 nodo" >>../Resultados/resultados-totales-sinc.txt
mpirun -mca btl openib,vader,self,tcp -mca mpi_cuda_support 0 -hostfile /export/home/nodos/hostfinal -map-by node -bind-to none -n 2 heat_pAsinc ../Cards/card >>../Resultados/resultados-totales-sinc.txt

echo "4 nodo" >>../Resultados/resultados-totales-sinc.txt
mpirun -mca btl openib,vader,self,tcp -mca mpi_cuda_support 0 -hostfile /export/home/nodos/hostfinal -map-by node -bind-to none -n 4 heat_pAsinc ../Cards/card >>../Resultados/resultados-totales-sinc.txt

echo "8 nodo" >>../Resultados/resultados-totales-sinc.txt
mpirun -mca btl openib,vader,self,tcp -mca mpi_cuda_support 0 -hostfile /export/home/nodos/hostfinal -map-by node -bind-to none -n 8 heat_pAsinc ../Cards/card >>../Resultados/resultados-totales-sinc.txt

echo "12 nodo" >>../Resultados/resultados-totales-sinc.txt
mpirun -mca btl openib,vader,self,tcp -mca mpi_cuda_support 0 -hostfile /export/home/nodos/hostfinal -map-by node -bind-to none -n 12 heat_pAsinc ../Cards/card >>../Resultados/resultados-totales-sinc.txt

echo "16 nodo" >>../Resultados/resultados-totales-sinc.txt
mpirun -mca btl openib,vader,self,tcp -mca mpi_cuda_support 0 -hostfile /export/home/nodos/hostfinal -map-by node -bind-to none -n 16 heat_pAsinc ../Cards/card >>../Resultados/resultados-totales-sinc.txt

echo "20 nodo" >>../Resultados/resultados-totales-sinc.txt
mpirun -mca btl openib,vader,self,tcp -mca mpi_cuda_support 0 -hostfile /export/home/nodos/hostfinal -map-by node -bind-to none -n 20 heat_pAsinc ../Cards/card >>../Resultados/resultados-totales-sinc.txt

echo "24 nodo" >>../Resultados/resultados-totales-sinc.txt
mpirun -mca btl openib,vader,self,tcp -mca mpi_cuda_support 0 -hostfile /export/home/nodos/hostfinal -map-by node -bind-to none -n 24 heat_pAsinc ../Cards/card >>../Resultados/resultados-totales-sinc.txt

echo "32 nodo" >>../Resultados/resultados-totales-sinc.txt
mpirun -mca btl openib,vader,self,tcp -mca mpi_cuda_support 0 -hostfile /export/home/nodos/hostfinal -map-by node -bind-to none -n 32 heat_pAsinc ../Cards/card >>../Resultados/resultados-totales-sinc.txt

echo "48 nodo" >>../Resultados/resultados-totales-sinc.txt
mpirun -mca btl openib,vader,self,tcp -mca mpi_cuda_support 0 -hostfile /export/home/nodos/hostfinal -map-by node -bind-to none -n 48 heat_pAsinc ../Cards/card >>../Resultados/resultados-totales-sinc.txt
