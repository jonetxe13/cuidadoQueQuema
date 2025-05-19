#!/bin/bash

mpicc -o heat_pAsinc heat_pAsinc.c diffusion_pAsinc.c faux_pAsinc.c
mpirun -mca btl openib,vader,self,tcp -mca mpi_cuda_support 0 -hostfile /export/home/nodos/hostfinal -map-by node -bind-to none -n 33 heat_pAsinc ../Cards/card
