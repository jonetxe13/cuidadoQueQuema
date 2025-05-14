mpicc -o heat_pAsinc heat_pAsinc.c diffusion_pAsinc.c faux_pAsinc.c
echo "2 nodo" >> resultadosd-totales-asinc.txt
mpirun -mca btl openib,vader,self,tcp -mca mpi_cuda_support 0 -hostfile /export/home/nodos/hostfinal -map-by node -bind-to none -n 2 heat_pAsinc card>> resultadosd-totales-asinc.txt
echo "4 nodo" >> resultadosd-totales-asinc.txt
mpirun -mca btl openib,vader,self,tcp -mca mpi_cuda_support 0 -hostfile /export/home/nodos/hostfinal -map-by node -bind-to none -n 4 heat_pAsinc card>> resultadosd-totales-asinc.txt
echo "8 nodo" >> resultadosd-totales-asinc.txt
mpirun -mca btl openib,vader,self,tcp -mca mpi_cuda_support 0 -hostfile /export/home/nodos/hostfinal -map-by node -bind-to none -n 8 heat_pAsinc card>> resultadosd-totales-asinc.txt
echo "16 nodo" >> resultadosd-totales-asinc.txt
mpirun -mca btl openib,vader,self,tcp -mca mpi_cuda_support 0 -hostfile /export/home/nodos/hostfinal -map-by node -bind-to none -n 16 heat_pAsinc card>> resultadosd-totales-asinc.txt
echo "24 nodo" >> resultadosd-totales-asinc.txt
mpirun -mca btl openib,vader,self,tcp -mca mpi_cuda_support 0 -hostfile /export/home/nodos/hostfinal -map-by node -bind-to none -n 24 heat_pAsinc card>> resultadosd-totales-asinc.txt
echo "32 nodo" >> resultadosd-totales-asinc.txt
mpirun -mca btl openib,vader,self,tcp -mca mpi_cuda_support 0 -hostfile /export/home/nodos/hostfinal -map-by node -bind-to none -n 32 heat_pAsinc card>> resultadosd-totales-asinc.txt
echo "64 nodo" >> resultadosd-totales-asinc.txt
mpirun -mca btl openib,vader,self,tcp -mca mpi_cuda_support 0 -hostfile /export/home/nodos/hostfinal -map-by node -bind-to none -n 64 heat_pAsinc card>> resultadosd-totales-asinc.txt
# ./heat_p card0
# vfinder card0_par.res


## poner ./script.sh "num_nodos"
