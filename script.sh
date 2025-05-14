mpicc -o heat_pAsinc heat_pAsinc.c diffusion_pAsinc.c faux_pAsinc.c
mpirun -mca btl openib,vader,self,tcp -mca mpi_cuda_support 0 -hostfile /export/home/nodos/hostfinal -map-by node -bind-to none -n $1 heat_pAsinc card0
# ./heat_p card0
# vfinder card0_par.res


## poner ./script.sh "num_nodos"
