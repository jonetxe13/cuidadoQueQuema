mpicc -o heat_p heat_p.c diffusion_p.c faux_p.c
mpirun -mca btl openib,vader,self,tcp -mca mpi_cuda_support 0 -hostfile /export/home/nodos/hostfinal -map-by node -bind-to none -n $1 heat_p card0
# ./heat_p card0
# vfinder card0_par.res
