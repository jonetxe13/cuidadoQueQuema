mpicc -o heat_p heat_p.c diffusion_p.c faux_p.c
./heat_p card0
vfinder card0_par.res
