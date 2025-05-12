/* File: diffusion_s.c */ 

#include <mpi.h>
#include "defines.h"

extern int pid;
extern int npr;

/************************************************************************************/
void thermal_update (struct info_param param, float *grid, float *grid_chips, int NROW_loc)
{
  int i, j, a, b;

  // heat injection at chip positions
  for (i=1; i<NROW_loc-1; i++)
  for (j=1; j<NCOL-1; j++) 
    if (grid_chips[i*NCOL+j] > grid[i*NCOL+j])
      grid[i*NCOL+j] += 0.05 * (grid_chips[i*NCOL+j] - grid[i*NCOL+j]);

  // air cooling at the middle of the card
  a = 0.44*(NCOL-2) + 1;
  b = 0.56*(NCOL-2) + 1;

  for (i=1; i<NROW_loc-1; i++)
  for (j=a; j<b; j++)
      grid[i*NCOL+j] -= 0.01 * (grid[i*NCOL+j] - param.t_ext);
}

/************************************************************************************/
double thermal_diffusion (struct info_param param, float *grid, float *grid_aux, int NROW_loc)
{
  int    i, j;
  double  T;
  double Tfull = 0.0;

  for (i=1; i<NROW_loc-1; i++)
    for (j=1; j<NCOL-1; j++)
    {
      T = grid[i*NCOL+j] + 
        0.10 * (grid[(i+1)*NCOL+j]   + grid[(i-1)*NCOL+j]   + grid[i*NCOL+(j+1)]     + grid[i*NCOL+(j-1)] + 
        grid[(i+1)*NCOL+j+1] + grid[(i-1)*NCOL+j+1] + grid[(i+1)*NCOL+(j-1)] + grid[(i-1)*NCOL+(j-1)] 
        - 8*grid[i*NCOL+j]);

      grid_aux[i*NCOL+j] = T;
      Tfull += T;
    }

  //new values for the grid
  for (i=1; i<NROW_loc-1; i++)
    for (j=1; j<NCOL-1; j++)
      grid[i*NCOL+j] = grid_aux[i*NCOL+j]; 

  return (Tfull);
}

/************************************************************************************/
double calculate_Tmean (struct info_param param, float *grid, float *grid_chips, float *grid_aux, int NROW_loc)
{
  int    i, j, end, niter, pid;
  double  Tfull, Tfull_tot;
  double Tmean, Tmean0 = param.t_ext;

  end = 0; niter = 0;
  MPI_Comm_rank( MPI_COMM_WORLD, &pid);

  while (end == 0)
  {
    niter++;
    Tmean = 0.0;

    // heat injection and air cooling 
    thermal_update (param, grid, grid_chips, NROW_loc);

    /*
    if (pid==0){
      MPI_Ssend(&grid[(NROW_loc*NCOL)-(2*NCOL)], NCOL, MPI_FLOAT, pid+1, 0, MPI_COMM_WORLD); //cada pid le envía su última fila(real) al siguiente pid
      MPI_Recv(&grid[(NROW_loc*NCOL)-NCOL], NCOL, MPI_FLOAT, pid+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //Todos los pid menos el 0 esperan a recibir la última fila del anterior pid para poder hacer los cálculos
    }
    else if(pid == npr-1){
      MPI_Recv(&grid[0], NCOL, MPI_FLOAT, pid-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //Todos los pid menos el 0 esperan a recibir la última fila del anterior pid para poder hacer los cálculos
      MPI_Ssend(&grid[NCOL], NCOL, MPI_FLOAT, pid-1, 0, MPI_COMM_WORLD); //cada pid le envía su última fila(real) al siguiente pid
    }
    else{
      MPI_Recv(grid, NCOL, MPI_FLOAT, pid-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //Todos los pid menos el 0 esperan a recibir la última fila del anterior pid para poder hacer los cálculos
      MPI_Ssend(&grid[NCOL], NCOL, MPI_FLOAT, pid-1, 0, MPI_COMM_WORLD); //cada pid le envía su última fila(real) al siguiente pid
      MPI_Ssend(&grid[(NROW_loc*NCOL)-2*NCOL], NCOL, MPI_FLOAT, pid+1, 0, MPI_COMM_WORLD); //cada pid le envía su última fila(real) al siguiente pid
      MPI_Recv(&grid[(NROW_loc*NCOL)-NCOL], NCOL, MPI_FLOAT, pid+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //Todos los pid menos el 0 esperan a recibir la última fila del anterior pid para poder hacer los cálculos
    }*/

    if (pid==0){
      MPI_Ssend(&grid[(NROW_loc-2)*NCOL+0], NCOL, MPI_FLOAT, pid+1, 0, MPI_COMM_WORLD); //cada pid le envía su última fila(real) al siguiente pid
      MPI_Recv(&grid[(NROW_loc-1)*NCOL+0], NCOL, MPI_FLOAT, pid+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //Todos los pid menos el 0 esperan a recibir la última fila del anterior pid para poder hacer los cálculos
    }
    else if(pid == npr-1){
      MPI_Recv(grid, NCOL, MPI_FLOAT, pid-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //Todos los pid menos el 0 esperan a recibir la última fila del anterior pid para poder hacer los cálculos
      MPI_Ssend(&grid[1*NCOL+0], NCOL, MPI_FLOAT, pid-1, 0, MPI_COMM_WORLD); //cada pid le envía su última fila(real) al siguiente pid
    }
    else{
      MPI_Recv(grid, NCOL, MPI_FLOAT, pid-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //Todos los pid menos el 0 esperan a recibir la última fila del anterior pid para poder hacer los cálculos
      MPI_Ssend(&grid[1*NCOL+0], NCOL, MPI_FLOAT, pid-1, 0, MPI_COMM_WORLD); //cada pid le envía su última fila(real) al siguiente pid
      MPI_Ssend(&grid[(NROW_loc-2)*NCOL+0], NCOL, MPI_FLOAT, pid+1, 0, MPI_COMM_WORLD); //cada pid le envía su última fila(real) al siguiente pid
      MPI_Recv(&grid[(NROW_loc-1)*NCOL+0], NCOL, MPI_FLOAT, pid+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //Todos los pid menos el 0 esperan a recibir la última fila del anterior pid para poder hacer los cálculos
    }

    // thermal diffusion
    Tfull = thermal_diffusion(param, grid, grid_aux, NROW_loc);


    // convergence every 10 iterations
    if (niter % 10 == 0)
    {
      MPI_Reduce(&Tfull, &Tfull_tot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      if (pid==0)
      {
        Tmean = Tfull_tot / ((NCOL-2)*(NROW-2));
        if ((fabs(Tmean - Tmean0) < param.t_delta) || (niter > param.max_iter))
          end = 1;
        else Tmean0 = Tmean;
      }
      MPI_Bcast(&end, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }
  } // end while 

  //sacar el total de las iteraciones y mandarlas a pid=0
 // int niter_total;
  //MPI_Reduce(&niter, &niter_total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  if(pid==0)
    //printf ("Iter (par): %d\t", niter_total);
    printf ("Iter (par): %d\t", niter);
  return (Tmean);
}

