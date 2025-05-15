/* File: diffusion_s.c */ 

#include <mpi.h>
#include "defines.h"

/************************************************************************************/
void thermal_update (float *grid, float *grid_chips, float t_ext, int NROW_loc , int NCOL_glob)
{
  int i, j, a, b;
  int pid, npr;

  MPI_Comm_rank( MPI_COMM_WORLD, &pid);
  MPI_Comm_size( MPI_COMM_WORLD, &npr);



  // heat injection at chip positions
  for (i=1; i<=NROW_loc; i++)
    for (j=1; j<NCOL_glob-1; j++) 
      if (grid_chips[i*NCOL_glob+j] > grid[i*NCOL_glob+j])
        grid[i*NCOL_glob+j] += 0.05 * (grid_chips[i*NCOL_glob+j] - grid[i*NCOL_glob+j]);

  // air cooling at the middle of the card
  a = 0.44*(NCOL_glob-2) + 1;
  b = 0.56*(NCOL_glob-2) + 1;

  for (i=1; i<=NROW_loc; i++)
    for (j=a; j<b; j++)
      grid[i*NCOL_glob+j] -= 0.01 * (grid[i*NCOL_glob+j] - t_ext);
  }

/************************************************************************************/
double thermal_diffusion (float *grid, float *grid_aux, int NROW_loc, int NCOL_glob, MPI_Request envio_ant, MPI_Request envio_sig, MPI_Request recv_ant, MPI_Request recv_sig)
{
  int    i, j;
  double  T;
  double Tfull = 0.0;

  //Se procesan todo el trozo menos la primera fila útil y la última
  for (i=2; i<=NROW_loc-1; i++)
    for (j=1; j<NCOL_glob-1; j++)
    {
      T = grid[i*NCOL_glob+j] + 
        0.10 * (grid[(i+1)*NCOL_glob+j]   + grid[(i-1)*NCOL_glob+j]   + grid[i*NCOL_glob+(j+1)]     + grid[i*NCOL_glob+(j-1)] + 
        grid[(i+1)*NCOL_glob+j+1] + grid[(i-1)*NCOL_glob+j+1] + grid[(i+1)*NCOL_glob+(j-1)] + grid[(i-1)*NCOL_glob+(j-1)] 
        - 8*grid[i*NCOL_glob+j]);

      grid_aux[i*NCOL_glob+j] = T;
      Tfull += T;
    }

  //Esperamos a que se reciba el marco superior y despues procesamos la primera fila útil
  MPI_Wait(&recv_ant, MPI_STATUS_IGNORE);
  for (j=1; j<NCOL_glob-1; j++)
  {
    T = grid[1*NCOL_glob+j] + 
      0.10 * (grid[(1+1)*NCOL_glob+j]   + grid[(1-1)*NCOL_glob+j]   + grid[1*NCOL_glob+(j+1)]     + grid[1*NCOL_glob+(j-1)] + 
      grid[(1+1)*NCOL_glob+j+1] + grid[(1-1)*NCOL_glob+j+1] + grid[(1+1)*NCOL_glob+(j-1)] + grid[(1-1)*NCOL_glob+(j-1)] 
      - 8*grid[1*NCOL_glob+j]);

    grid_aux[1*NCOL_glob+j] = T;
    Tfull += T;
  }

  //Esperamos a que se reciba el marco inferior y después procesamos la fila de abajo
  MPI_Wait(&recv_sig, MPI_STATUS_IGNORE);
  i= NROW_loc;
  for (j=1; j<NCOL_glob-1; j++)
  {
    T = grid[i*NCOL_glob+j] + 
      0.10 * (grid[(i+1)*NCOL_glob+j]   + grid[(i-1)*NCOL_glob+j]   + grid[i*NCOL_glob+(j+1)]     + grid[i*NCOL_glob+(j-1)] + 
      grid[(i+1)*NCOL_glob+j+1] + grid[(i-1)*NCOL_glob+j+1] + grid[(i+1)*NCOL_glob+(j-1)] + grid[(i-1)*NCOL_glob+(j-1)] 
      - 8*grid[i*NCOL_glob+j]);

    grid_aux[i*NCOL_glob+j] = T;
    Tfull += T;
  }


  //Le damos los nuevos valores al grid menos a la primera y la última fila útil
  for (i=2; i<=NROW_loc-1; i++)
    for (j=1; j<NCOL_glob-1; j++)
      grid[i*NCOL_glob+j] = grid_aux[i*NCOL_glob+j]; 

  //Esperamos a que se haya enviado la primera fila útil antes de cambiarle el valor
  MPI_Wait(&envio_ant, MPI_STATUS_IGNORE);
  for (j=1; j<NCOL_glob-1; j++)
    grid[1*NCOL_glob+j] = grid_aux[1*NCOL_glob+j]; 

  //Esperamos a que se haya enviado la última fila útil antes de cambiarle el valor
  MPI_Wait(&envio_sig, MPI_STATUS_IGNORE);
  for (j=1; j<NCOL_glob-1; j++)
    grid[NROW_loc*NCOL_glob+j] = grid_aux[NROW_loc*NCOL_glob+j]; 


  return (Tfull);
}

/************************************************************************************/
double calculate_Tmean (float *grid, float *grid_chips, float *grid_aux, float t_delta, int max_iter, float t_ext, int NROW_loc, int NROW_glob, int NCOL_glob)
{
  int    i, j, end, niter, pid, npr;
  double  Tfull, Tfull_tot;
  double Tmean, Tmean0 = t_ext;

  end = 0; niter = 0;
  MPI_Comm_rank( MPI_COMM_WORLD, &pid);
  MPI_Comm_size( MPI_COMM_WORLD, &npr);

  MPI_Request envio_ant, envio_sig, recv_ant, recv_sig;
  envio_ant = envio_sig = recv_ant = recv_sig = MPI_REQUEST_NULL;

  while (end == 0)
  {
    niter++;
    Tmean = 0.0;
    int nreq = 0;

    // heat injection and air cooling 
    thermal_update (grid, grid_chips, t_ext, NROW_loc, NCOL_glob);

    if (pid != 0) MPI_Isend(&grid[NCOL_glob], NCOL_glob, MPI_FLOAT, pid-1, 0, MPI_COMM_WORLD, &envio_ant); //Cada proceso envia su primera fila útil al anterior
    if (pid != npr-1) MPI_Isend(&grid[NROW_loc*NCOL_glob], NCOL_glob, MPI_FLOAT, pid+1, 0, MPI_COMM_WORLD, &envio_sig); //Cada proceso envia su última fila útil al siguiente


    if (pid != 0) MPI_Irecv(&grid[0], NCOL_glob, MPI_FLOAT, pid-1, 0, MPI_COMM_WORLD,&recv_ant); //Cada proceso recibe el marco de arriba del proceso anterior
    if (pid != npr-1) MPI_Irecv(&grid[(NROW_loc+1)*NCOL_glob], NCOL_glob, MPI_FLOAT, pid+1, 0, MPI_COMM_WORLD, &recv_sig); //Cada proceso recibe el marco de abajo del proceso siguiente

    // thermal diffusion
    Tfull = thermal_diffusion(grid, grid_aux, NROW_loc, NCOL_glob, envio_ant, envio_sig, recv_ant, recv_sig);


    // convergence every 10 iterations
    if (niter % 10 == 0)
    {
      MPI_Allreduce(&Tfull, &Tfull_tot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      Tmean = Tfull_tot / ((NCOL_glob-2)*(NROW_glob-2));
      if ((fabs(Tmean - Tmean0) < t_delta) || (niter > max_iter))
        end = 1;
      else Tmean0 = Tmean;
    }
  } // end while 

  if(pid==0)
    printf ("Iter (par): %d\t", niter);
  return (Tmean);
}

