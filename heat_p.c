/* heat_s.c 

     Difusion del calor en 2 dimensiones      Version en serie

     Se analizan las posiciones de los chips en una tarjeta, para conseguir la temperatura minima
     de la tarjeta. Se utiliza el metodo de tipo Poisson, y la tarjeta se discretiza en una rejilla 
     de puntos 2D.
     
     Entrada: card > la definicion de la tarjeta y las configuraciones a simular
     Salida: la mejor configuracion y la temperatura media
	      card_s.chips: situacion termica inicial
        card_s.res: la situacion termica final
 
     defines.h: definiciones de ciertas variables y estructuras de datos

     Compilar con estos dos ficheros: 
       diffusion.c: insertar calor, difundir y calcular la temperatura media hasta que se estabilice
       faux.c: ciertas funciones auxiliares

************************************************************************************************/

#include <stdio.h>
#include <values.h>
#include <time.h>
#include <mpi.h>

#include "defines.h"
#include "faux_p.h"
#include "diffusion_p.h"



/************************************************************************************/
void init_grid_chips (int conf, struct info_param param, struct info_chips *chips, int **chip_coord, float *grid_chips)
{
  int i, j, n;

  for (i=0; i<NROW; i++)
  for (j=0; j<NCOL; j++)  
    grid_chips[i*NCOL+j] = param.t_ext;

  for (n=0; n<param.nchip; n++)
  for (i = chip_coord[conf][2*n]   * param.scale; i < (chip_coord[conf][2*n] + chips[n].h) * param.scale; i++)
  for (j = chip_coord[conf][2*n+1] * param.scale; j < (chip_coord[conf][2*n+1]+chips[n].w) * param.scale; j++) 
    grid_chips[(i+1)*NCOL+(j+1)] = chips[n].tchip;
}

/************************************************************************************/
void init_grids (struct info_param param, float *grid, float *grid_aux)
{
  int i, j;

  for (i=0; i<NROW; i++)
  for (j=0; j<NCOL; j++) 
    grid[i*NCOL+j] = grid_aux[i*NCOL+j] = param.t_ext;
}

//Declaradas fuera para que se puedan acceder desde las funciones
int pid, npr;

/************************************************************************************/
/************************************************************************************/
int main (int argc, char *argv[])
{
  struct info_param param;
  struct info_chips *chips;
  int	 **chip_coord;

  float *grid, *grid_chips, *grid_aux;  
  struct info_results BT;
  
  int    conf, i, j, k;
  struct timespec t0, t1;
  double tej, Tmean;

  //Arrays para el reparto de filas
  int N;
  int *tam, *dis, *NROW_loc;
  //Para que cada proceso reciba sus trozos en el scatterv
  float *trozo, *trozo_chips, *trozo_aux;


  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &pid);
  MPI_Comm_size (MPI_COMM_WORLD, &npr);


  tam = malloc(npr * sizeof(int));
  dis = malloc(npr * sizeof(int));
  NROW_loc = malloc(npr * sizeof(int));

 // reading initial data file
  if (argc != 2) {
    printf ("\n\nERROR: needs a card description file \n\n");
    exit (-1);
  } 

  read_data (argv[1], &param, &chips, &chip_coord);

  //Numero de filas a repartir
  //int N = param.scale * 100;

  printf ("\n  ===================================================================");
  printf ("\n    Thermal diffusion - SERIAL version ");
  printf ("\n    %d x %d points, %d chips", RSIZE*param.scale, CSIZE*param.scale, param.nchip);
  printf ("\n    T_ext = %1.1f, Tmax_chip = %1.1f, T_delta: %1.3f, Max_iter: %d", param.t_ext, param.tmax_chip, param.t_delta, param.max_iter);
  printf ("\n  ===================================================================\n\n");
  
  if (pid==0) clock_gettime (CLOCK_REALTIME, &t0);

  grid = malloc(NROW*NCOL * sizeof(float));
  grid_chips = malloc(NROW*NCOL * sizeof(float));
  grid_aux = malloc(NROW*NCOL * sizeof(float));

  BT.bgrid = malloc(NROW*NCOL * sizeof(float));
  BT.cgrid = malloc(NROW*NCOL * sizeof(float));
  BT.Tmean = MAXDOUBLE;

  // loop to process chip configurations
  for (conf=0; conf<param.nconf; conf++)
  {
    // printf("holaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa");
    if (pid == 0){
      // inintial values for grids
      init_grid_chips (conf, param, chips, chip_coord, grid_chips);
      init_grids (param, grid, grid_aux);
    }

    //Repartimos las filas
    for (i=0; i<npr;i++){
      NROW_loc[i]= (i+1)*(NROW/npr) - (i*NROW/npr) + 1;
      if (i!=0 && i!=npr-1) NROW_loc[i]+=1;
      tam[i] = NROW_loc[i]*NCOL; //En el Scatterv no se puede hacer porque es un array
      if (i == 0) dis[i] = 0;
      else{
        dis[i] = dis[i-1] + tam[i-1] - 1;
      }
    }

    trozo = malloc(tam[pid]*NCOL*sizeof(float));
    trozo_chips = malloc(tam[pid]*NCOL*sizeof(float));
    trozo_aux = malloc(tam[pid]*NCOL*sizeof(float));

    for (i=0; i< tam[pid]* NCOL; i++){
      trozo[i] = trozo_aux[i] = param.t_ext;
    } 

    MPI_Scatterv(grid_chips, tam, dis, MPI_FLOAT, trozo_chips, tam[pid], MPI_FLOAT, 0, MPI_COMM_WORLD);


    // main loop: thermal injection/disipation until convergence (t_delta or max_iter)
    Tmean = calculate_Tmean (param, trozo, trozo_chips, trozo_aux, NROW_loc);
    if (pid==0) printf ("  Config: %2d    Tmean: %1.2f\n", conf + 1, Tmean);

    // //Transformamos tamaños y distancias para recibir
    // for (i=0;i<npr;i++){
    //   if (i==0){
    //     tam[i] = tam[i]-NCOL;
    //   }
    //   else if (i==npr){
    //     tam[i] = tam[i]-NCOL;
    //     dis[i] = dis[i-1]-NCOL;
    //   }
    //   else{
    //     tam[i] = tam[i]-(2*NCOL);
    //     dis[i]= dis[i-1] + tam[i-1];
    //   }
    // }

    //Recibimos cada trozo a grid y grid_chips
    // MPI_Gatherv(trozo, tam[pid], MPI_FLOAT, grid, tam, dis, MPI_FLOAT, 0, MPI_COMM_WORLD);
    // MPI_Gatherv(trozo_chips, tam[pid], MPI_FLOAT, grid_chips, tam, dis, MPI_FLOAT, 0, MPI_COMM_WORLD);

    // processing configuration results 
    if (pid==0) results_conf (conf, Tmean, param, grid, grid_chips, &BT);
  }

  if (pid==0){
    clock_gettime (CLOCK_REALTIME, &t1);
    tej = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec)/(double)1e9;
    printf ("\n\n >>> Best configuration: %2d    Tmean: %1.2f\n", BT.conf + 1, BT.Tmean); 
    printf ("   > Time (serial): %1.3f s \n\n", tej);
    // writing best configuration results
    results (param, &BT, argv[1]);
  }

  free (grid);free (grid_chips);free (grid_aux);
  free (BT.bgrid);free (BT.cgrid);
  free (chips);
  for (i=0; i<param.nconf; i++) free (chip_coord[i]);
  free (chip_coord);

  MPI_Finalize ();
  return (0);
}

