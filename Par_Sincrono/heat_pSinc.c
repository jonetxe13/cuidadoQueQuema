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
#include "faux_pSinc.h"
#include "diffusion_pSinc.h"



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
void init_grids (float t_ext, float *grid, float *grid_aux, int NROW_glob, int NCOL_glob)
{
  int i, j;

  for (i=0; i<NROW_glob; i++)
    for (j=0; j<NCOL_glob; j++) 
      grid[i*NCOL_glob+j] = grid_aux[i*NCOL_glob+j] = t_ext;
}


/************************************************************************************/
/************************************************************************************/
int main (int argc, char *argv[])
{
  struct info_param param;
  struct info_chips *chips;
  int	 **chip_coord;
  int pid, npr;

  float *grid, *grid_chips;  
  struct info_results BT;

  int    conf, i, j, k;
  double tej, Tmean, t0, t1;

  //Arrays para el reparto de filas
  int N;
  int *tam, *dis, *tam_marcos, *dis_marcos;

  char *pack_read_data;
  int pack_size, pos, max_iter, scale, nconf;
  float t_delta, t_ext;

  //Para que cada proceso reciba sus trozos en el scatterv
  float *trozo, *trozo_chips, *trozo_aux;

  int NROW_glob, NCOL_glob;


  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &pid);
  MPI_Comm_size (MPI_COMM_WORLD, &npr);

  // reading initial data file
  if (argc != 2) {
    printf ("\n\nERROR: needs a card description file \n\n");
    exit (-1);
  } 

  int pack1_size,pack2_size;
  MPI_Pack_size(3,MPI_INT,MPI_COMM_WORLD,&pack1_size);
  MPI_Pack_size(2,MPI_FLOAT,MPI_COMM_WORLD,&pack2_size);
  pack_size = pack1_size + pack2_size;
  pack_read_data = (char*) malloc (pack_size*sizeof(char));
  pos = 0;

  if (pid == 0) 
  {
    read_data (argv[1], &param, &chips, &chip_coord);
    MPI_Pack(&param.nconf,    1,MPI_INT,pack_read_data,pack_size,&pos,MPI_COMM_WORLD);
    MPI_Pack(&param.max_iter, 1,MPI_INT,pack_read_data,pack_size,&pos,MPI_COMM_WORLD);
    MPI_Pack(&param.scale,    1,MPI_INT,pack_read_data,pack_size,&pos,MPI_COMM_WORLD);
    MPI_Pack(&param.t_ext,    1,MPI_FLOAT,pack_read_data,pack_size,&pos,MPI_COMM_WORLD);
    MPI_Pack(&param.t_delta,  1,MPI_FLOAT,pack_read_data,pack_size,&pos,MPI_COMM_WORLD);
    nconf = param.nconf; 
    max_iter = param.max_iter; 
    t_ext = param.t_ext; 
    t_delta = param.t_delta; 
    scale = param.scale;
  }

  MPI_Bcast(pack_read_data,pack_size,MPI_PACKED,0,MPI_COMM_WORLD);

  if (pid != 0)
  {
    MPI_Unpack(pack_read_data,pack_size,&pos,&nconf,1,MPI_INT,MPI_COMM_WORLD);
    MPI_Unpack(pack_read_data,pack_size,&pos,&max_iter,1,MPI_INT,MPI_COMM_WORLD);
    MPI_Unpack(pack_read_data,pack_size,&pos,&scale,1,MPI_INT,MPI_COMM_WORLD);
    MPI_Unpack(pack_read_data,pack_size,&pos,&t_ext,1,MPI_FLOAT,MPI_COMM_WORLD);
    MPI_Unpack(pack_read_data,pack_size,&pos,&t_delta,1,MPI_FLOAT,MPI_COMM_WORLD);
  }


  NROW_glob = (200*scale + 2);  // extended row number
  NCOL_glob = (100*scale + 2);  // extended column number


  if(pid==0){
    printf ("\n  ===================================================================");
    printf ("\n    Thermal diffusion - PAR version ");
    printf ("\n    %d x %d points, %d chips", RSIZE*param.scale, CSIZE*param.scale, param.nchip);
    printf ("\n    T_ext = %1.1f, Tmax_chip = %1.1f, T_delta: %1.3f, Max_iter: %d", param.t_ext, param.tmax_chip, param.t_delta, param.max_iter);
    printf ("\n  ===================================================================\n\n");
    t0 = MPI_Wtime();
  }


  tam = (int*) malloc(npr * sizeof(int));
  dis = (int*) malloc(npr * sizeof(int));
  tam_marcos = (int*) malloc(npr * sizeof(int));
  dis_marcos = (int*) malloc(npr * sizeof(int));



  int resto = (NROW_glob-2) % npr;
  int cociente = (NROW_glob-2) / npr;

  for (i=0; i<npr; i++){
    tam[i] = cociente;
    if (i<resto) tam[i]++;
    tam_marcos[i] = tam[i] + 2;
    tam[i] *= NCOL_glob;
    tam_marcos[i] *= NCOL_glob;
    if (i==0) dis[i] = NCOL_glob;
    else dis[i] = dis[i-1] + tam[i-1];
    dis_marcos[i] = dis[i] - NCOL_glob;
  }

  if (pid==0){
    grid = malloc(NROW_glob*NCOL_glob * sizeof(float));
    grid_chips = malloc(NROW_glob*NCOL_glob * sizeof(float));

    BT.bgrid = malloc(NROW_glob*NCOL_glob * sizeof(float));
    BT.cgrid = malloc(NROW_glob*NCOL_glob * sizeof(float));
    BT.Tmean = MAXDOUBLE;
  }


  trozo = (float*) malloc(tam_marcos[pid]*sizeof(float));
  trozo_chips = (float*) malloc(tam_marcos[pid]*sizeof(float));
  trozo_aux = (float*) malloc(tam_marcos[pid]*sizeof(float));



  // loop to process chip configurations
  for (conf=0; conf<nconf; conf++)
  {

    if (pid == 0){
      // inintial values for grids
      init_grid_chips (conf, param, chips, chip_coord, grid_chips);
    }
    MPI_Scatterv(&grid_chips[NCOL_glob], tam, dis_marcos, MPI_FLOAT, &trozo_chips[NCOL_glob], tam[pid], MPI_FLOAT, 0, MPI_COMM_WORLD);

    init_grids (t_ext, trozo, trozo_aux, tam_marcos[pid]/NCOL_glob, NCOL_glob);

    // main loop: thermal injection/disipation until convergence (t_delta or max_iter)
    Tmean = calculate_Tmean (trozo, trozo_chips, trozo_aux, t_delta, max_iter, t_ext, tam[pid]/NCOL_glob, NROW_glob, NCOL_glob);


    //Recibimos cada trozo a grid y grid_chips
    MPI_Gatherv(&trozo[NCOL_glob], tam[pid], MPI_FLOAT, grid, tam, dis, MPI_FLOAT, 0, MPI_COMM_WORLD);


    if (pid==0) {
      printf ("  Config: %2d    Tmean: %1.2f\n", conf + 1, Tmean);
    // processing configuration results 
      results_conf (conf, Tmean, param, grid, grid_chips, &BT);
    }
  }

  if (pid==0){
    t1 = MPI_Wtime();
    tej = t1-t0;
    printf ("\n\n >>> Best configuration: %2d    Tmean: %1.2f\n", BT.conf + 1, BT.Tmean); 
    printf ("   > Time (par): %1.3f s \n\n", tej);
    // writing best configuration results
    results (param, &BT, argv[1]);
  }

  if (pid==0){
    free (grid);free (grid_chips);
    free (BT.bgrid);free (BT.cgrid);
    free (chips);
    for (i=0; i<param.nconf; i++) free (chip_coord[i]);
    free (chip_coord);
  }

  free(trozo); free(trozo_chips); free(trozo_aux);
  free(tam); free(dis); free(tam_marcos); free(dis_marcos);
  free(pack_read_data);

  MPI_Finalize ();
  return (0);
}

