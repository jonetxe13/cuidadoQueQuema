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
#include "faux_pAsinc.h"
#include "diffusion_pAsinc.h"



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

  MPI_Status status;

  MPI_Comm worker_comm = MPI_COMM_NULL; // Comm for the worker group
  int group_rank = MPI_UNDEFINED;       // Rank within the worker group
  int group_size = 0;                   // Size of the worker group
  int group_id = -1;                    // ID of the worker group

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

  int color = MPI_UNDEFINED; // Manager doesn't get a color
  int key = pid;             // Use world rank for ordering within groups

  // Assign colors to workers (ranks 1 to npr-1) based on 3 groups of 17 (51 workers / 3)
  int num_worker_groups = (npr - 1)/8; 
  //int workers_per_group = (npr > 1) ? (npr - 1) / num_worker_groups : 0; // 17 for npr=52
  int workers_per_group = 8; // 17 for npr=52
  // int worker_remainder = (npr > 1) ? (npr - 1) % num_worker_groups : 0; // 0 for npr=52


  if (pid > 0) { // If it's a worker process
    group_id = (pid - 1) / workers_per_group; // World rank 1-17 -> group 0, 18-34 -> group 1, 35-51 -> group 2
    color = group_id;
  }

  MPI_Comm_split(MPI_COMM_WORLD, color, key, &worker_comm);


  tam = (int*) malloc(npr * sizeof(int));
  dis = (int*) malloc(npr * sizeof(int));
  tam_marcos = (int*) malloc(npr * sizeof(int));
  dis_marcos = (int*) malloc(npr * sizeof(int));



  if (pid==0){
    grid = malloc(NROW_glob*NCOL_glob * sizeof(float));
    grid_chips = malloc(NROW_glob*NCOL_glob * sizeof(float));

    BT.bgrid = malloc(NROW_glob*NCOL_glob * sizeof(float));
    BT.cgrid = malloc(NROW_glob*NCOL_glob * sizeof(float));
    BT.Tmean = MAXDOUBLE;
  }


  
  conf = 0;
  int nconf_restantes = param.nconf;
  int tam_grupos = 15;

  if(pid == 0){
    while(nconf_restantes > 0){

      /*
      int resto = (NROW_glob-2) % tam_grupos;
      int cociente = (NROW_glob-2) / tam_grupos;

      for (i=0; i<tam_grupos; i++){
        tam[i] = cociente;
        if (i<resto) tam[i]++;
        tam_marcos[i] = tam[i] + 2;
        tam[i] *= NCOL_glob;
        tam_marcos[i] *= NCOL_glob;
        if (i==0) dis[i] = NCOL_glob;
        else dis[i] = dis[i-1] + tam[i-1];
        dis_marcos[i] = dis[i] - NCOL_glob;
      }
      */
      
    init_grid_chips (conf, param, chips, chip_coord, grid_chips);

      int req = 0;
      MPI_Status status;
      MPI_Recv(&req, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status); //esperar a recivir una solicitud para mandar bloques
      MPI_Send(&conf, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD );
      MPI_Send(grid_chips, NCOL_glob*NROW_glob, MPI_FLOAT, status.MPI_SOURCE, 1, MPI_COMM_WORLD );
      // int grupo = status.source;
      nconf_restantes--;

      conf++;
    }
  }

  if(pid != 0){
    MPI_Comm_rank(worker_comm, &group_rank);
    MPI_Comm_size(worker_comm, &group_size); // Should be 17 for all worker groups

    int conf = 0;
    int solic = 1;

    if(group_rank == 0){
      grid = malloc(NROW_glob*NCOL_glob * sizeof(float));
      grid_chips = malloc(NROW_glob*NCOL_glob * sizeof(float));

      MPI_Send(&solic, 1, MPI_INT, 0, 0, MPI_COMM_WORLD); //request
      MPI_Recv(&conf, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status); //esperar a recivir una solicitud para mandar bloques
      MPI_Recv(grid_chips, NCOL_glob*NROW_glob, MPI_FLOAT, 0, 1, MPI_COMM_WORLD, &status); //esperar a recivir grid_chips 
    }
    MPI_Bcast(&conf, 1, MPI_INT, 0, worker_comm);

      int resto = (NROW_glob-2) % group_size;
      int cociente = (NROW_glob-2) / group_size;

      for (i=0; i<group_size; i++){
        tam[i] = cociente;
        if (i<resto) tam[i]++;
        tam_marcos[i] = tam[i] + 2;
        tam[i] *= NCOL_glob;
        tam_marcos[i] *= NCOL_glob;
        if (i==0) dis[i] = NCOL_glob;
        else dis[i] = dis[i-1] + tam[i-1];
        dis_marcos[i] = dis[i] - NCOL_glob;
      }

  trozo = (float*) malloc(tam_marcos[group_rank]*sizeof(float));
  trozo_chips = (float*) malloc(tam_marcos[group_rank]*sizeof(float));
  trozo_aux = (float*) malloc(tam_marcos[group_rank]*sizeof(float));



    MPI_Scatterv(&grid_chips[NCOL_glob], tam, dis_marcos, MPI_FLOAT, &trozo_chips[NCOL_glob], tam[group_rank], MPI_FLOAT, 0, worker_comm);

    //init_grid_chips (conf, param, chips, chip_coord, grid_chips);
    init_grids (t_ext, trozo, trozo_aux, tam_marcos[group_rank]/NCOL_glob, NCOL_glob);

    // main loop: thermal injection/disipation until convergence (t_delta or max_iter)
    Tmean = calculate_Tmean (trozo, trozo_chips, trozo_aux, t_delta, max_iter, t_ext, tam[group_rank]/NCOL_glob, NROW_glob, NCOL_glob, group_rank, group_size, worker_comm);


    //Recibimos cada trozo a grid y grid_chips
    MPI_Gatherv(&trozo[NCOL_glob], tam[group_rank], MPI_FLOAT, grid, tam, dis, MPI_FLOAT, 0, worker_comm);

    if (group_rank == 0){
      MPI_Send(&conf, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
      MPI_Send(grid, NCOL_glob*NROW_glob, MPI_FLOAT, 0, 3, MPI_COMM_WORLD);
      MPI_Send(&Tmean, 1, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD);
    }

    //printf ("  Config: %2d    Tmean: %1.2f\n", conf + 1, Tmean);
    // processing configuration results 
    //results_conf (conf, Tmean, param, grid, grid_chips, &BT);
  }

  if (pid == 0){
    for (i=0; i<nconf; i++){
      MPI_Recv(&conf, 1, MPI_INT, MPI_ANY_SOURCE, 2, MPI_COMM_WORLD, &status);
      MPI_Recv(grid, NCOL_glob*NROW_glob, MPI_FLOAT, status.MPI_SOURCE, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&Tmean, 1, MPI_DOUBLE, status.MPI_SOURCE, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      init_grid_chips (conf, param, chips, chip_coord, grid_chips);
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

  if (group_rank==0){
    free (grid);free (grid_chips);
  }

  free(trozo); free(trozo_chips); free(trozo_aux);
  free(tam); free(dis); free(tam_marcos); free(dis_marcos);
  free(pack_read_data);

  MPI_Finalize ();
  return (0);
}

