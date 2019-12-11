
#include <mpi.h>
#include <math.h>
#include <iostream>
#include <fstream>

#include "function.hpp"
#include "matrices.hpp"
#include "DataFile.hpp"

using namespace std;


int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);

  if (argc < 2)
  {
    cout << "Please, enter the name of your data file." << endl;
    abort();
  }

  DataFile data_file(argv[1]);
  // Lecture du fichier de données
  data_file.ReadDataFile();


  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &me);
  // Get the number of processes
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &Np);
  // Get the rank of the process

  double t = 0;
  Init(&data_file);
  // !!$  ! On mesure le temps de calcul par proc, premier compteur t1
  // !!$  call CPU_TIME(t1)

  // ! Répartition des procs
  Charge2(Ny,Np,me,&i1,&iN,&q);
  int nb_per_proc = iN-i1+1;

  double* F = (double*)calloc(sizeof(double), nb_per_proc);
  double* D1 = (double*)calloc(sizeof(double), nb_per_proc);
  double* D2_m = (double*)calloc(sizeof(double), nb_per_proc);
  double* D2_p = (double*)calloc(sizeof(double), nb_per_proc);
  double* D3_m = (double*)calloc(sizeof(double), nb_per_proc);
  double* D3_p = (double*)calloc(sizeof(double), nb_per_proc);

  MatriceDF(D1,D2_m,D2_p,D3_m,D3_p,sx,sy,i1,iN, &data_file);

  int k = 0;
   /////////////////////////// BOUCLE EN TEMPS /////////////////////////
  while (t<tf)
  {
    sec_membre(dx, dy, F,t,i1,iN,&data_file);
    int toto = 0;
    for(toto = 0; toto < nb_per_proc; toto++){
      F[toto] += U0[toto];
    }
    grad_conj(D1,D2_m,D2_p,D3_m,D3_p,U, F,i1,iN);
    for (i = 0; i < nb_per_proc; i++)
      U0[i] = U[i];
    t += dt;
    k += 1;
    if (k%50 == 0){
      printf("%d %f\n",k, U0[4]);
    }
  }


  // Ecriture de la solution sur des fichiers .dat pour chaque proc

  ofstream file;
  ofstream plot;

  char file_name [50];
  char place_holder [50];
  char plot_name [50] = "Result/sol.plot";
  sprintf(file_name, "Result/sol%d.dat", me);

  file.open(file_name, ios::out);

  for(int k = 0; k < iN - i1 + 1; k++){
    file << Reste(k,Nx)*dx << " " << dy*(1 + (i1 + k)/Nx)  << " " << U[k] << "\n";
  }
  file.close();

  int cond_init = data_file.Get_cond_init();
  if (me == 0)
  {
    plot.open(plot_name, ios::out);
    plot << "set terminal png" << endl;
    plot << "set output 'Sol.png'" << endl;
    plot << "set ticslevel 0" << endl;
    plot << "splot 'Result/sol0.dat'";
    for (i = 1; i < Np; i++)
    {
      sprintf(place_holder, "'Result/sol%d.dat'", i);
      plot << ", " << place_holder;
    }
    if (cond_init == 1)
      plot << ", x*(1-x)*y*(1-y)" << endl;
    else if (cond_init == 2)
      plot << ", sin(x)+cos(y)" << endl;
    plot.close();
  }

  // !!$  ! Deuxième compteur en temps
  // !!$  call CPU_TIME(t2)
  // !!$
  // !!$  ! Temps de calcul final
  // !!$  temps=t2-t1
  // !!$  print*,'Le temps de calcul du proc',me,' pour résoudre le problème est :',temps,'s'

  MPI_Finalize();

  return 0;
}
