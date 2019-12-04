
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

  if (argc < 2)
  {
    cout << "Please, enter the name of your data file." << endl;
    abort();
  }

  DataFile data_file(argv[1]);
  // Lecture du fichier de données
  data_file.ReadDataFile();



  // ! Definition des variables
  // Real(PR), dimension(:), Allocatable :: Fx,D1,D2_m,D2_p,D3_m,D3_p
  // character(len=50) :: file_name, Me_string
  // Real(PR):: t1,t2,temps
  // ! Mise en place de l'environnement parallele
  // Call MPI_INIT(statinfo)
  // Call MPI_COMM_RANK(MPI_COMM_WORLD,me,statinfo)
  // Call MPI_COMM_SIZE(MPI_COMM_WORLD,Np,statinfo)
  MPI_Init(NULL, NULL);
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &me);
  // Get the number of processes
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &Np);
  // Get the rank of the process




  double t = 0;
  int cond_init = 0;
  // Choix condition initiale

  Init(&data_file);

  if (me == 0)
  {
    cout << " " << endl;
    cout << "\n Choisissez le jeu de conditions initiales (1,2 ou 3)" << endl;
    cout << "1 : f = 2(y - y^2 + x - x^2) ; g = 0 : h = 0" << endl;
    cout << "2 : f = sin(x) + cos(y) ; g = sin(x) + cos(y) : h = sin(x) + cos(y)" << endl;
    cout << "3 : f = exp(-(x - Lx/2)^2).exp(-(y - Ly/2)^2).cos(t.pi/2) ; g = 0 : h = 1" << endl;

    //cin >> cond_init;
    cond_init = 1;
    if (cond_init > 3 || cond_init < 1)
      cout <<  "Attention ! Le set (1) a ete pris par default" << endl;

  }

  // Partage du choix des conditions initiales aux autres processeurs
  // call MPI_BCAST(cond_init,1,MPI_INTEGER,0,MPI_COMM_WORLD,statinfo)
  MPI_Bcast(&cond_init, 1, MPI_INT, 0, MPI_COMM_WORLD);

  // !!$  ! On mesure le temps de calcul par proc, premier compteur t1
  // !!$  call CPU_TIME(t1)

  // ! Répartition des procs
  Charge2(Ny,Np,me,&i1,&iN,&q);
  cout << me << " " << i1 << " " << iN << " " << endl;
  int nb_per_proc = iN-i1+1;

  // Allocate(Fx(i1:iN),D1(i1:iN),D2_m(i1:iN),D2_p(i1:iN),D3_m(i1:iN),D3_p(i1:iN))
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
    if (k%100 == 0)
      printf("%d\n",k);
  }



  // Ecriture de la solution sur des fichiers .dat pour chaque proc

  // A voir plus tard
  // Write( Me_string, '(i10)' )  Me
  // ofstream s;
  ofstream file;
  //
  // s << "Result/res";
  // for(int boucle = 0; boucle < (log(Np) - log(me))/log(10); boucle++){
  //     s << '0';
  // }
  // s << me << ".data";
  char file_name [50];
  sprintf(file_name, "Result/sol%d.dat", me);
  //string name(s);

  if(me == 0){

    file.open(file_name, ios::out);
  }else{
    file.open(file_name, ios::out);
  }
  for(int i = 0; i < iN - i1+1; i++){
    file << dx << " " << dy << " " << U[i] << "\n";
  }
  file.close();
  // file_name = 'sol00' // trim(adjustl(Me_string)) // '.dat'
  // Open(10+Me, File=trim(file_name))
  // Do k=i1,iN
  // Write(10+Me,*) Reste(k,Nx)*dx, ((k-1)/Nx+1)*dy, U(k), k
  // End Do
  // close(10+Me)

  // !!$  ! Deuxième compteur en temps
  // !!$  call CPU_TIME(t2)
  // !!$
  // !!$  ! Temps de calcul final
  // !!$  temps=t2-t1
  // !!$  print*,'Le temps de calcul du proc',me,' pour résoudre le problème est :',temps,'s'

  MPI_Finalize();

  return 0;
}
