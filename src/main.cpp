
#include <mpi.h>
#include <math.h>
#include <iostream>
#include <fstream>

#include "function_2.hpp"
#include "matrices_2.hpp"
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




  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &me);
  // Get the number of processes
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &Np);
  // Get the rank of the process
  int choice = 0;
  choice = 1;
  if (choice == 0){

    DataFile data_file(argv[1]);
    // Lecture du fichier de données
    data_file.ReadDataFile();

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
  }

  else if (choice == 1)
  {
    int k = 0, l = 0;
    int i12, iN2;



    DataFile data_file(argv[1]);
    // Lecture du fichier de données
    data_file.ReadDataFile();

    double t = 0;
    Init(&data_file);
    // !!$  ! On mesure le temps de calcul par proc, premier compteur t1
    // !!$  call CPU_TIME(t1)

    int h = 2;
    // ! Répartition des procs
    Charge_part_domaine(Ny,Np,me,&i1,&i12,&iN,&iN2,h);
    int nb_per_proc = iN2-i12+1;
    printf("0 me : %d i12 : %d iN2 : %d ; l = %d\n ",me,i12,iN2,l);

    double* F = (double*)calloc(nb_per_proc, sizeof(double) );
    double* D1 = (double*)calloc(nb_per_proc, sizeof(double) );
    double* D2_m = (double*)calloc(nb_per_proc, sizeof(double) );
    double* D2_p = (double*)calloc(nb_per_proc, sizeof(double) );
    double* D3_m = (double*)calloc(nb_per_proc, sizeof(double) );
    double* D3_p = (double*)calloc(nb_per_proc, sizeof(double) );
    double* comp_1 = (double*)calloc(Nx, sizeof(double) );
    double* comp_2 = (double*)calloc(Nx, sizeof(double) );
    double* place_holder_1 = (double*)calloc(Nx, sizeof(double) );
    double* place_holder_2 = (double*)calloc(Nx, sizeof(double) );

    printf("01 me : %d i12 : %d iN2 : %d ; l = %d\n ",me,i12,iN2,l);
    MatriceDF(D1,D2_m,D2_p,D3_m,D3_p,sx,sy,i12,iN2, &data_file);
    printf("02 me : %d i12 : %d iN2 : %d ; l = %d\n ",me,i12,iN2,l);

    printf("03 me : %d i12 : %d iN2 : %d ; l = %d\n ",me,i12,iN2,l);
    // abort();
    // abort();


    /////////////////////////// BOUCLE EN TEMPS /////////////////////////
    while (t<tf)
    {
      printf("1 me : %d i12 : %d iN2 : %d ; function = %d ; line = %d : file = %d\n ",me,i12,iN2,__FUNC__,__LINE__);

      int test1 = 1;
      int test2 = 1;
      while (test1 == 1 || test2 == 1)
      {
        printf("2 me : %d i12 : %d iN2 : %d ; l = %d\n ",me,i12,iN2,l);

        sec_membre(dx, dy, F,t,i12,iN2,&data_file);
        printf("2.5 me : %d i12 : %d iN2 : %d ; l = %d\n ",me,i12,iN2,l);

        int toto = 0;
        for(toto = 0; toto < nb_per_proc; toto++){
          F[toto] += U0[toto];
        }
        grad_conj(D1,D2_m,D2_p,D3_m,D3_p,U, F,i12,iN2);
        for (i = 0; i < nb_per_proc; i++)
        U0[i] = U[i];

        printf(" 3 me : %d i12 : %d iN2 : %d ; l = %d\n ",me,i12,iN2,l);

        if (i12 == 0){
          MPI_Recv(comp_2, Nx, MPI_DOUBLE, me+1, tag, MPI_COMM_WORLD, &status);
          MPI_Send(U + iN2 - i12 - Nx*(h + 1) , Nx, MPI_DOUBLE, me+1, tag, MPI_COMM_WORLD);
          test1 = 0;
          for (l = iN2 - i12 - Nx*(h+1) + 1 ; l < iN2 - i12 - Nx*h ; l++)
            place_holder_2[l] = U[l];
          cmp_vect( place_holder_2, comp_2, Nx, &test2);

        }else if (iN2 == Nx*Ny - 1){
          MPI_Send(U + Nx*(h-1), Nx, MPI_DOUBLE, me-1, tag, MPI_COMM_WORLD);
          MPI_Recv(comp_1, Nx, MPI_DOUBLE, me-1, tag, MPI_COMM_WORLD, &status);
          for (l = Nx*h ; l < Nx*(h+1) - 1 ; l++)
            place_holder_1[l] = U[l];
          cmp_vect( place_holder_1, comp_1, Nx, &test1);
          test2 = 0;

        }else{
          MPI_Send(U + Nx*(h-1), Nx, MPI_DOUBLE, me-1, tag, MPI_COMM_WORLD);
          MPI_Recv(comp_1, Nx, MPI_DOUBLE, me-1, tag, MPI_COMM_WORLD, &status);
          MPI_Recv(comp_2, Nx, MPI_DOUBLE, me+1, tag, MPI_COMM_WORLD, &status);
          MPI_Send(U + iN2 - i12 - Nx*(h + 1) , Nx, MPI_DOUBLE, me+1, tag, MPI_COMM_WORLD);

          for (l = Nx*h ; l < Nx*(h+1) - 1 ; l++)
            place_holder_1[l] = U[l];
          cmp_vect( place_holder_1, comp_1, Nx, &test1);

          for (l = iN2 - i12 - Nx*(h+1) + 1 ; l < iN2 - i12 - Nx*h ; l++)
            place_holder_2[l] = U[l];
          cmp_vect( place_holder_2, comp_2, Nx, &test2);
        }

      }
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


  }

  MPI_Finalize();

  return 0;
}
