
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

  int choice = 1;
  int height = 1;
  if (argc < 2)
  {
    cerr << "Please, enter the name of your data file." << endl;
    abort();
  }
  if (argc > 3)
  {
    choice = atoi(argv[2]);
    height = atoi(argv[3]);
  }

  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &me);
  // Get the number of processes
  int world_size;

  // Get the rank of the process
  MPI_Comm_size(MPI_COMM_WORLD, &Np);

  // Lecture du fichier de données
  DataFile data_file(argv[1]);
  data_file.ReadDataFile();
  double t = 0;
  Init(&data_file);


  // Répartition des procs
  int nb_per_proc = -1;
  int i12, iN2;
  if(choice == 1 || choice == 2){
    Charge_part_domaine(Ny,Np,me,&i1,&i12,&iN,&iN2,height);
    nb_per_proc = iN2-i12+1;
  } else {
    Charge2(Ny,Np,me,&i1,&iN,&q);
    nb_per_proc = iN-i1+1;
  }

  U = (double*) calloc(nb_per_proc, sizeof(double));
  U0 = (double*) calloc(nb_per_proc, sizeof(double));

  double* F = (double*)calloc(sizeof(double), nb_per_proc);
  double* D1 = (double*)calloc(sizeof(double), nb_per_proc);
  double* D2_m = (double*)calloc(sizeof(double), nb_per_proc);
  double* D2_p = (double*)calloc(sizeof(double), nb_per_proc);
  double* D3_m = (double*)calloc(sizeof(double), nb_per_proc);
  double* D3_p = (double*)calloc(sizeof(double), nb_per_proc);

  if (choice == 0){

    MatriceDF(D1,D2_m,D2_p,D3_m,D3_p,sx,sy,i1,iN, &data_file);

    int k = 0;
    /////////////////////////// BOUCLE EN TEMPS /////////////////////////
    struct timespec t1, t2;
    clock_gettime(CLOCK_MONOTONIC, &t1);

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

    }

    clock_gettime(CLOCK_MONOTONIC, &t2);
      double error = 0, error_all;
      int i,j;
      for(i = 1; i <= Nx; i++)
        for(j = 1; j <= Ny; j++){
          double d = U[i + j * (Nx-1)] - i*dx*(1-i*dx)*j*dy*(1-j*dy);
          error += d*d;
        }
      MPI_Allreduce(&error , &error_all, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      error_all = sqrt(error_all);
    if(me == 0){
      cout << (t2.tv_sec - t1.tv_sec) * 1000000.0 + (t2.tv_nsec - t1.tv_nsec) / 1000.0 << "," << error_all << endl;
    }

  } else if (choice == 1) {

    int k = 0, l = 0;
    double* comp_1 = (double*)calloc(Nx, sizeof(double) );
    double* comp_2 = (double*)calloc(Nx, sizeof(double) );

    MatriceDF(D1,D2_m,D2_p,D3_m,D3_p,sx,sy,i12,iN2, &data_file);

    double error;

    struct timespec t1, t2;
    clock_gettime(CLOCK_MONOTONIC, &t1);

    /////////////////////////// BOUCLE EN TEMPS /////////////////////////
    while (t<tf) {

      error = 1e307;
      double norm[2] = {error, error};
      int current_norm = 0;
      int prev_norm = 1;
      int offset1 = Nx*height;
      int offset2 = iN-i1+1 - Nx;

      if (i12 == 0){
        MPI_Recv(comp_2, Nx, MPI_DOUBLE, me+1, tag, MPI_COMM_WORLD, &status);
        MPI_Send(U + iN-i1+1 - Nx, Nx, MPI_DOUBLE, me+1, tag, MPI_COMM_WORLD);

      }else if (iN2 == Nx*Ny - 1){
        MPI_Send(U + offset1, Nx, MPI_DOUBLE, me-1, tag, MPI_COMM_WORLD);
        MPI_Recv(comp_1, Nx, MPI_DOUBLE, me-1, tag, MPI_COMM_WORLD, &status);

      }else{
        MPI_Send(U + offset1, Nx, MPI_DOUBLE, me-1, tag, MPI_COMM_WORLD);
        MPI_Recv(comp_1, Nx, MPI_DOUBLE, me-1, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(comp_2, Nx, MPI_DOUBLE, me+1, tag, MPI_COMM_WORLD, &status);
        MPI_Send(U + offset2 , Nx, MPI_DOUBLE, me+1, tag, MPI_COMM_WORLD);

      }

      while (error > 0.00001) {
        int toto = 0;
        for(toto = 0; toto < nb_per_proc; toto++){
          F[toto] += U0[toto];
        }
        grad_conj(D1,D2_m,D2_p,D3_m,D3_p,U, F,i12,iN2);

        if (i12 == 0){
          MPI_Recv(comp_2, Nx, MPI_DOUBLE, me+1, tag, MPI_COMM_WORLD, &status);
          MPI_Send(U + offset2, Nx, MPI_DOUBLE, me+1, tag, MPI_COMM_WORLD);

        }else if (iN2 == Nx*Ny - 1){
          MPI_Send(U + offset1, Nx, MPI_DOUBLE, me-1, tag, MPI_COMM_WORLD);
          MPI_Recv(comp_1, Nx, MPI_DOUBLE, me-1, tag, MPI_COMM_WORLD, &status);

        }else{
          MPI_Send(U + offset1, Nx, MPI_DOUBLE, me-1, tag, MPI_COMM_WORLD);
          MPI_Recv(comp_1, Nx, MPI_DOUBLE, me-1, tag, MPI_COMM_WORLD, &status);
          MPI_Recv(comp_2, Nx, MPI_DOUBLE, me+1, tag, MPI_COMM_WORLD, &status);
          MPI_Send(U + offset2, Nx, MPI_DOUBLE, me+1, tag, MPI_COMM_WORLD);
        }

        update_sec_membre(dx, dy, F,t,i12,iN2,&data_file,comp_1,comp_2);

        double res = 0;
        for(i = 0; i < nb_per_proc; i++) res += U[i] * U[i];
        MPI_Allreduce(&res , norm + current_norm, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        norm[current_norm] = sqrt(norm[current_norm]);
        error = fabs(norm[current_norm] - norm[prev_norm]);
        prev_norm = current_norm;
        current_norm = (current_norm + 1) % 2;
      }

      for (i = 0; i < nb_per_proc; i++){
        U0[i] = U[i];
      }
      t += dt;
      k += 1;

    }

    clock_gettime(CLOCK_MONOTONIC, &t2);
      double error_all;
      error = 0;
      int i,j;
      for(i = 1; i <= Nx; i++)
        for(j = 1; j <= Ny; j++){
          double d = U[i + j * (Nx-1)] - i*dx*(1-i*dx)*j*dy*(1-j*dy);
          error += d*d;
        }
      MPI_Allreduce(&error , &error_all, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      error_all = sqrt(error_all);
    if(me == 0){
      cout << (t2.tv_sec - t1.tv_sec) * 1000000.0 + (t2.tv_nsec - t1.tv_nsec) / 1000.0 << "," << error_all << endl;
    }

    delete[] comp_1;
    delete[] comp_2;
  }

  ////////////////////////////////////////////////////////////////////
  ////////////////////////////DIRCHLET-NEUMANN////////////////////////
  ////////////////////////////////////////////////////////////////////
  else if (choice == 2)
  {

    int k = 0, l = 0;
    int alpha = 1;
    int beta = 1;

    double* comp_1 = (double*)calloc(2*Nx, sizeof(double) );
    double* comp_2 = (double*)calloc(2*Nx, sizeof(double) );

    MatriceDF3(D1,D2_m,D2_p,D3_m,D3_p,sx,sy,alpha, beta, i12,iN2, &data_file);

    double error;

    struct timespec t1, t2;
    clock_gettime(CLOCK_MONOTONIC, &t1);
    /////////////////////////// BOUCLE EN TEMPS /////////////////////////
    while (t<tf) {

      error = 1e307;
      double norm[2] = {error, error};
      int current_norm = 0;
      int prev_norm = 1;
      int offset1 = Nx*height;
      int offset2 = iN-i1+1 - Nx;

      if (i12 == 0){
        MPI_Recv(comp_2, 2*Nx, MPI_DOUBLE, me+1, tag, MPI_COMM_WORLD, &status);
        MPI_Send(U + offset2, 2*Nx, MPI_DOUBLE, me+1, tag, MPI_COMM_WORLD);


      }else if (iN2 == Nx*Ny - 1){
        MPI_Send(U + offset1 - Nx, 2*Nx, MPI_DOUBLE, me-1, tag, MPI_COMM_WORLD);
        MPI_Recv(comp_1, 2*Nx, MPI_DOUBLE, me-1, tag, MPI_COMM_WORLD, &status);

      }else{
        MPI_Send(U + offset1 - Nx, 2*Nx, MPI_DOUBLE, me-1, tag, MPI_COMM_WORLD);
        MPI_Recv(comp_1, 2*Nx, MPI_DOUBLE, me-1, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(comp_2, 2*Nx, MPI_DOUBLE, me+1, tag, MPI_COMM_WORLD, &status);
        MPI_Send(U + offset2, 2*Nx, MPI_DOUBLE, me+1, tag, MPI_COMM_WORLD);

      }


      while (error > 0.00001) {
        int toto = 0;
        for(toto = 0; toto < nb_per_proc; toto++){
          F[toto] += U0[toto];
        }
        grad_conj(D1,D2_m,D2_p,D3_m,D3_p,U, F,i12,iN2);

        if (i12 == 0){
          MPI_Recv(comp_2, 2*Nx, MPI_DOUBLE, me+1, tag, MPI_COMM_WORLD, &status);
          MPI_Send(U + offset2, 2*Nx, MPI_DOUBLE, me+1, tag, MPI_COMM_WORLD);


        }else if (iN2 == Nx*Ny - 1){
          MPI_Send(U + offset1 - Nx, 2*Nx, MPI_DOUBLE, me-1, tag, MPI_COMM_WORLD);
          MPI_Recv(comp_1, 2*Nx, MPI_DOUBLE, me-1, tag, MPI_COMM_WORLD, &status);

        }else{
          MPI_Send(U + offset1 - Nx, 2*Nx, MPI_DOUBLE, me-1, tag, MPI_COMM_WORLD);
          MPI_Recv(comp_1, 2*Nx, MPI_DOUBLE, me-1, tag, MPI_COMM_WORLD, &status);
          MPI_Recv(comp_2, 2*Nx, MPI_DOUBLE, me+1, tag, MPI_COMM_WORLD, &status);
          MPI_Send(U + offset2, 2*Nx, MPI_DOUBLE, me+1, tag, MPI_COMM_WORLD);

        }

        update_sec_membre3(dx, dy, F,t,alpha, beta, i12,iN2,&data_file,comp_1,comp_2);

        double res = 0;
        for(i = 0; i < nb_per_proc; i++) res += U[i] * U[i];
        MPI_Allreduce(&res , norm + current_norm, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        norm[current_norm] = sqrt(norm[current_norm]);
        error = fabs(norm[current_norm] - norm[prev_norm]);
        prev_norm = current_norm;
        current_norm = (current_norm + 1) % 2;
      }

      for (i = 0; i < nb_per_proc; i++){
        U0[i] = U[i];
      }
      t += dt;
      k += 1;

    }

    clock_gettime(CLOCK_MONOTONIC, &t2);
      double error_all;
      error = 0;
      int i,j;
      for(i = 1; i <= Nx; i++)
        for(j = 1; j <= Ny; j++){
          double d = U[i + j * (Nx-1)] - i*dx*(1-i*dx)*j*dy*(1-j*dy);
          error += d*d;
        }
      MPI_Allreduce(&error , &error_all, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      error_all = sqrt(error_all);
    if(me == 0){
      cout << (t2.tv_sec - t1.tv_sec) * 1000000.0 + (t2.tv_nsec - t1.tv_nsec) / 1000.0 << "," << error_all << endl;
    }
    delete[] comp_1;
    delete[] comp_2;
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
  if (me == 0) {
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

  delete[] U0;
  delete[] U;
  delete[] F;
  delete[] D1;
  delete[] D2_m;
  delete[] D2_p;
  delete[] D3_m;
  delete[] D3_p;

  MPI_Finalize();

  return 0;
}
