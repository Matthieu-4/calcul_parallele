#include <mpi.h>
#include <math.h>
#include <iostream>
#include <fstream>

#include "Matrices.hpp"
#include "DataFile.hpp"


Solver::Solver(std::string file_name) {
  data_file = new DataFile(file_name);
  double q;
  MPI_Comm_rank(MPI_COMM_WORLD, &me);
  MPI_Comm_size(MPI_COMM_WORLD, &Np);
  Nx = data_file->Get_Nx();
  Ny = data_file->Get_Ny();
  Lx = data_file->Get_Lx();
  Ly = data_file->Get_Ly();
  D = data_file->Get_D();
  dt = data_file->Get_dt();
  eps1 = data_file->Get_epsilon();
  kmax = data_file->Get_kmax();
  tf = data_file->Get_tf();
  cond_init = data_file->Get_cond_init();

  dx = Lx/(Nx+1);
  dy = Ly/(Ny+1);
  sx = D*dt/(dx*dx);
  sy = D*dt/(dy*dy);

  Charge2(&q);

  U = (double*) calloc((iN-i1+1), sizeof(double));
  U0 = (double*) calloc((iN-i1+1), sizeof(double));
  for(int i = 0; i < iN-i1+1; i++){
    U0[i] = 1.;
  }

  nb_per_proc = iN-i1+1;

  F = (double*)calloc(sizeof(double), nb_per_proc);
  D1 = (double*)calloc(sizeof(double), nb_per_proc);
  D2_m = (double*)calloc(sizeof(double), nb_per_proc);
  D2_p = (double*)calloc(sizeof(double), nb_per_proc);
  D3_m = (double*)calloc(sizeof(double), nb_per_proc);
  D3_p = (double*)calloc(sizeof(double), nb_per_proc);

  this->MatriceDF();
}

Solver::~Solver(void){
  delete[] U;
  delete[] U0;
  delete[] F;
  delete[] D1;
  delete[] D2_m;
  delete[] D2_p;
  delete[] D3_m;
  delete[] D3_p;
}

void Solver::Charge2(double* q2){

  int q = Ny / Np;
  int r = Ny % Np;

  if(r == 0){
    i1 = q * me;
    iN = q * (me + 1) - 1;
  }else{
    if( me < r){
      i1 = me * (q + 1);
      iN = (me + 1) * (q + 1) - 1;
    }else{
      i1 = me * q + r;
      iN = (me + 1) * q + r - 1;
    }
  }

  // q est le nombre de lignes connues par le proc
  *q2 = iN - i1 + 1;

  // On travaille sur n = nx*ny éléments
  i1 = i1 * Nx;
  iN = (iN + 1) * Nx - 1;
}


// Construction de la matrice DF (cas linéaire)
void Solver::MatriceDF(void) {
  int k = 0;
  for (k = 0; k < nb_per_proc; k++){
    D1[k] = 1.0 + 2.0*D*(sx+sy);
    D2_m[k] = -D*sx;
    D2_p[k] = -D*sx;
    D3_m[k] = -D*sy;
    D3_p[k] = -D*sy;
  }

  if (iN == Nx*Ny-1){
    for (k = nb_per_proc - Nx ; k < nb_per_proc; k++)
      D3_p[k] = 0.0 ;
  }

  if (i1 == 0){
    for (k = 0 ; k < Nx; k++)
      D3_m[k] = 0.0;
  }

  for( k = 0; k < nb_per_proc ; k += Nx){
    D2_p[k+Nx-1] = 0.0;
    D2_m[k] = 0.0;
  }
}

  void Solver::ProdMatVect(double x[], double y[]){



    int i = 0;
    if(i1 == 0){

      double* x1 = (double*)malloc(sizeof(double) * (iN + Nx + 1));
      for(i = 0; i < iN + 1; i++){
        x1[i] = x[i];
      }

      //double q = 0;
      MPI_Recv(x1 + iN + 1, Nx, MPI_DOUBLE, 1, tag, MPI_COMM_WORLD, &status);

      MPI_Send(x1 + iN - Nx + 1, Nx, MPI_DOUBLE, 1, tag, MPI_COMM_WORLD);

      y[0] = D1[0] * x1[0] + D2_p[0] * x1[1] + D3_p[0] * x1[Nx];
      for(i = 1; i < Nx; i++){
        y[i] = D2_m[i]*x1[i-1]+D1[i]*x1[i]+D2_p[i]*x1[i+1]+D3_p[i]*x1[i+Nx];
      }
      for(i = Nx; i < iN + 1; i++){
        y[i] = D3_m[i] * x1[i-Nx] + D2_m[i] * x1[i-1] + D1[i] * x1[i] + D2_p[i] * x1[i+1] + D3_p[i] * x1[i+Nx];
      }


      free(x1);

    } else if (iN == Nx*Ny -1){

      double* x3 = (double*)calloc(sizeof(double),iN + Nx - i1 + 1);
      for(i = Nx; i < Nx + iN - i1 + 1; i++){
        x3[i] = x[i-Nx];
      }

      MPI_Send(x3 + Nx, Nx, MPI_DOUBLE, me-1, tag, MPI_COMM_WORLD);
      MPI_Recv(x3, Nx, MPI_DOUBLE, me-1, tag, MPI_COMM_WORLD, &status);


      for(i = Nx; i < iN - i1 + 1; i++){
        y[i - Nx] = D3_m[i-Nx]*x3[i-Nx]+D2_m[i-Nx]*x3[i-1]+D1[i-Nx]*x3[i]+D2_p[i-Nx]*x3[i+1]+D3_p[i-Nx]*x3[i+Nx];
      }

      for(i = iN - i1 + 1; i < iN + Nx - i1 ; i++){
        y[i - Nx] = D3_m[i-Nx]*x3[i-Nx]+D2_m[i-Nx]*x3[i-1]+D1[i-Nx]*x3[i]+D2_p[i-Nx]*x3[i+1];
      }

      y[iN - i1] = D3_m[iN - i1]*x3[iN - i1]+D2_m[iN - i1]*x3[iN + Nx - i1 - 1]+D1[iN - i1]*x3[iN + Nx - i1];


      free(x3);

    } else {

      double* x2 = (double*)calloc(sizeof(double),iN + 2 * Nx - i1 + 1);
      for(i = Nx; i < Nx + iN - i1 + 1; i++){
        x2[i] = x[i-Nx];
      }

      MPI_Send(x2 + Nx, Nx, MPI_DOUBLE, me-1, tag, MPI_COMM_WORLD);

      MPI_Recv(x2, Nx, MPI_DOUBLE, me-1, tag, MPI_COMM_WORLD, &status);



      MPI_Recv(x2 + iN - i1 + 1 + Nx, Nx, MPI_DOUBLE, me+1, tag, MPI_COMM_WORLD, &status);

      MPI_Send(x2 + iN - i1 + 1 , Nx, MPI_DOUBLE, me+1, tag, MPI_COMM_WORLD);

      for(i = Nx; i < Nx + iN - i1 + 1; i++){
        y[i - Nx] = D3_m[i-Nx]*x2[i-Nx] + D2_m[i-Nx]*x2[i-1] + D1[i-Nx]*x2[i]+D2_p[i-Nx]*x2[i+1] + D3_p[i-Nx]*x2[i+Nx];
      }

      free(x2);
    }

  }


void Solver::computeSecondMember(double t) {

      int i = 0, j = 0, k = 0;
      double B, C;

      //A = 1.0+2.0*D*(sx+sy);
      B = -D*sx;
      C = -D*sy;
      double* Fx = F;
      // Premier proc
      if (i1 == 0)
      {
        Fx[0] = dt*f(dx,dy,t) - C*g(dx,0.0,t)-B*h(0.0,dy,t);
        for (k = 1; k < Nx-1; k++)
          Fx[k] = dt*f(Reste(k)*dx,dy,t) - C*g(Reste(k)*dx,0.0,t);


        Fx[Nx-1] = dt*f(Reste(k)*dx,dy,t) - C*g(Lx-dx,0.0,t)-B*h(Lx,dy,t);
        for (k = Nx; k < nb_per_proc; k++)                // k(i,j) = i + Nx*(j-1)
        {
          i = Reste(k);                // i(k) = reste de k/Nx (+ voir fonction Reste)
          j = (i1+k)/Nx + 1.0;           // j(k) = (quotient de k-1 divis� par Nx) + 1
          Fx[k] = dt*f(i*dx,j*dy,t);


          if (k%Nx == 0)
            Fx[k] = Fx[k]-B*h(0.0,j*dy,t);
          else if (k%Nx == Nx - 1)
            Fx[k] = Fx[k]-B*h(Lx,j*dy,t);
        }
      } else if (iN == Nx*Ny-1) { // Dernier proc

        for ( k = 0; k < nb_per_proc - Nx; k++)
        {
          i = Reste(k);
          j = (i1+k)/Nx + 1.0;
          Fx[k] = dt*f(i*dx,j*dy,t);


          if (k%Nx == 0)
            Fx[k] = Fx[k]-B*h(0.0,j*dy,t);
          else if (k%Nx == Nx - 1)
            Fx[k] = Fx[k]-B*h(Lx,j*dy,t);
        }

        Fx[nb_per_proc - Nx] = dt*f(dx,Ly-dy,t)-C*g(dx,Ly,t)-B*h(0.0,Ly-dy,t);


        for (k= nb_per_proc - Nx + 1; k < nb_per_proc - 1; k++)
          Fx[k] = dt*f(Reste(k)*dx,Ly-dy,t) - C*g(Reste(k)*dx,Ly,t);

        Fx[nb_per_proc-1] = dt*f(Lx-dx,Ly-dy,t)-C*g(Lx-dx,Ly,t)-B*h(Lx,Ly-dy,t);

      } else {// Les autres procs
        for (k = 0; k < nb_per_proc; k++)
        {
          i = Reste(k);
          j = (i1+k)/Nx + 1;

          Fx[k] = dt*f(i*dx,j*dy,t);

          if (k%Nx == 0)
            Fx[k] = Fx[k]-B*h(0.0,j*dy,t);
          else if (k%Nx == Nx-1){
            Fx[k] = Fx[k]-B*h(Lx,j*dy,t);
        }
      }
    }
    int toto = 0;
    for(toto = 0; toto < nb_per_proc; toto++){
      F[toto] += U0[toto];
    }
}

double Solver::prodscal(const double X[], const double Y[]){

    int i;
    double result = 0, somme = 0;
    for (i = 0; i < iN - i1 + 1; i++) {
      result += X[i] * Y[i];
    }

    MPI_Allreduce(&result, &somme, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return somme;
  }

int Solver::Reste(int k) const {
  int i = (k+1) % Nx;
  if(i == 0){
    return Nx;
  }
  return i;
}


void Solver::computeConjugateGradient(void){


  double* x = U;
  double* b = F;
  int i = 0;
  double* r = (double*)malloc(sizeof(double) * (iN - i1 + 1));
  double* p = (double*)malloc(sizeof(double) * (iN - i1 + 1));
  double* r2 = (double*)malloc(sizeof(double) * (iN - i1 + 1));
  double* z = (double*)malloc(sizeof(double) * (iN - i1 + 1));
  double* y = (double*)malloc(sizeof(double) * (iN - i1 + 1));

  for(i = 0; i < iN - i1 + 1; i++)
    x[i] = 293.0;
  ProdMatVect(x, y);

  for(i = 0; i < iN - i1 + 1; i++){

    r[i] = b[i] - y[i];
    p[i] = r[i];

  }
  double beta = sqrt(prodscal(r, r));
  int k = 0;
  double alpha = 0, gamma = 0;
  while (beta > eps1 && k <= kmax){
    ProdMatVect(p, z);
    alpha = prodscal(r, r)/(prodscal(z, p));
    for(i = 0; i < iN - i1 + 1; i++){
      x[i] = x[i] + alpha * p[i];
      r2[i] = r[i];
      r[i] = r[i] - alpha * z[i];
    }
    gamma = prodscal(r, r)/prodscal(r2, r2);

    for(i = 0; i < iN - i1 + 1; i++){
      p[i] = r[i] + gamma * p[i];
    }
    beta = sqrt(prodscal(r2, r2));
    k++;
  }

  for (i = 0; i < nb_per_proc; i++) U0[i] = U[i];

  delete[] r;
  delete[] p;
  delete[] r2;
  delete[] z;
  delete[] y;
}


void Solver::saveSolution(const char* file_name) const {
  std::ofstream file;
  file.open(file_name, std::ios::out);

  for(int k = 0; k < iN - i1 + 1; k++){
    file << Reste(k)*dx << " " << dy*(1 + (i1 + k)/Nx)  << " " << U[k] << "\n";
  }
  file.close();
}


void Solver::makePlot(const char* plot_name) const {
  if (me == 0){
    std::ofstream plot;
    char place_holder [50];
    plot.open(plot_name, std::ios::out);
    plot << "set terminal png" << std::endl;
    plot << "set output 'Sol.png'" << std::endl;
    plot << "set ticslevel 0" << std::endl;
    plot << "splot 'Result/sol0.dat'";
    for (int i = 1; i < Np; i++)
    {
      sprintf(place_holder, "'Result/sol%d.dat'", i);
      plot << ", " << place_holder;
    }
    if (cond_init == 1)
      plot << ", x*(1-x)*y*(1-y)" << std::endl;
    else if (cond_init == 2)
      plot << ", sin(x)+cos(y)" << std::endl;
    plot.close();
  }
}


double Solver::f(double x, double y, double t){
  if(cond_init == 3){
    return exp(-((x - Lx/2) * (x - Lx/2))) * exp(-((y - Ly/2) * (y - Ly/2))) * cos(t*M_PI/2);
  }else if(cond_init == 2){
    return sin(x) + cos(y);
  }
  return 2 * (y - y * y + x - x * x);

}


double Solver::g(double x, double y, double t){
  if(cond_init == 2){
    return sin(x) + cos(y);
  }
  return 0;

}

double Solver::h(double x, double y, double t){
  if(cond_init == 3){
    return 1.;
  }else if(cond_init == 2){
    return sin(x) + cos(y);
  }
  return 0;

}
