#include <mpi.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <malloc.h>

#include "function.hpp"
using namespace std;

MPI_Status status;
int PR = 4;
int tag = 100;
double *U, *U0;
int cond_init;
double sx,sy,dx,dy,dt,Lx,Ly,D,eps1,t,tf, q;
int i, me, Np, statinfo, i1, i12, iN, iN2, j_1, jN, Ny, Nx, n, k, kmax;


void Charge2(int n,
            int Np,
            int me,
            int* i1,
            int* iN,
            double* q2){

  int q = n / Np;
  int r = n % Np;

  if(r == 0){
    *i1 = q * me;
    *iN = q * (me + 1) - 1;
  }else{
    if( me < r){
      *i1 = me * (q + 1);
      *iN = (me + 1) * (q + 1) - 1;
    }else{
      *i1 = me * q + r;
      *iN = (me + 1) * q + r - 1;
    }
  }

  // q est le nombre de lignes connues par le proc
  *q2 = (*iN) - (*i1) + 1;

  // On travaille sur n = nx*ny éléments
  *i1 = (*i1) * Nx;
  *iN = (*iN + 1) * Nx - 1;
}

void Charge_part_domaine(int n,
            int Np,
            int me,
            int* i1,
            int* i12,
            int* iN,
            int* iN2,
            int h){

  int q = n / Np;
  int r = n % Np;

  if(r == 0){
    *i1 = q * me;
    *iN = q * (me + 1) - 1;
  }else{
    if( me < r){
      *i1 = me * (q + 1);
      *iN = (me + 1) * (q + 1) - 1;
    }else{
      *i1 = me * q + r;
      *iN = (me + 1) * q + r - 1;
    }
  }

  // q est le nombre de lignes connues par le proc
  // *q2 = (*iN) - (*i1) + 1;

  // On travaille sur n = nx*ny éléments
  *i1 = (*i1) * Nx;
  *iN = (*iN + 1) * Nx - 1;

  if (*i1 == 0){
    i12 = 0;
    iN2 = iN + Nx*h;
  }else if (*iN == Nx*Ny-1){
    i12 = i1 - Nx*h;
    iN2 = iN;
  }else{
    i12 = i1 - Nx*h;
    iN2 = iN + Nx*h;
  }
}


void Init(DataFile* dataFile){
  int i1, iN, i12, iN2;
  int h = 2;
  double q;

  Nx = dataFile->Get_Nx();
  Ny = dataFile->Get_Ny();
  Lx = dataFile->Get_Lx();
  Ly = dataFile->Get_Ly();
  D = dataFile->Get_D();
  dt = dataFile->Get_dt();
  eps1 = dataFile->Get_epsilon();
  kmax = dataFile->Get_kmax();
  tf = dataFile->Get_tf();
  cond_init = dataFile->Get_cond_init();

  dx = Lx/(Nx+1);
  dy = Ly/(Ny+1);
  sx = D*dt/(dx*dx);
  sy = D*dt/(dy*dy);

  Charge_part_domaine(Ny,Np,me,&i1,&i12,&iN,&iN2,h);

  U = (double*) calloc((iN2-i12+1), sizeof(double));
  U0 = (double*) calloc((iN2-i12+1), sizeof(double));
  for(i = 0; i < iN2-i12+1; i++){
    U0[i] = 1.;
  }
}



double f(double x, double y, double t){
  if(cond_init == 3){
    return exp(-((x - Lx/2) * (x - Lx/2))) * exp(-((y - Ly/2) * (y - Ly/2))) * cos(t*M_PI/2);
  }else if(cond_init == 2){
    return sin(x) + cos(y);
  }else{
    return 2 * (y - y * y + x - x * x);
  }
}


double g(double x, double y, double t){
  if(cond_init == 2){
    return sin(x) + cos(y);
  }else{
    return 0;
  }
}

double h(double x, double y, double t){
  if(cond_init == 3){
    return 1.;
  }else if(cond_init == 2){
    return sin(x) + cos(y);
  }else{
    return 0;
  }
}



void ProdMatVect(double D1[],
                 double D2_m[],
                 double D2_p[],
                 double D3_m[],
                 double D3_p[],
                 double x[],
                 double y[],
                 int i1,
                 int iN){



  int i = 0;



  if(i1 == 0){

    double* x1 = (double*)malloc(sizeof(double) * (iN + Nx + 1));
    for(i = 0; i < iN + 1; i++){
      x1[i] = x[i];
    }

    double q = 0;
    printf("[%s:%d] %d %d\n", __FILE__, __LINE__, i1, iN);
    MPI_Recv(x1 + iN + 1, Nx, MPI_DOUBLE, 1, tag, MPI_COMM_WORLD, &status);
    printf("%s : %d\n", __FILE__, __LINE__);

    MPI_Send(x1 + iN - Nx + 1, Nx, MPI_DOUBLE, 1, tag, MPI_COMM_WORLD);
    printf("%s : %d\n", __FILE__, __LINE__);

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
    printf("%s : %d\n", __FILE__, __LINE__);

    MPI_Send(x3 + Nx, Nx, MPI_DOUBLE, me-1, tag, MPI_COMM_WORLD);
    printf("%s : %d\n", __FILE__, __LINE__);

    MPI_Recv(x3, Nx, MPI_DOUBLE, me-1, tag, MPI_COMM_WORLD, &status);
    printf("%s : %d\n", __FILE__, __LINE__);


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
    printf("%s : %d\n", __FILE__, __LINE__);

    MPI_Send(x2 + Nx, Nx, MPI_DOUBLE, me-1, tag, MPI_COMM_WORLD);
    printf("%s : %d\n", __FILE__, __LINE__);

    MPI_Recv(x2, Nx, MPI_DOUBLE, me-1, tag, MPI_COMM_WORLD, &status);

    printf("%s : %d\n", __FILE__, __LINE__);


    MPI_Recv(x2 + iN - i1 + 1 + Nx, Nx, MPI_DOUBLE, me+1, tag, MPI_COMM_WORLD, &status);
    printf("%s : %d\n", __FILE__, __LINE__);

    MPI_Send(x2 + iN - i1 + 1 , Nx, MPI_DOUBLE, me+1, tag, MPI_COMM_WORLD);
    printf("%s : %d\n", __FILE__, __LINE__);

    for(i = Nx; i < Nx + iN - i1 + 1; i++){
      y[i - Nx] = D3_m[i-Nx]*x2[i-Nx] + D2_m[i-Nx]*x2[i-1] + D1[i-Nx]*x2[i]+D2_p[i-Nx]*x2[i+1] + D3_p[i-Nx]*x2[i+Nx];
    }

    free(x2);
  }

}

double prodscal(const double X[],
                const double Y[],
                const int i1,
                const int iN){

  int i;
  double result = 0, somme = 0;
  for (i = 0; i < iN - i1 + 1; i++) {
    result += X[i] * Y[i];
  }
  printf("%s : %d\n", __FILE__, __LINE__);

  MPI_Allreduce(&result, &somme, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  printf("%s : %d\n", __FILE__, __LINE__);

  return somme;
}

int Reste(int k, int Nx){
  int i = (k+1) % Nx;
  if(i == 0){
    return Nx  ;
  }
  return i;
}

void grad_conj(double D1[],
                 double D2_m[],
                 double D2_p[],
                 double D3_m[],
                 double D3_p[],
                 double x[],
                 double b[],
                 int i1,
                 int iN){



  int i = 0;
  double* r = (double*)malloc(sizeof(double) * (iN - i1 + 1));
  double* p = (double*)malloc(sizeof(double) * (iN - i1 + 1));
  double* r2 = (double*)malloc(sizeof(double) * (iN - i1 + 1));
  double* z = (double*)malloc(sizeof(double) * (iN - i1 + 1));
  double* y = (double*)malloc(sizeof(double) * (iN - i1 + 1));

  for(i = 0; i < iN - i1 + 1; i++)
    x[i] = 293.0;
  ProdMatVect(D1,D2_m,D2_p,D3_m,D3_p,x,y,i1,iN);

  for(i = 0; i < iN - i1 + 1; i++){

    r[i] = b[i] - y[i];
    p[i] = r[i];

  }
  double beta = sqrt(prodscal(r,r,i1,iN));
  int k = 0;
  double alpha = 0, gamma = 0;
  while (beta > eps1 && k <= kmax){
    ProdMatVect(D1,D2_m,D2_p,D3_m,D3_p,p,z,i1,iN);
    alpha = prodscal(r,r,i1,iN)/(prodscal(z,p,i1,iN));
    for(i = 0; i < iN - i1 + 1; i++){
      x[i] = x[i] + alpha * p[i];
      r2[i] = r[i];
      r[i] = r[i] - alpha * z[i];
    }
    gamma = prodscal(r,r,i1,iN)/prodscal(r2,r2,i1,iN);

    for(i = 0; i < iN - i1 + 1; i++){
      p[i] = r[i] + gamma * p[i];
    }
    beta = sqrt(prodscal(r2,r2,i1,iN));
    k++;
  }

  delete[] r;
  delete[] p;
  delete[] r2;
  delete[] z;
  delete[] y;
}

int cmp_vect(const double x[],
             const double y[],
             const int Nx,
             int *test)    // ajouter epsilon (précision)
{
  *test = 1;
  int i = 0;
  while (i < Nx && fabs(x[i] - y[i]) > 0.0001 )  // epsilon = 000.00001
    i += 1;

  if ( i == Nx)
    *test = 0;
  return *test;
}