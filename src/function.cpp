//#include <mpi>
#include <math.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <malloc.h>

#include "DataFile.hpp"

int *Status; // MPI_STATUS_SIZE
int PR = 4;
int tag = 100;
double *U, *U0;
int cond_init;
double sx,sy,dx,dy,dt,Lx,Ly,D,eps1,t,tf, q;
int i, me, Np, statinfo, i1, iN, j_1, jN, Ny, Nx, n, k, kmax;


void Charge2(int n,
            int Np,
            int me,
            int* i1,
            int* iN,
            double* q2){


  int q = n / Np;
  int r = n % Np;

  if(r == 0){
    *i1 = 1 + q * me;
    *iN = q * (me + 1);
  }else{
    if( me < r){
      *i1 = me * (q + 1) + 1;
      *iN = (me + 1) * (q + 1);
    }else{
      *i1 = 1 + me * q + r;
      *iN = (*i1) + q - 1;
    }
  }


  // q est le nombre de lignes connues par le proc
  *q2 = (*iN) - (*i1) + 1;

  // On travaille sur n = nx*ny éléments
  *i1 = ((*i1) - 1) * Nx + 1;
  *iN = (*iN) * Nx;
}


void Init(void){
  int i1, iN;
  double q;
  //&Nx, &Ny, &Lx, &Ly, &D, &dt, &eps1, &kmax, &tf
  DataFile dataFile("data.txt");

  dx = Lx/(Nx+1);
  dy = Ly/(Ny+1);
  sx = D*dt/(dx*dx);
  sy = D*dt/(dy*dy);

  Charge2(Ny,Np,me,&i1,&iN,&q);

  U = (double*) calloc((iN-i1), sizeof(double));
  U0 = (double*) calloc((iN-i1), sizeof(double));
  for(i = 0; i < iN-i1; i++){
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
  (void) t;
  if(cond_init == 2){
    return sin(x) + cos(y);
  }else{
    return 2 * (y - y * y + x - x * x);
  }
}

double h(double x, double y, double t){
  if(cond_init == 3){
    return 1.;
  }else if(cond_init == 2){
    return sin(x) + cos(y);
  }else{
    return 2 * (y - y * y + x - x * x);
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
  if(i1 == 1){

    double* x1 = (double*)malloc(sizeof(double) * (iN + Nx - i1));
    for(i = 0; i < iN - i1; i++){
      x1[i] = x[i];
    }
    double q = 0;
    Charge2(Ny, Np, me+1, &j_1, &jN, &q);
    //MPI_RECV(x1(j_1:j_1+Nx-1),Nx,MPI_REAL,1,tag,MPI_COMM_WORLD,status,statinfo)
    //MPI_SEND(x1(iN-Nx+1:iN),Nx,MPI_REAL,1,tag,MPI_COMM_WORLD,statinfo)

    y[0] = D1[0] * x1[0] + D2_p[0] * x1[1] + D3_p[0] * x1[Nx];
    for(i = 1; i < iN; i++){
      y[i] = D2_m[i]*x1[i-1]+D1[i]*x1[i]+D2_p[i]*x1[i+1]+D3_p[i]*x1[i+Nx];
    }
    for(i = Nx; i < iN; i++){
      y[i] = D3_m[i] * x1[i-Nx] + D2_m[i] * x1[i-1] + D1[i] * x1[i] + D2_p[i] * x1[i+1] + D3_p[i] * x1[i+Nx];
    }
    delete[] x1;
  } else if (iN == Nx*Ny){

    double* x3 = (double*)malloc(sizeof(double) * (iN + Nx - i1));
    for(i = 0; i < iN - i1; i++){
      x3[i] = x[i];
    }
    //MPI_SEND(x3(i1:i1+Nx-1),Nx,MPI_REAL,me-1,tag,MPI_COMM_WORLD,statinfo)
    Charge2(Ny, Np, me-1, &j_1, &jN, &q);
    //MPI_RECV(x3(jN-Nx+1:jN),Nx,MPI_REAL,me-1,tag,MPI_COMM_WORLD,status,statinfo)
    for(i = 0; i < iN - Nx - i1; i++){
      y[i] = D3_m[i]*x3[i-Nx]+D2_m[i]*x3[i-1]+D1[i]*x3[i]+D2_p[i]*x3[i+1]+D3_p[i]*x3[i+Nx];
    }
    for(i = iN - Nx + 1 - i1; i < iN - i1 - 1; i++){
      y[i] = D3_m[i]*x3[i-Nx]+D2_m[i]*x3[i-1]+D1[i]*x3[i]+D2_p[i]*x3[i+1];
    }
    y[Nx*Ny] = D3_m[Nx*Ny]*x3[Nx*Ny-Nx]+D2_m[Nx*Ny]*x3[Nx*Ny-1]+D1[Nx*Ny]*x3[Nx*Ny];
    delete[] x3;
  } else {

    double* x2 = (double*)malloc(sizeof(double) * (iN + 2 * Nx - i1));
    for(i = 0; i < iN - i1; i++){
      x2[i] = x[i];
    }

    //MPI_SEND(x2(i1:i1+Nx-1),Nx,MPI_REAL,me-1,tag,MPI_COMM_WORLD,statinfo)
    Charge2(Ny, Np, me-1, &j_1, &jN, &q);
    //MPI_RECV(x2(jN-Nx+1:jN),Nx,MPI_REAL,me-1,tag,MPI_COMM_WORLD,status,statinfo)
    Charge2(Ny, Np, me+1, &j_1, &jN, &q);
    //MPI_RECV(x2(j_1:j_1+Nx-1),Nx,MPI_REAL,me+1,tag,MPI_COMM_WORLD,status,statinfo)
    //MPI_SEND(x2(iN-Nx+1:iN),Nx,MPI_REAL,me+1,tag,MPI_COMM_WORLD,statinfo)

    for(i = 0; i < iN - i1; i++){
      y[i] = D3_m[i]*x2[i-Nx]+D2_m[i]*x2[i-1]+D1[i]*x2[i]+D2_p[i]*x2[i+1]+D3_p[i]*x2[i+Nx];
    }
    delete[] x2;
  }
}

double prodscal(const double X[],
                const double Y[],
                const int i1,
                const int iN){

  int i;
  double result = 0, Somme = 0;
  for (i = 0; i < iN - i1; i++) {
    result += X[i] * Y[i];
  }

  //MPI_ALLREDUCE(result,Somme,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,statinfo)

  return Somme;
}

int Reste(int k, int Nx){
  int r = k % Nx;
  if(r != 0){
    return r;
  }
  return Nx;
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
  double* r = (double*)malloc(sizeof(double) * (iN - i1));
  double* p = (double*)malloc(sizeof(double) * (iN - i1));
  double* r2 = (double*)malloc(sizeof(double) * (iN - i1));
  double* z = (double*)malloc(sizeof(double) * (iN - i1));
  double* y = (double*)malloc(sizeof(double) * (iN - i1));

  ProdMatVect(D1,D2_m,D2_p,D3_m,D3_p,x,y,i1,iN);
  for(i = 0; i < iN - i1; i++){
    x[i] = 293;
    r[i] = b[i] - y[i];
    p[i] = r[i];
  }
  double beta = sqrt(prodscal(r,r,i1,iN));
  int k = 0;
  double alpha = 0, gamma = 0;
  while (beta > eps1 && k <= kmax){
    ProdMatVect(D1,D2_m,D2_p,D3_m,D3_p,p,z,i1,iN);
    alpha = prodscal(r,r,i1,iN)/(prodscal(z,p,i1,iN));
    for(i = 0; i < iN - i1; i++){
      x[i] = x[i] + alpha * p[i];
      r2[i] = r[i];
      r[i] = r[i] - alpha * z[i];
    }
    gamma = prodscal(r,r,i1,iN)/prodscal(r2,r2,i1,iN);
    for(i = 0; i < iN - i1; i++){
      p[i] = r[i] + gamma * p[i];
    }
    beta = sqrt(prodscal(r2,r2,i1,iN));
    k++;
  }

  delete[] r;
  delete[] p;
  delete[] r2;
  delete[] z;
}
