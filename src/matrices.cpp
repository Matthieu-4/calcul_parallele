//#include <mpi>
#include <math.h>
#include <iostream>
#include <fstream>
#include <malloc.h>

#include "function.hpp"
#include "matrices.hpp"
using namespace std;

using namespace std;

// Construction de la matrice DF (cas linéaire)
void MatriceDF(double *D1,
               double *D2_m,
               double *D2_p,
               double *D3_m,
               double *D3_p,
               double sx,
               double sy,
               int i1,
               int iN,
               DataFile* data_file)
  {
    int k = 0;
    int nb_per_proc = iN-i1+1;

    int Nx = data_file->Get_Nx();
    int Ny = data_file->Get_Ny();
    double D = data_file->Get_D();

    for (k = 0; k < nb_per_proc; k++)
    {
      D1[k] = 1.0 + 2.0*D*(sx+sy);
      D2_m[k] = -D*sx;
      D2_p[k] = -D*sx;
      D3_m[k] = -D*sy;
      D3_p[k] = -D*sy;
    }

    if (iN == Nx*Ny-1)
    {
      for (k = nb_per_proc - Nx ; k < nb_per_proc; k++)
      D3_p[k] = 0.0 ;
    }

    if (i1 == 0)
    {
      for (k = 0 ; k < Nx; k++)
      D3_m[k] = 0.0;
    }

    for( k = 0; k < nb_per_proc ; k += Nx)
    {
      D2_p[k+Nx-1] = 0.0;
      D2_m[k] = 0.0;
    }
  }


void sec_membre(double dx,
                double dy,
                double *Fx,
                double t,
                int i1,
                int iN,
                DataFile* data_file)
  {
    int Nx = data_file->Get_Nx();
    int Ny = data_file->Get_Ny();
    double D = data_file->Get_D();
    int Lx = data_file->Get_Lx();
    int Ly = data_file->Get_Ly();
    double dt = data_file->Get_dt();

    int nb_per_proc = iN-i1+1;

    int i = 0, j = 0, k = 0;
    double sx, sy, A, B, C;


    sx = dt/(dx*dx);
    sy = dt/(dy*dy);
    A = 1.0+2.0*D*(sx+sy);
    B = -D*sx;
    C = -D*sy;

    // Premier proc
    if (i1 == 0)
    {
      Fx[0] = dt*f(dx,dy,t) - C*g(dx,0.0,t)-B*h(0.0,dy,t);

      for (k = 1; k < Nx-1; k++)
        Fx[k] = dt*f(Reste(k,Nx)*dx,dy,t) - C*g(Reste(k,Nx)*dx,0.0,t);


      Fx[Nx-1] = dt*f(Reste(k,Nx)*dx,dy,t) - C*g(Lx-dx,0.0,t)-B*h(Lx,dy,t);



      for (k = Nx; k < nb_per_proc; k++)                // k(i,j) = i + Nx*(j-1)
      {
        i = Reste(k,Nx);                // i(k) = reste de k/Nx (+ voir fonction Reste)
        j = (i1+k)/Nx + 1.0;           // j(k) = (quotient de k-1 divis� par Nx) + 1
        Fx[k] = dt*f(i*dx,j*dy,t);


        if (k%Nx == 0)
          Fx[k] = Fx[k]-B*h(0.0,j*dy,t);
        else if (k%Nx == Nx - 1)
          Fx[k] = Fx[k]-B*h(Lx,j*dy,t);
      }
    }

    // Dernier proc
    else if (iN == Nx*Ny-1)
    {

      for ( k = 0; k < nb_per_proc - Nx; k++)
      {
        i = Reste(k,Nx);
        j = (i1+k)/Nx + 1.0;
        Fx[k] = dt*f(i*dx,j*dy,t);


        if (k%Nx == 0)
          Fx[k] = Fx[k]-B*h(0.0,j*dy,t);
        else if (k%Nx == Nx - 1)
          Fx[k] = Fx[k]-B*h(Lx,j*dy,t);
      }

      Fx[nb_per_proc - Nx] = dt*f(dx,Ly-dy,t)-C*g(dx,Ly,t)-B*h(0.0,Ly-dy,t);


      for (k= nb_per_proc - Nx + 1; k < nb_per_proc - 1; k++)
        Fx[k] = dt*f(Reste(k,Nx)*dx,Ly-dy,t) - C*g(Reste(k,Nx)*dx,Ly,t);

      Fx[nb_per_proc] = dt*f(Lx-dx,Ly-dy,t)-C*g(Lx-dx,Ly,t)-B*h(Lx,Ly-dy,t);

    }

    // Les autres procs
    else
    {
      for (k = 0; k < nb_per_proc; k++)
      {
        i = Reste(k,Nx);
        j = (i1+k)/Nx + 1;

        Fx[k] = dt*f(i*dx,j*dy,t);

        if (k%Nx == 0)
          Fx[k] = Fx[k]-B*h(0.0,j*dy,t);
        else if (k%Nx == Nx-1){
          Fx[k] = Fx[k]-B*h(Lx,j*dy,t);
      }
    }
  }
