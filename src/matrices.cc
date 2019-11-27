//#include <mpi>
#include <math.h>
#include <iostream>
#include <fstream>
#include <malloc.h>
#include "function.h"

// Construction de la matrice DF (cas linéaire)
void MatriceDF(double *D1,
               double *D2_m,
               double *D2_p,
               double *D3_m,
               double *D3_p,
               double D,
               double sx,
               double sy,
               int *i1,
               int *iN)
  {
    int k = 0;
    int nb_per_proc = iN-i1+1;

    for (k = 0; k < nb_per_proc, k++)
    {
      D1[k] = 1.0 + 2.0*D*(sx+sy);
      D2_m[k] = -D*sx;
      D2_p[k] = -D*sx;
      D3_m[k] = -D*sy;
      D3_p[k] = -D*sy;
    }

    If (iN == Nx*Ny-1)
    {
      for (k = nb_per_proc - Nx ; k < nb_per_proc; k++)
      D3_p[k] = 0.0 ;
    }

    If (i1 == 0)
    {
      for (k = 0 ; k < Nx; k++)
      D3_m[k] = 0.0;
    }

    for( k = 0; k < nb_per_proc ; k += Nx)
    {
      D2_p[k+Nx-1] = 0.0;
      D2_m[k] = 0.0;
    }


void sec_membre(int Nx,
                int Ny,
                double dx,
                double dy,
                double dt,
                double Lx,
                double Ly,
                double D,
                double *Fx,
                double t,
                double *i1,
                double *iN)

  {

    int i = 0, j = 0, k = 0;
    double sx, sy, A, B, C;

    sx = dt/(dx**2);
    sy = dt/(dy**2);
    A = 1.0+2.0*D*(sx+sy);
    B = -D*sx;
    C = -D*sy;


    // Premier proc
    if (i1 == 0)
    {
      Fx[0] = dt*f(dx,dy,t) - C*g(dx,0.0,t)-B*h(0.0,dy,t);

      for (k = 1; k < Nx-2; k++)
        Fx[k] = dt*f(Reste(k,Nx)*dx,dy,t) - C*g(Reste(k,Nx)*dx,0.0,t);

      Fx[Nx-1] = dt*f(Reste(k,Nx)*dx,dy,t) - C*g(Lx-dx,0.0,t)-B*h(Lx,dy,t);

      for (k = Nx, k < nb_per_proc; k++)                // k(i,j) = i + Nx*(j-1)
      {
        i = Reste(k,Nx);        // i(k) = reste de k/Nx (+ voir fonction Reste)
        j = (k-1.0)/Nx + 1.0;           !// j(k) = (quotient de k-1 divis� par Nx) + 1
        Fx[k] = dt*f(i*dx,j*dy,t);
        if (k%Nx == 1)
          Fx[k] = Fx[k]-B*h(0.0,j*dy,t);
        else if (k%Nx == 0)
          Fx[k] = Fx[k]-B*h(Lx,j*dy,t);
      }
    }

    // Dernier proc
    if (iN == Nx*Ny-1)
    {
      for ( k = 0; k < nb_per_proc - Nx; k++)
      {
        i = Reste(k,Nx);
        j = (k-1.0)/Nx + 1;
        Fx[k] = dt*f(i*dx,j*dy,t);
        if (k%Nx == 1)
          Fx[k] = Fx[k]-B*h(0._PR,j*dy,t);
        else if (k%Nx == 0)
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
        j = (k-1.0)/Nx + 1;

        Fx[k] = dt*f(i*dx,j*dy,t);
        If (k%Nx == 1)
        Fx[k] = Fx[k]-B*h(0.0,j*dy,t);
        Else If (Mod(k,Nx) == 0)
        Fx[k] = Fx[k]-B*h(Lx,j*dy,t);
      }
    }
