#ifndef _FONCTION_H
#define _FONCTION_H

#include <mpi.h>

#include "DataFile.hpp"


extern MPI_Status status;
extern int PR;
extern int tag;
extern double *U, *U0;
extern int cond_init;
extern double sx,sy,dx,dy,dt,Lx,Ly,D,eps1,t,tf, q;
extern int i, me, Np, statinfo, i1, iN, j_1, jN, Ny, Nx, n, k, kmax;


void Charge2(int n,
  int Np,
  int me,
  int* i1,
  int* iN,
  double* q2);


void Charge_part_domaine(int n,
  int Np,
  int me,
  int* i1,
  int* i12,
  int* iN,
  int* iN2,
  int h);    // h hauteur de la partie partagée

void Init(DataFile*);
double f(double x, double y, double t);
double g(double x, double y, double t);
double h(double x, double y, double t);

void ProdMatVect(double D1[],
  double D2_m[],
  double D2_p[],
  double D3_m[],
  double D3_p[],
  double x[],
  double y[],
  int i1,
  int iN);

double prodscal(const double X[],
  const double Y[],
  const int i1,
  const int iN);


double prodscal(const double X[],
  const double Y[],
  const int i1,
  const int iN);


int Reste(int k, int Nx);

void grad_conj(double D1[],
  double D2_m[],
  double D2_p[],
  double D3_m[],
  double D3_p[],
  double x[],
  double b[],
  int i1,
  int iN);


  int cmp_vect(const double x[],
               const double y[],
               const int Nx);    // ajouter epsilon (précision)

#endif
