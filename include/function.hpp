#ifndef _MATRICES_H
#define _MATRICES_H

#include "DataFile.hpp"

void Charge2(int n,
  int Np,
  int me,
  int* i1,
  int* iN,
  double* q2);

void Init(void);
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


#endif
