#ifndef MATRICE_H
#define MATRICE_H

#include <mpi.h>
#include <math.h>
#include <iostream>
#include <fstream>

#include "DataFile.hpp"

class Solver{
private:
  MPI_Status status;
  int tag = 100;
  DataFile* data_file;
  double sx, sy, D, dx, dy, dt, eps1, tf, Lx, Ly;
  int i1, iN, Nx, Ny, kmax, nb_per_proc, me, Np, cond_init;
  double *D1, *D2_m, *D2_p, *D3_m, *D3_p, *U, *U0, *F;


  int Reste(int k) const;
  void Charge2(double* q2);
  void MatriceDF(void);
  void ProdMatVect(double x[], double y[]);
  double prodscal(const double X[], const double Y[]);
  double f(double x, double y, double t);
  double g(double x, double y, double t);
  double h(double x, double y, double t);

public:
  Solver(std::string file_name);
  ~Solver(void);

  int getMe(void) const { return me; }
  double getTf(void) const { return tf; }
  double getDt(void) const { return dt; }
  int getNp(void) const { return Np; }
  int getNbPerProc(void) const { return nb_per_proc; }
  double* getSolution(void) const { return U; }

  void computeConjugateGradient(void);
  void computeSecondMember(double t);
  void saveSolution(const char* fileName) const;
  void makePlot(const char* fileName) const;
};

#endif
