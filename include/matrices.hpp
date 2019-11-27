#ifndef _MATRICES_H
#define _MATRICES_H

#include "function.hpp"

void MatriceDF(double *D1,
  double *D2_m,
  double *D2_p,
  double *D3_m,
  double *D3_p,
  double sx,
  double sy,
  int i1,
  int iN,
  DataFile* data_file);



void sec_membre(double dx,
  double dy,
  double *Fx,
  double t,
  int i1,
  int iN,
  DataFile* data_file);




#endif
