#ifndef _DATA_FILE_H
#define _DATA_FILE_H

#include <string>
#include <vector>
#include <iostream>
// Définition de la classe

class DataFile {
 private:

  std::string _file_name;
  int _Nx, _Ny, _kmax, _height, _cond_init;
  double _Lx, _Ly, _D, _dt, _epsilon, _tf;
  std::string _results;

  bool _if_Nx;
  bool _if_Ny;
  bool _if_Lx;
  bool _if_Ly;
  bool _if_D;
  bool _if_dt;
  bool _if_tf;
  bool _if_kmax;
  bool _if_epsilon;
  bool _if_results;
  bool _if_cond_init;
  bool _if_height;

 public:
  DataFile(std::string file_name);
  std::string Get_file_name() const {return _file_name;}
  void ReadDataFile();
  double Get_Nx() const {return _Nx;};
  double Get_Ny() const {return _Ny;};
  double Get_Lx() const {return _Lx;};
  double Get_Ly() const {return _Ly;};
  double Get_D() const {return _D;};
  double Get_dt() const { return _dt;}
  double Get_tf() const { return _tf;}
  double Get_kmax() const {return _kmax;};
  double Get_epsilon() const {return _epsilon;};
  double Get_cond_init() const {return _cond_init;};
  int Get_height() const {return _height;};
  std::string Get_results() const {return _results;};
  //void SaveResult();


};

#endif
