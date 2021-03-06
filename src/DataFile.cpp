#ifndef _DATA_FILE_CPP

#include <fstream>
#include <iostream>
#include <cmath>
#include <assert.h>

#include "DataFile.hpp"

using namespace std;


// Classe pour la lecture du fichier data.txt
DataFile::DataFile(std::string file_name)
: _if_Nx(false), _if_Ny(false), _if_Lx(false), _if_Ly(false), _if_D(false), _if_dt(false), _if_tf(false), _if_kmax(false), _if_epsilon(false), _if_results(false), _if_cond_init(false), _file_name(file_name), _if_height(false)
{}



  void DataFile::ReadDataFile()
  {
    ifstream data_file(_file_name.data());
    if (!data_file.is_open())
    {
      cout << "Unable to open file " << _file_name << endl;
      abort();
    }
    else
    {
      //cout << "-------------------------------------------------" << endl;
      //cout << "Reading data file " << _file_name << endl;
    }

    string file_line;

    while (!data_file.eof())
    {
      getline(data_file, file_line);
      if (file_line.size() > 0)
        if (file_line[0] == '#')
          continue;

      if (file_line.find("Nx") != std::string::npos)
      {
        data_file >> _Nx; _if_Nx = true;
        //continue;
      }

      if (file_line.find("Ny") != std::string::npos)
      {
        data_file >> _Ny; _if_Ny = true;
      }

      if (file_line.find("Lx") != std::string::npos)
      {
        data_file >> _Lx; _if_Lx = true;
      }

      if (file_line.find("Ly") != std::string::npos)
      {
        data_file >> _Ly; _if_Ly = true;
      }

      if (file_line.find("D") != std::string::npos)
      {
        data_file >> _D; _if_D = true;
      }

      if (file_line.find("dt") != std::string::npos)
      {
        data_file >> _dt; _if_dt = true;
      }

      if (file_line.find("tf") != std::string::npos)
      {
        data_file >> _tf; _if_tf = true;
      }

      if (file_line.find("kmax") != std::string::npos)
      {
        data_file >> _kmax; _if_kmax = true;
      }

      if (file_line.find("epsilon") != std::string::npos)
      {
        data_file >> _epsilon; _if_epsilon = true;
      }

      if (file_line.find("height") != std::string::npos)
      {
        data_file >> _height; _if_height = true;
      }

      if (file_line.find("set de conditions initiales") != std::string::npos)
      {
        data_file >> _cond_init; _if_cond_init = true;
      }

      if (file_line.find("results") != std::string::npos)
      {
        data_file >> _results; _if_results = true;
      }
    }

    if (!_if_results)
    {

      //cout << "-------------------------------------------------" << endl;
      cout << "Beware - The default results folder name (results) is used." << endl;
      _results = "results";
    }

    if (!_if_results)
    {

      //cout << "-------------------------------------------------" << endl;
      cout << "Beware - The default initial conditions (set 1) are used." << endl;
      _cond_init = 1;
    }

  }

  #define _DATA_FILE_CPP
  #endif
