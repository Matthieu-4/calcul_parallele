
#include <mpi.h>
#include <math.h>
#include <iostream>
#include <fstream>

#include "Matrices.hpp"
#include "DataFile.hpp"

using namespace std;


int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);

  if (argc < 2)
  {
    cout << "Please, enter the name of your data file." << endl;
    abort();
  }

  Solver* solver = new Solver(argv[1]);
  const double* solution = solver->getSolution();


   /////////////////////////// BOUCLE EN TEMPS /////////////////////////

  int k = 0;
  double t = 0;
  double tf = solver->getTf();
  double dt = solver->getDt();
  while (t < tf)
  {
    solver->computeSecondMember(t);
    solver->computeConjugateGradient();

    t += dt;
    k += 1;
    if (k%50 == 0){
      printf("%d %f\n",k, solution[4]);
    }
  }


  // Ecriture de la solution sur des fichiers .dat pour chaque proc


  char file_name [50];
  char plot_name [50] = "Result/sol.plot";
  sprintf(file_name, "Result/sol%d.dat", solver->getMe());

  solver->saveSolution(file_name);
  solver->makePlot(plot_name);


  MPI_Finalize();

  return 0;
}
