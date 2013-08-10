#include <iostream>
#include <fstream>
#include <cmath>

#include <TROOT.h>
#include <TRint.h>

#include "Plotting.hh"

int main(int argc, char * argv[]) {

  std::cout << "Welcome to Garfield++\n";
  Garfield::plottingEngine.SetDefaultStyle();
  
  TRint* app = new TRint("Garfield", &argc, argv, 0, 0);
  
  app->Run(kTRUE);

}
