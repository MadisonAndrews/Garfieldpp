#include <iostream>
// #include <fstream>
// #include <sstream>
#include <string>
#include <cmath>

#include <Riostream.h>

// #include <TApplication.h>
#include <TCanvas.h>
// #include <TH1F.h>
// #include <TGeoManager.h>
// #include <TGeoMaterial.h>
// #include <TGeoMedium.h>
// #include <TGeoVolume.h>
// #include <TGeoBBox.h>
// #include <TGeoTube.h>
// #include <TGeoPcon.h>
// #include <TGeoHalfSpace.h>
// #include <TGeoMatrix.h>
// #include <TGeoCompositeShape.h>
#include <TTree.h>
#include <TFile.h>
#include <TROOT.h>
// #include <TDatime.h>
#include <stdio.h>
#include <stdlib.h>

#include "../Include/ComponentAnsys123.hh"
#include "../Include/ViewField.hh"
#include "../Include/MediumMagboltz.hh"
#include "../Include/Sensor.hh"
#include "../Include/AvalancheMicroscopic.hh"
#include "../Include/AvalancheMC.hh"
#include "../Include/Random.hh"
#include "../Include/Plotting.hh"
#include "../Include/ViewSignal.hh"

// #include <Python.h>
// #include <boost/python/def.hpp>                                                                                                    
// #include <boost/python/module.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
                                                                             
using namespace Garfield;
using namespace std;


int main(int argc, char * argv[]) 
{
  boost::property_tree::ptree pt;
  boost::property_tree::ini_parser::read_ini("config.ini", pt);

  // Load the field map.
  ComponentAnsys123* fm = new ComponentAnsys123();
  const std::string path( pt.get<std::string>("fieldmap.path") );
  const std::string efile( path + pt.get<std::string>("fieldmap.efile") );
  const std::string nfile( path + pt.get<std::string>("fieldmap.nfile") );
  const std::string mfile( path + pt.get<std::string>("fieldmap.mfile") );
  const std::string sfile( path + pt.get<std::string>("fieldmap.mfile") );
  const std::string wfile( path + pt.get<std::string>("fieldmap.wfile") );
  const std::string unit( pt.get<std::string>("fieldmap.unit") );
  fm->Initialise(efile, nfile, mfile, sfile, unit);
  fm->EnableMirrorPeriodicityX();
  fm->EnableMirrorPeriodicityY();
  fm->PrintRange();

  // Dimensions of the GEM
  const double pitch      = pt.get<double>("gem.pitch");
  const double kapton     = pt.get<double>("gem.kapton");
  const double metal      = pt.get<double>("gem.metal");
  const double outdia     = pt.get<double>("gem.outdia");
  const double middia     = pt.get<double>("gem.middia");
  // const double drift      = pt.get<double>("gem.drift");
  const double transfer1  = pt.get<double>("gem.transfer1");
  const double transfer2  = pt.get<double>("gem.transfer2");
  const double induction  = pt.get<double>("gem.induction");

  // Setup the gas.
  MediumMagboltz* gas = new MediumMagboltz();
  gas->SetComposition(pt.get<std::string>("gas.gas1"),
		      pt.get<double>("gas.gas1_mix"),
		      pt.get<std::string>("gas.gas2"),
		      pt.get<double>("gas.gas2_mix") );		     
  gas->SetTemperature(pt.get<double>("gas.temperature") );
  gas->SetPressure(pt.get<double>("gas.pressure") );
  gas->EnableDebugging();
  gas->Initialise();  
  gas->DisableDebugging();
  gas->EnablePenningTransfer( pt.get<double>("gas.r_penning"),
			      pt.get<double>("gas.lambda_penning"),
			      pt.get<std::string>("gas.gas1") );
  gas->LoadIonMobility( pt.get<std::string>("gas.ion_mobility") );
  
  // Associate the gas with the corresponding field map material. 
  const int nMaterials = fm->GetNumberOfMaterials();
  for (int i = 0; i < nMaterials; ++i) {
    const double eps = fm->GetPermittivity(i);
    if (fabs(eps - 1.) < 1.e-3) fm->SetMedium(i, gas);
  }
  fm->PrintMaterials();

  // Create the sensor.
  Sensor* sensor = new Sensor();
  sensor->AddComponent(fm);
  sensor->SetArea(pt.get<double>("sensor.x_min") * pitch,
		  pt.get<double>("sensor.x_max") * pitch,
		  pt.get<double>("sensor.y_min") * pitch,
		  pt.get<double>("sensor.y_max") * pitch,
		  pt.get<double>("sensor.z_min"),
		  pt.get<double>("sensor.z_max") );

  AvalancheMicroscopic* aval = new AvalancheMicroscopic();
  aval->SetSensor(sensor);

  AvalancheMC* drift = new AvalancheMC();
  drift->SetSensor(sensor);
  drift->SetDistanceSteps(pt.get<double>("drift.distance_steps") );

  const int nEvents( pt.get<int>("sim.number_of_events") );
  for (int i = nEvents; i--;) { 
    std::cout << i << "/" << nEvents << "\n";
    // Randomize the initial position.
    const double smear = pitch / 2.; 
    double x0 = -smear + RndmUniform() * smear;
    double y0 = -smear + RndmUniform() * smear;
    double z0 = 0.025; 
    double t0 = 0.;
    double e0 = 0.1;
    aval->AvalancheElectron(x0, y0, z0, t0, e0, 0., 0., 0.);
    int ne = 0, ni = 0;
    aval->GetAvalancheSize(ne, ni);
    const int np = aval->GetNumberOfElectronEndpoints();
    double xe1, ye1, ze1, te1, e1;
    double xe2, ye2, ze2, te2, e2;
    double xi1, yi1, zi1, ti1;
    double xi2, yi2, zi2, ti2;
    int status;
    for (int j = np; j--;) {
      aval->GetElectronEndpoint(j, xe1, ye1, ze1, te1, e1, 
                                   xe2, ye2, ze2, te2, e2, status);
    }
  }

}

  
