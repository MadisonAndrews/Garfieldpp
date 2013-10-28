#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>

#include <TApplication.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TFile.h>
#include <TPolyMarker.h>
#include <TPaveText.h>

#include "../../Include/ComponentAnsys123.hh"
#include "../../Include/ViewField.hh"
#include "../../Include/ViewFEMesh.hh"
#include "../../Include/MediumMagboltz.hh"
#include "../../Include/Sensor.hh"
#include "../../Include/AvalancheMicroscopic.hh"
#include "../../Include/AvalancheMC.hh"
#include "../../Include/Random.hh"
#include "../../Include/Plotting.hh"
#include "../../Include/GarfieldConstants.hh"

using namespace Garfield;
using namespace std;

struct charge{
  double x,y,z;
  double t;
  double e;
  double dx,dy,dz;
};

//
// Makes an animation of the passage of electrons through GEMs
//
int main(int argc, char * argv[]) 
{
  TApplication app("app", &argc, argv);
  plottingEngine.SetDefaultStyle();
 
  ComponentAnsys123* fm = new ComponentAnsys123();
  
  // string path = "/afs/cern.ch/user/d/dildick/FEM/lis_files/std_1000_3000_500/";
  string path = "";
  const double rPenning = 0.55;
  const double lambdaPenning = 0.0;
  string fmname[4];
  string lisfile[] = {"ELIST.lis", "NLIST.lis", "MPLIST.lis", "PRNSOL.lis"};
  char buf[100];

  for (int i=0; i<sizeof(lisfile)/sizeof(lisfile[0]); i++){
    fmname[i] = path + lisfile[i];
  }
  fm->DisableDeleteBackgroundElements();
  fm->Initialise(fmname[0], fmname[1], fmname[2], fmname[3], "mm");
  fm->EnableMirrorPeriodicityX();
  fm->EnableMirrorPeriodicityY();
  fm->PrintRange();
  
  const double pitch = 0.014; 
  const double kapton = 50.e-4;
  const double metal = 5.e-4;
  const double outdia = 70.e-4;
  const double middia = 50.e-4;

  // Setup the gas.
  MediumMagboltz* gas = new MediumMagboltz();
  gas->SetComposition("ar", 70., "co2", 30.);
  gas->SetTemperature(293.15);
  gas->SetPressure(760.);
  gas->EnableDebugging();
  gas->Initialise();  
  gas->DisableDebugging();
  gas->EnablePenningTransfer(rPenning, lambdaPenning, "ar");
  // gas->LoadIonMobility("/afs/cern.ch/user/d/dildick/Garfield++/Data/IonMobility_Ar+_Ar.txt");
  gas->LoadIonMobility("../../Data/IonMobility_Ar+_Ar.txt");

  const int nMaterials = fm->GetNumberOfMaterials();
  for (int i = 0; i < nMaterials; ++i) {
    const double eps = fm->GetPermittivity(i);
    if (fabs(eps - 1.) < 1.e-3) fm->SetMedium(i, gas);
  }
  fm->PrintMaterials();

  Sensor* sensor = new Sensor();
  sensor->AddComponent(fm);
  sensor->SetArea(-5 * pitch, -5 * pitch, -0.03,
                   5 * pitch,  5 * pitch,  0.03);

  AvalancheMicroscopic* aval = new AvalancheMicroscopic();
  aval->SetSensor(sensor);

  AvalancheMC* drift = new AvalancheMC();
  drift->SetSensor(sensor);

  const bool plotDrift = false;
  ViewDrift* driftView = new ViewDrift();
  if (plotDrift) {
    driftView->SetArea(-2 * pitch, -2 * pitch, -0.02,
		       2 * pitch,  2 * pitch,  0.02);
    aval->SetCollisionSteps(10); 
    aval->EnablePlotting(driftView);
    drift->EnablePlotting(driftView);
  }

  std::vector<charge> oldElectrons;
  std::vector<charge> newElectrons;
  std::vector<charge> oldIons;
  std::vector<charge> newIons;
  
  newElectrons.clear();
  oldElectrons.clear();  
  newIons.clear();
  oldIons.clear();
  TCanvas* c1 = new TCanvas();
  
  // Visualize the FE map
  ViewFEMesh* vFE = new ViewFEMesh();
  vFE->SetCanvas(c1);
  vFE->SetComponent(fm);
  vFE->SetPlane(0, -1, 0, 0, 0, 0);
  vFE->SetFillMesh(false);
  vFE->SetColor(1, kBlue);
  vFE->SetColor(2, kOrange);
  vFE->SetColor(3, kGreen);
  if (plotDrift) {
    vFE->SetViewDrift(driftView);
  }
  vFE->SetArea(-0.01, -0.03, -0.02, 0.01, 0.01, 0.02);
  vFE->Plot();

  const double smear = 0; 
  double x0 = -smear + RndmUniform() * smear;
  double y0 = -smear + RndmUniform() * smear;
  double z0 = 0.005; 
  double e0 = 0.1;
  charge seed;
  seed.x = x0; seed.y = y0; seed.z = z0; seed.t = 0.;
  seed.e = e0;
  seed.dx = 0.; seed.dy = 0.; seed.dz = 0.;
  newElectrons.push_back(seed);
  newIons.push_back(seed);
  bool finished = false;
  double tStep = 0.05;
  // Put this at a better place:
  drift->SetTimeSteps(tStep);
  double tMin = 0.;
  double tMax = tMin + tStep;
  int counter = 0;
  while (!finished){
    std::cout << tMin << '\n';
    drift->SetTimeWindow(tMin,tMax);
    aval->SetTimeWindow(tMin,tMax);
    newElectrons.swap(oldElectrons);
    newElectrons.clear();
    newIons.swap(oldIons);
    newIons.clear();
    int ne = oldElectrons.size();
    for(int j = ne; j--;){
      double x0 = oldElectrons[j].x;
      double y0 = oldElectrons[j].y;
      double z0 = oldElectrons[j].z;
      double t0 = oldElectrons[j].t;
      double e0 = oldElectrons[j].e;
      double dx0 = oldElectrons[j].dx;
      double dy0 = oldElectrons[j].dy;
      double dz0 = oldElectrons[j].dz;
	
      aval->AvalancheElectron(x0, y0, z0, t0, e0, dx0, dy0, dz0);  
	
      const int np = aval->GetNumberOfElectronEndpoints();
      double xe1, ye1, ze1, te1, e1;
      double xe2, ye2, ze2, te2, e2;
      double dxe2, dye2,dze2;
      int status;
      for (int l = np; l--;) {
	aval->GetElectronEndpoint(l, xe1, ye1, ze1, te1, e1, 
				  xe2, ye2, ze2, te2, e2, dxe2,dye2, dze2,status);
	if (status == StatusOutsideTimeWindow){
	  charge newElectron;
	  newElectron.x  = xe2;
	  newElectron.y  = ye2;
	  newElectron.z  = ze2;
	  newElectron.dx  = dxe2;
	  newElectron.dy = dye2;
	  newElectron.dz = dze2;
	  newElectron.e = e2;
	  newElectron.t = te2;
	  newElectrons.push_back(newElectron);
	}
	double distance = sqrt(pow(xe1 - x0, 2) + 
                               pow(ye1 - y0, 2) + 
                               pow(ze1 - z0, 2));
	if(distance < 1.e-8) {
	  continue;
	}
	charge newIon;
	newIon.x  = xe1;
	newIon.y  = ye1;
	newIon.z  = ze1;
	newIon.t = te1;
	oldIons.push_back(newIon);
      }
    }
    const int nIons = oldIons.size();
    for(int j = nIons; j--;){
      double x0 = oldIons[j].x;
      double y0 = oldIons[j].y;
      double z0 = oldIons[j].z;
      double t0 = oldIons[j].t;
	
      drift->DriftIon(x0, y0, z0, t0);
	
      double xi1, yi1, zi1, ti1;
      double xi2, yi2, zi2, ti2;
      int status;
 
      drift->GetIonEndpoint(0, xi1, yi1, zi1, ti1, 
		        	  xi2, yi2, zi2, ti2, status);
      if (status == StatusOutsideTimeWindow) {
	charge newIon;
	newIon.x  = xi2;
	newIon.y  = yi2;
	newIon.z  = zi2;
	newIon.t = ti2;
	newIons.push_back(newIon);
      }
    
    }
    const int nNewElectrons = newElectrons.size();
    const int nNewIons = newIons.size();
    if (nNewElectrons == 0 && nNewIons > 0) tStep = 50.; 
    tMin = tMax;
    tMax += tStep;
    if (nNewElectrons < 1 && nNewIons < 1) { 
      finished = true;
    } else {
      //TPointSet3D p(nNew);
      TPolyMarker p(nNewElectrons);
      TPolyMarker q(nNewIons);
      p.SetMarkerStyle(20);
      p.SetMarkerSize(0.5);
      p.SetMarkerColor(kBlue); 
      q.SetMarkerStyle(20);
      q.SetMarkerSize(0.5);
      q.SetMarkerColor(kPink);
      for (int k = nNewElectrons; k--;){
	//p.SetNextPoint(newElectrons[k].x, 
	//       newElectrons[k].y,
	//       newElectrons[k].z);		
	p.SetNextPoint(newElectrons[k].x, newElectrons[k].z);
      } 
      for (int k = nNewIons; k--;){
	//p.SetNextPoint(newElectrons[k].x, 
	//       newElectrons[k].y,
	//       newElectrons[k].z);		
	q.SetNextPoint(newIons[k].x, newIons[k].z);
      }
      c1->cd();
      vFE->Plot();
      p.Draw("same");
      q.Draw("same");
      TPaveText pave(0.8, 0.85, 0.95, 0.9, "NDC");
      pave.SetTextFont(42);
      pave.SetFillColor(0);
      pave.SetBorderSize(0);
      char lbl[100];
      sprintf(lbl, "t = %04.2f ns", tMin);
      pave.AddText(lbl);
      pave.Draw("same");
      c1->Update();
      sprintf(buf,"/afs/cern.ch/user/d/dildick/work/GEM/StandaloneSimulations/Garfieldpp/Test/GEM/pics/movie_%d.png", counter);
      c1->Print(buf);
    }
    counter += 1;
  } 
  app.Run(true);
}
