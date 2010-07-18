#include <iostream>
#include <string>
#include <sstream>
#include <cmath>

#include "MediumView.hh"
#include "Plotting.hh"

namespace Garfield {

MediumView::MediumView() :
  debug(false),
  canvas(0), hasExternalCanvas(false),
  medium(0),
  eMin(0.), eMax(1000.), bMin(0.), bMax(1.e5) {
  
  functions.clear();

}

MediumView::~MediumView() {

  if (!hasExternalCanvas && canvas != 0) delete canvas;

}

void
MediumView::SetCanvas(TCanvas* c) {

  if (c == 0) return;
  if (!hasExternalCanvas && canvas != 0) {
    delete canvas;
    canvas = 0;
  }
  canvas = c;
  hasExternalCanvas = true;

}

void
MediumView::SetMedium(Medium* m) {

  if (m == 0) {
    std::cerr << "MediumView::SetMedium:" << std::endl;
    std::cerr << "    Medium is not defined." << std::endl;
    return;
  }

  medium = m;

}

void
MediumView::SetElectricFieldRange(const double emin, const double emax) {

  if (emin >= emax || emin < 0.) {
    std::cerr << "MediumView::SetElectricFieldRange:" << std::endl;
    std::cerr << "    Incorrect field range." << std::endl;
    return;
  }

  eMin = emin; eMax = emax;

}

void
MediumView::SetMagneticFieldRange(const double bmin, const double bmax) {

  if (bmin >= bmax || bmin < 0.) { 
    std::cerr << "MediumView::SetMagneticFieldRange:" << std::endl;
    std::cerr << "    Incorrect field range." << std::endl;
    return;
  }

  bMin = bmin; bMax = bmax;

}

void 
MediumView::PlotElectronVelocity(const bool keep) {

  SetupCanvas();
  AddFunction(eMin, eMax, keep,
              "electric field [V/cm]", "drift velocity [cm/ns]", 0);
  canvas->Update();
  
}

void 
MediumView::PlotHoleVelocity(const bool keep) {

  SetupCanvas();
  AddFunction(eMin, eMax, keep,
              "electric field [V/cm]", "drift velocity [cm/ns]", 10);
  canvas->Update();

}

void 
MediumView::PlotIonVelocity(const bool keep) {

  SetupCanvas();
  AddFunction(eMin, eMax, keep,
              "electric field [V/cm]", "drift velocity [cm/ns]", 20);
  canvas->Update();

}

void 
MediumView::PlotElectronTownsend(const bool keep) {

  SetupCanvas();
  AddFunction(eMin, eMax, keep,
              "electric field [V/cm]", "Townsend coefficient [1/cm]", 3);
  canvas->Update();

}

void 
MediumView::PlotHoleTownsend(const bool keep) {

  SetupCanvas();
  AddFunction(eMin, eMax, keep,
              "electric field [V/cm]", "Townsend coefficient [1/cm]", 13);
  canvas->Update();

}

void 
MediumView::PlotElectronAttachment(const bool keep) {

  SetupCanvas();
  AddFunction(eMin, eMax, keep,
              "electric field [V/cm]", "Attachment coefficient [1/cm]", 4);
  canvas->Update();

}

void 
MediumView::PlotHoleAttachment(const bool keep) {

  SetupCanvas();
  AddFunction(eMin, eMax, keep,
              "electric field [V/cm]", "Attachment coefficient [1/cm]", 14);
  canvas->Update();
  
}

void
MediumView::SetupCanvas() {

  if (canvas == 0) {
    canvas = new TCanvas();
    canvas->SetTitle("Medium View");
    if (hasExternalCanvas) hasExternalCanvas = false;
  }
  canvas->cd();

}

void
MediumView::AddFunction(const double xmin, const double xmax, const bool keep,
                        const std::string xlabel, const std::string ylabel,
                        const int type) {

  if (medium == 0) {
    std::cerr << "MediumView::AddFunction:" << std::endl;
    std::cerr << "    Medium is not defined." << std::endl;
    return;
  }

  int idx = 0;
  std::string fname = "fMediumView_0";
  while (gROOT->GetListOfFunctions()->FindObject(fname.c_str())) {
    ++idx;
    std::stringstream ss;
    ss << "fMediumView_";
    ss  << idx;
    fname = ss.str();
  }

  if (!keep) functions.clear();

  TF1 fNew(fname.c_str(), this, &MediumView::EvaluateFunction, 
            xmin, xmax, 1, "MediumView", "EvaluateFunction");
  functions.push_back(fNew);

  const std::string title = medium->GetName() + ";" + 
                            xlabel + ";" + ylabel;
  functions.back().SetRange(xmin, xmax);
  functions.back().GetXaxis()->SetTitle(xlabel.c_str());
  functions.back().GetYaxis()->SetTitle(ylabel.c_str());
  functions.back().SetTitle(title.c_str());
  functions.back().SetParameter(0, type);
  if (type < 10) {
    functions.back().SetLineColor(kOrange);
  } else if (type < 20) {
    functions.back().SetLineColor(kGreen + 2);
  } else {
    functions.back().SetLineColor(kRed);
  }
  if (keep) {
    functions.back().Draw("LSAME");
  } else {
    functions.back().DrawCopy("");
  }
 
}

double
MediumView::EvaluateFunction(double* pos, double* par) {

  if (medium == 0) return 0.;
  
  int type = int(par[0]);
  const double x = pos[0];
  double y = 0.;

  // Auxiliary variables
  double a = 0., b = 0., c = 0.;

  switch (type) {
    case 0:
      // Electron drift velocity
      if (!medium->ElectronVelocity(x, 0, 0, 0, 0, 0, a, b, c)) return 0.;
      y = fabs(a);
      break;
    case 1:
      // Electron transverse diffusion
      if (!medium->ElectronDiffusion(x, 0, 0, 0, 0, 0, a, b)) return 0.;
      y = b;
      break;
    case 2:
      // Electron longitudinal diffusion
      if (!medium->ElectronDiffusion(x, 0, 0, 0, 0, 0, a, b)) return 0.;
      y = a;
      break;
    case 3:
      // Electron Townsend coefficient
      if (!medium->ElectronTownsend(x, 0, 0, 0, 0, 0, a)) return 0.;
      y = a;
      break;
    case 4:
      // Electron attachment coefficient
      if (!medium->ElectronAttachment(x, 0, 0, 0, 0, 0, a)) return 0.;
      y = a;
      break;
    case 10:
      // Hole drift velocity
      if (!medium->HoleVelocity(x, 0, 0, 0, 0, 0, a, b, c)) return 0.;
      y = a;
      break;
    case 11:
      // Hole transverse diffusion
      if (!medium->HoleDiffusion(x, 0, 0, 0, 0, 0, a, b)) return 0.;
      y = b;
      break;
    case 12:
      // Hole longitudinal diffusion
      if (!medium->HoleDiffusion(x, 0, 0, 0, 0, 0, a, b)) return 0.;
      y = a;
      break;
    case 13:
      // Hole Townsend coefficient
      if (!medium->HoleTownsend(x, 0, 0, 0, 0, 0, a)) return 0.;
      y = a;
      break;
    case 14:
      // Hole attachment coefficient
      if (!medium->HoleAttachment(x, 0, 0, 0, 0, 0, a)) return 0.;
      y = a;
      break;
    case 20:
      // Ion drift velocity
      if (!medium->IonVelocity(x, 0, 0, 0, 0, 0, a, b, c)) return 0.;
      y = fabs(a);
      break;
    case 21:
      // Ion transverse diffusion
      if (!medium->IonDiffusion(x, 0, 0, 0, 0, 0, a, b)) return 0.;
      y = b;
      break;
    case 22:
      // Ion longitudinal diffusion
      if (!medium->IonDiffusion(x, 0, 0, 0, 0, 0, a, b)) return 0.;
      y = a;
      break;
    default:
      std::cerr << "MediumView::EvaluateFunction:" << std::endl;
      std::cerr << "    Unknown type of transport coefficient requested." 
                << std::endl;
      std::cerr << "    Program bug!" << std::endl;
      return 0.;
  }

  return y;
    
}

}
