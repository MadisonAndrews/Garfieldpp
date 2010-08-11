#include <iostream>
#include "ComponentBase.hh"

namespace Garfield {

ComponentBase::ComponentBase() :
  className("ComponentBase"),
  theGeometry(0),
  ready(false), 
  xPeriodic(false),         yPeriodic(false),         zPeriodic(false),
  xMirrorPeriodic(false),   yMirrorPeriodic(false),   zMirrorPeriodic(false),
  xAxiallyPeriodic(false),  yAxiallyPeriodic(false),  zAxiallyPeriodic(false),
  xRotationSymmetry(false), yRotationSymmetry(false), zRotationSymmetry(false),
  bx0(0.), by0(0.), bz0(0.),
  debug(false) {

}

void 
ComponentBase::SetGeometry(GeometryBase* geo) {

  // Make sure the geometry is defined
  if (geo == 0) {
    std::cerr << "ComponentBase::SetGeometry:\n";
    std::cerr << "    Geometry pointer is null.\n";
    return;
  }
 
  theGeometry = geo;

}

bool 
ComponentBase::GetMedium(const double x, const double y, const double z, 
                         Medium*& m) {
 
  if (theGeometry == 0) return false;
  return theGeometry->GetMedium(x, y, z, m);

}

void 
ComponentBase::Clear() {

  theGeometry = 0; 
  Reset();

}


void 
ComponentBase::WeightingField(const double x, const double y, const double z,
                              double& wx, double& wy, double& wz,
                              const std::string label) {
                              
  std::cerr << className << "::WeightingField:\n";
  std::cerr << "    This function is not implemented.\n";
  wx = wy = wz = 0.;

}

double
ComponentBase::WeightingPotential(
                    const double x, const double y, const double z,
                    const std::string label) {

  std::cerr << className << "::WeightingPotential:\n";
  std::cerr << "    This function is not implemented.\n";
  return 0.;

}

void 
ComponentBase::MagneticField(const double x, const double y, const double z,
    	                     double& bx, double& by, double& bz, int& status) {

  bx = bx0; by = by0; bz = bz0;
  status = 0;

}

void 
ComponentBase::SetMagneticField(const double bx, 
                                const double by, const double bz) {

  bx0 = bx; by0 = by; bz0 = bz;

}

bool
ComponentBase::GetBoundingBox(double& xmin, double& ymin, double& zmin,
                              double& xmax, double& ymax, double& zmax) {

  if (theGeometry == 0) return false;
  return (theGeometry->GetBoundingBox(xmin, ymin, zmin, xmax, ymax, zmax));

}

}
