#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>

#include "MediumSilicon.hh"
#include "Random.hh"
#include "GarfieldConstants.hh"
#include "FundamentalConstants.hh"

namespace Garfield {

MediumSilicon::MediumSilicon() :
  Medium(), 
  bandGap(1.12), 
  dopingType('i'), dopingConcentration(0.),
  mLongX(0.916), mTransX(0.191),
  eLatticeMobility(1.35e-6), hLatticeMobility(0.48e-6),
  eMobility(1.43e-6), hMobility(0.46e-6),
  eBetaCanali(1.109), hBetaCanali(1.213),
  eSatVel(1.02e-2), hSatVel(0.72e-2),
  eHallFactor(1.15), hHallFactor(0.7),
  eTrapCs(1.e-15), hTrapCs(1.e-15),
  eTrapDensity(1.e13), hTrapDensity(1.e13),
  eTrapTime(0.), hTrapTime(0.),
  trappingModel(0),  
  eImpactA0(3.318e5), eImpactA1(0.703e6), eImpactA2(0.),
  eImpactB0(1.135e6), eImpactB1(1.231e6), eImpactB2(0.),
  hImpactA0(1.582e6), hImpactA1(0.671e6), hImpactA2(0.),
  hImpactB0(2.036e6), hImpactB1(1.693e6), hImpactB2(0.),
  hasUserMobility(false), hasUserSaturationVelocity(false),
  latticeMobilityModel(0),
  dopingMobilityModel(0),
  saturationVelocityModel(0),
  highFieldMobilityModel(0),
  impactIonisationModel(0),
  useCfOutput(false),
  useNonParabolicity(false), useAnisotropy(false),
  eFinal(2.), eStep(0.001),
  nLevelsX(0),
  cfNullElectrons(0.),
  hasOpticalData(false), opticalDataFile("OpticalData_Si_V1.txt") {

  className = "MediumSilicon";
  name = "Si";

  SetTemperature(300.);
  SetDielectricConstant(11.9);
  SetAtomicNumber(14.);
  SetAtomicWeight(28.0855);
  SetMassDensity(2.329);
  
  EnableDrift();
  EnablePrimaryIonisation();
  microscopic = true;

  wValue = 3.6;
  fanoFactor = 0.11;  
  
  cfTotElectronsX.clear();
  cfElectronsX.clear();
  energyLossElectronsX.clear();
  scatTypeElectronsX.clear();

}

void 
MediumSilicon::SetDoping(const char type, const double c) {

  if (toupper(type) == 'N') {
    dopingType = 'n';
    if (c > Small) {
      dopingConcentration = c;
    } else {
      std::cerr << className << "::SetDoping:\n";
      std::cerr << "    Doping concentration must be greater than zero.\n"; 
      std::cerr << "    Using default value for n-type silicon "
                << "(10^12 cm-3) instead.\n";
      dopingConcentration = 1.e12;
    }
  } else if (toupper(type) == 'P') {
    dopingType = 'p';
    if (c > Small) {
      dopingConcentration = c;
    } else {
      std::cerr << className << "::SetDoping:\n";
      std::cerr << "    Doping concentration must be greater than zero.\n"; 
      std::cerr << "    Using default value for p-type silicon "
                << "(10^18 cm-3) instead.\n";
      dopingConcentration = 1.e18;
    }
  } else if (toupper(type) == 'I') {
    dopingType = 'i';
    dopingConcentration = 0.;
  } else {
    std::cerr << className << "::SetDoping:\n";
    std::cerr << "    Unknown dopant type (" << type << ").\n";
    std::cerr << "    Available types are n, p and i (intrinsic).\n";
    return;
  }
  
  isChanged = true;
  
}

void 
MediumSilicon::GetDoping(char& type, double& c) const {

  type = dopingType; c = dopingConcentration;
  
}

void
MediumSilicon::SetTrapCrossSection(const double ecs, const double hcs) {

  if (ecs < 0.) {
    std::cerr << className << "::SetTrapCrossSection:\n";
    std::cerr << "    Capture cross-section [cm2] must positive.\n"; 
  } else {
    eTrapCs = ecs;
  }
  
  if (hcs < 0.) {
    std::cerr << className << "::SetTrapCrossSection:\n";
    std::cerr << "    Capture cross-section [cm2] must be positive.n"; 
  } else {
    hTrapCs = hcs;
  }
  
  trappingModel = 0;

}

void
MediumSilicon::SetTrapDensity(const double n) {

  if (n < 0.) {
    std::cerr << className << "::SetTrapDensity:\n";
    std::cerr << "    Trap density [cm-3] must be greater than zero.\n"; 
  } else {
    eTrapDensity = n;
    hTrapDensity = n;
  }
  
  trappingModel = 0;

}

void
MediumSilicon::SetTrappingTime(const double etau, const double htau) {

  if (etau <= 0.) {
    std::cerr << className << "::SetTrappingTime:\n";
    std::cerr << "    Trapping time [ns-1] must be positive.\n"; 
  } else {
    eTrapTime = etau;
  }
  
  if (htau <= 0.) {
    std::cerr << className << "::SetTrappingTime:\n";
    std::cerr << "    Trapping time [ns-1] must be positive.\n"; 
  } else {
    hTrapTime = htau;
  }
  
  trappingModel = 1;

}

bool 
MediumSilicon::ElectronVelocity(
            const double ex, const double ey, const double ez, 
            const double bx, const double by, const double bz, 
            double& vx, double& vy, double& vz) {

  if (isChanged) {
    UpdateTransportParameters();
    isChanged = false;
  }

  if (hasElectronVelocityE) {
    // Interpolation in user table.
    return Medium::ElectronVelocity(ex, ey, ez, bx, by, bz, vx, vy, vz);
  }

  const double e = sqrt(ex * ex + ey * ey + ez * ez);

  // Calculate the mobility
  double mu;
  switch (highFieldMobilityModel) {
    case HighFieldMobilityModelMinimos:
      ElectronMobilityMinimos(e, mu);
      break;
    case HighFieldMobilityModelCanali:
      ElectronMobilityCanali(e, mu);
      break;
    case HighFieldMobilityModelReggiani:
      ElectronMobilityReggiani(e, mu);
      break;
    default:
      mu = eMobility;
      break;
  }
  mu = -mu;

  // Hall mobility
  const double muH = eHallFactor * mu;
  // Compute drift velocity using the Langevin equation
  const double c1 = mu / (1. + muH * muH * (bx * bx + by * by + bz * bz));
  const double c2 = muH * muH * (bx * ex + by * ey + bz * ez);
  vx = c1 * (ex + muH * (ey * bz - ez * by) + c2 * bx);
  vy = c1 * (ey + muH * (ez * bx - ex * bz) + c2 * by);
  vz = c1 * (ez + muH * (ex * by - ey * bx) + c2 * bz);
  return true;

}

bool 
MediumSilicon::ElectronTownsend(
            const double ex, const double ey, const double ez,
            const double bx, const double by, const double bz,
            double& alpha) {
                         
  if (isChanged) {
    UpdateTransportParameters();
    isChanged = false;
  }
 
  if (hasElectronTownsend) {
    // Interpolation in user table.
    return Medium::ElectronTownsend(ex, ey, ez, bx, by, bz, alpha);
  }
 
  const double e = sqrt(ex * ex + ey * ey + ez * ez);
  
  switch (impactIonisationModel) {
    case ImpactIonisationModelVanOverstraeten:
      return ElectronImpactIonisationVanOverstraetenDeMan(e, alpha);
      break;
    case ImpactIonisationModelGrant:
      return ElectronImpactIonisationGrant(e, alpha);
      break;
    default:
      std::cerr << className << "::ElectronTownsend:\n";
      std::cerr << "    Unknown model activated. Program bug!\n";    
      break;
  }  
  return false;
  
}  

bool 
MediumSilicon::ElectronAttachment(
            const double ex, const double ey, const double ez,
            const double bx, const double by, const double bz,
            double& eta) {

  if (isChanged) {
    UpdateTransportParameters();
    isChanged = false;
  }

  if (hasElectronAttachment) {
    // Interpolation in user table.
    return Medium::ElectronAttachment(ex, ey, ez, bx, by, bz, eta);
  }
  
  switch (trappingModel) {
    case 0:
      eta = eTrapCs * eTrapDensity;
      break;
    case 1:
      double vx, vy, vz;
      ElectronVelocity(ex, ey, ez, bx, by, bz, vx, vy, vz);
      eta = eTrapTime * sqrt(vx * vx + vy * vy + vz * vz);
      if (eta > 0.) eta = 1. / eta;
      break;
    default:
      std::cerr << className << "::ElectronAttachment:\n";
      std::cerr << "    Unknown model activated. Program bug!\n";
      return false;
      break;      
  }
    
  return true;

}            

bool 
MediumSilicon::HoleVelocity(
            const double ex, const double ey, const double ez, 
            const double bx, const double by, const double bz, 
            double& vx, double& vy, double& vz) {

  if (isChanged) {
    UpdateTransportParameters();
    isChanged = false;
  }

  if (hasHoleVelocityE) {
    // Interpolation in user table.
    return Medium::HoleVelocity(ex, ey, ez, bx, by, bz, vx, vy, vz);
  }

  const double e = sqrt(ex * ex + ey * ey + ez * ez);

  // Calculate the mobility
  double mu;
  switch (highFieldMobilityModel) {
    case HighFieldMobilityModelMinimos:
      HoleMobilityMinimos(e, mu);
      break;
    case HighFieldMobilityModelCanali:
      HoleMobilityCanali(e, mu);
      break;
    case HighFieldMobilityModelReggiani:
      HoleMobilityReggiani(e, mu);
      break;
    default:
      mu = hMobility;
  }

  // Hall mobility
  const double muH = hHallFactor * mu;
  // Compute drift velocity using the Langevin equation
  const double c1 = mu / (1. + muH * muH * (bx * bx + by * by + bz * bz));
  const double c2 = muH * muH * (bx * ex + by * ey + bz * ez);
  vx = c1 * (ex + muH * (ey * bz - ez * by) + c2 * bx);
  vy = c1 * (ey + muH * (ez * bx - ex * bz) + c2 * by);
  vz = c1 * (ez + muH * (ex * by - ey * bx) + c2 * bz);
  return true;

}

bool 
MediumSilicon::HoleTownsend(
            const double ex, const double ey, const double ez,
            const double bx, const double by, const double bz,
            double& alpha) {
                         
  if (isChanged) {
    UpdateTransportParameters();
    isChanged = false;
  }
 
  if (hasHoleTownsend) {
    // Interpolation in user table.
    return Medium::HoleTownsend(ex, ey, ez, bx, by, bz, alpha);
  }
 
  const double e = sqrt(ex * ex + ey * ey + ez * ez);
  
  switch (impactIonisationModel) {
    case ImpactIonisationModelVanOverstraeten:
      return HoleImpactIonisationVanOverstraetenDeMan(e, alpha);
      break;
    case ImpactIonisationModelGrant:
      return HoleImpactIonisationGrant(e, alpha);
      break;
    default:
      std::cerr << className << "::HoleTownsend:\n";
      std::cerr << "    Unknown model activated. Program bug!\n";      
      break;
  }  
  return false;

}

bool 
MediumSilicon::HoleAttachment(
            const double ex, const double ey, const double ez,
            const double bx, const double by, const double bz,
            double& eta) {

  if (isChanged) {
    UpdateTransportParameters();
    isChanged = false;
  }

  if (hasHoleAttachment) {
    // Interpolation in user table.
    return Medium::HoleAttachment(ex, ey, ez, bx, by, bz, eta);
  }
  
  switch (trappingModel) {
    case 0:
      eta = hTrapCs * hTrapDensity;
      break;
    case 1:
      double vx, vy, vz;
      HoleVelocity(ex, ey, ez, bx, by, bz, vx, vy, vz);
      eta = hTrapTime * sqrt(vx * vx + vy * vy + vz * vz);
      if (eta > 0.) eta = 1. / eta;
      break;
    default:
      std::cerr << className << "::HoleAttachment:\n";
      std::cerr << "    Unknown model activated. Program bug!\n";
      return false;
      break;      
  }
  
  return true;

}

void
MediumSilicon::SetLowFieldMobility(const double mue, const double muh) {

  if (mue <= 0. || muh <= 0.) {
    std::cerr << className << "::SetLowFieldMobility:\n";
    std::cerr << "    Mobility must be greater than zero.\n";
    return;
  }
  
  eMobility = mue;
  hMobility = muh;
  hasUserMobility = true;
  isChanged = true;

}

void 
MediumSilicon::SetLatticeMobilityModelMinimos() {

  latticeMobilityModel = LatticeMobilityModelMinimos;  
  hasUserMobility = false;
  isChanged = true;
  
}

void 
MediumSilicon::SetLatticeMobilityModelSentaurus() {

  latticeMobilityModel = LatticeMobilityModelSentaurus;
  hasUserMobility = false;
  isChanged = true;
  
}

void
MediumSilicon::SetLatticeMobilityModelReggiani() {

  latticeMobilityModel = LatticeMobilityModelReggiani;
  hasUserMobility = false;
  isChanged = true;

}

void 
MediumSilicon::SetDopingMobilityModelMinimos() {

  dopingMobilityModel = DopingMobilityModelMinimos;
  hasUserMobility = false;
  isChanged = true;
  
}

void 
MediumSilicon::SetDopingMobilityModelMasetti() {

  dopingMobilityModel = DopingMobilityModelMasetti;
  hasUserMobility = false;
  isChanged = true;
  
}

void
MediumSilicon::SetSaturationVelocity(const double vsate, const double vsath) {

  if (vsate <= 0. || vsath <= 0.) {
    std::cout << className << "::SetSaturationVelocity:\n";
    std::cout << "    Restoring default values.\n";
    hasUserSaturationVelocity = false;
  } else {
    eSatVel = vsate; hSatVel = vsath;
    hasUserSaturationVelocity = true;
  }
  
  isChanged = true;

}

void
MediumSilicon::SetSaturationVelocityModelMinimos() {

  saturationVelocityModel = SaturationVelocityModelMinimos;
  hasUserSaturationVelocity = false;
  isChanged = true;

}

void
MediumSilicon::SetSaturationVelocityModelCanali() {

  saturationVelocityModel = SaturationVelocityModelCanali;
  hasUserSaturationVelocity = false;
  isChanged = true;

}

void
MediumSilicon::SetSaturationVelocityModelReggiani() {

  saturationVelocityModel = SaturationVelocityModelReggiani;
  hasUserSaturationVelocity = false;
  isChanged = true;

}

void 
MediumSilicon::SetHighFieldMobilityModelMinimos() {

  highFieldMobilityModel = HighFieldMobilityModelMinimos;
  isChanged = true;
  
}

void 
MediumSilicon::SetHighFieldMobilityModelCanali() {

  highFieldMobilityModel = HighFieldMobilityModelCanali;
  isChanged = true;
  
}

void
MediumSilicon::SetHighFieldMobilityModelReggiani() {

  highFieldMobilityModel = HighFieldMobilityModelReggiani;
  isChanged = true;

}

void
MediumSilicon::SetHighFieldMobilityModelConstant() {

  highFieldMobilityModel = HighFieldMobilityModelConstant;
  
}

void 
MediumSilicon::SetImpactIonisationModelVanOverstraetenDeMan() {

  impactIonisationModel = ImpactIonisationModelVanOverstraeten;
  isChanged = true;  

}

void
MediumSilicon::SetImpactIonisationModelGrant() {

  impactIonisationModel = ImpactIonisationModelGrant;
  isChanged = true;
 
}

bool
MediumSilicon::SetMaxElectronEnergy(const double e) {

  if (e <= Small) {
    std::cerr << className << "::SetMaxElectronEnergy:\n";
    std::cerr << "    Provided upper electron energy limit (" << e
              << " eV) is too small.\n";
    return false;
  }

  eFinal = e;
  // Determine the energy interval size
  eStep = eFinal / nEnergySteps;

  isChanged = true;

  return true;

}

double
MediumSilicon::GetElectronEnergy(
               const double px, const double py, const double pz, 
               double& vx, double& vy, double& vz, const int band) {

  double e = 0.;
  double a = 0.;

  // Non-parabolicity parameter
  double alpha = 0.;
  if (useNonParabolicity) {
    alpha = 0.5;
  }
  
  switch (band) {
    case 0:
    case 1:
      // X 100, -100  
      a = 0.5 * inverseElectronMass * (px * px / mLongX + 
                                      (py * py + pz * pz) / mTransX);
      if (useNonParabolicity) {
        e = 0.5 * (sqrt(1. + 4 * alpha * a) - 1.) / alpha;
      } else {
        e = a;
      }
      a = SpeedOfLight * inverseElectronMass / (1. + 2 * alpha * e); 
      vx = a * px / mLongX; vy = a * py / mTransX; vz = a * pz / mTransX;
      return e;
      break;
    case 2:
    case 3:
      // X 010, 0-10
      a = 0.5 * inverseElectronMass * (py * py / mLongX + 
                                      (px * px + pz * pz) / mTransX);
      if (useNonParabolicity) {
        e = 0.5 * (sqrt(1. + 4 * alpha * a) - 1.) / alpha;
      } else {
        e = a;
      }
      a = SpeedOfLight * inverseElectronMass / (1. + 2 * alpha * e); 
      vx = a * px / mTransX; vy = a * py / mLongX; vz = a * pz / mTransX;
      return e;
      break;
    case 4:
    case 5:
      // X 001, 00-1
      a = 0.5 * inverseElectronMass * (pz * pz / mLongX + 
                                      (px * px + py * py) / mTransX);
      if (useNonParabolicity) {
        e = 0.5 * (sqrt(1. + 4 * alpha * a) - 1.) / alpha;
      } else {
        e = a;
      } 
      a = SpeedOfLight * inverseElectronMass / (1. + 2 * alpha * e); 
      vx = a * px / mTransX; vy = a * py / mTransX; vz = a * pz / mLongX;
      return e;
      break;     
    default:
      std::cerr << className << "::GetElectronEnergy:\n";
      std::cerr << "    Unknown band index " << band << "!\n";
      break;
  }
  
  vx = SpeedOfLight * inverseElectronMass * px;
  vy = SpeedOfLight * inverseElectronMass * py;
  vz = SpeedOfLight * inverseElectronMass * pz;
  return 0.5 * inverseElectronMass * (px * px + py * py + pz * pz);

}

void
MediumSilicon::GetElectronMomentum(const double e,
                                   double& px, double& py, double& pz,
                                   const int band) {

  double alpha = 0.;
  if (useNonParabolicity) alpha = 0.5;

  const double pstar = sqrt(2. * ElectronMass * e * (1. + alpha * e));
  const double ctheta = 1. - 2. * RndmUniform();
  const double stheta = sqrt(1. - ctheta * ctheta);
  const double phi = TwoPi * RndmUniform();
  
  if (band >= 0 && band < 6) {
    const double pl = pstar * sqrt(mLongX);
    const double pt = pstar * sqrt(mTransX);
    switch (band) {
      case 0:
      case 1:
        // 100
        px = pl * cos(phi) * stheta;
        py = pt * sin(phi) * stheta;
        pz = pt * ctheta;
        return;
        break;
      case 2:
      case 3:
        // 010
        px = pt * cos(phi) * stheta;
        py = pl * sin(phi) * stheta;
        pz = pt * ctheta;
        return;
        break;
      case 4:
      case 5:
        // 001
        px = pt * cos(phi) * stheta;
        py = pt * sin(phi) * stheta;
        pz = pl * ctheta;
        return;
        break;
      default:
        break;
    }
  } 
  
  std::cerr << className << "::GetElectronMomentum:\n";
  std::cerr << "    Unknown band index " << band << "!\n";
  
  px = pstar * stheta * cos(phi);
  py = pstar * stheta * sin(phi);
  pz = pstar * ctheta;

}

double
MediumSilicon::GetElectronNullCollisionRate() {

  if (isChanged) {
    if (!UpdateTransportParameters()) {
      std::cerr << className << "::GetElectronNullCollisionRate:\n";
      std::cerr << "    Error calculating the collision rates table.\n";
      return 0.;
    }
    isChanged = false;
  }

  return cfNullElectrons;

}

double
MediumSilicon::GetElectronCollisionRate(const double e, const int band) {

  if (e <= 0.) {
    std::cerr << className << "::GetElectronCollisionRate:\n";
    std::cerr << "    Electron energy must be positive.\n";
    return cfTotElectronsX[0];
  }

  if (e > eFinal) {
    std::cerr << className << "::GetElectronCollisionRate:\n";
    std::cerr << "    Collision rate at " << e
              << " eV is not included in the current table.\n";
    std::cerr << "    Increasing energy range to " << 1.05 * e
              << " eV.\n";
    SetMaxElectronEnergy(1.05 * e);
  }

  if (isChanged) {
    if (!UpdateTransportParameters()) {
      std::cerr << className << "::GetElectronCollisionRate:\n";
      std::cerr << "    Error calculating the collision rates table.\n";
      return 0.;
    }
    isChanged = false;
  }

  if (band >= 0 && band < 6) {
    return cfTotElectronsX[int(e / eStep)];
  }
  
  std::cerr << className << "::GetElectronCollisionRate:\n";
  std::cerr << "    Unknown band!\n";
  return 0.;
  
}

bool
MediumSilicon::GetElectronCollision(const double e, 
                                    int& type, int& level, double& e1, 
                                    double& px, double& py, double& pz, 
                                    int& nsec, double& esec,
                                    int& band) {

  if (e > eFinal) {
    std::cerr << className << "::GetElectronCollision:\n";
    std::cerr << "    Requested electron energy (" << e;
    std::cerr << " eV) exceeds current energy range (" << eFinal;
    std::cerr << " eV).\n";
    std::cerr << "    Increasing energy range to " << 1.05 * e
              << " eV.\n";
    SetMaxElectronEnergy(1.05 * e);
  } else if (e <= 0.) {
    std::cerr << className << "::GetElectronCollision:\n";
    std::cerr << "    Electron energy must be greater than zero.\n";
    return false;
  }

  if (isChanged) {
    if (!UpdateTransportParameters()) {
      std::cerr << className << "::GetElectronCollision:\n";
      std::cerr << "    Error calculating the collision rates table.\n";
      return false;
    }
    isChanged = false;
  }

  // Get the energy interval.
  int iE = int(e / eStep);
  if (iE >= nEnergySteps) iE = nEnergySteps - 1;

  // Sample the scattering process.
  if (band >= 0 && band < 6) {
    const double r = RndmUniform();
    int iLow = 0; 
    int iUp = nLevelsX - 1;
    if (r <= cfElectronsX[iE][iLow]) {
      level = iLow;
    } else if (r >= cfElectronsX[iE][nLevelsX - 1]) {
      level = iUp;
    } else {
      int iMid;
      while (iUp - iLow > 1) {
        iMid = (iLow + iUp) >> 1;
        if (r < cfElectronsX[iE][iMid]) {
          iUp = iMid;
        } else {
          iLow = iMid;
        }
      }
      level = iUp;
    }

    // Get the collision type.
    type = scatTypeElectronsX[level];
    if (type == 11 || type == 12) {
      // Intervalley scattering
      switch (band) {
        case 0: 
          band = 1; break;
        case 1:
          band = 0; break;
        case 2:
          band = 3; break;
        case 3:
          band = 2; break;
        case 4:
          band = 5; break;
        case 5:
          band = 4; break;
        default:
          break;
      }
    } else if (type == 13 || type == 14) {
      // Intervalley scattering
      switch (band) {
        case 0:
        case 1:
          band = int(RndmUniform() * 4) + 2;
          break;
        case 2:
        case 3:
          band = int(RndmUniform() * 4);
          if (band > 1) band += 2;
          break;
        case 4:
        case 5:
          band = int(RndmUniform() * 4);          
          break;
      }
    }
      
    // Get the energy loss.
    double loss = energyLossElectronsX[level];
    // Secondary electron energy (none by default)
    esec = 0.;
    nsec = 0;
    // Ionising collision
    if (type == 1) {
      esec = RndmUniform() * (e - loss);
      loss += esec;
      if (esec < Small) esec = Small;
      nsec = 1;
    }

    if (e < loss) loss = e - 0.0001;
    // Update the energy.
    e1 = e - loss;
    if (e1 < Small) e1 = Small;
  
    // Update the momentum.
    double alpha = 0.;
    if (useNonParabolicity) alpha = 0.5;

    const double pstar = sqrt(2. * ElectronMass * e1 * (1. + alpha * e1));
    const double ctheta = 1. - 2. * RndmUniform();
    const double stheta = sqrt(1. - ctheta * ctheta);
    const double phi = TwoPi * RndmUniform();

    const double pl = pstar * sqrt(mLongX);
    const double pt = pstar * sqrt(mTransX);
    switch (band) {
      case 0:
      case 1:
        // 100
        px = pl * cos(phi) * stheta;
        py = pt * sin(phi) * stheta;
        pz = pt * ctheta;
        break;
      case 2:
      case 3:
        // 010
        px = pt * cos(phi) * stheta;
        py = pl * sin(phi) * stheta;
        pz = pt * ctheta;
        break;
      case 4:
      case 5:
        // 001
        px = pt * cos(phi) * stheta;
        py = pt * sin(phi) * stheta;
        pz = pl * ctheta;
        break;
      default:
        return false;
        break;
    }
    
    return true;

  } 
  
  std::cerr << className << "::GetElectronCollision:\n";
  std::cerr << "   Unknown band!\n";
  e1 = e;
  type = 0;
  return false;

}

bool 
MediumSilicon::GetOpticalDataRange(double& emin, double& emax, const int i) {

  if (i != 0) {
    std::cerr << className << "::GetOpticalDataRange:\n";
    std::cerr << "    Medium has only one component.\n";
  }

  // Make sure the optical data table has been loaded.
  if (!hasOpticalData) {
    if (!LoadOpticalData(opticalDataFile)) {
      std::cerr << className << "::GetOpticalDataRange:\n";
      std::cerr << "    Optical data table could not be loaded.\n";
      return false;
    }
    hasOpticalData = true;
  }
   
  emin = opticalDataTable[0].energy;
  emax = opticalDataTable.back().energy;
  if (debug) {
    std::cout << className << "::GetOpticalDataRange:\n";
    std::cout << "    " << emin << " < E [eV] < " << emax << "\n";
  }
  return true;  
  
}

bool 
MediumSilicon::GetDielectricFunction(const double e, 
                                     double& eps1, double& eps2, 
                                     const int i) {
                        
  if (i != 0) {
    std::cerr << className << "::GetDielectricFunction:\n";
    std::cerr << "    Medium has only one component.\n";
    return false;
  }
                        
  // Make sure the optical data table has been loaded.
  if (!hasOpticalData) {
    if (!LoadOpticalData(opticalDataFile)) {
      std::cerr << className << "::GetDielectricFunction:\n";
      std::cerr << "    Optical data table could not be loaded.\n";
      return false;
    }
    hasOpticalData = true;
  }
  
  // Make sure the requested energy is within the range of the table.
  const double emin = opticalDataTable[0].energy;
  const double emax = opticalDataTable.back().energy;    
  if (e < emin || e > emax) {
    std::cerr << className << "::GetDielectricFunction:\n";
    std::cerr << "    Requested energy (" << e << " eV) " 
              << " is outside the range of the optical data table.\n";
    std::cerr << "    " << emin << " < E [eV] < " << emax << "\n";
    eps1 = eps2 = 0.;
    return false;
  }

  // Locate the requested energy in the table.
  int iLow = 0;
  int iUp = opticalDataTable.size() - 1;
  int iM;
  while (iUp - iLow > 1) {
    iM = (iUp + iLow) >> 1;
    if (e >= opticalDataTable[iM].energy) {
      iLow = iM;
    } else {
      iUp = iM;
    }
  }
  
  // Interpolate the real part of dielectric function.
  // Use linear interpolation if one of the values is negative,
  // Otherwise use log-log interpolation.
  const double logX0 = log(opticalDataTable[iLow].energy);
  const double logX1 = log(opticalDataTable[iUp].energy);
  const double logX = log(e);
  if (opticalDataTable[iLow].eps1 <= 0. || opticalDataTable[iUp].eps1 <= 0.) {
    eps1 = opticalDataTable[iLow].eps1 + (e - opticalDataTable[iLow].energy) * 
           (opticalDataTable[iUp].eps1 - opticalDataTable[iLow].eps1) / 
          (opticalDataTable[iUp].energy - opticalDataTable[iLow].energy);  
  } else {
    const double logY0 = log(opticalDataTable[iLow].eps1);
    const double logY1 = log(opticalDataTable[iUp].eps1);
    eps1 = logY0 + (logX - logX0) * (logY1 - logY0) / (logX1 - logX0);
    eps1 = exp(eps1);
  }
      
  // Interpolate the imaginary part of dielectric function,
  // using log-log interpolation.
  const double logY0 = log(opticalDataTable[iLow].eps2);
  const double logY1 = log(opticalDataTable[iUp].eps2);  
  eps2 = logY0 + (log(e) - logX0) * (logY1 - logY0) / (logX1 - logX0);
  eps2 = exp(eps2);
  return true;
  
}

bool 
MediumSilicon::UpdateTransportParameters() {

  // Calculate impact ionisation coefficients
  switch (impactIonisationModel) {
    case ImpactIonisationModelVanOverstraeten:
      UpdateImpactIonisationVanOverstraetenDeMan();
      break;
    case ImpactIonisationModelGrant:
      UpdateImpactIonisationGrant();
      break;
    default:
      std::cerr << className << "::UpdateTransportParameters:\n";
      std::cerr << "    Unknown impact ionisation model. Bug!\n";
      break;
  }
  
  if (!hasUserMobility) {        
    // Calculate lattice mobility
    switch (latticeMobilityModel) {
      case LatticeMobilityModelMinimos:
        UpdateLatticeMobilityMinimos();
        break;
      case LatticeMobilityModelSentaurus:
        UpdateLatticeMobilitySentaurus();
        break;   
      case LatticeMobilityModelReggiani:
        UpdateLatticeMobilityReggiani();
        break; 
      default:
        std::cerr << className << "::UpdateTransportParameters:\n"; 
        std::cerr << "    Unknown lattice mobility model.\n";
        std::cerr << "    Program bug!\n"; 
        break;
    }

  
    // Calculate doping mobility
    switch (dopingMobilityModel) {
      case DopingMobilityModelMinimos:
        UpdateDopingMobilityMinimos();
        break;
      case DopingMobilityModelMasetti:
        UpdateDopingMobilityMasetti();
        break;
      default:
        std::cerr << className << "::UpdateTransportParameters:\n"; 
        std::cerr << "    Unknown doping mobility model.\n";
        std::cerr << "    Program bug!\n"; 
        break;
    }
  }
      
  // Calculate saturation velocity
  if (!hasUserSaturationVelocity) {
    switch (saturationVelocityModel) {
      case SaturationVelocityModelMinimos:
        UpdateSaturationVelocityMinimos();
        break;
      case SaturationVelocityModelCanali:
        UpdateSaturationVelocityCanali();
        break;
      case SaturationVelocityModelReggiani:
        UpdateSaturationVelocityReggiani();
        break;
    }
  }

  // Calculate high field saturation parameters
  switch (highFieldMobilityModel) {
    case HighFieldMobilityModelCanali:
      UpdateHighFieldMobilityCanali();
      break;
  }
  
  if (debug) {
    std::cout << className << "::UpdateTransportParameters:\n";
    std::cout << "    Low-field mobility [cm2 V-1 ns-1]\n";
    std::cout << "      Electrons: " << eMobility << "\n";
    std::cout << "      Holes:     " << hMobility << "\n";
    if (highFieldMobilityModel > 2) {
      std::cout << "    Mobility is not field-dependent.\n";
    } else {
      std::cout << "    Saturation velocity [cm / ns]\n";
      std::cout << "      Electrons: " << eSatVel << "\n";
      std::cout << "      Holes:     " << hSatVel << "\n";
    }
  }

  if (!ElectronScatteringRates()) return false;
  return true;

}

void 
MediumSilicon::UpdateLatticeMobilityMinimos() {

  // References:
  // - S. Selberherr, W. Haensch, M. Seavey, J. Slotboom,
  //   Solid State Electronics 33 (1990), 1425-1436
  // - Minimos 6.1 User's Guide (1999)

  // Lattice mobilities at 300 K [cm2 / (V ns)]
  const double eMu0 = 1.43e-6; 
  const double hMu0 = 0.46e-6;
  // Temperature normalized to 300 K
  const double t = temperature / 300.;
  // Temperature dependence of lattice mobility
  eLatticeMobility = eMu0 * pow(t, -2.);
  hLatticeMobility = hMu0 * pow(t, -2.18);
  
}

void 
MediumSilicon::UpdateLatticeMobilitySentaurus() {

  // References:
  // - C. Lombardi et al.,
  //   IEEE Trans. CAD 7 (1988), 1164-1171
  // - Sentaurus Device User Guide (2007)

  // Lattice mobilities at 300 K [cm2 / (V ns)]
  const double eMu0 = 1.417e-6; 
  const double hMu0 = 0.4705e-6;
  // Temperature normalized to 300 K
  const double t = temperature / 300.;
  // Temperature dependence of lattice mobility
  eLatticeMobility = eMu0 * pow(t, -2.5);
  hLatticeMobility = hMu0 * pow(t, -2.2);
  
}

void
MediumSilicon::UpdateLatticeMobilityReggiani() {

  // Reference:
  // - M. A. Omar, L. Reggiani
  //   Solid State Electronics 30 (1987), 693-697

  // Lattice mobilities at 300 K [cm2 / (V ns)]
  const double eMu0 = 1.320e-6; 
  const double hMu0 = 0.460e-6;
  // Temperature normalized to 300 K
  const double t = temperature / 300.;
  // Temperature dependence of lattice mobility
  eLatticeMobility = eMu0 * pow(t, -2.);
  hLatticeMobility = hMu0 * pow(t, -2.2);

}

void
MediumSilicon::UpdateDopingMobilityMinimos() {

  // References:
  // - S. Selberherr, W. Haensch, M. Seavey, J. Slotboom,
  //   Solid State Electronics 33 (1990), 1425-1436
  // - Minimos 6.1 User's Guide (1999)

  // Mobility reduction due to ionised impurity scattering
  // Surface term not taken into account
  double eMuMin = 0.080e-6;
  double hMuMin = 0.045e-6;
  if (temperature > 200.) {
    eMuMin *= pow(temperature / 300., -0.45);
    hMuMin *= pow(temperature / 300., -0.45);
  } else {
    eMuMin *= pow(2. / 3., -0.45) * pow(temperature / 200., -0.15);
    hMuMin *= pow(2. / 3., -0.45) * pow(temperature / 200., -0.15);
  }
  const double eRefC = 1.12e17 * pow(temperature / 300., 3.2);
  const double hRefC = 2.23e17 * pow(temperature / 300., 3.2);
  const double alpha = 0.72 * pow(temperature / 300., 0.065);
  // Assume impurity concentration equal to doping concentration
  eMobility =  eMuMin + (eLatticeMobility - eMuMin) / 
               (1. + pow(dopingConcentration / eRefC, alpha));
  hMobility =  hMuMin + (hLatticeMobility - hMuMin) / 
               (1. + pow(dopingConcentration / hRefC, alpha));

}

void 
MediumSilicon::UpdateDopingMobilityMasetti() {

  // Reference:
  // - G. Masetti, M. Severi, S. Solmi,
  //   IEEE Trans. Electron Devices 30 (1983), 764-769
  // - Sentaurus Device User Guide (2007)
  // - Minimos NT User Guide (2004)

  if (dopingConcentration < 1.e13) {
    eMobility = eLatticeMobility;
    hMobility = hLatticeMobility;
    return;
  }
  
  // Parameters adopted from Minimos NT User Guide
  const double eMuMin1 = 0.0522e-6;
  const double eMuMin2 = 0.0522e-6;
  const double eMu1    = 0.0434e-6;
  const double hMuMin1 = 0.0449e-6;
  const double hMuMin2 = 0.;
  const double hMu1    = 0.029e-6;
  const double eCr = 9.68e16;
  const double eCs = 3.42e20;
  const double hCr = 2.23e17;
  const double hCs = 6.10e20;
  const double hPc = 9.23e16;
  const double eAlpha = 0.68;
  const double eBeta = 2.;
  const double hAlpha = 0.719;
  const double hBeta = 2.;
  
  eMobility = eMuMin1 + (eLatticeMobility - eMuMin2) / 
                        (1. + pow(dopingConcentration / eCr, eAlpha)) - eMu1 / 
                        (1. + pow(eCs / dopingConcentration, eBeta));

  hMobility = hMuMin1 * exp(-hPc / dopingConcentration) + 
                        (hLatticeMobility - hMuMin2) / 
                        (1. + pow(dopingConcentration / hCr, hAlpha)) - hMu1 / 
                        (1. + pow(hCs / dopingConcentration, hBeta));

}
  

void 
MediumSilicon::UpdateSaturationVelocityMinimos() {

  // References:
  // - R. Quay, C. Moglestue, V. Palankovski, S. Selberherr,
  //   Materials Science in Semiconductor Processing 3 (2000), 149-155
  // - Minimos NT User Guide (2004)
  
  // Temperature-dependence of saturation velocities [cm / ns]
  eSatVel = 1.e-2    / (1. + 0.74 * (temperature / 300. - 1.));
  hSatVel = 0.704e-2 / (1. + 0.37 * (temperature / 300. - 1.));
  
}

void 
MediumSilicon::UpdateSaturationVelocityCanali() {

  // References:
  // - C. Canali, G. Majni, R. Minder, G. Ottaviani,
  //   IEEE Transactions on Electron Devices 22 (1975), 1045-1047
  // - Sentaurus Device User Guide (2007)

  eSatVel = 1.07e-2 * pow(300. / temperature, 0.87);
  hSatVel = 8.37e-3 * pow(300. / temperature, 0.52);

}

void
MediumSilicon::UpdateSaturationVelocityReggiani() {

  // Reference:
  // - M. A. Omar, L. Reggiani
  //   Solid State Electronics 30 (1987), 693-697

  eSatVel = 1.470e-2 * sqrt(tanh(150. / temperature));
  hSatVel = 0.916e-2 * sqrt(tanh(300. / temperature));

}

void
MediumSilicon::UpdateHighFieldMobilityCanali() {

  // References:
  // - C. Canali, G. Majni, R. Minder, G. Ottaviani,
  //   IEEE Transactions on Electron Devices 22 (1975), 1045-1047
  // - Sentaurus Device User Guide (2007)

  // Temperature dependent exponent in high-field mobility formula
  eBetaCanali = 1.109 * pow(temperature / 300., 0.66);
  hBetaCanali = 1.213 * pow(temperature / 300., 0.17);

}

void 
MediumSilicon::UpdateImpactIonisationVanOverstraetenDeMan() {

  // References:
  //  - R. van Overstraeten and H. de Man, 
  //    Solid State Electronics 13 (1970), 583-608
  //  - W. Maes, K. de Meyer and R. van Overstraeten, 
  //    Solid State Electronics 33 (1990), 705-718
  // - Sentaurus Device User Guide (2007)

  // Temperature dependence as in Sentaurus Device
  // Optical phonon energy
  const double hbarOmega = 0.063;
  // Temperature scaling coefficient
  const double gamma = tanh(hbarOmega / (2. * BoltzmannConstant * 300.)) /
                       tanh(hbarOmega / (2. * BoltzmannConstant * temperature));
  
  // Low field coefficients taken from Maes, de Meyer, van Overstraeten
  eImpactA0 = gamma * 3.318e5; eImpactB0 = gamma * 1.135e6;
  eImpactA1 = gamma * 7.03e5;  eImpactB1 = gamma * 1.231e6;
  
  hImpactA0 = gamma * 1.582e6; hImpactB0 = gamma * 2.036e6;
  hImpactA1 = gamma * 6.71e5;  hImpactB1 = gamma * 1.693e6;

}

void 
MediumSilicon::UpdateImpactIonisationGrant() {

  // References:
  // - W. N. Grant, 
  //   Solid State Electronics 16 (1973), 1189 - 1203
  // - Sentaurus Device User Guide (2007)

  // Temperature dependence as in Sentaurus Device
  // Optical phonon energy
  double hbarOmega = 0.063;
  // Temperature scaling coefficient
  double gamma = tanh(hbarOmega / (2. * BoltzmannConstant * 300.)) /
                 tanh(hbarOmega / (2. * BoltzmannConstant * temperature));
                 
  eImpactA0 = 2.60e6 * gamma; eImpactB0 = 1.43e6 * gamma;
  eImpactA1 = 6.20e5 * gamma; eImpactB1 = 1.08e6 * gamma;
  eImpactA2 = 5.05e5 * gamma; eImpactB2 = 9.90e5 * gamma;
  
  hImpactA0 = 2.00e6 * gamma; hImpactB0 = 1.97e6 * gamma;
  hImpactA1 = 5.60e5 * gamma; hImpactB1 = 1.32e6 * gamma;
  
}

bool 
MediumSilicon::ElectronMobilityMinimos(const double e, double& mu) const {

  // Reference:
  // - Minimos User's Guide (1999)
  
  if (e < Small) {
    mu = 0.;
  } else {
    mu = 2. * eMobility / 
         (1. + sqrt(1. + pow(2. * eMobility * e / eSatVel, 2.)));
  }  
  return true;
  
}

bool 
MediumSilicon::ElectronMobilityCanali(const double e, double& mu) const {

  // Reference:
  // - Sentaurus Device User Guide (2007)
  
  if (e < Small) {
    mu = 0.;
  } else {
    mu = eMobility / 
         pow(1. + pow(eMobility * e / eSatVel, eBetaCanali), 1. / eBetaCanali);
  }  
  return true;
  
}

bool
MediumSilicon::ElectronMobilityReggiani(const double e, double& mu) const {

  // Reference:
  // - M. A. Omar, L. Reggiani
  //   Solid State Electronics 30 (1987), 693-697

  if (e < Small) {
    mu = 0.;
  } else {
    mu = eMobility / pow(1 + pow(eMobility * e / eSatVel, 1.5), 1. / 1.5);
  }
  return true;

}

bool 
MediumSilicon::ElectronImpactIonisationVanOverstraetenDeMan(const double e, 
                                                         double& alpha) const {

  // References:
  //  - R. van Overstraeten and H. de Man, 
  //    Solid State Electronics 13 (1970), 583-608
  //  - W. Maes, K. de Meyer and R. van Overstraeten, 
  //    Solid State Electronics 33 (1990), 705-718

  if (e < Small) {
    alpha = 0.;
  } else if (e < 1.2786e5) {
    alpha = eImpactA0 * exp(-eImpactB0 / e);
  } else {
    alpha = eImpactA1  * exp(-eImpactB1 / e);
  }
  return true;

}

bool 
MediumSilicon::ElectronImpactIonisationGrant(const double e, 
                                             double& alpha) const {

  // Reference:
  //  - W. N. Grant, Solid State Electronics 16 (1973), 1189 - 1203

  if (e < Small) {
    alpha = 0.;
  } else if (e < 2.4e5) {
    alpha = eImpactA0 * exp(-eImpactB0 / e);
  } else if (e < 5.3e5) {
    alpha = eImpactA1 * exp(-eImpactB1 / e);
  } else {
    alpha = eImpactA2 * exp(-eImpactB2 / e);
  }
  return true;

}

bool 
MediumSilicon::HoleMobilityMinimos(const double e, double& mu) const {

  // Reference:
  // - Minimos User's Guide (1999)

  if (e < Small) {
    mu = 0.;
  } else {
    mu = hMobility / (1. + hMobility * e / eSatVel);
  }  
  return true;
  
}

bool 
MediumSilicon::HoleMobilityCanali(const double e, double& mu) const {

  // Reference:
  // - Sentaurus Device User Guide (2007)
  
  if (e < Small) {
    mu = 0.;
  } else {
    mu = hMobility / 
         pow(1. + pow(hMobility * e / hSatVel, hBetaCanali), 1. / hBetaCanali);
  }  
  return true;
  
}

bool 
MediumSilicon::HoleMobilityReggiani(const double e, double& mu) const {

  // Reference:
  // - M. A. Omar, L. Reggiani
  //   Solid State Electronics 30 (1987), 693-697

  if (e < Small) {
    mu = 0.;
  } else {
    mu = hMobility / 
         pow(1. + pow(hMobility * e / hSatVel, 2.), 0.5);
  }
  return true;

}

bool 
MediumSilicon::HoleImpactIonisationVanOverstraetenDeMan(const double e, 
                                                       double& alpha) const {

  // Reference:
  //  - R. van Overstraeten and H. de Man, 
  //    Solid State Electronics 13 (1970), 583-608

  if (e < Small) {
    alpha = 0.;
  } else {
    alpha = hImpactA1 * exp(-hImpactB1 / e);
  }
  return true;

}

bool 
MediumSilicon::HoleImpactIonisationGrant(const double e, double& alpha) const {

  // Reference:
  //  - W. N. Grant, Solid State Electronics 16 (1973), 1189 - 1203

  if (e < Small) {
    alpha = 0.;
  } else if (e < 5.3e5) {
    alpha = hImpactA0 * exp(-hImpactB0 / e);
  } else {
    alpha = hImpactA1 * exp(-hImpactB1 / e);
  }
  return true;

}


bool 
MediumSilicon::LoadOpticalData(const std::string filename) {

  // Get the path to the data directory.
  char* pPath = getenv("GARFIELD_HOME");
  if (pPath == 0) {
    std::cerr << className << "::LoadOpticalData:\n";
    std::cerr << "    Environment variable GARFIELD_HOME is not set.\n"; 
    return false;
  }
  std::string filepath = pPath;
  filepath = filepath + "/Data/" + filename;

  // Open the file.
  std::ifstream infile;
  infile.open(filepath.c_str(), std::ios::in);
  // Make sure the file could actually be opened.
  if (!infile) {
    std::cerr << className << "::LoadOpticalData:\n";
    std::cerr << "    Error opening file " << filename << ".\n";
    return false;
  }
  
  // Clear the optical data table.
  opticalDataTable.clear();
  
  double lastEnergy = -1.;
  double energy, eps1, eps2, loss;  
  opticalData data;
  // Read the file line by line.
  std::string line;
  std::istringstream dataStream;  
  int i = 0;
  while (!infile.eof()) {
    ++i;
    // Read the next line.
    std::getline(infile, line);
    // Strip white space from the beginning of the line.
    line.erase(line.begin(), std::find_if(line.begin(), line.end(), 
               not1(std::ptr_fun<int, int>(isspace))));
    // Skip comments.
    if (line[0] == '#' ||
        (line[0] == '/' && line[1] == '/')) continue;
    // Extract the values.
    dataStream.str(line);
    dataStream >> energy >> eps1 >> eps2 >> loss;
    if (dataStream.eof()) break;
    // Check if the data has been read correctly.
    if (infile.fail()) {
      std::cerr << className << "::LoadOpticalData:\n";
      std::cerr << "    Error reading file "
                << filename << " (line " << i << ").\n";
      return false;
    }
    // Reset the stringstream.
    dataStream.str("");
    dataStream.clear();
    // Make sure the values make sense.
    // The table has to be in ascending order
    //  with respect to the photon energy.
    if (energy <= lastEnergy) {
      std::cerr << className << "::LoadOpticalData:\n";
      std::cerr << "    Table is not in monotonically " 
                << "increasing order (line " << i << ").\n";
      std::cerr << "    " << lastEnergy << "  " << energy << "  " << eps1 << "  " << eps2 << "\n";
      return false;
    }
    // The imaginary part of the dielectric function has to be positive.
    if (eps2 < 0.) {
      std::cerr << className << "::LoadOpticalData:\n";
      std::cerr << "    Negative value of the loss function "
                << "(line " << i << ").\n";
      return false;
    }
    // Ignore negative photon energies.
    if (energy <= 0.) continue;
    // Add the values to the list.
    data.energy = energy;
    data.eps1 = eps1;
    data.eps2 = eps2;
    opticalDataTable.push_back(data);
    lastEnergy = energy;
  }
  
  const int nEntries = opticalDataTable.size();
  if (nEntries <= 0) {
    std::cerr << className << "::LoadOpticalData:\n";
    std::cerr << "    Import of data from file " << filepath << "failed.\n";
    std::cerr << "    No valid data found.\n";
    return false;
  }
  
  if (debug) {
    std::cout << className << "::LoadOpticalData:\n";
    std::cout << "    Read " << nEntries << " values from file " 
              << filepath << "\n";
  }
  return true;

}


bool
MediumSilicon::ElectronScatteringRates() {

  // Reset the scattering rates
  cfTotElectronsX.resize(nEnergySteps);
  cfElectronsX.resize(nEnergySteps);
  for (int i = nEnergySteps; i--;) {
    cfTotElectronsX[i] = 0.;
    cfElectronsX[i].clear();
  }
  energyLossElectronsX.clear();
  scatTypeElectronsX.clear();
  cfNullElectrons = 0.;
  
  nLevelsX = 0;
  // Fill the scattering rates table
  ElectronAcousticScatteringRatesX();
  // ElectronImpurityScatteringRates();
  ElectronIntervalleyScatteringRatesXX();
  ElectronIonisationRates();

  std::ofstream outfile;  
  if (useCfOutput) outfile.open("rates.txt", std::ios::out);

  for (int i = 0; i < nEnergySteps; ++i) {
    // Sum up the scattering rates of all processes 
    for (int j = nLevelsX; j--;) cfTotElectronsX[i] += cfElectronsX[i][j];
    
    if (useCfOutput) {
      outfile << i * eStep << " " << cfTotElectronsX[i] << " ";
      for (int j = 0; j < nLevelsX; ++j) outfile << cfElectronsX[i][j] << " ";
      outfile << "\n";
    }

    if (cfTotElectronsX[i] > cfNullElectrons) {
      cfNullElectrons = cfTotElectronsX[i];
    }
   
    // Make sure the total scattering rate is positive
    if (cfTotElectronsX[i] <= 0.) { 
      std::cerr << className << "::ElectronScatteringRates:\n";
      std::cerr << "    Scattering rate at " << i * eStep << " <= 0.\n"; 
      return false;
    }

    // Normalise the rates
    for (int j = 0; j < nLevelsX; ++j) {
      cfElectronsX[i][j] /= cfTotElectronsX[i];
      if (j > 0) cfElectronsX[i][j] += cfElectronsX[i][j - 1];
    }  

  }
  if (useCfOutput) outfile.close();

  return true;

}

bool 
MediumSilicon::ElectronAcousticScatteringRatesX() {

  // Reference:
  //  - C. Jacoboni and L. Reggiani,
  //    Rev. Mod. Phys. 55, 645-705

  // Mass density [(eV/c2)/cm3]
  const double rho = density * atomicWeight * AtomicMassUnitElectronVolt;
  // Lattice temperature [eV]
  const double kbt = BoltzmannConstant * temperature;  

  // Acoustic phonon intraband scattering  
  // Acoustic deformation potential [eV]
  const double defpot = 9.;
  // Longitudinal and transverse velocity of sound [cm/ns]
  const double ut = 9.0e-4;
  const double ul = 5.3e-4;
  // Average velocity of sound [cm/ns]
  const double u = (ul + 2. * ut) / 3.;    
  // Prefactor for acoustic deformation potential scattering
  const double cIntra = TwoPi * SpeedOfLight * SpeedOfLight * 
                        kbt * defpot * defpot /
                        (Hbar * u * u * rho); 
  
  double en = Small;
  for (int i = 0; i < nEnergySteps; ++i) {
    cfElectronsX[i].push_back(cIntra * 
                              GetConductionBandDensityOfStates(en, 0));
    en += eStep;
  }
  
  energyLossElectronsX.push_back(0.);
  scatTypeElectronsX.push_back(10);
  ++nLevelsX;

  return true;

}

bool 
MediumSilicon::ElectronIntervalleyScatteringRatesXX() {

  // Reference:
  //  - C. Jacoboni and L. Reggiani,
  //    Rev. Mod. Phys. 55, 645-705

  // Mass density [(eV/c)/cm3]
  const double rho = density * atomicWeight * AtomicMassUnitElectronVolt;
  // Lattice temperature [eV]
  const double kbt = BoltzmannConstant * temperature;  

  const int nPhonons = 6;
  // f-type scattering: transition between orthogonal axes (multiplicity 4)
  // g-type scattering: transition between opposite axes (multiplicity 1)
  // Sequence of transitions in the table:
  // TA (g) - LA (g) - LO (g) - TA (f) - LA (f) - TO (f)
  // Coupling constants [eV/cm]
  const double dtk[nPhonons] = {0.5e8, 0.8e8, 1.1e9, 
                                0.3e8, 2.0e8, 2.0e8};
  // Phonon energies [eV]
  const double eph[nPhonons] = {12.06e-3, 18.53e-3, 62.04e-3, 
                                18.86e-3, 47.39e-3, 59.03e-3};
  // Phonon cccupation numbers
  double nocc[nPhonons];
  // Prefactors
  const double c0 = HbarC * SpeedOfLight * Pi / rho;
  double c[nPhonons];

  for (int j = 0; j < 6; ++j) {
    nocc[j] = 1. / (exp(eph[j] / kbt) - 1);
    c[j] = c0 * dtk[j] * dtk[j] / eph[j];
    if (j > 2) c[j] *= 4;
  }

  double en = 0.;
  double dos = 0.;
  for (int i = 0; i < nEnergySteps; ++i) {
    for (int j = 0; j < nPhonons; ++j) {
     // Absorption
      dos = GetConductionBandDensityOfStates(en + eph[j], 0);
      cfElectronsX[i].push_back(c[j] * nocc[j] * dos);      
      // Emission
      if (en > eph[j]) {
        dos = GetConductionBandDensityOfStates(en - eph[j], 0);
        cfElectronsX[i].push_back(c[j] * (nocc[j] + 1) * dos);
      } else {
        cfElectronsX[i].push_back(0.);
      }
    }
    en += eStep;
  }

  for (int j = 0; j < nPhonons; ++j) {
    // Absorption
    energyLossElectronsX.push_back(-eph[j]);
    // Emission
    energyLossElectronsX.push_back(eph[j]);
    if (j <= 2) {
      scatTypeElectronsX.push_back(11);
      scatTypeElectronsX.push_back(12);
    } else {
      scatTypeElectronsX.push_back(13);
      scatTypeElectronsX.push_back(14);
    }
  }

  nLevelsX += 2 * nPhonons;

  return true;

}

bool
MediumSilicon::ElectronIonisationRates() {

  // Reference:
  // - E. Cartier, M. V. Fischetti, E. A. Eklund and F. R. McFeely,
  //   Appl. Phys. Lett 62, 3339-3341

  // Coefficients [ns-1]
  const double p[3] = {6.25e1, 3.e3, 6.8e5};
  // Threshold energies [eV]
  const double eth[3] = {1.1, 1.8, 3.45};

  double en = 0.;
  double fIon = 0.;
  for (int i = 0; i < nEnergySteps; ++i) {
    fIon = 0.;
    if (en > eth[0]) {
      fIon += p[0] * (en - eth[0]) * (en - eth[0]);
    }
    if (en > eth[1]) {
      fIon += p[1] * (en - eth[1]) * (en - eth[1]);
    }
    if (en > eth[2]) {
      fIon += p[2] * (en - eth[2]) * (en - eth[2]);
    }
    cfElectronsX[i].push_back(fIon);
    en += eStep;
  }

  energyLossElectronsX.push_back(eth[0]);
  scatTypeElectronsX.push_back(1);
  ++nLevelsX;

  return true;

}

bool
MediumSilicon::ElectronImpurityScatteringRates() {

  // Lattice temperature [eV]
  const double kbt = BoltzmannConstant * temperature;

  // Band parameters
  // Density of states effective mass
  const double md = ElectronMass * 
                    pow(mLongX * mTransX * mTransX, 1. / 3.);
  // Non-parabolicity coefficient [eV-1]
  // const double alpha = 0.5;
  const double alpha = 0.;

  // Dielectric constant
  const double eps = GetDielectricConstant();
  // Impurity concentration
  const double impurityConcentration = dopingConcentration;
  if (impurityConcentration < Small) return true;

  // Screening length
  const double ls = sqrt(eps * kbt / (4 * Pi * 
                                      FineStructureConstant * HbarC * 
                                      impurityConcentration));
  const double eb = 0.5 * HbarC * HbarC / (md * ls * ls);

  // Prefactor
  // const double c = pow(2., 2.5) * Pi * impurityConcentration * 
  //                 pow(FineStructureConstant * HbarC, 2) *
  //                 SpeedOfLight / (eps * eps * sqrt(md) * eb * eb);
  // Use momentum-transfer cross-section
  const double c = impurityConcentration * Pi * 
                   pow(FineStructureConstant * HbarC, 2) * 
                   SpeedOfLight / 
                   (sqrt(2 * md) * eps * eps);
  
  double en = 0.;
  for (int i = 0; i < nEnergySteps; ++i) {
    const double gamma = en * (1 + alpha * en);
    // cfElectrons[i][iLevel] = c * sqrt(gamma) * (1. + 2 * alpha * en) /
    //                         (1. + 4. * gamma / eb);
    if (gamma <= 0.) {
      cfElectronsX[i].push_back(0.);
      en += eStep;
      continue;
    }
    const double b = 4 * gamma / eb;
    cfElectronsX[i].push_back((c / pow(gamma, 1.5)) * 
                              (log(1. + b) - b / (1. + b)));
    en += eStep;
  }

  energyLossElectronsX.push_back(0.);
  scatTypeElectronsX.push_back(15); 
  ++nLevelsX;

  return true;

}

double
MediumSilicon::GetConductionBandDensityOfStates(const double e, 
                                                const int band) {

  if (band >= 0 && band < 6) {
    // X valleys
    // Density-of-states effective mass (cube)
    const double md3 = mLongX * mTransX * mTransX;
  
    // Non-parabolicity parameter
    double alpha = 0.;
    if (useNonParabolicity) alpha = 0.5;

    return ElectronMass * 
           sqrt(ElectronMass * md3 * e * (1. + alpha * e) / 2.) * 
           (1. + 2. * alpha * e) / (Pi2 * pow(HbarC, 3.));
  }
  
  return ElectronMass * sqrt(ElectronMass * e / 2.) / 
         (Pi2 * pow(HbarC, 3.));

}
  
}
