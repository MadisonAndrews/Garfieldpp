#include <iostream>

#include "heed++/code/HeedPhoton.h"

#include "TrackHeed.hh"
#include "TrackHeedGlobals.hh"
#include "FundamentalConstants.hh"
#include "Random.hh"

namespace Garfield {

TrackHeed::TrackHeed() : 
  ready(false), hasActiveTrack(false),
  mediumDensity(-1.), mediumName(""),
  databasePath(""), isPathSet(false),
  particle(), node(0),
  matter(0), gas(0), material(0),
  atPacs(0), molPacs(0),
  energyMesh(0), transferCs(0),
  elScat(0), lowSigma(0), pairProd(0), deltaCs(0) {
  
  Garfield::HeedInterface::sensor = 0;
  Garfield::HeedInterface::useEfield = false;
  Garfield::HeedInterface::useBfield = false;
  
}

TrackHeed::~TrackHeed() {

  if (node       != 0) delete node;
  if (matter     != 0) delete matter;
  if (gas        != 0) delete gas;
  if (material   != 0) delete material;
  if (atPacs     != 0) delete atPacs;
  if (molPacs    != 0) delete molPacs;
  if (energyMesh != 0) delete energyMesh;
  if (transferCs != 0) delete transferCs;
  if (elScat     != 0) delete elScat;
  if (lowSigma   != 0) delete lowSigma;
  if (pairProd   != 0) delete pairProd;
  if (deltaCs    != 0) delete deltaCs;
  
  Garfield::HeedInterface::sensor = 0;

}

void
TrackHeed::NewTrack(
            const double x0, const double y0, const double z0, const double t0,
            const double dx0, const double dy0, const double dz0) {

  node = 0;
  hasActiveTrack = false;
  
  // Make sure the sensor has been set.
  if (sensor == 0) {
    std::cerr << "TrackHeed::NewTrack:\n";
    std::cerr << "    Sensor is not defined.\n";
    ready = false;
    return;
  }
  
  HeedInterface::sensor = sensor;
  
  // Make sure the initial position is inside an ionisable medium.
  Medium* medium;
  if (!sensor->GetMedium(x0, y0, z0, medium)) {
    std::cerr << "TrackHeed::NewTrack:\n";
    std::cerr << "    No medium at initial position.\n";
    ready = false;
    return;
  } else if (!medium->IsIonisable()) {
    std::cerr << "TrackHeed:NewTrack:\n";
    std::cerr << "    Medium at initial position is not ionisable.\n";
    ready = false;
    return;
  }

  // Check if the medium has changed since the last call.
  if (medium->GetName()        != mediumName || 
      medium->GetMassDensity() != mediumDensity) {
    isChanged = true;
  }
  
  // If medium or particle have changed, update the cross-sections.
  if (isChanged) {
    ready = false;
    if (!Setup(medium)) return;
    ready = true;
    isChanged = false;
    mediumName    = medium->GetName();
    mediumDensity = medium->GetMassDensity();

  }
  
  particle_bank.clear();
  cluster_bank.allocate_block(100);
  chamber.conduction_electron_bank.allocate_block(1000);
  
  // Check the direction vector.
  double dx = dx0, dy = dy0, dz = dz0;
  const double d = sqrt(dx * dx + dy * dy + dz * dz);
  if (d <= 0.) {
    // Null vector. Sample the direction isotropically.
    const double phi = TwoPi * RndmUniform();
    const double ctheta = 1. - 2. * RndmUniform();
    const double stheta = sqrt(1. - ctheta * ctheta);
    dx = cos(phi) * stheta;
    dy = sin(phi) * stheta;
    dz = ctheta;
  } else {
    // Normalise the direction vector.
    dx /= d; dy /= d; dz /= d;
  }
  vec velocity(dx, dy, dz);
  velocity = velocity * mparticle::speed_of_light * GetBeta();
  
  if (debug) {
    std::cout << "TrackHeed::NewTrack:\n";
    std::cout << "    Track starts at (" 
              << x0 << ", " << y0 << ", " << z0 << ") at time " 
              << t0 << "\n";
    std::cout << "    Initial direction: ("
              << dx << ", " << dy << ", " << dz << ")\n";
  }

  // Initial position (convert from cm to mm).
  point p0(x0 * 10., y0 * 10., z0 * 10.);
  // Setup the particle.
  last_particle_number = 0;
  particle = HeedParticle(&chamber, 
                          p0, velocity, t0, 
                          &proton_def);
  // Transport the particle.
  particle.fly();
  node = particle_bank.get_first_node();
  hasActiveTrack = true;

}

bool
TrackHeed::GetCluster(double& xcls, double& ycls, double& zcls, double& tcls,
                      int& n, double& e, double& extra) {

  // Initial settings.
  xcls = ycls = zcls = tcls = 0.;
  extra = 0.;
  n = 0;
  e = 0.;
  
  // Make sure NewTrack has successfully been called. 
  if (!ready) {
    std::cerr << "TrackHeed::GetCluster:\n";
    std::cerr << "    Track has not been initialized.\n";
    std::cerr << "    Call NewTrack first.\n";
    return false;
  }
  
  if (!hasActiveTrack) {
    std::cerr << "TrackHeed::GetCluster:\n";
    std::cerr << "    There are no more clusters.\n";
    return false;
  }
  
  // Make sure the particle bank is not empty.
  if (node == 0) {
    hasActiveTrack = false;
    return false;
  }
  
  // Convert the particle to a (virtual) photon.
  HeedPhoton* virtualPhoton = dynamic_cast<HeedPhoton*>(node->el.get());
  if (virtualPhoton == 0) {
    std::cerr << "TrackHeed::GetCluster:\n";
    std::cerr << "    Particle is not a virtual photon.\n";
    std::cerr << "    Program bug!\n";
    return false;
  }
  if (virtualPhoton->parent_particle_number != 0) {
    std::cerr << "TrackHeed::GetCluster:\n";
    std::cerr << "    Virtual photon has an unexpected parent particle.\n";
    return false;
  }
  virtualPhoton->fly();
  // Get the location of the interaction (convert from mm to cm).
  xcls = virtualPhoton->currpos.pt.v.x * 0.1;
  ycls = virtualPhoton->currpos.pt.v.y * 0.1;
  zcls = virtualPhoton->currpos.pt.v.z * 0.1;
  // Get the transferred energy (convert from MeV to eV).
  e = virtualPhoton->energy * 1.e6;
  
  // Make a list of parent particle id numbers. 
  std::vector<int> ids;
  ids.clear();
  // At the beginning, ther is only the virtual photon.
  ids.push_back(virtualPhoton->particle_number);
  int nIds = 1;
  
  // Look for daughter particles.
  chamber.conduction_electron_bank.allocate_block(1000);
  bool deleteNode = false;
  HeedDeltaElectron* delta = 0;
  HeedPhoton* photon = 0;
  AbsListNode<ActivePtr<gparticle> >* nextNode = node->get_next_node();
  AbsListNode<ActivePtr<gparticle> >* tempNode = 0;
  // Loop over the particle bank.
  while (nextNode != 0) {
    deleteNode = false;
    // Check if it is a delta electron.
    delta = dynamic_cast<HeedDeltaElectron*>(nextNode->el.get());
    if (delta != 0) {
      // Check if the delta electron was produced by one of the photons
      // belonging to this cluster.
      for (int i = nIds; i--;) {
        if (delta->parent_particle_number == ids[i]) {
          // Transport the delta electron.
          delta->fly();
          deleteNode = true;
          break;
        }
      }
    } else {
      // Check if it is a real photon.
      photon = dynamic_cast<HeedPhoton*>(nextNode->el.get());
      if (photon == 0) {
        std::cerr << "TrackHeed::GetCluster:\n";
        std::cerr << "    Particle is neither an electron nor a photon.\n";
        return false;
      }
      for (int i = nIds; i--;) {
        if (photon->parent_particle_number == ids[i]) {
          // Transport the photon and add its number to the list of ids.
          photon->fly();
          deleteNode = true;
          ids.push_back(photon->particle_number);
          ++nIds;
          break;
        }
      }
    }
    // Proceed with the next node in the particle bank.
    if (deleteNode) {
      tempNode = nextNode->get_next_node();
      particle_bank.erase(nextNode);
      nextNode = tempNode;
    } else {
      nextNode = nextNode->get_next_node();
    }
  }
  
  // Get the total number of conduction electrons produced in this step.
  n = chamber.conduction_electron_bank.get_qel();
  
  // Next virtual photon.
  nextNode = node->get_next_node();
  particle_bank.erase(node);
  node = nextNode;

  return true;

}

bool
TrackHeed::GetElectron(const int i, double& x, double& y, double& z) {

  // Make sure NewTrack has successfully been called.
  if (!ready) {
    std::cerr << "TrackHeed::GetElectron:\n";
    std::cerr << "    Track has not been initialized.\n";
    std::cerr << "    Call NewTrack first.\n";
    return false;
  }

  // Make sure an electron with this number exists.
  const int n = chamber.conduction_electron_bank.get_qel();
  if (i < 0 || i >= n) {
    std::cerr << "TrackHeed::GetElectron:\n";
    std::cerr << "    Electron number out of range.\n";
    return false;
  }
  
  x = chamber.conduction_electron_bank[i].ptloc.v.x * 0.1;
  y = chamber.conduction_electron_bank[i].ptloc.v.y * 0.1;
  z = chamber.conduction_electron_bank[i].ptloc.v.z * 0.1;

  return true;

}

void
TrackHeed::TransportDeltaElectron(
      const double x0, const double y0, const double z0, 
      const double t0, const double e0, 
      const double dx0, const double dy0, const double dz0,
      int& nel) {

  // Make sure the kinetic energy is positive.
  if (e0 <= 0.) {
    std::cerr << "TrackHeed::TransportDeltaElectron:\n";
    std::cerr << "    Kinetic energy must be positive.\n";
    return;
  }
  
  // Make sure the sensor has been set.
  if (sensor == 0) {
    std::cerr << "TrackHeed::TransportDeltaElectron:\n";
    std::cerr << "    Sensor is not defined.\n";
    ready = false;
    return;
  }
  
  HeedInterface::sensor = sensor;
  
  // Make sure the initial position is inside an ionisable medium.
  Medium* medium;
  if (!sensor->GetMedium(x0, y0, z0, medium)) {
    std::cerr << "TrackHeed::TransportDeltaElectron:\n";
    std::cerr << "    No medium at initial position.\n";
    return;
  } else if (!medium->IsIonisable()) {
    std::cerr << "TrackHeed:TransportDeltaElectron:\n";
    std::cerr << "    Medium at initial position is not ionisable.\n";
    ready = false;
    return;
  }

  // Check if the medium has changed since the last call.
  if (medium->GetName()        != mediumName || 
      medium->GetMassDensity() != mediumDensity) {
    isChanged = true;
    ready = false;
    hasActiveTrack = false;
    if (!Setup(medium)) return;
    ready = true;
    mediumName    = medium->GetName();
    mediumDensity = medium->GetMassDensity();
  }
    
  chamber.conduction_electron_bank.allocate_block(1000);
  
  // Check the direction vector.
  double dx = dx0, dy = dy0, dz = dz0;
  const double d = sqrt(dx * dx + dy * dy + dz * dz);
  if (d <= 0.) {
    // Null vector. Sample the direction isotropically.
    const double phi = TwoPi * RndmUniform();
    const double ctheta = 1. - 2. * RndmUniform();
    const double stheta = sqrt(1. - ctheta * ctheta);
    dx = cos(phi) * stheta;
    dy = sin(phi) * stheta;
    dz = ctheta;
  } else {
    // Normalise the direction vector.
    dx /= d; dy /= d; dz /= d;
  }
  vec velocity(dx, dy, dz);
  
  // Calculate the speed for the given kinetic energy.
  const double gamma = 1. + e0 / ElectronMass;
  const double beta = sqrt(1. - 1. / (gamma * gamma)); 
  double speed = mparticle::speed_of_light * beta;
  velocity = velocity * speed;
  
  // Initial position (convert from cm to mm).
  point p0(x0 * 10., y0 * 10., z0 * 10.);
 
  // Transport the electron.
  HeedDeltaElectron delta(&chamber, p0, velocity, t0, 0);
  delta.fly();
  
  nel = chamber.conduction_electron_bank.get_qel();
  
}
      
void
TrackHeed::SetDatabasePath(const std::string dbpath) {

  if (dbpath == "") {
    std::cerr << "TrackHeed::SetDatabasePath:\n";
    std::cerr << "    String is empty.\n";
    return;
  }
  
  databasePath = dbpath;
  // Append '/' if necessary.
  if (dbpath[dbpath.length() - 1] != '/') {
    databasePath.append("/");
  }
  std::cout << "TrackHeed::SetDatabasePath:\n";
  std::cout << "    Database path set to " << databasePath << ".\n";
  isPathSet = true;
     
}

void
TrackHeed::EnableElectricField() {

  HeedInterface::useEfield = true;
  
}

void
TrackHeed::DisableElectricField() {

  HeedInterface::useEfield = false;
  
}

void
TrackHeed::EnableMagneticField() {

  HeedInterface::useBfield = true;
  
}

void
TrackHeed::DisableMagneticField() {

  HeedInterface::useBfield = false;

}

bool
TrackHeed::Setup(Medium* medium) {

  // Make sure the path to the Heed database is known.
  if (!isPathSet) {
    char* dbPath = getenv("HEED_DATABASE");
    if (dbPath == 0) {
      std::cerr << "TrackHeed::Setup:\n";
      std::cerr << "    Database path is not defined.\n";
      std::cerr << "    Environment variable HEED_DATABASE is not set.\n";
      std::cerr << "    Cannot proceed with initialization.\n";
      return false;
    }
    databasePath = std::string(dbPath) + "/"; 
  }
  
  // Check once more that the medium exists.
  if (medium == 0) {
    std::cerr << "TrackHeed::Setup:\n";
    std::cerr << "    Medium pointer is null.\n";
    return false;
  }
    
  // Setup the energy mesh.
  if (energyMesh != 0) {
    delete energyMesh; energyMesh = 0;
  }
  energyMesh = new EnergyMesh(2.0e-6, 2.0e-1, 200);
  
  if (medium->IsGas()) {
    if (!SetupGas(medium)) return false;
  } else {
    if (!SetupMaterial(medium)) return false;
  }
    
  // Energy transfer cross-section
  // Set a flag indicating whether the primary particle is an electron.
  int sel = 0;
  if (isElectron) sel = 1;
  const double gamma = GetGamma();

  if (transferCs != 0) {
    delete transferCs;
    transferCs = 0;
  }
  transferCs = new EnTransfCS(mass / 1.e6, gamma - 1, sel, matter, long(q));
  
  if (!SetupDelta()) return false;  

  if (debug) {
    const double nc = transferCs->quanC;
    const double dedx = transferCs->meanC1 * 1.e3;
    const double w = matter->W * 1.e6;
    const double f = matter->F;
    std::cout << "TrackHeed::Setup:\n";
    std::cout << "    Cluster density: " << nc << " cm-1\n";
    std::cout << "    Stopping power:  " << dedx << " keV/cm\n";
    std::cout << "    W value:         " << w << " eV\n";
    std::cout << "    Fano factor:     " << f << "\n";
  }
  
  fixsyscoor primSys(point(0., 0., 0.), basis("primary"), "primary");
  chamber = Chamber(primSys, transferCs, deltaCs);

  return true;
  
}

bool
TrackHeed::SetupGas(Medium* medium) {

  // Get temperature and pressure.
  double pressure = medium->GetPressure();
  pressure = (pressure / AtmosphericPressure) * atmosphere;
  double temperature = medium->GetTemperature();
  
  const int nComponents = medium->GetNumberOfComponents();
  if (nComponents < 1) {
    std::cerr << "TrackHeed::SetupGas:\n";
    std::cerr << "    Gas " << medium->GetName() 
              << " has zero constituents.\n";
    return false;
  }

  if (molPacs != 0) {
    delete molPacs;
    molPacs = 0;
  }
  molPacs = new MolecPhotoAbsCS*[nComponents];
  DynLinArr<std::string> notations; notations.clear();
  DynLinArr<double> fractions; fractions.clear();
  
  for (int i = 0; i < nComponents; ++i) {
    std::string gasname;
    double frac;
    medium->GetComponent(i, gasname, frac);
    // If necessary, change the Magboltz name to the Heed internal name.
    if (gasname == "He-3") gasname = "He";
    if (gasname == "CD4") gasname = "CH4";
    if (gasname == "iC4H10" || gasname == "nC4H10") gasname = "C4H10";
    if (gasname == "neoC5H12" || gasname == "nC5H12") gasname = "C5H12";
    if (gasname == "H2O") gasname = "Water";
    if (gasname == "D2") gasname = "H2";
    if (gasname == "cC3H6") gasname = "C3H6";
    // Find the corresponding photoabsorption cross-section.
    if (gasname == "CF4")         molPacs[i] = &CF4_MPACS;
    else if (gasname == "Ar")     molPacs[i] = &Ar_MPACS;
    else if (gasname == "He")     molPacs[i] = &He_MPACS;
    else if (gasname == "Ne")     molPacs[i] = &Ne_MPACS;
    else if (gasname == "Kr")     molPacs[i] = &Kr_MPACS;
    else if (gasname == "Xe")     molPacs[i] = &Xe_MPACS;
    else if (gasname == "CH4")    molPacs[i] = &CH4_MPACS;
    else if (gasname == "C2H6")   molPacs[i] = &C2H6_MPACS;
    else if (gasname == "C3H8")   molPacs[i] = &C3H8_MPACS;
    else if (gasname == "C4H10")  molPacs[i] = &C4H10_MPACS;
    else if (gasname == "CO2")    molPacs[i] = &CO2_MPACS;
    else if (gasname == "C5H12")  molPacs[i] = &C5H12_MPACS;
    else if (gasname == "Water")  molPacs[i] = &H2O_MPACS;
    else if (gasname == "O2")     molPacs[i] = &O2_MPACS;
    else if (gasname == "N2")     molPacs[i] = &N2_MPACS;
    else if (gasname == "NO")     molPacs[i] = &NO_MPACS;
    else if (gasname == "N2O")    molPacs[i] = &N2O_MPACS;
    else if (gasname == "C2H4")   molPacs[i] = &C2H4_MPACS;
    else if (gasname == "C2H2")   molPacs[i] = &C2H2_MPACS;
    else if (gasname == "H2")     molPacs[i] = &H2_MPACS;
    else if (gasname == "CO")     molPacs[i] = &CO_MPACS;
    else if (gasname == "Methylal") molPacs[i] = &Methylal_MPACS;
    else if (gasname == "DME")    molPacs[i] = &DME_MPACS;
    else if (gasname == "C2F6")   molPacs[i] = &C2F6_MPACS;
    else if (gasname == "SF6")    molPacs[i] = &SF6_MPACS;
    else if (gasname == "NH3")    molPacs[i] = &NH3_MPACS;
    else if (gasname == "C3H6")   molPacs[i] = &C3H6_MPACS;
    else if (gasname == "CH3OH")  molPacs[i] = &CH3OH_MPACS;
    else if (gasname == "C2H5OH") molPacs[i] = &C2H5OH_MPACS;
    else if (gasname == "C3H7OH") molPacs[i] = &C3H7OH_MPACS;
    else if (gasname == "Cs")     molPacs[i] = &Cs_MPACS;
    else if (gasname == "F2")     molPacs[i] = &F2_MPACS;
    else if (gasname == "CS2")    molPacs[i] = &CS2_MPACS;
    else if (gasname == "COS")    molPacs[i] = &COS_MPACS;
    else if (gasname == "CD4")    molPacs[i] = &CH4_MPACS;
    else if (gasname == "BF3")    molPacs[i] = &BF3_MPACS;
    else if (gasname == "C2HF5")  molPacs[i] = &C2HF5_MPACS;
    else if (gasname == "CHF3")   molPacs[i] = &CHF3_MPACS;
    else if (gasname == "CF3Br")  molPacs[i] = &CF3Br_MPACS;
    else if (gasname == "C3F8")   molPacs[i] = &C3F8_MPACS;
    else if (gasname == "O3")     molPacs[i] = &O3_MPACS;
    else if (gasname == "Hg")     molPacs[i] = &Hg_MPACS;
    else if (gasname == "H2S")    molPacs[i] = &H2S_MPACS;
    else if (gasname == "GeH4")   molPacs[i] = &GeH4_MPACS;
    else if (gasname == "SiH4")   molPacs[i] = &SiH4_MPACS;
    else {
      std::cerr << "TrackHeed::SetupGas:\n";
      std::cerr << "    Photoabsorption cross-section data for " << gasname 
                << " are not available.\n";
      return false;
    }
    notations.increment(gasname);
    fractions.increment(frac);
  }
  std::string gasname = medium->GetName();
  if (gas != 0) {
    delete gas; gas = 0;
  }

  gas = new GasDef(gasname, gasname, nComponents,
                   notations, fractions, pressure, temperature);
    
  double w = medium->GetW() * 1.e-6;
  if (w < 0.) w = 0.;
  double f = medium->GetFanoFactor();
  if (f <= 0.) f = standard_factor_Fano;

  if (matter != 0) {
    delete matter;
    matter = 0;
  }
  matter = new HeedMatterDef(energyMesh, gas, molPacs, w, f);
  
  return true;

}

bool
TrackHeed::SetupMaterial(Medium* medium) {

  // Get temperature and density.
  double temperature = medium->GetTemperature();
  double density = medium->GetNumberDensity() / cm3;
  
  const int nComponents = medium->GetNumberOfComponents();
  if (atPacs != 0) {
    delete atPacs;
    atPacs = 0;
  }
  atPacs = new AtomPhotoAbsCS*[nComponents];
  
  DynLinArr<std::string> notations; notations.clear();
  DynLinArr<double> fractions; fractions.clear();
  for (int i = 0; i < nComponents; ++i) {
    std::string materialName;
    double frac;
    medium->GetComponent(i, materialName, frac);
    if (materialName == "C")       atPacs[i] = &Carbon_PACS;
    else if (materialName == "Si") atPacs[i] = &Silicon_crystal_PACS;
    else if (materialName == "Ge") atPacs[i] = &Germanium_crystal_PACS;
    else {
      std::cerr << "TrackHeed::SetupMaterial:\n";
      std::cerr << "    Photoabsorption cross-section data for " << materialName 
                << " are not implemented.\n";
      return false;
    }
    notations.increment(materialName);
    fractions.increment(frac);
  }
  if (material != 0) {
    delete material; material = 0;
  }
  std::string materialName = medium->GetName();
  material = new MatterDef(materialName, materialName, nComponents,
                           notations, fractions, density, temperature);


  double w = medium->GetW() * 1.e-6;
  if (w < 0.) w = 0.;
  double f = medium->GetFanoFactor();
  if (f <= 0.) f = standard_factor_Fano;
    
  if (matter != 0) {
    delete matter;
    matter = 0;
  }
  matter = new HeedMatterDef(energyMesh, material, atPacs, w, f);
  
  return true;
 
}

bool
TrackHeed::SetupDelta() {

  // Load elastic scattering data.
  std::string filename = databasePath + "cbdel.dat";
  if (elScat != 0) {
    delete elScat; elScat = 0;
  }
  elScat = new ElElasticScat(filename);
  
  filename = databasePath + "elastic_disp.dat";
  if (lowSigma != 0) {
    delete lowSigma; lowSigma = 0;
  }
  lowSigma = new ElElasticScatLowSigma(elScat, filename);
  
  // Load data for calculation of ionization.
  // Get W value and Fano factor.
  const double w = matter->W * 1.e6;
  const double f = matter->F;
  filename = databasePath + "delta_path.dat";
  if (pairProd != 0) {
    delete pairProd; pairProd = 0;
  }
  pairProd = new PairProd(filename, w, f);
  
  if (deltaCs != 0) {
    delete deltaCs; deltaCs = 0;
  }
  deltaCs = new HeedDeltaElectronCS(matter, elScat, lowSigma, pairProd);

}

TrackHeed::Chamber::Chamber(const abssyscoor& fcsys,
                            const EnTransfCSType etcst,
                            const HeedDeltaElectronCSType hdecst) :
    sh_manip_absvol(fcsys),
    box(1000. * mm, 1000. * mm, 15. * mm, "chamber"),
    EnTransfCSType(etcst), HeedDeltaElectronCSType(hdecst) {
    
}

}