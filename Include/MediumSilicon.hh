// Solid crystalline silicon

#ifndef G_MEDIUM_SILICON_H
#define G_MEDIUM_SILICON_H

#include <string>
#include <vector>

#include "Medium.hh"

namespace Garfield {

class MediumSilicon : public Medium {

  public:
    // Constructor
    MediumSilicon();
    // Destructor
    ~MediumSilicon() {}

    // Doping concentration [cm-3] and type ('i', 'n', 'p')
    void SetDoping(const char type, const double c);
    void GetDoping(char& type, double& c) const;
    // Trapping cross-section
    void SetTrapCrossSection(const double ecs, const double hcs);
    void SetTrapDensity(const double n);
    void SetTrappingTime(const double etau, const double htau);
  
    // Electron transport parameters
    bool ElectronVelocity(const double ex, const double ey, const double ez,
                          const double bx, const double by, const double bz, 
                          double& vx, double& vy, double& vz);
    bool ElectronTownsend(const double ex, const double ey, const double ez,
                          const double bx, const double by, const double bz,
                          double& alpha);
    bool ElectronAttachment(const double ex, const double ey, const double ez,
                            const double bx, const double by, const double bz,
                            double& eta);
    // Hole transport parameters
    bool HoleVelocity(const double ex, const double ey, const double ez,
                      const double bx, const double by, const double bz, 
                      double& vx, double& vy, double& vz);
    bool HoleTownsend(const double ex, const double ey, const double ez,
                      const double bx, const double by, const double bz,
                      double& alpha);
    bool HoleAttachment(const double ex, const double ey, const double ez,
                        const double bx, const double by, const double bz,
                        double& eta);

    void SetLatticeMobilityModelMinimos();
    void SetLatticeMobilityModelSentaurus();
    void SetLatticeMobilityModelSah();
    
    void SetDopingMobilityModelMinimos();
    void SetDopingMobilityModelMasetti();
    
    void SetHighFieldMobilityModelMinimos();
    void SetHighFieldMobilityModelCanali();
    void SetHighFieldMobilityModelConstant();
    
    void SetImpactIonisationModelVanOverstraetenDeMan();
    void SetImpactIonisationModelGrant();

    bool GetOpticalDataRange(double& emin, double& emax, const int i = 0);
    bool GetDielectricFunction(const double e, double& eps1, double& eps2, const int i = 0);

  private:

    double bandGap;
    // Doping
    char   dopingType;
    double dopingConcentration;

    // Effective masses
    double eEffMass, hEffMass;
    // Lattice mobility
    double eLatticeMobility, hLatticeMobility;
    // Low-field mobility
    double eMobility, hMobility;
    // High-field mobility parameters
    double eBetaCanali, hBetaCanali;
    // Saturation velocity
    double eSatVel, hSatVel;
    // Hall factor
    double eHallFactor, hHallFactor;
    
    // Trapping parameters
    double eTrapCs, hTrapCs;
    double eTrapDensity, hTrapDensity;
    double eTrapTime, hTrapTime;
    int trappingModel;
    
    // Impact ionisation parameters
    double eImpactA0, eImpactA1, eImpactA2;
    double eImpactB0, eImpactB1, eImpactB2;
    double hImpactA0, hImpactA1, hImpactA2;
    double hImpactB0, hImpactB1, hImpactB2;    
    
    // Models
    int latticeMobilityModel;
    int dopingMobilityModel;
    int highFieldMobilityModel;
    int impactIonisationModel;

    void UpdateTransportParameters();
    void UpdateLatticeMobilityMinimos();
    void UpdateLatticeMobilitySentaurus();
    void UpdateLatticeMobilitySah();
    
    void UpdateDopingMobilityMinimos();
    void UpdateDopingMobilityMasetti();
    
    void UpdateSaturationVelocityMinimos();
    void UpdateSaturationVelocityCanali();
    
    void UpdateImpactIonisationVanOverstraetenDeMan();
    void UpdateImpactIonisationGrant();

    bool ElectronMobilityMinimos(const double e, double& mu) const;
    bool ElectronMobilityCanali(const double e, double& mu) const;
    bool ElectronImpactIonisationVanOverstraetenDeMan(
                                       const double e, double& alpha) const;
    bool ElectronImpactIonisationGrant(const double e, double& alpha) const;
    bool HoleMobilityMinimos(const double e, double& mu) const;
    bool HoleMobilityCanali(const double e, double& mu) const;
    bool HoleImpactIonisationVanOverstraetenDeMan(
                                   const double e, double& alpha) const;
    bool HoleImpactIonisationGrant(const double e, double& alpha) const;
        
    // Optical data
    bool hasOpticalData;
    std::string opticalDataFile;
    struct opticalData {
      // Energy [eV]
      double energy;
      // Dielectric function
      double eps1, eps2;
    };
    std::vector<opticalData> opticalDataTable;    
    bool LoadOpticalData(const std::string filename);

    void ElectronScatteringRates();
    void HoleScatteringRates();
    void ComputeHoleIntegrals(const double x, double& f3, double& f4, double& g3, double& g4);

};

}

#endif
