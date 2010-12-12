#ifndef G_COMPONENT_ANSYS123_H
#define G_COMPONENT_ANSYS123_H

#include <RQ_OBJECT.h>

#include "ComponentFieldMap.hh"

namespace Garfield {

class ComponentAnsys123: public ComponentFieldMap { 

  RQ_OBJECT("ComponentAnsys123")
    
  public:
    // Constructor
    ComponentAnsys123();
    // Destructor
    ~ComponentAnsys123() {}
    
    void ElectricField(const double x, const double y, const double z,
                       double& ex, double& ey, double& ez, 
                       Medium*& m, int& status);
    void ElectricField(const double x, const double y, const double z, 
                       double& ex, double& ey, double& ez, double& v, 
                       Medium*& m, int& status);

    void WeightingField(const double x, const double y, const double z,
                        double& wx, double& wy, double& wz,
                        const std::string label);

    double WeightingPotential(const double x, const double y, const double z,
                              const std::string label);
    
    bool GetMedium(const double x, const double y, const double z,
                   Medium*& medium);
           
    bool IsInBoundingBox(const double x, const double y, const double z) {
      return x >= xMinBoundingBox && x <= xMaxBoundingBox && 
             y >= yMinBoundingBox && y <= yMaxBoundingBox && 
             z >= zMinBoundingBox && y <= zMaxBoundingBox;
    }
        
    bool Initialise(std::string elist = "ELIST.lis", 
                    std::string nlist = "NLIST.lis",
                    std::string mplist = "MPLIST.lis", 
                    std::string prnsol = "PRNSOL.lis",
                    std::string unit = "cm");

    bool SetWeightingField(std::string prnsol, std::string label);

  protected:
  
    // Verify periodicities
    void UpdatePeriodicity() {UpdatePeriodicityCommon();}
    
  };
  
}
#endif
