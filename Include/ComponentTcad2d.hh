// Interpolation in a two-dimensional field map created by Sentaurus Device

#ifndef G_COMPONENT_TCAD_2D_H
#define G_COMPONENT_TCAD_2D_H

#include <string>
#include <vector>

#include "ComponentBase.hh"

namespace Garfield {

class ComponentTcad2d : public ComponentBase {

  public:
    // Constructor
    ComponentTcad2d();
    // Destructor
    ~ComponentTcad2d() {}
    
    void ElectricField(const double x, const double y, const double z,
                       double& ex, double& ey, double& ez, double& v,
                       Medium*& m, int& status);
    void ElectricField(const double x, const double y, const double z,
                       double& ex, double& ey, double& ez,
                       Medium*& m, int& status);    
    
    bool GetMedium(const double x, const double y, const double z,
                   Medium*& medium);

    bool GetVoltageRange(double& vmin, double& vmax);
    bool GetBoundingBox(double& xmin, double& ymin, double& zmin,
                        double& xmax, double& ymax, double& zmax); 
    void SetRangeZ(const double zmin, const double zmax);
                        
    // Import mesh and field map from files.
    bool Initialise(const std::string gridfilename, 
                    const std::string datafilename);

    // Get the number of regions in the device.
    int  GetNumberOfRegions() const {return nRegions;}
    void GetRegion(const int i, std::string& name, bool& active);
    void SetDriftRegion(const int ireg);
    void UnsetDriftRegion(const int ireg);
    // Set/get the medium for a given region.
    void SetMedium(const int ireg, Medium* m);
    bool GetMedium(const int ireg, Medium*& m) const;

    int GetNumberOfElements() const {return nElements;}
    bool GetElement(const int i, double& vol, 
                    double& dmin, double& dmax, int& type);

  private:

    // Max. number of vertices per element
    static const int nMaxVertices = 4;
  
    // Regions
    int nRegions;
    struct region {
      // Name of region (from Tcad)
      std::string name;
      // Flag indicating if the region is active (i. e. a drift medium)
      bool drift;
      Medium* medium;
    };
    std::vector<region> regions;

    // Vertices
    int nVertices;
    struct vertex {
      // Coordinates [cm]
      double x, y;
      // Potential [V] and electric field [V / cm]
      double p, ex, ey;
      // Flag indicating if vertex belongs to more than one region
      bool   isShared;
    };
    std::vector<vertex> vertices;

    // Elements
    int nElements;
    struct element {
      // Indices of vertices
      int vertex[nMaxVertices];
      // Type of element
      // 0: Point
      // 1: Segment (line)
      // 2: Triangle
      // 3: Rectangle
      // 4: Polygon
      // Types 1 - 3 are supported by this class.
      int type;
      // Associated region
      int region; 
    };
    std::vector<element> elements;

    // Voltage range
    double pMin, pMax;

    // Bounding box
    bool hasBoundingBox;
    bool hasRangeZ;
    double xMinBoundingBox, yMinBoundingBox, zMinBoundingBox;
    double xMaxBoundingBox, yMaxBoundingBox, zMaxBoundingBox;

    // Element from the previous call
    int lastElement;
    // Shape functions for interpolation 
    // (local coordinates)
    double w[nMaxVertices];
    
    // Reset the component
    void Reset();
    // Periodicities
    void UpdatePeriodicity();    

    bool CheckRectangle(const double x, const double y, const int i);
    bool CheckTriangle(const double x, const double y, const int i);
    bool CheckLine(const double x, const double y, const int i);

    bool LoadGrid(const std::string gridfilename);
    bool LoadData(const std::string datafilename);
    void Cleanup();

};

}
#endif
