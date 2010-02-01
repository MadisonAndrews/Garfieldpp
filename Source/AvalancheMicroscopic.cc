#include <iostream>
#include <cmath>
#include <string>

#include "AvalancheMicroscopic.hh"
#include "FundamentalConstants.hh"
#include "Random.hh"

namespace Garfield {

// Numerical prefactors
double AvalancheMicroscopic::c1 = SpeedOfLight * sqrt(2. / ElectronMass);
double AvalancheMicroscopic::c2 = c1 * c1 / 4.;

AvalancheMicroscopic::AvalancheMicroscopic() :
  sensor(0), 
  nElectrons(0), nIons(0),  
  histogram(0), hasEnergyHistogram(false),
  useSignal(false), useDriftLines(false),
  deltaCut(0.),
  nCollSkip(100),
  hasUserHandleAttachment(false), 
  hasUserHandleInelastic(false),
  hasUserHandleIonisation(false),
  userHandleAttachment(0), userHandleInelastic(0),
  userHandleIonisation(0), 
  debug(false), warning(false) {

}

void 
AvalancheMicroscopic::SetSensor(Sensor* s) {

  if (s == 0) {
    std::cerr << "AvalancheMicroscopic::SetSensor:" << std::endl;
    std::cerr << "    Sensor is not defined." << std::endl;
    return;
  }
  sensor = s;
}

void 
AvalancheMicroscopic::EnableEnergyHistogramming(TH1F* histo) {

  if (histo == 0) {
    std::cerr << "AvalancheMicroscopic::EnableEnergyHistogramming:" 
              << std::endl;
    std::cerr << "    Histogram is not defined." << std::endl;
    return;
  }
  
  histogram = histo;
  hasEnergyHistogram = true;
  
}

void 
AvalancheMicroscopic::DisableEnergyHistogramming() {

  hasEnergyHistogram = false;
  
}

void 
AvalancheMicroscopic::SetCollisionSteps(const int n) {

  if (n <= 0) {
    std::cerr << "AvalancheMicroscopic::SetCollisionSteps:" << std::endl;
    std::cerr << "    Number of collisions to be skipped set to" 
              << " default value (100)." << std::endl;
    nCollSkip = 100;
    return;
  }
  
  nCollSkip = n;
  
}

void 
AvalancheMicroscopic::GetEndpoint(const int i, 
  double& x0, double& y0, double& z0, double& t0, double& e0,
  double& x1, double& y1, double& z1, double& t1, double& e1,
  int& status) const {
  
  if (i < 0 || i >= nEndpoints) {
    std::cerr << "AvalancheMicroscopic::GetEndpoint:" << std::endl;
    std::cerr << "    Endpoint " << i << " does not exist." << std::endl;
    return;
  }

  x0 = endpoints[i].x0; y0 = endpoints[i].y0; z0 = endpoints[i].z0;
  t0 = endpoints[i].t0; e0 = endpoints[i].e0;
  x1 = endpoints[i].x;  y1 = endpoints[i].y;  z1 = endpoints[i].z;
  t1 = endpoints[i].t;  e1 = endpoints[i].energy;  
  status = -1;  
  if (e1 < 0) {
    // Electron stopped because of attachment
    e1 = -e1;
    status = -7;
  }

}

int 
AvalancheMicroscopic::GetNumberOfDriftLinePoints(const int i) const {

  if (i < 0 || i >= nEndpoints) {
    std::cerr << "AvalancheMicroscopic::GetNumberOfDriftLinePoints:" << std::endl;
    std::cerr << "    Endpoint " << i << " does not exist." << std::endl;
    return 0;
  }
  
  if (!useDriftLines) return 2;

  return endpoints[i].driftLine.size() + 2;

}

void 
AvalancheMicroscopic::GetDriftLinePoint(
  double& x, double& y, double& z, double &t,
  const int ip, const int iel) const {
  
  if (iel < 0 || iel >= nEndpoints) {
    std::cerr << "AvalancheMicroscopic::GetDriftLinePoint:" << std::endl;
    std::cerr << "    Endpoint " << iel << " does not exist." << std::endl;
    return;
  }

  if (ip <= 0) {
    x = endpoints[iel].x0;
    y = endpoints[iel].y0;
    z = endpoints[iel].z0;
    t = endpoints[iel].t0;
    return;
  }

  if (ip > endpoints[iel].driftLine.size()) {
    x = endpoints[iel].x;
    y = endpoints[iel].y;
    z = endpoints[iel].z;
    t = endpoints[iel].t;
    return;
  }

  x = endpoints[iel].driftLine[ip - 1].x;
  y = endpoints[iel].driftLine[ip - 1].y;
  z = endpoints[iel].driftLine[ip - 1].z;
  t = endpoints[iel].driftLine[ip - 1].t;

}

void 
AvalancheMicroscopic::SetUserHandleAttachment(
    void (*f)(double x, double y, double z, double t, 
              int type, int level, Medium* m)) {
         
  userHandleAttachment = f;
  hasUserHandleAttachment = true;
  
}

void 
AvalancheMicroscopic::UnsetUserHandleAttachment() {
  
  userHandleAttachment = 0;
  hasUserHandleAttachment = false;
  
}

void 
AvalancheMicroscopic::SetUserHandleInelastic(
    void (*f)(double x, double y, double z, double t, 
              int type, int level, Medium* m)) {
         
  userHandleInelastic = f;
  hasUserHandleInelastic = true;
  
}

void 
AvalancheMicroscopic::UnsetUserHandleInelastic() {
  
  userHandleInelastic = 0;
  hasUserHandleInelastic = false;
  
}

void 
AvalancheMicroscopic::SetUserHandleIonisation(
    void (*f)(double x, double y, double z, double t, 
              int type, int level, Medium* m)) {
         
  userHandleIonisation = f;
  hasUserHandleIonisation = true;
  
}

void 
AvalancheMicroscopic::UnsetUserHandleIonisation() {
  
  userHandleIonisation = 0;
  hasUserHandleIonisation = false;
  
}

bool 
AvalancheMicroscopic::AvalancheElectron(
    const double x0, const double y0, const double z0, const double t0, 
    const double e0, const double dx0, const double dy0, const double dz0) {
  
  // Make sure that the sensor is defined
  if (sensor == 0) {
    std::cerr << "AvalancheMicroscopic::AvalancheElectron:" << std::endl;
    std::cerr << "    Sensor is not defined." << std::endl;
    return false;
  }

  // Make sure that the starting point is inside a medium
  Medium* medium;
  if (!sensor->GetMedium(x0, y0, z0, medium)) {
    std::cerr << "AvalancheMicroscopic::AvalancheElectron:" << std::endl;
    std::cerr << "    No medium at initial position." << std::endl;
    return false;
  }
  
  // Make sure that the medium is "driftable" and microscopic
  if (!medium->IsDriftable() || !medium->IsMicroscopic()) {
    std::cerr << "AvalancheMicroscopic::AvalancheElectron:" << std::endl;
    std::cerr << "    Medium at initial position does not provide " 
              << " microscopic tracking data." << std::endl;
    return false;
  }
  
  if (debug) {
    std::cout << "AvalancheMicroscopic::AvalancheElectron:" << std::endl;
    std::cout << "    Starting to drift in medium " 
              << medium->GetName() << "." << std::endl;
  }
  
  // Get the id number of the drift medium
  int id = medium->GetId();    
   
  // Clear the stack
  stack.clear();
  endpoints.clear();
  
  // Reset the particle counters
  nElectrons = 1;
  nIons = 0;
  nEndpoints = 0;
  
  // Check if the initial energy is higher than the transport cut
  if (e0 <= deltaCut) {
    std::cerr << "AvalancheMicroscopic::AvalancheElectron:" << std::endl;
    std::cerr << "    Initial electron energy (" << e0 << " eV)"
              << " is lower than the transport threshold" 
              << " (" << deltaCut << " eV)." << std::endl;
    return false;
  }
    
  // Null-collision rate
  double fLim = medium->GetNullCollisionRate();
  if (fLim <= 0.) {
    std::cerr << "AvalancheMicroscopic::AvalancheElectron:" << std::endl;
    std::cerr << "    Got null-collision rate <= 0." << std::endl;
    return false;
  }
  // Real collision rate
  double fReal;
   
  // Count collisions between updates
  int nCollTemp = 0;  
  
  // Electric field
  double ex, ey, ez;
  int status;
         
  // Current position, direction and energy
  double x, y, z, t;
  double dx, dy, dz;
  double energy;
  // Timestep (squared)
  double dt, dt2;
  // Direction and energy after a step
  double newDx, newDy, newDz, newEnergy;
  
  // Collision type (elastic, ionisation, attachment, inelastic)
  int cstype;
  // Cross-section term
  int level;
  // Energy loss in inelastic collisions
  double eloss;
  // Scattering parameters
  double s1, s2;
  // Scattering angles
  double ctheta0, stheta0;
  double phi0, cphi0, sphi0;
  double phi, theta;
  double ctheta, stheta;
  double arg, q, d;
  // Secondary electron energy
  double esec;
  
  // Random number
  double r;
  // Prefactors
  double a, b;
    
  // Add the first electron
  electron newElectron;
  newElectron.x0 = x0;  newElectron.y0 = y0;  newElectron.z0 = z0;  
  newElectron.x  = x0;  newElectron.y  = y0;  newElectron.z = z0;  
  newElectron.t0 = t0;  newElectron.t  = t0;
  newElectron.dx = dx0; newElectron.dy = dy0; newElectron.dz = dz0;
  newElectron.e0 = Max(e0, Small); newElectron.energy = newElectron.e0;
  newElectron.driftLine.clear();
  stack.push_back(newElectron);

  // Check the given initial direction
  d = sqrt(dx0 * dx0 + dy0 * dy0 + dz0 * dz0);
  if (fabs(d) < Small) {
    // Direction has zero norm, draw a random direction
    phi = TwoPi * RndmUniform();
    ctheta = 1. - 2. * RndmUniform();
    stheta = sin(acos(ctheta));
    stack[0].dx = cos(phi) * stheta;
    stack[0].dy = sin(phi) * stheta;
    stack[0].dz = ctheta;
  } else if (fabs(d - 1.) > Small) {
    // Renormalise direction to 1
    stack[0].dx /= d; stack[0].dy /= d; stack[0].dz /= d;
  }

  // Status flag
  bool ok = true;
  // Stack size
  int nSize = 1;
  // Index of the electron in the stack
  int iEl;
  while (1) {
    nSize = stack.size();
    if (nSize <= 0) break;
    // Loop over all electrons in the avalanche
    for (iEl = nSize; iEl--;) {
      // Get the electron from the stack
      x = stack[iEl].x; y = stack[iEl].y; z = stack[iEl].z;
      energy = stack[iEl].energy; t = stack[iEl].t;      
      dx = stack[iEl].dx; dy = stack[iEl].dy; dz = stack[iEl].dz;

      ok = true;
      nCollTemp = 0;

      while (1) {

        if (energy < deltaCut) {
          stack[iEl].x = x; stack[iEl].y = y; stack[iEl].z = z;
          stack[iEl].t = t; stack[iEl].energy = energy;
          stack[iEl].dx = dx; stack[iEl].dy = dy; stack[iEl].dz = dz;
          endpoints.push_back(stack[iEl]);
          stack.erase(stack.begin() + iEl);
          ok = false;
          break;
        }

        if (hasEnergyHistogram) {
          histogram->Fill(energy);
        }

        // Get the local electric field and medium
        sensor->ElectricField(x, y, z, ex, ey, ez, medium, status);
        // Sign change
        ex = -ex; ey = -ey; ez = -ez;
        
        if (status != 0) {
          // Electron has left all drift media
          stack[iEl].x = x; stack[iEl].y = y; stack[iEl].z = z;
          stack[iEl].t = t; stack[iEl].energy = energy;
          stack[iEl].dx = dx; stack[iEl].dy = dy; stack[iEl].dz = dz;          
          endpoints.push_back(stack[iEl]);
          stack.erase(stack.begin() + iEl);
          ok = false;
          break;
        }
        
        if (medium->GetId() != id) {
          // Medium has changed
          if (!medium->IsMicroscopic()) {
            // Electron has left the microscopic drift medium
            stack[iEl].x = x; stack[iEl].y = y; stack[iEl].z = z;
            stack[iEl].t = t; stack[iEl].energy = energy;
            stack[iEl].dx = dx; stack[iEl].dy = dy; stack[iEl].dz = dz;
            endpoints.push_back(stack[iEl]);
            stack.erase(stack.begin() + iEl);
            ok = false;
            break;
          }
          id = medium->GetId();
          // Update the null-collision rate
          fLim = medium->GetNullCollisionRate();
          if (fLim <= 0.) {
            std::cerr << "AvalancheMicroscopic::AvalancheElectron:" << std::endl;
            std::cerr << "    Got null-collision rate <= 0." << std::endl;
            return false;
          }          
        }

        a = c1 * (dx * ex + dy * ey + dz * ez) * sqrt(energy);
        b = c2 * (ex * ex + ey * ey + ez * ez);

        // Determine the timestep
        dt = 0.;
        while (1) {
          // Determine the flight time
          r = RndmUniformPos();
          dt += - log(r) / fLim;
          // Update the energy
          newEnergy = Max(energy + (a + b * dt) * dt, Small);
          // Get the real collision rate at the updated energy
          fReal = medium->GetCollisionRate(newEnergy);
          if (fReal > fLim) {
            // Real collision rate is higher than null-collision rate
            dt += log(r) / fLim;
            // Increase the null collision rate and try again
            std::cerr << "AvalancheMicroscopic::AvalancheElectron:" << std::endl;
            std::cerr << "    Increasing the null-collision rate by 5%." 
                      << std::endl;
            fLim *= 1.05;
            continue;
          }
          // Check for real or null collision
          if (RndmUniform() > fReal / fLim) continue;
          break;
        }
      
        ++nCollTemp;

        // Update the directions (at instant before collision)
        dt2 = dt * dt;
        a = sqrt(energy / newEnergy);
        b = 0.5 * c1 * dt / sqrt(newEnergy);
        newDx = dx * a + ex * b; 
        newDy = dy * a + ey * b; 
        newDz = dz * a + ez * b;
        // Update the position
        a = c1 * dt * sqrt(energy);
        b = dt2 * c2;            
        x += dx * a + ex * b;
        y += dy * a + ey * b;
        z += dz * a + ez * b;
        t += dt;
      
        // Verify the new position
        if (!sensor->IsInArea(x, y, z)) {
          stack[iEl].x = x; stack[iEl].y = y; stack[iEl].z = z;
          stack[iEl].t = t; stack[iEl].energy = energy;
          stack[iEl].dx = dx; stack[iEl].dy = dy; stack[iEl].dz = dz;          
          endpoints.push_back(stack[iEl]);
          stack.erase(stack.begin() + iEl);
          ok = false;
          break;
        }
        
        // Get the collision type and parameters
        medium->GetCollision(newEnergy, cstype, level, s1, ctheta0, eloss, esec);
        
        switch (cstype) {
          // Elastic collision
          case 0:
            break;
          // Ionising collision
          case 1:
            if (hasUserHandleIonisation) {
              userHandleIonisation(x, y, z, t, cstype, level, medium);
            }          
            // Randomise secondary electron direction
            phi = TwoPi * RndmUniform();
            ctheta = 1. - 2. * RndmUniform();
            stheta = sin(acos(ctheta));
            // Add the secondary electron to the stack
            newElectron = stack[iEl];
            newElectron.x0 = x; newElectron.x = x;
            newElectron.y0 = y; newElectron.y = y;
            newElectron.z0 = z; newElectron.z = z;
            newElectron.t0 = t; newElectron.t = t; 
            newElectron.energy = Max(esec, Small);
            newElectron.e0 = newElectron.energy;
            newElectron.dx = cos(phi) * stheta;
            newElectron.dy = sin(phi) * stheta;
            newElectron.dz = ctheta;
            newElectron.driftLine.clear();
            stack.push_back(newElectron);
            // Increment the electron and ion counters         
            ++nElectrons; ++nIons;
            break;          
          // Attachment
          case 2:          
            if (hasUserHandleAttachment) {
              userHandleAttachment(x, y, z, t, cstype, level, medium);
            }
            // Decrement the electron counter
            --nElectrons;
            stack[iEl].x = x; stack[iEl].y = y; stack[iEl].z = z; 
            stack[iEl].t = t; stack[iEl].energy = - energy;
            endpoints.push_back(stack[iEl]);
            stack.erase(stack.begin() + iEl);
            ok = false;
            break;
          // Inelastic collision
          case 3:
            if (hasUserHandleInelastic) {
              userHandleInelastic(x, y, z, t, cstype, level, medium);
            }          
            break;
          // Super-elastic collision
          case 4:
            break;
          default:
            std::cerr << "AvalancheMicroscopic::AvalancheElectron:" << std::endl;
            std::cerr << "    Unknown collision type." << std::endl;
            ok = false;
            break;
        }

        if (!ok) break;

        // Determine scattering angles
        s2 = (s1 * s1) / (s1 - 1.);
        stheta0 = sin(acos(ctheta0));
        phi0 = TwoPi * RndmUniform();
        sphi0 = sin(phi0);
        cphi0 = cos(phi0);
        
        arg = Max(1. - s1 * eloss / newEnergy, Small);
        d = 1. - ctheta0 * sqrt(arg);
        // Update the energy
        energy = Max(newEnergy * (1. - eloss / (s1 * newEnergy) - 2. * d / s2), Small);
        q = Min(sqrt((newEnergy / energy) * arg) / s1, 1.);
        theta = asin(q * stheta0);
        ctheta = cos(theta);
        if (ctheta0 < 0.) {
          double u = (s1 - 1.) * (s1 - 1.) / arg;
          if (ctheta0 * ctheta0 > u) ctheta *= -1.;
        }
        stheta = sin(theta);
        
        newDz = Min(newDz, 1.);       
        arg = sqrt(newDx * newDx + newDy * newDy);
        if (arg == 0.) {
          dz = ctheta;
          dx = cphi0 * stheta;
          dy = sphi0 * stheta;
        } else {
          a = stheta / arg;
          dz = newDz * ctheta + arg * stheta * sphi0;
          dy = newDy * ctheta + a * (newDx * cphi0 - newDy * newDz * sphi0);
          dx = newDx * ctheta - a * (newDy * cphi0 + newDx * newDz * sphi0);
        }
        
        // Continue with the next electron in the stack?
        if (nCollTemp > nCollSkip) break;
        
      }
      
      if (!ok) continue;
      
      // Normalise the direction vector
      d = sqrt(dx * dx + dy * dy + dz * dz);
      dx = dx / d; dy = dy / d; dz = dz / d;
      // Update the stack
      stack[iEl].energy = energy; stack[iEl].t = t;
      stack[iEl].x = x; stack[iEl].y = y; stack[iEl].z = z;
      stack[iEl].dx = dx; stack[iEl].dy = dy; stack[iEl].dz = dz;
      // Add a new point to the drift line (if enabled)
      if (useDriftLines) {
        point newPoint;
        newPoint.x = x; newPoint.y = y; newPoint.z = z; newPoint.t = t;
        stack[iEl].driftLine.push_back(newPoint);
      }
    }
  }
  nEndpoints = endpoints.size();
  
  return true;
    
}

}
