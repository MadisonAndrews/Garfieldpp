#ifndef HEEDPHOTON_H
#define HEEDPHOTON_H

#include "heed++/code/HeedMatterDef.h"
#include "wcpplib/geometry/gparticle.h"
/*
Definition of the photon which can be emitted at atomic relaxation cascades
and traced through the geometry.

2003, I. Smirnov
*/

//#define SFER_PHOTOEL  // make direction of photoelectron absolutely random

class HeedPhoton: public gparticle
{public:
  //PassivePtr< EnTransfCS > encs;
  long particle_number;
  long parent_particle_number;
  int s_print_listing;  // convenient to print internal algorithms 
                        // of a selected event
  
  //PassivePtr< gparticle > parent_part;
  double energy;  // MeV
  int s_photon_absorbed;  // used in physics_after_new_speed
  long na_absorbing;  // number of absorbing atom
  long ns_absorbing;  // number of absorbing shell
#ifdef SFER_PHOTOEL
  int s_sfer_photoel;
#endif
  // sign that delta-electrons are already generated (or cannot be created) 
  int s_delta_generated;  

  HeedPhoton(manip_absvol* primvol, const point& pt, 
	     const vec& vel, vfloat time, 
	     long fparent_particle_number,
	     //PassivePtr< gparticle > fparent_part,
	     double fenergy,
	     int fs_print_listing = 0);
  //EnTransfCS* encs);
  void physics_after_new_speed(void);
  virtual void physics(void);
  HeedPhoton(void): gparticle() {;}
  virtual void print(std::ostream& file, int l) const;
  macro_copy_total(HeedPhoton);
  virtual ~HeedPhoton() {;}
  
};

#endif
