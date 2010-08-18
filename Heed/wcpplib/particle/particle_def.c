#include <stdlib.h>
#include <string.h>
#include "wcpplib/particle/particle_def.h"
#include "wcpplib/clhep_units/WPhysicalConstants.h"
#include "wcpplib/stream/prstream.h"
#include "wcpplib/util/FunNameStack.h"
/*
1998 - 2004,   I. Smirnov
*/
spin_def::spin_def(float ftotal, float fprojection):
  total(ftotal), projection(fprojection) 
{
  mfunname("spin_def::spin_def(float ftotal, float fprojection)");
  check_econd11(total , < 0, mcerr); 
  check_econd12(total , < , projection, mcerr); 
}

 ostream & operator << (ostream & file, const spin_def & f)
{
  Ifile<<"spin_def: total="<<f.total<<" projection="<<f.projection;
  return file;
}


//AbsList< particle_def* > particle_def::logbook;
// This call should be before all particle_def
// otherwise they are all lost.


particle_def electron_def("electron", "e-", electron_mass_c2/c_squared, 
		      electron_charge, 1, 0, 0.5, spin_def(0.0, 0.0));
particle_def positron_def("positron", "e+", electron_def);
particle_def proton_def("proton", "p+", proton_mass_c2/c_squared, eplus, 
			0, 1, 0.5, spin_def(0.5, 0.5));
particle_def anti_proton_def("", "p-", proton_def);
particle_def neutron_def("neutron","n", neutron_mass_c2/c_squared, 
			 0, 0, 0, 0.5, spin_def(0.5, -0.5));
particle_def anti_neutron_def("","",neutron_def);

particle_def P11_def("P11", "P11", 
		       1440.0*MeV/c_squared, 1*eplus, 
		       0, 1, 0.5, spin_def(0.5, 0.5));
particle_def D13_def("D13", "D13", 
		       1520.0*MeV/c_squared, 1*eplus, 
		       0, 1, 1.5, spin_def(0.5, 0.5));
particle_def S11_def("S11", "S11", 
		       1535.0*MeV/c_squared, 1*eplus, 
		       0, 1, 0.5, spin_def(0.5, 0.5));

//particle_def P_1_1_def("roper", "roper", 
//		       1440.0*MeV/c_squared, 1*eplus, 
//		       0, 1, 0.5, spin_def(0.5, 0.5));

// light unflavored mesons
particle_def pi_plus_meson_def("pi_plus_meson","pi+",139.56755*MeV/c_squared, 
			       eplus, 0, 0, 1.0, spin_def(1.0, 1.0));
particle_def pi_minus_meson_def("pi_minus_meson",
				"pi-",139.56755*MeV/c_squared, 
				-eplus, 0, 0, 1.0, spin_def(1.0, -1.0));
particle_def pi_0_meson_def("pi_0_meson","pi0",134.9734*MeV/c_squared, 
			    0, 0, 0, 0, spin_def(1.0, 0.0));
particle_def eta_meson_def("eta_meson_def","eta", 548.8*MeV/c_squared, 
			       0, 0, 0, 1.0, spin_def(0.0, 0.0));
particle_def K_plus_meson_def("K_plus_meson_def","K+", 493.677*MeV/c_squared, 
			       1, 
			      0, 0, 0.0, spin_def(0.0, 0.0)); // this I
// don't know for the moment.


particle_def deuteron_def("deuteron","dtr",
				1875.613*MeV/c_squared, 
				eplus, 0, 0, 0.0, spin_def(0.0, 0.0));
particle_def alpha_particle_def("alpha_particle","alpha",
				3727.417*MeV/c_squared, 
				2*eplus, 0, 0, 0.0, spin_def(0.0, 0.0));
/*
particle_def electron("electron", 0.51099906, -1, 1, 0, 0.5, 0.0);
particle_def positron("positron",electron);
particle_def proton("proton", 938.27231, 1, 0, 1, 0.5, 0.5);
particle_def anti_proton("",proton);
particle_def neutron("neutron",939.56563, 0, 0, 0, 0.5, -0.5);
particle_def anti_neutron("",neutron);
particle_def pi_plus_meson("pi_plus_meson",139.56755, 1, 0, 0, 1.0, 1.0);
particle_def pi_minus_meson("pi_plus_meson",139.56755, -1, 0, 0, -1.0, -1.0);
particle_def pi_0_meson("pi_0_meson",134.9734, 0, 0, 0, 0, 0);
particle_def alpha_particle("alpha_particle",3727.417, 2, 0, 0, 0.0, 0.0);
*/
/*
particle_def* allapardef[pqallapardef]={
  &electron_def,
  &positron_def,
  &proton_def,
  &anti_proton_def,
  &neutron_def,
  &anti_neutron_def,
  &pi_plus_meson_def,
  &pi_0_meson_def,
  &pi_minus_meson_def,
  &alpha_particle_def,
  NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 
  NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 
  NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 
  NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL}; 
int qallapardef=10; // user can define more particles
*/
particle_def::particle_def(const String& fname, const String& fnotation,
			   double fmass, double fcharge, 
			   int flepton_n, int fbarion_n, 
			   float fspin, const spin_def& fisospin)
{
  //strcpy(name,fname);
  name=fname;
  notation=fnotation;
  //mcout<<"particle_def::particle_def: name="<<name<<'\n';
  mass=fmass;
  charge=fcharge;
  barion_n=fbarion_n;
  lepton_n=flepton_n;
  spin=fspin; 
  isospin=fisospin;
  verify();
  particle_def::get_logbook().append(this);
  //printall(mcout);
}

particle_def::particle_def(const String& fname, 
			   const String& fnotation, particle_def& p) 
    // creates anti-particle through the call of anti_particle(p)
{
  *this = anti_particle(p);
  //if(strlen(fname) > 0)
  //strcpy(name,fname);
  if(!(fname=="" || fname==" "))
    name=fname;
  if(!(fnotation=="" || fnotation==" "))
    notation=fnotation;
  verify(); 
  particle_def::get_logbook().append(this);
}

/*
void particle_def::operator=(const particle_def& f)
{
  name=f.name;
  notation=f.notation;
  mass=f.mass;
  charge=f.charge;
  barion_n=f.barion_n;
  lepton_n=f.lepton_n;
  spin=f.spin; 
  isospin=f.isospin;
}
*/
particle_def particle_def::anti_particle(const particle_def& p)
{
  String aname=String("anti-") + p.name;
  String anot=String("anti-") + p.notation;
  //char s[100];
  //strcpy(s,"anti-");
  //strcat(s,p.name);
  return particle_def(aname, anot, p.mass, -p.charge, 
		      -p.lepton_n, -p.barion_n, -p.spin, p.isospin);
}  
AbsList< particle_def* >& particle_def::get_logbook(void)
{
  static AbsList< particle_def* > logbook;
  return logbook;
}

const AbsList< particle_def* >& particle_def::get_const_logbook(void)
{
  return particle_def::get_logbook();
}

particle_def* particle_def::get_particle_def(const String& fnotation)
{
  AbsList< particle_def* >& logbook = particle_def::get_logbook();
  AbsListNode<particle_def*>* an=NULL;
  while( (an = logbook.get_next_node(an)) != NULL)
  { 
    if(an->el->notation == fnotation)
    {
      return an->el;
    }
  }
  return NULL;
}

void particle_def::print(ostream & file)
{
  file<<(*this);
  /*
  file<<name<<" mass="<<mass<<" mass/(GeV/c_squared)="<<mass/(GeV/c_squared)
      <<" charge="<<charge<<" charge/eplus="<<charge/eplus
      <<" lepton_n="<<lepton_n<<" barion_n="<<barion_n
      <<" spin="<<spin<<" isospin="<<isospin<<'\n';
  */
}
void particle_def::printall(ostream & file)
{
  Ifile<<"particle_def::printall:\n";
  particle_def* apd=NULL;
  int n;
  AbsList< particle_def* >& logbook = particle_def::get_logbook();
  AbsListNode<particle_def*>* an=NULL;
  while( (an = logbook.get_next_node(an)) != NULL)
  { 
    file<<(*(an->el));
  }
  /*
#ifdef USE_STLLIST
  list< particle_def* >::const_iterator i = cont.begin();
  while(i != cont.end() )
  {
    file<<(**i);
    i++;
  }
#else
  for( n=0; (apd = particle_def::cont[n]) != NULL ; n++)
  {
    file<<(*apd);
  }
#endif
  */
}


/*
void particle_def::verify(void)
{
  mfunname("void particle_def::verify(void)");
  if(name != "none")
  {
    particle_def* apd=NULL;
    int n;
    for( n=0; (apd = cont[n]) != NULL ; n++)
    {
      if( name == apd->name)
      {
	funnw.ehdr(mcerr);
	mcerr<<"another registered particle definition with the same name found\n";
	spexit(mcerr);
      }
    }
  }
} 
*/
ostream & operator << (ostream & file, const particle_def & f)
{
  Ifile<<"particle_def: name="<<f.name<<" notation="<<f.notation<<'\n';
  Ifile<<"mass="<<f.mass<<" mass/(GeV/c_squared)="<<f.mass/(GeV/c_squared)
       <<" charge="<<f.charge<<" charge/eplus="<<f.charge/eplus<<'\n';
  Ifile<<"lepton_n="<<f.lepton_n<<" barion_n="<<f.barion_n<<'\n';
  Ifile<<"spin="<<f.spin<<" isospin="<<f.isospin<<'\n';
  return file;
}

particle_type::particle_type(const char* name, int s)
{
  mfunname("particle_type::particle_type(const char* name, int s)");
  particle_def* apd=NULL;
  int n;
  //mcout<<"particle_type::particle_type(char* name):\n";
  //particle_def::printall(mcout);
  AbsListNode<particle_def*>* an=NULL;
  AbsList< particle_def* >& logbook = particle_def::get_logbook();
  while( (an = logbook.get_next_node(an)) != NULL)
  { 
   if( name == an->el->notation)
   { 
     pardef = an->el;
     return; 
   }
  }
  an = NULL; // to start from beginning
  while( (an = logbook.get_next_node(an)) != NULL)
  { 
   if( name == an->el->name)
   { 
     pardef = an->el;
     return; 
   }
  }
  /*
#ifdef USE_STLLIST
  list< particle_def* >::const_iterator i = particle_def::cont.begin();
  while(i != particle_def::cont.end() )
  {
    apd = *i;
#else
  for( n=0; (apd = particle_def::cont[n]) != NULL ; n++)
  {
#endif
    //mcout<<"apd->name="<<apd->name<<'\n';
    if( name == apd->notation)
      //if( ! strcmp(name, apd->name) )
    {
      pardef = apd;
      return;
    }
#ifdef USE_STLLIST
    i++;
#endif
  }
#ifdef USE_STLLIST
  i = particle_def::cont.begin();
  while(i != particle_def::cont.end() )
  {
    apd = *i;
#else
  for( n=0; (apd = particle_def::cont[n]) != NULL ; n++)
  {
#endif
    //mcout<<"apd->name="<<apd->name<<'\n';
    if( name == apd->name)
      //if( ! strcmp(name, apd->name) )
    {
      pardef = apd;
      return;
    }
#ifdef USE_STLLIST
    i++;
#endif
  }
  */
  /*
  for( n=0; n<qallapardef; n++)
  {
    //mcout<<"name="<<name<<" allapardef[n]->name="<<allapardef[n]->name<<'\n';
    if(!strcmp(name, allapardef[n]->name ))
    {
      pardef = allapardef[n];
      return;
    }
  }
  */
  if(s==0)
  {
    mcerr<<"this type of particle is absent, name="<<name<<'\n';
    spexit(mcerr);
  }
  pardef = NULL;
}

void particle_type::print_notation(ostream& file) const
{
  if(pardef.get()==NULL)
    file<<"none";
  else
    file<<pardef->notation;
}

ostream& operator<< (ostream& file, const particle_type& f)
{
  if(f.pardef.get()==NULL)
    file<<"type is not initialized";
  else
    file<<(f.pardef->name);
  return file;
}

