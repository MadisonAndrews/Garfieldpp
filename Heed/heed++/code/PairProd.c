#include <fstream>
#include <cmath>
#include <cfloat>
#include <climits>
#include "heed++/code/PairProd.h"
#include "wcpplib/random/ranluxint.h"
/*
2003, I. Smirnov
*/
//#define USE_GET_ELOSS_CUT
#ifdef USE_GET_ELOSS_CUT
const double w_cut_ratio = 0.2;
#else
const double V_ratio = 0.5;
#endif

PairProd::PairProd(const String& file_name, double fwa, double ffactorFano):
  wa(fwa), factorFano(ffactorFano)
{
  mfunnamep("PairProd::PairProd(const String& file_name, double fwa, double ffactorFano)");
  
#ifdef USE_STLSTRING
  ifstream file(file_name.c_str());
#else
  ifstream file(file_name);
#endif
  if( !file )
  {
    funnw.ehdr(mcerr);
    mcerr<<"cannot open file "<<file_name<<endl;
    spexit(mcerr);
  }
  long q;
  file>>wa_table>>I_table>>J_table>>factorFano_table>>q;
  if( !file.good() )
  {
    funnw.ehdr(mcerr);
    mcerr<<"error at reading file"<<endl;
    spexit(mcerr);
  }
  DynLinArr< double > xx(q);
  DynLinArr< double > yy(q);
  long n;
  for( n=0; n<q; n++)
  {
    file>>xx[n]>>yy[n];
  }
  pran = PointsRan(xx, yy, I_table, J_table);
  k = sqrt( factorFano * wa * wa / (factorFano_table * wa_table * wa_table) );
  s = wa  - k * wa_table;
  //s = wa - wa / factorFano * ( 2.0 * factorFano - factorFano_table);
  //k = wa / (factorFano * wa_table) * ( 2.0 * factorFano - factorFano_table);
}

double PairProd::get_eloss(void) const
{
  mfunname("double PairProd::get_eloss(void) const");
  //mcout<<"PairProd::get_eloss started\n";
  double e_loss = k * pran.ran(SRANLUX()) + s; // working wariant
  //double e_loss = pran.ran(SRANLUX()); // for debug
  //mcout<<"PairProd::get_eloss finished\n";
  return e_loss;
}

#ifdef USE_GET_ELOSS_CUT

double PairProd::get_eloss(double e_cur) const
{
  mfunname("double PairProd::get_eloss(double ecur) const");
  //mcout<<"PairProd::get_eloss started\n";
  double e_loss = k * pran.ran(SRANLUX()) + s;
  if(e_cur - e_loss <  w_cut_ratio * wa)
    e_loss = 1.0e20;  // to stop tracing

  //mcout<<"PairProd::get_eloss finished\n";
  return e_loss;
}

#else

double PairProd::get_eloss(double e_cur) const
{
  mfunname("double PairProd::get_eloss(double ecur) const");
  //mcout<<"PairProd::get_eloss started\n";
  double e_loss = k * pran.ran(SRANLUX()) + s;
  double c;
  if( e_cur <= V_ratio * wa)
  {
    c = DBL_MAX;
  }
  else
  {
    //c = 1.0 / (1.0 - V_ratio * wa / e_cur);
    c = 1.0 / (1.0 - pow( V_ratio * wa / e_cur, 2.0 ));
  }
  e_loss = e_loss * c;

  //mcout<<"PairProd::get_eloss finished\n";
  return e_loss;
}

#endif

void PairProd::print(ostream& file) const 
{
  Ifile<<"PairProd:\n";
  indn.n+=2;
  Ifile<<"wa="<<wa<<" factorFano="<<factorFano<<'\n';
  Ifile<<"wa_table="<<wa_table<<" factorFano_table="<<factorFano_table<<'\n';
  Ifile<<" I_table="<<I_table<<" J_table="<<J_table<<" k="<<k<<" s="<<s<<'\n';
  pran.print(file);
  indn.n-=2;
}
  
 

