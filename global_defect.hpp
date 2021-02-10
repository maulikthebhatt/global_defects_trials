#ifndef GLOBAL_DEFECT_HPP
#define GLOBAL_DEFECT_HPP

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <sys/time.h>
#include<thread>
#include<chrono>
#include <fstream>
#include <bits/stdc++.h>

#include "LATfield2.hpp"
using namespace LATfield2;
#include "powerSpectra.hpp"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <bits/stdc++.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "metadata.hpp"
#include "parser.hpp"
#include "background.hpp"

class GlobalDefect
{
private:

    double dx_;
    Lattice * lat_;
    Lattice * klat_;
    metadata * sim_;
    defects_metadata * defects_sim_;
    
    Field<double> phi_defect_;
    Field<double> pi_defect_;

public:
  void initialize(Lattice * lat, Lattice * klat, double dx, metadata * sim, defects_metadata * defects_sim);
  void generate_init_cond();
  void update_phi(double dt);
  void update_pi(double dt, double a, double adot_overa);
  
  unsigned long int random_seed();
  double potentialprime(Site & x, int comp);
};


void GlobalDefect::initialize(Lattice * lat, Lattice * klat, double dx, metadata * sim, defects_metadata * defects_sim)
{
  dx_ = dx;
  lat_ = lat;
  klat_ = klat;
  defects_sim_ = defects_sim;
  sim_ = sim;

  phi_defect_.initialize(*lat_, defects_sim_->nComponents);
  phi_defect_.alloc();

  pi_defect_.initialize(*lat_, defects_sim_->nComponents);
  pi_defect_.alloc();

  generate_init_cond();
}


unsigned long int GlobalDefect::random_seed()
{
  struct timeval tv;
  gettimeofday(&tv,0);
  return(tv.tv_sec + tv.tv_usec);
}


void GlobalDefect::generate_init_cond()
{
  Site x(phi_defect_.lattice());

  const gsl_rng_type * T;
  gsl_rng * r;

  gsl_rng_env_setup();
  
  gsl_rng_default_seed = random_seed();

  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  for(x.first();x.test();x.next())
  {
    double phiNorm2 = 0;
    for(int c=0;c< defects_sim_->nComponents;c++)
    {
      phi_defect_(x,c) = gsl_ran_gaussian (r,1);
      phiNorm2 += phi_defect_(x,c)*phi_defect_(x,c);
    }
    double ratio =  sqrt( defects_sim_->eta2/phiNorm2);
    for(int c=0; c< defects_sim_->nComponents;c++)
    {
      phi_defect_(x,c) *= ratio;
      pi_defect_(x,c) = 0;
    }
  }
  
  gsl_rng_free (r);

  COUT<< "Initial Condition for defect is generated"<<endl<<endl;
}


void GlobalDefect::update_phi(double dt)
{ 
  double dt_ = dt;
  Site x(phi_defect_.lattice());
  for(x.first();x.test();x.next())
  {
    for(int c = 0; c < defects_sim_->nComponents;c++)
    {
      phi_defect_(x,c) += dt_ * pi_defect_(x,c);
    }
  }

  phi_defect_.updateHalo(); 
}

double GlobalDefect::potentialprime(Site & x, 
                                        int comp)
{
  double phiNorm2 = 0;
  for(int i =0; i < defects_sim_->nComponents; i++)phiNorm2 += phi_defect_(x,i)*phi_defect_(x,i);
  return 2.0 * defects_sim_->lambda * ( phiNorm2 - defects_sim_->eta2) *  phi_defect_(x,comp);
}


void GlobalDefect::update_pi(double dt,
                                double a, 
                                double adot_overa)
{
  double dt_ = dt;
  double a_ = a;
  double adot_overa_ = adot_overa;
  
  Site x(pi_defect_.lattice()); 

  double c1 = (1.0 - dt_ * (adot_overa_)) / (1.0 + dt_ * (adot_overa_));
  double c2 = dt_ / (1.0 + dt_ * adot_overa_);
  double a2 = a_*a_;


  for(x.first();x.test();x.next())
  {
    for(int c = 0; c < defects_sim_->nComponents; c++)
    {
      double lapPhi = -6.0 * phi_defect_(x,c) ;
      for(int i = 0 ; i<3 ; i++)lapPhi += phi_defect_(x+i,c) + phi_defect_(x-i,c);
      lapPhi /= dx_*dx_;
      pi_defect_(x,c) = c1 * pi_defect_(x,c) + c2 * ( lapPhi -  a2 * potentialprime(x,c) );
    }
 }

}
#endif
