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

class Global_defect
{
private:
    string runID_;
    string path_;
    int step = 0;
    
    Lattice lat_;
    Lattice klat_;
    Field<double> phi_defect_;
    Field<double> pi_defect_;

    Field<double> pi_defect_prev_;
    Field<double> rho_;
    Field<double> P_;
    
    Field<Imag> rho_k_;
    
    PlanFFT<Imag> planrho_;

    double dt_;
    double dx_;
    
    double a_;
    double adot_overa_;
    
    metadata sim;
    global_defects gdefects;

public:
    Global_defect(){;}
    Global_defect(const Lattice& latx, const Lattice& latk, const double& D_x, const double& D_t, const global_defects & gdefects)
    {
      initialize(latx, latk, D_x, D_t, gdefects);
    }
    ~Global_defect()
    {

    }
    unsigned long int random_seed();
    void initialize(const Lattice& latx, const Lattice& latk, const double& D_x, const double& D_t, const global_defects & gdefects);

    void generate_initCond(const global_defects & gdefects);
    
    void evolve(const double& dt, const double& sf, const double& adota);
    void field_leapfrog_update_phi();
    void field_leapfrog_update_pi();
    double potential(Site & x);
    double potentialprime(Site & x, int comp);
    
    void averagephidefect();
    void averagerhodefect();
    
    double modsqphi(Site & x);
    double modsqrho(Site & x);
    
    void compute_rho_P_();
    void compute_pk_();
    
    bool output_now();
    void output();

};

////////////////////////////////////////////////////////////////
////    Initialize
//
////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
//  This function initialises the value 
//
// Parameters:
//     dx_ = the distance between the lattice points
//     dt_ = the time step taken
//     lat_ = the real space lattice
//     klat_ = the fourier space lattice
//     phi_defect = the scalar field of the defect
//     pi_defect = the phi dot term used for the leapfrog algorithm
//     rho_ = the T00 term of the defect
//     rho_k_ = the T00 term of the defect in real space
//     P_ = the Tii term of the defect
//     rho_c0 = the current critical density
//     
///////////////////////////////////////////////////////////////

void Global_defect::initialize(const Lattice& latx, const Lattice& latk, const double& D_x, const double& D_t, const global_defects & gdefects)
{
    //loadSettings(settings_filename);
    path_ = "/media/maulik/DATA/UNIGE/Thesis/gevolution-defect/defect/";
    runID_ = "001";
    
    lat_ = latx;
    
    klat_ = latk;
    
    dx_ = D_x;
    
    dt_ = D_t;
    
    //dx_ = physicalSize_/latticeSize_;
    //dt_ = courantFactor_ * dx_;

    //lat_.initialize(3,latticeSize_,1);
    //klat_.initializeRealFFT(lat_,0);
    
    phi_defect_.initialize(lat_,gdefects.nComponents);
    phi_defect_.alloc();

    pi_defect_.initialize(lat_,gdefects.nComponents);
    pi_defect_.alloc();
	
    
    cout << endl << endl << " NO OF COMPONENTS ARE = " << gdefects.nComponents << endl << endl;
     
    pi_defect_prev_.initialize(lat_,gdefects.nComponents);
    pi_defect_prev_.alloc();

    rho_.initialize(lat_);
    rho_.alloc();
    
    rho_k_.initialize(klat_);
    rho_k_.alloc();
    planrho_.initialize(&rho_,&rho_k_);
  //  planrho_.alloc();  
    
    P_.initialize(lat_);
    P_.alloc();

    generate_initCond(gdefects);
    
    //COUT<<"Saving the files to:"<<path_<<endl<<endl;
    COUT<< "The initial conditions set are:"<<endl;
    //COUT<< "start time is = "<<t_start_<<endl;
    //COUT<< " End time is = "<<t_end_<<endl;
    COUT<< "Time interval is = "<<dt_<<endl;
    COUT<< "Lattice interval is = "<< dx_<<endl<<endl;
    //COUT<< "The value set for univ is: "<<univ<<endl;
    //COUT<< "Omega matter is = "<< omega_m<<endl;
    //COUT<< "Omega radiation is = "<< omega_r<<endl;
    //COUT<< "Omega lambda is = "<< omega_lambda<<endl;
    //COUT<< "H0 is = "<< H0 <<endl;
    //COUT<< "initial scale factor is = "<< a_i<<endl;
    //COUT<< "The G is ="<<G<<endl;
    //COUT<<"The current critical density is ="<<rho_c0<<endl;
    //COUT<<"The val of pi is="<<M_PI<<endl<<endl;
}

/////////////////////////////////////////////////////////////
//  random_seed
//
/////////////////////////////////////////////////////////////
// This function is used to generate the initial seed to 
// generate the phi_defect and pi_defect 
//
////////////////////////////////////////////////////////////

unsigned long int Global_defect::random_seed()
    {
     struct timeval tv;
     gettimeofday(&tv,0);
     return(tv.tv_sec + tv.tv_usec);
    }

/////////////////////////////////////////////////////////////
//  generate_initCond
//
/////////////////////////////////////////////////////////////
// This function generates the initial conditions of phi_defect
// and pi_defect.In this the initial conditions are 
// taken to be a gaussian field
//
////////////////////////////////////////////////////////////

void Global_defect::generate_initCond(const global_defects & gdefects)
{

    // TODO: understand how to change seed with gsl....
    Site x(lat_);

    const gsl_rng_type * T;
    gsl_rng * r;

    gsl_rng_env_setup();
    
    gsl_rng_default_seed = random_seed();

    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

    for(x.first();x.test();x.next())
    {
      double phiNorm2 = 0;
      for(int c=0;c<gdefects.nComponents;c++)
      {
        phi_defect_(x,c) = gsl_ran_gaussian (r,1);
        phiNorm2 += phi_defect_(x,c)*phi_defect_(x,c);
      }
      double ratio =  sqrt(gdefects.eta2/phiNorm2);
      for(int c=0;c<gdefects.nComponents;c++)
      {
        phi_defect_(x,c) *= ratio;
        pi_defect_(x,c) = 0;
      }
    }
    
    gsl_rng_free (r);

    phi_defect_.saveHDF5(path_ + runID_ + "_phi_defect_initCond.h5");
    pi_defect_.saveHDF5(path_  + runID_ + "_pi_defect_initCond.h5");

    COUT<< "Initial Condition for defect is generated"<<endl<<endl;

}

/////////////////////////////////////////////////////////////
//  evolve
//
/////////////////////////////////////////////////////////////
// This function is the main function which has steps to 
// evolve the fields using the leapfrog algorithm
//
////////////////////////////////////////////////////////////

void Global_defect::evolve(const double& dt, const double& sf, const double& adota)
{
    dt_ = dt;
    a_ = sf;
    adot_overa_ = adota;
	     
    COUT<< "Evolving defects after dt "<<dt_<<" and at scale factor "<<a_<< " and adota = " << adot_overa_ <<endl;
    if(step ==0)
    {
        averagephidefect();
        averagerhodefect();
    }
    else
    {
        field_leapfrog_update_phi();
        field_leapfrog_update_pi();
        compute_rho_P_();
        averagephidefect();
        averagerhodefect();
        
//        if(output_now())
//        {
//          output();
//        }
    }
    step++;
}

/////////////////////////////////////////////////////////////
//  field_leapfrog
//
/////////////////////////////////////////////////////////////
// This function has the leapfrog algorithm implemented on 
// the equation of motion of the global defect.
//
////////////////////////////////////////////////////////////
// Parameters:
//    
//    fric_term_1 = the coefficient used to change the friction term set using the settings file
//    fric_term_2 = the coefficient used to change the friction term set using the settings file
/////////////////////////////////////////////////////////////

void Global_defect::field_leapfrog_update_phi()
{

    Site x(lat_);

    for(x.first();x.test();x.next())
    {
      for(int c = 0;c<gdefects.nComponents;c++)
      {

        phi_defect_(x,c) += dt_ * pi_defect_(x,c);
      }
    }

    phi_defect_.updateHalo(); //update the value of phi in the halo
}

void Global_defect::field_leapfrog_update_pi()
{
    Site x(lat_);
    double fric_term;
    
    if(gdefects.dissipation)
    {
        if (step<gdefects.diss_end)
        {
          fric_term = gdefects.friction_coeff;
        }
        else
        {
          fric_term = 1;
        }
    }
    else
    {
        fric_term = 1;
    }

    double c1 = (1.0 - dt_ * (fric_term*adot_overa_)) / (1.0 + dt_ * (fric_term*adot_overa_));
    double c2 = dt_ / (1.0 + dt_ * adot_overa_);
    double a2 = a_*a_;

    // put what is in pi in pi_defect_prev
    //.... we then switch the data between pi and pi_defect_prev:
    double * temp = pi_defect_prev_.data_;
    pi_defect_prev_.data_ = pi_defect_.data_;
    pi_defect_.data_ = temp;

    for(x.first();x.test();x.next())
    {
      for(int c = 0;c<gdefects.nComponents;c++)
      {

        double lapPhi = -6.0 * phi_defect_(x,c) ;
        for(int i = 0 ; i<3 ; i++)lapPhi += phi_defect_(x+i,c) + phi_defect_(x-i,c);
        lapPhi /= dx_*dx_;
        pi_defect_(x,c) = c1 * pi_defect_prev_(x,c) + c2 * ( lapPhi -  a2 * potentialprime(x,c) );
      }
    }

}

/////////////////////////////////////////////////////////////
//  modsqphi
//
/////////////////////////////////////////////////////////////
// This function computes.the modulus of phi_defect square
//
////////////////////////////////////////////////////////////

double Global_defect::modsqphi(Site &x)
{
    double phiNorm2 = 0;
    for(int i =0;i<gdefects.nComponents;i++)phiNorm2 += phi_defect_(x,i)*phi_defect_(x,i);
    return pow(phiNorm2,0.5);
}

/////////////////////////////////////////////////////////////
//  averagephidefect
//
/////////////////////////////////////////////////////////////
// This function computes the <phi_defect> and writes it to a file
//
////////////////////////////////////////////////////////////

void Global_defect::averagephidefect()
{
    Site x(lat_);
    double phisum_ = 0;
    double phiavg_ = 0;
    for(x.first();x.test();x.next())
    {
      phisum_ += modsqphi(x);
    }
    parallel.sum(phisum_);
    phiavg_ = phisum_/pow(sim.numpts,3);
    
    if(parallel.rank() == 0)
	  {
      ofstream phifile;
      phifile.open (path_ + runID_+"_average_phi_defect.txt",std::ios_base::app);
      phifile << a_ <<" "<<phiavg_<<endl;
      phifile.close();
    }
    //averageField(phi_defect_,"/media/vilasini/DATA/UNIGE/Thesis/plots/"+ runID_+"_average_phi.txt",0);
}

/////////////////////////////////////////////////////////////
//  averagerhodefect
//
/////////////////////////////////////////////////////////////
// This function computes the <rho_> and writes it to a file
//
////////////////////////////////////////////////////////////

void Global_defect::averagerhodefect()
{
    Site x(lat_);
    double rhosum_ = 0;
    double rhoavg_ = 0;
    for(x.first();x.test();x.next())
    {
      rhosum_ += rho_(x);
    }
    parallel.sum(rhosum_);
    rhoavg_ = rhosum_/pow(sim.numpts,3);
    if(parallel.rank() == 0)
    {
      ofstream rhofile;
      rhofile.open (path_ + runID_ + "_average_rho_defect.txt",std::ios_base::app);
      rhofile << a_ <<" "<<rhoavg_<<endl;
      rhofile.close();
    }
}

////////////////////////////////////////////////////
//  potential
///////////////////////////////////////////////////
// This function computes the potential of the scalar field
//
//////////////////////////////////////////////////

double Global_defect::potential(Site & x)
{
    double phiNorm2 = 0;
    for(int i =0;i<gdefects.nComponents;i++)phiNorm2 += phi_defect_(x,i)*phi_defect_(x,i);
    return gdefects.lambda * ( phiNorm2 - gdefects.eta2) * ( phiNorm2 - gdefects.eta2) / 2.0;
}

////////////////////////////////////////////////
//  potentialprime
//
////////////////////////////////////////////////
// This function computes the derivative of potential
//
///////////////////////////////////////////////

double Global_defect::potentialprime(Site & x, int comp)
{
    double phiNorm2 = 0;
    for(int i =0;i<gdefects.nComponents;i++)phiNorm2 += phi_defect_(x,i)*phi_defect_(x,i);
    return 2.0 * gdefects.lambda * ( phiNorm2 - gdefects.eta2) *  phi_defect_(x,comp);
}


bool Global_defect::output_now()
{
    return step%10==0?true:false;
}

////////////////////////////////////////////////
//  output
//
////////////////////////////////////////////////
// This function writws the output 
//
///////////////////////////////////////////////

void Global_defect::output()
{
    COUT<<"outputing field at a="<< a_ <<endl;
    string filename_end= int2string(a_,9999)+".h5";
    

    phi_defect_.saveHDF5(path_ + runID_ + "_phi_defect_" + filename_end);
    pi_defect_.saveHDF5(path_ + runID_ + "_pi_defect_" + filename_end);

    //computeT();
    rho_.saveHDF5(path_ + runID_ + "_rho_" + filename_end);
    compute_pk_();
}

////////////////////////////////////////////////
//  compute_rho_P
//
////////////////////////////////////////////////
// This function computes the T00 and Tii 
//
///////////////////////////////////////////////

void Global_defect::compute_rho_P_()
{
    Site x(lat_);

    double a2 = a_*a_;
    for(x.first();x.test();x.next())
    {
      double mpidot = 0;
      double temp;
      double gradPhi2 = 0;
      for(int c=0;c<gdefects.nComponents;c++)
      {
        temp = (pi_defect_prev_(x,c)+pi_defect_(x,c))/2.0;
        mpidot = temp*temp;
        for(int i = 0;i<3;i++)
        {
          temp = ( phi_defect_(x+i,c) - phi_defect_(x-i,c) ) / 2.0 / dx_;
          gradPhi2 += temp*temp;
        }
      }  
      rho_(x) = mpidot / 2.0 / a2 + potential(x) + gradPhi2  / 2.0 / a2;
      P_(x) = mpidot / 2.0 / a2 - potential(x) - gradPhi2 / 6.0 / a2;
    }
}

////////////////////////////////////////////////
//  compute_pk
//
////////////////////////////////////////////////
// This function computes the powerspectrum 
//
///////////////////////////////////////////////
//  Parameters:
//
//    rho_k_ = the field used to compute the pk
//    binnos = contains the no of bins in k: taken from settings file
//    physicalSize_ = the box size: taken from settings file 
///////////////////////////////////////////////

void Global_defect::compute_pk_()
{
    planrho_.execute(FFT_FORWARD); 
    string filename_end= path_ + runID_ + "_powerspectrum" + int2string(a_,99999)+".txt";
    output_powerSpectrum(rho_k_,
                            filename_end,
                            sim.numbins,
                            sim.boxsize,
                            false,
                            false,
                            true,
                            false);
}

#endif


