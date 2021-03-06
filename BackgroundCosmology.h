#ifndef _BACKGROUNDCOSMOLOGY_HEADER
#define _BACKGROUNDCOSMOLOGY_HEADER
#include <iostream>
#include <fstream>
#include "Utils.h"

using Vector = std::vector<double>;

class BackgroundCosmology{
  private:

    // Cosmological parameters
    double h;                       // Little h = H0/(100km/s/Mpc)
    double OmegaB;                  // Baryon density today
    double OmegaCDM;                // CDM density today
    double OmegaLambda;             // Dark energy density today
    double Neff;                    // Effective number of relativistic species (3.046 or 0 if ignoring neutrinos)
    double TCMB;                    // Temperature of the CMB today in Kelvin

    // Derived parameters
    double OmegaR;                  // Photon density today (follows from TCMB)
    double OmegaNu;                 // Neutrino density today (follows from TCMB and Neff)
    double OmegaK;                  // Curvature density = 1 - OmegaM - OmegaR - OmegaNu - OmegaLambda
    double H0;                      // The Hubble parameter today H0 = 100h km/s/Mpc

    //double OmegaR = 2*pow(M_PI,2)/30 * pow(Constants.k_b*TCMB,4)/ (pow(Constants.hbar,3)*pow(Constants.c,5))*8*M_PI*Constants.G/(3*pow(Constants.H0_over_h,2));
    //double OmegaLambda = 1 - (OmegaK+OmegaB+OmegaCDM+Omega_gamma+OmegaNu);

    // Start and end of x-integration (can be changed)
    double x_start = Constants.x_start;
    double x_end   = Constants.x_end;


    // Splines to be made
    Spline eta_of_x_spline{"eta"};
    Spline t_of_x_spline{"t"};



  public:

    // Constructors
    BackgroundCosmology() = delete;
    BackgroundCosmology(
        double h,
        double OmegaB,
        double OmegaCDM,
        double OmegaK,
        double Neff,
        double TCMB
        );

    // Print some useful info about the class
    void info() const;

    // Do all the solving
    void solve();
    void solve_t();
    // Output some results to file
    void output(const std::string filename) const;

    // Get functions that we must implement
    double eta_of_x(double x) const;
    double t_of_x(double x) const;
    double H_of_x(double x) const;
    double dHdx_of_x(double x) const;
    double ddHddx_of_x(double x) const;
    double Hp_of_x(double x) const;
    double dHpdx_of_x(double x) const;
    double ddHpddx_of_x(double x) const;

    double get_OmegaB(double x = 0.0) const;
    double get_OmegaM(double x = 0.0) const;
    double get_OmegaR(double x = 0.0) const;
    double get_OmegaRtot(double x = 0.0) const;
    double get_OmegaNu(double x = 0.0) const;
    double get_OmegaCDM(double x = 0.0) const;
    double get_OmegaLambda(double x = 0.0) const;
    double get_OmegaK(double x = 0.0) const;
    double get_OmegaMnu(double x = 0.0) const;

    double get_z(double x) const;
    double get_dA( double x=0) const;
    double get_dL( double x=0) const;
    double get_chi( double x=0) const;
    double get_H0() const;
    double get_h() const;
    double get_Neff() const;
    double get_TCMB(double x = 0.0) const;

    double eta_Hp_c(double x = 0.0) const;
    double t_of_0(double x = 0.0) const;
    double time_today(double x = 0.0) const;

};

#endif
